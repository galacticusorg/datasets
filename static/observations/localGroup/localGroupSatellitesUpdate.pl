#!/usr/bin/env perl
use strict;
use warnings;
use WWW::Curl::Easy;
use JSON::PP qw(encode_json decode_json);
use Unicode::Normalize;

# Read arguments.
die("Usage: localGroupSatellitesUpdate.pl <apiToken>")
    unless ( scalar(@ARGV) == 1 );
my $apiToken   = $ARGV[0];

# Construct a curl object.
my $curl = WWW::Curl::Easy->new();
$curl->setopt(CURLOPT_HEADER,1);
$curl->setopt(CURLOPT_HTTPHEADER, ['Authorization: Bearer:'.$apiToken]);

# Initialize a dictionary of known URL updates.
my %knownUpdatedURLs;

# Initialize a dictionary of bibcodes;
my %bibCodes;
{
    open(my $in ,"<","localGroupSatellites.xml"       );
    open(my $out,">","localGroupSatellites.xml.stage1");
    while ( my $line = <$in> ) {
	# Update URLs that point to an arXiv paper on NASA ADS.
	if ( $line =~ m/"(https:\/\/ui\.adsabs\.harvard\.edu\/abs\/\d+arXiv.*?)"/ ) {
	    my $oldURL = $1;
	    unless ( exists($knownUpdatedURLs{$oldURL}) ) {
		open(my $curl,"curl -Ls -o /dev/null -w %{url_effective} ".$oldURL."|");
		my $newURL = <$curl>;
		close($curl);
		chomp($newURL);
		$newURL =~ s/\/abstract//;
		$knownUpdatedURLs{$oldURL} = $newURL;
		print "Updating URL '".$oldURL."' to '".$newURL."'\n";
	    }
	    $line =~ s/$oldURL/$knownUpdatedURLs{$oldURL}/;
	}
	# Update URLs that point to an arXiv paper.
	if ( $line =~ m/"(https:\/\/arxiv\.org\/abs\/(\d+\.\d+))"/ ) {
	    my $originalURL = $1;
	    my $oldURL      = "https://ui.adsabs.harvard.edu/abs/arXiv:".$2;
	    unless ( exists($knownUpdatedURLs{$oldURL}) ) {
		open(my $curl,"curl -Ls -o /dev/null -w %{url_effective} ".$oldURL."|");
		my $newURL = <$curl>;
		close($curl);
		chomp($newURL);
		$newURL =~ s/\/abstract//;
		$knownUpdatedURLs{$oldURL} = $newURL;
		print "Updating URL '".$originalURL."' to '".$newURL."'\n";
	    }
	    $line =~ s/$originalURL/$knownUpdatedURLs{$oldURL}/;
	}
	# Update bibliographic references. 
	if ( $line =~ m/"https:\/\/ui\.adsabs\.harvard\.edu\/abs\/(.*?)"/ ) {
	    my $bibCode = $1;
	    $bibCodes{$bibCode} = "unknown";
	}
	# Output the updated line.
	print $out $line;
    }
    close($in);
    close($out);
}

# Get bibliographic records from NASA ADS.
my $records;
{
    my $countRecords = scalar(keys(%bibCodes));
    $curl->setopt(CURLOPT_URL, 'https://api.adsabs.harvard.edu/v1/search/bigquery?q=*:*&rows='.$countRecords.'&fl=bibcode,title,author,year,pub,volume,page');
    $curl->setopt(CURLOPT_HTTPHEADER, ['Authorization: Bearer:'.$apiToken,"Content-Type: big-query/csv"]);
    my $response_body;
    $curl->setopt(CURLOPT_WRITEDATA,\$response_body);
    $curl->setopt(CURLOPT_POST, 1);
    $curl->setopt(CURLOPT_POSTFIELDS, "bibcode\n".join("\n",keys(%bibCodes)));
    my $retcode = $curl->perform;
    if ($retcode == 0) {
	my $response_code = $curl->getinfo(CURLINFO_HTTP_CODE);
	if ( $response_code == 200 ) {
	    # Extract the JSON.
	    my $json;
	    my $startFound = 0;
	    open(my $response,"<",\$response_body);
	    while ( my $line = <$response> ) {
		$startFound = 1
		    if ( $line =~ m/^\{/ );
		$json .= $line
		    if ( $startFound );
	    }
	    close($response);
	    $records = decode_json($json);
	} else {
	    die("Failed to retrieve record identifiers: ".$response_code.$response_body);
	}
    } else {
	die("Failed to retrieve record identifiers: ".$retcode." ".$curl->strerror($retcode)." ".$curl->errbuf);
    }
}

# Journal abbreviations.
my %journalAbbr =
    (
     "The Astronomical Journal"                                => "AJ"          ,
     "The Astrophysical Journal"                               => "ApJ"         ,
     "Astronomy and Astrophysics"                              => "AA"          ,
     "Acta Astronomica"                                        => "Acta Astron.",
     "Monthly Notices of the Royal Astronomical Society"       => "MNRAS"       ,
     "arXiv e-prints"                                          => "arXiv"       ,
     "Research Notes of the American Astronomical Society"     => "RNAAS"       ,
     "Publications of the Astronomical Society of the Pacific" => "PASJ"        ,
     "Publications of the Astronomical Society of Japan"       => "PASP"
    );

# Construct bibliographic entries for each case.
foreach my $record ( @{$records->{'response'}->{'docs'}} ) {
    die("No journal abbreviation found for '".$record->{'pub'}."'")
	unless ( exists($journalAbbr{$record->{'pub'}}) );
    my @authors = map {$_ =~ s/([^,]+),.*/$1/; $_} @{$record->{'author'}};
    my $author;
    if      ( scalar(@authors) == 1 ) {
	$author = $authors[0];
    } elsif ( scalar(@authors) == 2 ) {
	$author = $authors[0]." &amp; ".$authors[1];
    } elsif ( scalar(@authors) == 3 ) {
	$author = $authors[0].", ".$authors[1]." &amp; ".$authors[2];
    } else {
	$author = $authors[0].", ".$authors[1].", ".$authors[2]." et al.";
    }
    $author .= "; ";
    $author  = NFKD($author);
    $author  =~ s/\p{NonspacingMark}//g;
    my $year    =                                            $record->{'year'  }      ."; "     ;
    my $volume  = exists($record->{'volume'}) ?              $record->{'volume'}      ."; " : "";
    my $journal =                               $journalAbbr{$record->{'pub'   }     }."; "     ;
    my $page    =                                            $record->{'page'  }->[0]           ;
    $bibCodes{$record->{'bibcode'}} =
	$author.$year.$journal.$volume.$page;
}

# Update references to standardized format.
{
    open(my $in ,"<","localGroupSatellites.xml.stage1");
    open(my $out,">","localGroupSatellites.xml.stage2");
    while ( my $line = <$in> ) {
	# Update references.
	if ( $line =~ m/"https:\/\/ui\.adsabs\.harvard\.edu\/abs\/(.*)"/ ) {
	    my $bibCode = $1;
	    if ( exists($bibCodes{$bibCode}) ) {
		$line =~ s/\sreference="[^"]+"/ reference="$bibCodes{$bibCode}"/g;
	    }
	}

	# Output the updated line.
	print $out $line;
    }
    close($in);
    close($out);
}

# Update the original file.
system("mv localGroupSatellites.xml.stage2 localGroupSatellites.xml");
unlink("localGroupSatellites.xml.stage1");

exit;
