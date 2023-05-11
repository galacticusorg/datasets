import xmlschema

my_schema = xmlschema.XMLSchema11('localGroupSatellites.xsd')
my_schema.validate('localGroupSatellites.xml')
