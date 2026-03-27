#!/usr/bin/env python3
"""Generate random points within the angular mask of the UKIDSS UDS survey.

Construct a set of random points that lie within the angular mask of the UKIDSS UDS sample
used by Caputi et al. (2011; http://adsabs.harvard.edu/abs/2011MNRAS.413..162C). The mask
is defined by a set of boundaries, plus a set of rectangular cut-outs. The boundaries are
hard-coded into this script. The rectangles required are found by randomly placing rectangles
of an initial minimum size into the field and growing them to be as large as possible while
not containing any galaxies. Once enough rectangles have been placed, random points are
generated and rejected if they lie outside of the boundaries or inside any cut-out rectangle.
"""

import os
from datetime import datetime, timezone

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np


def rectangle_is_acceptable(x, y, size_x, size_y, ra, dec, boundaries, rectangles):
    """Return True if the rectangle is within survey bounds and contains no galaxies."""
    # Reject if any galaxies fall inside the rectangle.
    if np.any((np.abs(ra - x) < 0.5 * size_x) & (np.abs(dec - y) < 0.5 * size_y)):
        return False

    # Check boundary constraints (using center point only, matching original behaviour).
    for b in boundaries:
        if b["type"] in ("above", "below"):
            y_limit = y + 0.5 * size_y if b["type"] == "above" else y - 0.5 * size_y
            boundary_y = (
                (b["y2"] - b["y1"]) * (x - b["x1"]) / (b["x2"] - b["x1"]) + b["y1"]
            )
            if b["type"] == "above" and boundary_y < y_limit:
                return False
            if b["type"] == "below" and boundary_y > y_limit:
                return False
        elif b["type"] in ("left", "right"):
            x_limit = x + 0.5 * size_x if b["type"] == "right" else x - 0.5 * size_x
            boundary_x = (
                (b["x2"] - b["x1"]) * (y - b["y1"]) / (b["y2"] - b["y1"]) + b["x1"]
            )
            if b["type"] == "right" and boundary_x < x_limit:
                return False
            if b["type"] == "left" and boundary_x > x_limit:
                return False

    # Reject if center of new rectangle falls inside an existing rectangle.
    for r in rectangles:
        if (
            r["x"] + 0.5 * r["size_x"] > x
            and r["x"] - 0.5 * r["size_x"] < x
            and r["y"] + 0.5 * r["size_y"] > y
            and r["y"] - 0.5 * r["size_y"] < y
        ):
            return False

    return True


def points_are_acceptable(p_ra, p_dec, boundaries, rectangles):
    """Return indices of points inside survey boundaries and outside exclusion rectangles."""
    mask = np.ones(len(p_ra), dtype=bool)

    for b in boundaries:
        dx = b["x2"] - b["x1"]
        dy = b["y2"] - b["y1"]
        if b["type"] in ("above", "below"):
            boundary_y = dy * (p_ra - b["x1"]) / dx + b["y1"]
            if b["type"] == "above":
                mask &= p_dec < boundary_y
            else:
                mask &= p_dec > boundary_y
        elif b["type"] in ("left", "right"):
            boundary_x = dx * (p_dec - b["y1"]) / dy + b["x1"]
            if b["type"] == "right":
                mask &= p_ra < boundary_x
            else:
                mask &= p_ra > boundary_x

    for r in rectangles:
        inside = (
            (p_ra >= r["x"] - 0.5 * r["size_x"])
            & (p_ra <= r["x"] + 0.5 * r["size_x"])
            & (p_dec >= r["y"] - 0.5 * r["size_y"])
            & (p_dec <= r["y"] + 0.5 * r["size_y"])
        )
        mask &= ~inside

    return np.where(mask)[0]


def main():
    work_directory = os.environ["GALACTICUS_DATA_PATH"] + "/surveys/UKIDSS_UDS/"
    data_file = (
        os.environ["GALACTICUS_DATA_PATH"]
        + "/surveyGeometry/UKIDDS_UDS/surveyGeometry.txt"
    )

    # Read galaxy positions (skip comment lines starting with #).
    data = np.loadtxt(data_file, comments="#")
    ra = data[:, 0]
    dec = data[:, 1]

    rectangle_size_minimum = 0.01

    # Survey boundary definitions: two points on a line plus a direction indicating
    # which side of the line is excluded.
    boundaries = [
        {"x1": 34.460, "y1": -4.648, "x2": 34.750, "y2": -4.648, "type": "above"},
        {"x1": 34.000, "y1": -4.908, "x2": 34.460, "y2": -4.648, "type": "above"},
        {"x1": 34.900, "y1": -4.920, "x2": 34.750, "y2": -4.648, "type": "above"},
        {"x1": 34.440, "y1": -5.518, "x2": 34.905, "y2": -5.255, "type": "below"},
        {"x1": 34.210, "y1": -5.518, "x2": 34.440, "y2": -5.518, "type": "below"},
        {"x1": 34.210, "y1": -5.518, "x2": 34.000, "y2": -5.150, "type": "below"},
        {"x1": 34.000, "y1": -4.908, "x2": 34.000, "y2": -5.150, "type": "left"},
        {"x1": 34.900, "y1": -4.920, "x2": 34.905, "y2": -5.255, "type": "right"},
    ]

    # Place exclusion rectangles until 100 consecutive placement attempts fail.
    rectangles = []
    fail_count = 0
    fail_maximum = 100

    while fail_count < fail_maximum:
        size_x = rectangle_size_minimum
        size_y = rectangle_size_minimum

        # Try up to 1000 random positions to find a valid initial placement.
        good = False
        for _ in range(1000):
            x = np.random.uniform(0, 1) * 0.9 + 34.0
            y = np.random.uniform(0, 1) * 1.0 - 5.6
            if rectangle_is_acceptable(x, y, size_x, size_y, ra, dec, boundaries, rectangles):
                good = True
                break

        if good:
            # Grow the rectangle on all four sides until none can expand further.
            max_size_reached = 0
            grow_factor = 0.01
            while max_size_reached < 4:
                max_size_reached = 0

                # Grow left edge.
                x -= 0.5 * grow_factor * size_x
                size_x *= 1.0 + grow_factor
                if not rectangle_is_acceptable(x, y, size_x, size_y, ra, dec, boundaries, rectangles):
                    max_size_reached += 1
                    size_x /= 1.0 + grow_factor
                    x += 0.5 * grow_factor * size_x

                # Grow right edge.
                x += 0.5 * grow_factor * size_x
                size_x *= 1.0 + grow_factor
                if not rectangle_is_acceptable(x, y, size_x, size_y, ra, dec, boundaries, rectangles):
                    max_size_reached += 1
                    size_x /= 1.0 + grow_factor
                    x -= 0.5 * grow_factor * size_x

                # Grow bottom edge.
                y -= 0.5 * grow_factor * size_y
                size_y *= 1.0 + grow_factor
                if not rectangle_is_acceptable(x, y, size_x, size_y, ra, dec, boundaries, rectangles):
                    max_size_reached += 1
                    size_y /= 1.0 + grow_factor
                    y += 0.5 * grow_factor * size_y

                # Grow top edge.
                y += 0.5 * grow_factor * size_y
                size_y *= 1.0 + grow_factor
                if not rectangle_is_acceptable(x, y, size_x, size_y, ra, dec, boundaries, rectangles):
                    max_size_reached += 1
                    size_y /= 1.0 + grow_factor
                    y -= 0.5 * grow_factor * size_y

            rectangles.append({"x": x, "y": y, "size_x": size_x, "size_y": size_y})
            print(
                f"Placed rectangle {len(rectangles)} at ({x},{y}) "
                f"with size {size_x}x{size_y}"
            )
        else:
            fail_count += 1

    # Generate 1 million random points uniformly distributed on the sphere within
    # the bounding box of the survey.
    pi = np.pi
    cos_theta_min = np.cos(1.6517967)
    cos_theta_max = np.cos(1.6675065)
    phi_min = 0.5931317
    phi_max = 0.60971425

    random_count = 1_000_000
    phi = np.random.uniform(0, 1, random_count) * (phi_max - phi_min) + phi_min
    cos_theta = (
        np.random.uniform(0, 1, random_count) * (cos_theta_max - cos_theta_min)
        + cos_theta_min
    )
    theta = np.arccos(cos_theta)
    random_ra = phi * 180.0 / pi
    random_dec = 90.0 - theta * 180.0 / pi

    # Filter to points within the survey mask.
    random_accept = points_are_acceptable(random_ra, random_dec, boundaries, rectangles)

    # Monte Carlo estimate of the survey solid angle.
    solid_angle_total = (phi_max - phi_min) * (cos_theta_min - cos_theta_max)
    solid_angle_survey = solid_angle_total * len(random_accept) / random_count
    print(f"Survey solid angle is: {solid_angle_survey}")

    # Rotate coordinates so that the survey center becomes the new pole.
    theta0 = 95.05 * pi / 180.0
    phi0 = 34.45 * pi / 180.0
    c1 = np.array([
        np.sin(theta0) * np.cos(phi0),
        np.sin(theta0) * np.sin(phi0),
        np.cos(theta0),
    ])
    c2 = np.array([np.sin(phi0), -np.cos(phi0), 0.0])
    c3 = np.cross(c1, c2)

    # Cartesian coordinates of all random points, shape (3, N).
    x_cart = np.array([
        np.sin(theta) * np.cos(phi),
        np.sin(theta) * np.sin(phi),
        np.cos(theta),
    ])

    x1 = c1 @ x_cart  # dot product → (N,)
    x2 = c2 @ x_cart
    x3 = c3 @ x_cart

    theta_prime = np.arccos(x1)
    phi_prime = np.arctan2(x3, x2)

    # Write accepted points to HDF5.
    os.makedirs(work_directory + "data", exist_ok=True)
    with h5py.File(work_directory + "data/surveyGeometryRandoms.hdf5", "w") as hdf:
        hdf.create_dataset("theta", data=theta_prime[random_accept])
        hdf.create_dataset("phi", data=phi_prime[random_accept])
        now = datetime.now(tz=timezone.utc).astimezone()
        ms = now.strftime("%f")[:3]
        tz_offset = now.strftime("%z")
        tz_formatted = tz_offset[:-2] + ":" + tz_offset[-2:]
        timestamp = now.strftime("%Y-%m-%dT%H:%M:%S") + "." + ms + tz_formatted
        hdf.attrs["createdBy"] = (
            "Galacticus; constraints/dataAnalysis/"
            "stellarMassFunctions_UKIDSS_UDS_z3_5/surveyGeometryRandoms.py"
        )
        hdf.attrs["description"] = (
            "Random points in survey geometry of UKIDSS UDS for Caputi et al. "
            "(2011; http://adsabs.harvard.edu/abs/2011MNRAS.413..162C)"
        )
        hdf.attrs["timestamp"] = timestamp

    # Plot the survey mask, exclusion rectangles, galaxies, and sample random points.
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_title("Survey Mask for Caputi et al. (2011) Sample")
    ax.set_xlabel("right ascension [°]")
    ax.set_ylabel("declination [°]")
    ax.set_xlim(33.98, 34.93)
    ax.set_ylim(-5.53, -4.63)

    for b in boundaries:
        ax.plot([b["x1"], b["x2"]], [b["y1"], b["y2"]], color="indianred", linewidth=2)

    for r in rectangles:
        rect_patch = patches.Rectangle(
            (r["x"] - 0.5 * r["size_x"], r["y"] - 0.5 * r["size_y"]),
            r["size_x"],
            r["size_y"],
            linewidth=1,
            edgecolor="cornflowerblue",
            facecolor="cornflowerblue",
            alpha=0.5,
        )
        ax.add_patch(rect_patch)

    ax.scatter(ra, dec, s=1, color="mediumseagreen", label="Galaxies", zorder=3)

    sample_mask = np.random.uniform(0, 1, len(random_accept)) <= 0.1
    sample_indices = random_accept[sample_mask]
    ax.scatter(
        random_ra[sample_indices],
        random_dec[sample_indices],
        s=1,
        color="orange",
        label="Random points",
        zorder=2,
    )

    ax.legend(loc="upper left", markerscale=5)
    plot_file = work_directory + "surveyMask.pdf"
    fig.savefig(plot_file, bbox_inches="tight")
    plt.close(fig)
    print(f"Plot saved to {plot_file}")


if __name__ == "__main__":
    main()
