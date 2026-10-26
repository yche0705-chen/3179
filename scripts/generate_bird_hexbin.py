#!/usr/bin/env python3
"""Regenerate half-card hexbin centres from the 706-species source data."""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data" / "bird_mass_longevity.csv"
OUTPUT = ROOT / "data" / "bird_hexbin_counts.csv"

PLOT_WIDTH = 400.0
PLOT_HEIGHT = 300.0
HEX_RADIUS = 8.5
HEX_HALF_WIDTH = math.sqrt(3.0) * HEX_RADIUS / 2.0
X_DOMAIN = (4.0, 3000.0)
Y_DOMAIN = (0.0, 90.0)

# Keep the full hexagon and its hover stroke inside both quantitative axes.
X_INSET = HEX_HALF_WIDTH + 1.0
Y_INSET = HEX_RADIUS + 1.0


def chart_x(mass_g: float) -> float:
    ratio = math.log(mass_g / X_DOMAIN[0]) / math.log(X_DOMAIN[1] / X_DOMAIN[0])
    return X_INSET + ratio * (PLOT_WIDTH - 2.0 * X_INSET)


def chart_y(longevity: float) -> float:
    ratio = (longevity - Y_DOMAIN[0]) / (Y_DOMAIN[1] - Y_DOMAIN[0])
    return PLOT_HEIGHT - Y_INSET - ratio * (PLOT_HEIGHT - 2.0 * Y_INSET)


def data_mass(pixel_x: float) -> float:
    ratio = pixel_x / PLOT_WIDTH
    return X_DOMAIN[0] * math.exp(ratio * math.log(X_DOMAIN[1] / X_DOMAIN[0]))


def data_longevity(pixel_y: float) -> float:
    ratio = 1.0 - pixel_y / PLOT_HEIGHT
    return Y_DOMAIN[0] + ratio * (Y_DOMAIN[1] - Y_DOMAIN[0])


def make_centres() -> list[tuple[float, float]]:
    centres: list[tuple[float, float]] = []
    row = 0
    pixel_y = Y_INSET
    horizontal_step = 2.0 * HEX_HALF_WIDTH
    vertical_step = 1.5 * HEX_RADIUS

    while pixel_y <= PLOT_HEIGHT - Y_INSET:
        pixel_x = X_INSET + (HEX_HALF_WIDTH if row % 2 else 0.0)
        while pixel_x <= PLOT_WIDTH - X_INSET:
            centres.append((pixel_x, pixel_y))
            pixel_x += horizontal_step
        row += 1
        pixel_y = Y_INSET + row * vertical_step

    return centres


def nearest_centre(point_x: float, point_y: float, centres: list[tuple[float, float]]) -> int:
    return min(
        range(len(centres)),
        key=lambda index: (centres[index][0] - point_x) ** 2 + (centres[index][1] - point_y) ** 2,
    )


def main() -> None:
    centres = make_centres()
    bins: dict[str, dict[int, list[str]]] = {
        "All birds": defaultdict(list),
        "Parrots": defaultdict(list),
        "Passerines": defaultdict(list),
    }

    with SOURCE.open(newline="", encoding="utf-8") as source_file:
        records = list(csv.DictReader(source_file))

    for record in records:
        centre_index = nearest_centre(
            chart_x(float(record["adult_mass_g"])),
            chart_y(float(record["maximum_longevity_years"])),
            centres,
        )
        bins["All birds"][centre_index].append(record["species"])
        bins[record["group"]][centre_index].append(record["species"])

    with OUTPUT.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(
            output_file,
            fieldnames=["group", "mass_g", "longevity_years", "count", "species_names"],
        )
        writer.writeheader()

        for group in ("All birds", "Parrots", "Passerines"):
            for centre_index in sorted(bins[group]):
                pixel_x, pixel_y = centres[centre_index]
                names = sorted(bins[group][centre_index])
                writer.writerow(
                    {
                        "group": group,
                        "mass_g": f"{data_mass(pixel_x):.12g}",
                        "longevity_years": f"{data_longevity(pixel_y):.12g}",
                        "count": len(names),
                        "species_names": "; ".join(names),
                    }
                )

    print(f"Wrote {OUTPUT.relative_to(ROOT)} from {len(records)} records using {len(centres)} centres")


if __name__ == "__main__":
    main()
