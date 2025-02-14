#!/usr/bin/env python3
import argparse
import os
import numpy as np
import re
import matplotlib.pyplot as plt
import pandas as pd
import sys
from laura import Laura
from jeremy import Jeremy


# Parser is used to input via terminal the required arguments for Emma
parser = argparse.ArgumentParser(
    description="""Takes the ground state geometry of a molecule and compares it to optimized TDDFT geometries run in ORCA.\
        \nIt compares the energy, RMSD and final root, plotting the results in an image. """,
    epilog="""It should run appropriately with ORCA 5.0.X""",
)

parser.add_argument(
    "-fc",
    type=str,
    required=True,
    default=42,
    help="Optimized ground state .xyz file",
)
parser.add_argument(
    "-ex",
    type=str,
    required=True,
    nargs="+",
    default=42,
    help='Excited state optimization .out "filenames" (it is important that the name ends with *_N.in.out)',
)
parser.add_argument(
    "-O",
    type=str,
    default=42,
    help='Output file extension (default=".in.out")',
)
parser.add_argument("-o", type=str, default=42, help="Output image filename")

parser.add_argument(
    "-r", type=str, default=42, help="Generate a report table in output file"
)

parser.add_argument(
    "-md",
    type=str,
    default=True,
    help="Generate report in MD format (default is True)",
)
parser.add_argument(
    "-i",
    type=bool,
    default=False,
    help="Interactive plot mode (default is False)",
)

# display help message if there is no arguments


class Emma:

    def __init__(
        self,
        fc_filename: str,
        excited_opt_templatess: list[str],
        WorkDir: str = ".",
    ):
        self._fc_filename = fc_filename
        self._excited_opt_templates = excited_opt_templatess
        self._maxroot = 0
        self._WorkDir = WorkDir
        self._geometry_files = []

        self._determine_targets()
        self._generate_geometries()
        self._generate_geometry_objects()

    def _determine_targets(self):
        self._targets = []
        for template in self._excited_opt_templates:
            startname = "_".join(template.strip().split("_")[:-1]).replace(
                self._WorkDir, ""
            )
            for file in os.listdir(self._WorkDir):
                if startname in file and ".xyz" not in file:
                    self._targets.append(self._WorkDir + file)

    def _generate_geometries(self):
        for target in self._targets:
            print(target)
            l = Laura.from_file(target)
            l.generate_xyzfile(target.replace(".in.out", ".xyz"))
            self._geometry_files.append(target.replace(".in.out", ".xyz"))

    def _generate_geometry_objects(self):
        self._fc_obj = Jeremy(self._fc_filename)

        self._target_objects = []

        for geom in self._geometry_files:
            J = Jeremy(geom)
            J.override_connectivity_matrix(self._fc_obj.connectivity_matrix)

            self._target_objects.append(J)

    def internal_deviation(self):
        for geom in self._target_objects:
            geom.rmsd = np.mean(
                (self._fc_obj.internal_values - geom.internal_values) ** 2
            )
        for geom in self._target_objects:
            print(geom.xyzfile, geom.rmsd)
            print("Max deviation in coordinate")


if __name__ == "__main__":
    # if len(sys.argv) == 1:
    #     parser.print_help(sys.stderr)
    #     sys.exit(1)

    # args = parser.parse_args()

    # a = emma(args.ex, args.fc, args.O, args.r, args.md, args.i, args.o)
    # print(a)

    e = Emma(
        "upgrades/xanthine_FC.xyz",
        [
            "upgrades/xanthine_opt_followiroot_N.in.out",
            "upgrades/xanthine_opt_iroot_N.in.out",
        ],
        WorkDir="upgrades/",
    )

