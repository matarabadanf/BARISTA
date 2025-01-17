#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import sys
import argparse
import pandas as pd

from misc.NXSpecReader import NXSpecReader

from typing import List, Dict, Optional
from numpy.typing import NDArray

# Parser is used to input via terminal the required arguments for Rodrigo
parser = argparse.ArgumentParser(
    description="""Generates the TDDFT absorption spectre from an ORCA TDDFT calculation. \
        It generates spectre .jpg image.\n \
        It can generate a peaks.dat and a spectra.dat that contains the peaks and spectre.""",
    epilog="""It should run appropriately with ORCA 5.0.X""",
)

parser.add_argument(
    "-F", type=str, required=True, default=42, help="Output ORCA file"
)

parser.add_argument(
    "-o",
    type=str,
    default=None,
    help="Output image filename. Default is output_file_spectre.jpg",
)

parser.add_argument(
    "-u", type=str, default="eV", help="Units, nm or eV. Default is eV."
)

parser.add_argument(
    "--interactive",
    type=bool,
    default=False,
    help="Interactive mode, will not save a .jpg file",
)

parser.add_argument(
    "-exp", type=str, default=None, help="Experimental data for comparison"
)

parser.add_argument(
    "-semi",
    type=str,
    default=None,
    help="Semiclassical nx log file",
)


# display help message if there is no arguments


class Rodrigo:
    """Class to plot and manage spectra"""

    def __init__(self, filename, units="eV"):
        self.filename = filename
        self.units = units

        self._initialize()

    def _initialize(self):
        self._load_file(self.filename)
        self._locate_peaks_in_file()
        self._extract_peaks()
        self._normalize_peaks()

    def _load_file(self, filename):

        with open(filename, "r") as fi:
            self.orca_output_content = fi.readlines()

    def _locate_peaks_in_file(self):

        for index, line in enumerate(self.orca_output_content):
            if (
                "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"
                in line
            ):
                self.init_index = index
            elif (
                "ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS"
                in line
            ):
                self.end_index = index

    def _extract_peaks(self):

        self.spectra_peaks = np.array(
            [
                np.array(
                    [
                        float(line.strip().split()[2]),
                        float(line.strip().split()[3]),
                    ]
                )
                for line in self.orca_output_content[
                    self.init_index + 5 : self.end_index - 2
                ]
            ]
        )

    def _normalize_peaks(self):
        self.spectra_peaks[:, 1] = self.spectra_peaks[:, 1] / max(
            self.spectra_peaks[:, 1]
        )

    def get_peaks_nm(self):
        return np.copy(self.spectra_peaks).T

    def get_peaks_eV(self):

        spectra_ev = np.copy(self.spectra_peaks).T
        spectra_ev[0] = 1239.8 / spectra_ev[0]

        return spectra_ev

    def _gaussian(
        self,
        X: NDArray[np.float64],
        height: float = 1,
        center: float = 15,
        spread: float = 0.05,
    ) -> NDArray[np.float64]:
        """
        Generate a Gaussian function.

        Parameters
        ----------
        X : NDArray[np.float64]
            X values for Gaussian calculation
        height : float, optional
            Peak height (default: 1)
        center : float, optional
            Peak center (default: 15)
        spread : float, optional
            Peak width (default: 0.05)

        Returns
        -------
        NDArray[np.float64]
            Gaussian function values
        """
        return height * np.exp(-((X - center) ** 2) / (2 * spread**2))

    def load_experimental(self, experimental_filename):
        experimental = np.loadtxt(experimental_filename, delimiter=",")
        experimental[:, 1] /= max(experimental[:, 1])

        if max(experimental[:, 0]) > 50:
            experimental[:, 0] = 1239.8 / experimental[:, 0]

        return experimental

    def _plot_vertical_nm(self):

        x, ymax = self.get_peaks_nm()
        ymin = np.zeros_like(x)

        plt.vlines(x, ymin, ymax, colors="black", label="Vertical spectrum")

    def _plot_vertical_eV(self):

        x, ymax = self.get_peaks_eV()
        ymin = np.zeros_like(x)

        plt.vlines(x, ymin, ymax, colors="black", label="Vertical spectrum")

    def plot_vertical(self):
        if self.units == "eV":
            self._plot_vertical_eV()

        elif self.units == "nm":
            self._plot_vertical_nm()
            plt.xlim(right=400)

    def plot_additional_spectra(
        self, data, units="eV", label="Imported Spectrum"
    ):
        if units != self.units:
            data[0] = 1239.8 / data[0]
        plt.plot(data[0], data[1], label=label)

    def plot(self):
        plt.legend()
        plt.show()

    def save_image(self, imagename):
        plt.legend()
        plt.savefig(imagename, dpi=300)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    # a = rodrigo(
    #     args.F,
    #     args.o,
    #     args.u,
    #     args.i,
    #     args.sp,
    #     args.cf,
    #     args.r,
    #     args.nopeaks,
    #     args.md,
    #     args.shift,
    #     args.exp,
    # )

    a = Rodrigo(args.F, args.u)

    a.plot_vertical()

    if args.semi is not None:
        nx_specturm = NXSpecReader(args.semi).return_semiclassical_spectrum()

        a.plot_additional_spectra(
            nx_specturm, units="eV", label="Semiclassical spectrum"
        )

    if args.exp is not None:
        a.plot_additional_spectra(
            a.load_experimental("beta-carboline_exp.csv").T,
            units="eV",
            label="Exp. spectrum",
        )

    if args.o is not None:
        a.save_image(args.o)
    else:
        a.save_image(args.F.replace(".in.out", ".jpg"))
