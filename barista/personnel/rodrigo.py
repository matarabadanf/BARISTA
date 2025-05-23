#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import sys
import argparse
import pandas as pd
from functools import cached_property

from barista.personnel.misc.NXSpecReader import NXSpecReader

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
    "--gaussian",
    action='store_true',
    required=False,
    help="Include gaussian wrap for peaks.",
)

parser.add_argument(
    "-gauss_disp", 
    type=float, 
    default=None, 
    required=False,
    help='Gaussian dispersion in the current energy unit.',
)

parser.add_argument(
    "-exp", type=str, default=None, help="Experimental data for comparison."
)

parser.add_argument(
    "-exp_units", 
    type=str, 
    default='eV', 
    required=False,
    help='Experimental data energy units. Default is eV.',
)

parser.add_argument(
    "-exp_shift", 
    type=float, 
    default=0.0, 
    required=False,
    help='Experimental spectrum shift in current energy units. Default is 0.',
)

parser.add_argument(
    "-semi",
    type=str,
    default=None,
    help="Semiclassical nx log file.",
)

parser.add_argument(
    "--interactive",
    action='store_true',
    help="Interactive mode, will not save a .jpg file.",
)


class Rodrigo:

    def __init__(self, filename, units="eV"):
        self.filename = filename
        self.units = units

        self._initialize()

    def _initialize(self):
        self._load_file(self.filename)
        self._normalize_peaks()

    def _load_file(self, filename):
        with open(filename, "r") as fi:
            self.orca_output_content = fi.readlines()

    @cached_property
    def spectra_peaks(self):
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
        s_p = np.array(
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

        return np.copy(s_p)
        


    def _normalize_peaks(self):
        self.spectra_peaks[:, 1] = self.spectra_peaks[:, 1] / max(
            self.spectra_peaks[:, 1]
        )

    @cached_property
    def get_peaks_nm(self):
        return np.copy(self.spectra_peaks).T

    @cached_property
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

    def _plot_gaussian_nm(self, dispersion: float = 7.5):

        x_array = np.linspace(
            min(self.get_peaks_nm[0]) - 20,
            max(self.get_peaks_nm[0]) + 20,
            1000,
        )


        y_array = np.zeros(len(x_array))

        for peak in self.get_peaks_nm.T:
            y_array += self._gaussian(x_array, peak[1], peak[0], dispersion)

        plt.plot(x_array, y_array, label='Convoluted spectrum')
    
    def export_gaussian_eV(self, dispersion: float = 0.2):
        x_array = np.linspace(0, 15, 1000)
        y_array = np.zeros(len(x_array))

        for peak in self.get_peaks_eV.T:
            y_array += self._gaussian(x_array, peak[1], peak[0], dispersion)

        array = np.array([x_array,y_array])
        return array

    def _plot_gaussian_eV(self, dispersion: float = 0.2):

        x_array = np.linspace(
            min(self.get_peaks_eV[0]) - 1,
            max(self.get_peaks_eV[0]) + 1,
            1000,
        )

        y_array = np.zeros(len(x_array))

        for peak in self.get_peaks_eV.T:
            y_array += self._gaussian(x_array, peak[1], peak[0], dispersion)

        plt.plot(x_array, y_array/max(y_array), label='Convoluted spectrum')

    def plot_gaussian(self, dispersion=None):
        if self.units == "eV":
            if dispersion is not None:
                self._plot_gaussian_eV(dispersion)
            else:
                self._plot_gaussian_eV()

        elif self.units == "nm":
            if dispersion is not None:
                self._plot_gaussian_nm(dispersion)
            else:
                self._plot_gaussian_nm()

    def load_experimental(self, experimental_filename, units='eV'):
        experimental = np.loadtxt(experimental_filename, delimiter=",")
        experimental[:, 1] /= max(experimental[:, 1])
        
        if self.units != units and units == 'eV':
            experimental[:, 0] = 1239.8 / experimental[:, 0]

        return experimental

    def _plot_vertical_nm(self):

        x, ymax = self.get_peaks_nm
        ymin = np.zeros_like(x)

        plt.vlines(x, ymin, ymax, colors="black", label="Vertical spectrum")

    def _plot_vertical_eV(self):

        x, ymax = self.get_peaks_eV
        ymin = np.zeros_like(x)

        plt.vlines(x, ymin, ymax, colors="black", label="Vertical spectrum")

    def plot_vertical(self):
        if self.units == "eV":
            self._plot_vertical_eV()

        elif self.units == "nm":
            self._plot_vertical_nm()

    def plot_additional_spectra(
        self, data, units="eV", label="Imported Spectrum", shift=0.0
    ):
        if units != self.units:
            data[0] = 1239.8 / data[0]
        data[1] /= max(data[1])
        plt.plot(data[0]+shift, data[1], label=label)

    def plot(self):
        plt.xlabel(f'Energy / {self.units}')
        plt.ylabel('Norm. osc. str. / arb. u.')
        plt.legend()
        plt.show()

    def save_image(self, imagename):
        plt.xlabel(f'Energy / {self.units}')
        plt.ylabel('Norm. osc. str. / arb. u.')
        plt.legend()
        plt.savefig(imagename, dpi=300)

    
    def save_data(self, name:str='None', dispersion: float = 0.2):
        np.savetxt(name.replace('.in.out', '_peaks.dat'), self.get_peaks_eV.T)

        np.savetxt(
            name.replace('.in.out', '_convolved.dat'), 
            self.export_gaussian_eV(dispersion = 0.2).T
        )


if __name__ == "__main__":
    
    # display help message if there is no arguments
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    a = Rodrigo(args.F, args.u)
    # a = Rodrigo("tests/xanthine_spectra_tda.in.out", "nm")

    a.save_data(name=args.F)
    print(a.get_peaks_eV)

    a.plot_vertical()
    
    if args.gaussian:
        a.plot_gaussian(args.gauss_disp)

    if args.semi is not None:
        nx_specturm = NXSpecReader(args.semi).return_semiclassical_spectrum()

        a.plot_additional_spectra(
            nx_specturm, units="eV", label="Semiclassical spectrum"
        )

    if args.exp is not None:
        label = 'Experimental spectrum'
        if args.exp_shift != 0.0:
            label = f'Exp. spectrum shifted {args.exp_shift:.2f} {args.u}'
        a.plot_additional_spectra(
            a.load_experimental(args.exp, units=args.exp_units).T,
            units=args.exp_units,
            label=label,
            shift=args.exp_shift
        )
    
    if args.interactive:
        a.plot()
    elif args.o is not None:
        a.save_image(args.o)
    else:
        a.save_image(args.F.replace(".in.out", ".jpg"))



