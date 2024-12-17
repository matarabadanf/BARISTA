#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import sys
import argparse
import pandas as pd


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
    default=False,
    help="Output image filename. Default is output_file_spectre.jpg",
)
parser.add_argument(
    "-r",
    type=str,
    default=42,
    help='Name for the _peaks.dat and _spectra.dat. If "D" is given,\
    it will choose the ORCA calculation file as default orca_filename_spectra.dat and orca_filename_peaks.dat',
)
parser.add_argument(
    "-u", type=str, default="eV", help="Units, nm or eV. Default is nm."
)
parser.add_argument(
    "-cf",
    type=float,
    default=-1,
    help="Max peak energy cuttoff. Default is nm.",
)
parser.add_argument(
    "-sp",
    type=float,
    default=-1,
    help="Gaussian envelope spread. Give in the requeted units (nm or eV)",
)
parser.add_argument(
    "-nopeaks",
    type=bool,
    default=False,
    help="Remove the peaks in the spectra .jpg file",
)
parser.add_argument(
    "-i",
    type=bool,
    default=False,
    help="Interactive mode, will not save a .jpg file",
)
parser.add_argument(
    "-md", type=bool, default=False, help="Request Markdown table format"
)
parser.add_argument(
    "-exp", type=str, default="", help="Experimental data for comparison"
)
parser.add_argument(
    "-shift", type=float, default=0, help="Shift in nm of the spectrum"
)




# display help message if there is no arguments


# Gaussian function to hug the peaks
def gaussian(X, a, b, c):
    return a * np.exp(-((X - b) ** 2) / (2 * c**2))


def rodrigo(
    filename: str,
    image_filename: bool = False,
    units: str = "eV",
    interactive: bool = False,
    envelope_spread: float = -1,
    energy_cutoff: float = -1,
    report_filename: str = 42,
    nopeaks: bool = False,
    request_markdown: bool = False,
    experimental_spectra_shift: float = 0,
    experimental_filename: str = "",
) -> pd.DataFrame:
    # Parsing the output file to obtain the peaks
    fi = open(filename, "r")
    cont = fi.readlines()
    
    print(units)

    for line in cont:
        if (
            "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"
            in line
        ):
            init_index = cont.index(line)
        elif (
            "ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS"
            in line
        ):
            end_index = cont.index(line)

    spectra_peaks = np.array(
        [
            np.array(
                [
                    float(line.strip().split()[2]),
                    float(line.strip().split()[3]),
                ]
            )
            for line in cont[init_index + 5 : end_index - 2]
        ]
    )
    spectra_peaks[:, 1] = spectra_peaks[:, 1] / max(spectra_peaks[:, 1])

    # Generate the grid depending on the units and set the spread for the gaussian function
    if units == "nm":
        X = np.linspace(100, 450, 1000)
    elif units == "eV":
        X = np.linspace(0, 15, 1000)

    if envelope_spread == -1:
        if units == "nm":
            spread = 10
        else:
            spread = 0.5
    else:
        spread = envelope_spread

    # gets the peaks and constrains the cutoff condition
    if energy_cutoff != -1:
        if str(units) == "nm":
            spectra_peaks = np.array(
                [
                    np.array(
                        [
                            float(line.strip().split()[2]),
                            float(line.strip().split()[3]),
                        ]
                    )
                    for line in cont[init_index + 5 : end_index - 2]
                    if float(line.strip().split()[2]) > energy_cutoff
                ]
            )
            print(spectra_peaks)
        elif str(units) == "eV":
            spectra_peaks = np.array(
                [
                    np.array(
                        [
                            1239.8 / float(line.strip().split()[2]),
                            float(line.strip().split()[3]),
                        ]
                    )
                    for line in cont[init_index + 5 : end_index - 2]
                    if 1239.8 / float(line.strip().split()[2]) < energy_cutoff
                ]
            )
    else:
        if units == "nm":
            spectra_peaks = np.array(
                [
                    np.array(
                        [
                            float(line.strip().split()[2]),
                            float(line.strip().split()[3]),
                        ]
                    )
                    for line in cont[init_index + 5 : end_index - 2]
                ]
            )
        elif units == "eV":
            spectra_peaks = np.array(
                [
                    np.array(
                        [
                            1239.8 / float(line.strip().split()[2]),
                            float(line.strip().split()[3]),
                        ]
                    )
                    for line in cont[init_index + 5 : end_index - 2]
                ]
            )

    spectra_peaks[:, 1] = spectra_peaks[:, 1] / max(spectra_peaks[:, 1])

    # Generates the spectra array
    final_spectra = np.zeros(1000)
    for peak in spectra_peaks:
        final_spectra += gaussian(X, peak[1], peak[0], spread)
        if nopeaks == False:
            plt.vlines(x=peak[0], ymin=0, ymax=peak[1], colors="plum", label='Oscillator strength')

    print(spectra_peaks)

    plt.plot(X, final_spectra / max(final_spectra), c="yellowgreen")

    name = filename.replace(".out", "").replace(".in", "")

    plt.title("Absorption Spectra\n" + name)
    plt.xlabel("Energy / %s" % units)
    plt.ylabel("Norm. Osc. strenght")

    if report_filename != 42:
        if report_filename == "D":
            name = filename.replace(".out", "").replace(".in", "")
        else:
            name = report_filename

        if request_markdown == True:
            peaksfile = open(name + "_peaks.md", "w")
            peaksfile.write(
                "|State|Wavelength|Energy / ev | Norm. Osc. Str.|\n|-----|-----|-----|-----|\n"
            )
            for i, peak in enumerate(spectra_peaks):
                peaksfile.write(
                    "|S%i|%.1f|%.5f|%.5f|\n"
                    % (i + 1, peak[0], 1239.8 / peak[0], peak[1])
                )
            peaksfile.close()

        else:
            peaksfile = open(name + "_peaks.dat", "w")
            for peak in spectra_peaks:
                peaksfile.write("%.1f %.5f\n" % (peak[0], peak[1]))
            peaksfile.close()
            spectre_file = open(name + "_spectre.dat", "w")
            for point in final_spectra:
                spectre_file.write("%.1f %.5f\n" % (peak[0], peak[1]))
            spectre_file.close()

    if experimental_filename != "":
        experimental = np.loadtxt(experimental_filename, delimiter=",")
        experimental[:, 1] /= max(experimental[:, 1])
        experimental[:, 0] -= experimental_spectra_shift
        if experimental_spectra_shift > 0 and experimental_spectra_shift != 0:
            lab = "Experimental, shifted -%.2f nm" % experimental_spectra_shift

        plt.plot(experimental[:, 0], experimental[:, 1], label=lab)

    plt.legend()
    

    if interactive == True:
        plt.show()
    elif image_filename:
        name = filename.replace(".out", "").replace(".in", "")
        plt.savefig(name + "_spectre.jpg", dpi=800)
        
    df = pd.DataFrame({
        'Final state' : [i for i in range(len(spectra_peaks))],
        'Energy' : spectra_peaks[:,0],
        'Osc. Strength' : spectra_peaks[:,1],
        }
    )
        
    return df



if __name__ == "__main__":
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    

    a = rodrigo(
        args.F,
        args.o,
        args.u,
        args.i,
        args.sp,
        args.cf,
        args.r,
        args.nopeaks,
        args.md,
        args.shift,
        args.exp,
    )
    
    print(a)
