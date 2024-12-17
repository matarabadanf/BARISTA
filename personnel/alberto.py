#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys

# Parser is used to input via terminal the required arguments for Emma
parser = argparse.ArgumentParser(
    description="""Takes a excited state optimization in ORCA and plots the energy of each state at each step of the opitmizaiton and the actual root.""",
    epilog="""It should run appropriately with ORCA 5.0.X""",
)

parser.add_argument(
    "-f", type=str, required=True, default=42, help="ORCA optimization file"
)
parser.add_argument("-o", type=str, default=True, help="Output image filename")
parser.add_argument(
    "-en",
    type=float,
    default=0,
    help="Reference energy of the ground state optimized energy",
)
parser.add_argument("-u", type=str, default="eV", help="Energy units")

# display help message if there is no arguments


def alberto(
    filename: str,
    reference_energy: float = 0,
    output_image: str = True,
    units: str = "eV",
) -> np.array:
    # get the file and content

    file = open(filename, "r")
    cont = file.readlines()
    # grep the lines of interest
    counter = 0
    n_os = []

    for line in cont:
        if "TD-DFT/TDA EXCITED STATES (SINGLETS)" in line:
            counter += 1
            if counter >= 2:
                break
        elif (
            "STATE" in line
            and "TD-DFT/TDA EXCITED STATES (SINGLETS)" not in line
            and "EXCITED STATE GRADIENT DONE" not in line
        ):
            n_os.append(int(line.strip().split()[1].replace(":", "")))

    total_list = [[] for i in range(0, max(n_os))]
    cis_energies = []

    for line in cont:
        if (
            "STATE" in line
            and "TD-DFT/TDA EXCITED STATES (SINGLETS)" not in line
            and "EXCITED STATE GRADIENT DONE" not in line
        ):
            index = int(line.strip().split()[1].replace(":", "")) - 1
            energy = float(line.strip().split()[3].replace(":", ""))
            total_list[index].append(energy)
        elif "E(SCF)" in line:
            cis_energies.append(float(line.strip().split()[2]))

    # print(total_list)

    cis_array = np.array(cis_energies) - reference_energy
    total_arrays = [np.array(l) for l in total_list]

    if units == "eV":
        cis_array *= 27.2114
        for i in total_arrays:
            i *= 27.2114
    else:
        units = "Hartree"

    x = np.arange(1, len(total_list[0]) + 1, 1)

    plt.plot(x, cis_array, label="Ground state")

    # obtain the state of the actual root
    curr_energy = []
    curr_energy_index = []
    for line in cont:
        if "DE(CIS) =" in line:
            curr_energy_index.append(
                float(line.strip().split()[-1].replace(")", ""))
            )

    for i, root in enumerate(curr_energy_index):
        curr_energy.append(total_arrays[int(root) - 1][i])

    if output_image is True:
        output_image = filename

    for i in range(len(total_arrays)):
        plt.plot(x, cis_array + total_arrays[i], label="Root %i" % (i + 1))

    plt.scatter(
        x,
        cis_array + curr_energy,
        label="Active root",
        marker="x",
        c="rebeccapurple",
    )
    plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    plt.xlabel("Step")
    plt.ylabel("Energy / %s" % units)
    plt.title(output_image.replace(".in.out", ""))
    plt.savefig(
        output_image.replace(".in.out", ".jpg"), dpi=250, bbox_inches="tight"
    )
    # plt.show()

    to_return = [np.array(cis_array + array) for array in total_arrays]

    to_return.append(cis_array + curr_energy)

    return to_return


if __name__ == "__main__":
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    a = alberto(args.f, args.en, args.o, args.u)
