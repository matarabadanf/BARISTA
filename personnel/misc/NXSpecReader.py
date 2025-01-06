import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from typing import List, Dict, Optional
from numpy.typing import NDArray


class NXSpecReader:
    """
    Reads and processes spectral data files to generate convoluted spectra.

    Parameters
    ----------
    filename : str
        Path to the input spectral data file
    npoints : int, optional
        Number of points for spectrum generation (default: 1000)
    """

    def __init__(self, filename: str, npoints: int = 1000) -> None:
        self.filename: str = filename
        self.npoints: int = npoints
        self.energy_array: Optional[NDArray[np.float64]] = None
        self.energy_dictionary: Dict[str, List[float]] = {}
        self.intensity_dictionary: Dict[str, List[float]] = {}
        self.file_content: List[str] = []
        self._max_state: int = 0
        self.energy_intensities_df: pd.DataFrame = pd.DataFrame()
        self.x_grid: NDArray[np.float64] = np.array([])
        self.x_grid_nm: NDArray[np.float64] = np.array([])
        self.full_spectrum: NDArray[np.float64] = np.array([])
        self._initialize()

    def _initialize(self) -> None:
        """Initialize all necessary data structures and perform initial calculations."""
        self._read_file()
        self._get_max_state()
        self._generate_dictionaries()
        self._extract_states()
        self._transform_dictionaries_to_df()
        self._generate_x_grid()
        self._generate_xgrid_nm()
        print(self._max_state)

    def _read_file(self) -> None:
        """Read the input file and store its contents."""
        with open(self.filename) as f:
            self.file_content = f.readlines()

    def _get_max_state(self) -> None:
        """Determine the highest excited state number in the file."""
        max_state = 0
        for line in self.file_content:
            if "Excited State" in line:
                curr_state = int(line.strip().split()[2].replace(":", ""))
                if curr_state > max_state:
                    max_state = curr_state
                else:
                    self._max_state = max_state
                    break

    def _generate_dictionaries(self) -> None:
        """Initialize dictionaries to store energy and intensity data for each state."""
        keys = [str(i) + ":" for i in range(1, self._max_state + 1)]
        for key in keys:
            self.energy_dictionary.update({key: []})
            self.intensity_dictionary.update({key: []})

    def _extract_states(self) -> None:
        """Extract energy and oscillator strength values for each excited state."""
        for line in self.file_content:
            if "Excited State" in line:
                split_line = line.strip().split()
                key = split_line[2]
                energy = split_line[4]
                osc_strenght = split_line[8].replace("f=", "")
                self.energy_dictionary[key].append(float(energy))
                self.intensity_dictionary[key].append(float(osc_strenght))

    def _transform_dictionaries_to_df(self) -> pd.DataFrame:
        """
        Convert energy and intensity dictionaries to a pandas DataFrame.

        Returns
        -------
        pd.DataFrame
            DataFrame containing state, energy, and oscillator strength data
        """
        df = (
            pd.DataFrame(
                {
                    "State": self.energy_dictionary.keys(),
                    "Energy": self.energy_dictionary.values(),
                    "Oscillator Strength": self.intensity_dictionary.values(),
                }
            )
            .explode(["Energy", "Oscillator Strength"])
            .reset_index(drop=True)
        )
        self.energy_intensities_df = df
        return df

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

    def _generate_x_grid(self, npoints: Optional[int] = None) -> None:
        """Generate energy grid for spectrum calculation."""
        if npoints is None:
            npoints = self.npoints
        minimum = min(self.energy_intensities_df["Energy"]) - 2
        maximum = max(self.energy_intensities_df["Energy"]) + 2
        self.x_grid = np.linspace(minimum, maximum, npoints)

    def _generate_xgrid_nm(self) -> None:
        """Generate wavelength grid in nanometers."""
        self.x_grid_nm = 1239.8 / self.x_grid

    def _convolute_spectrum(
        self, filtered_df: pd.DataFrame
    ) -> NDArray[np.float64]:
        """
        Generate convolved spectrum for given states.

        Parameters
        ----------
        filtered_df : pd.DataFrame
            DataFrame containing energy and oscillator strength values

        Returns
        -------
        NDArray[np.float64]
            Convolved spectrum
        """
        spectrum = np.zeros(len(self.x_grid))
        for index, row in filtered_df.iterrows():
            energy = row["Energy"]
            oscillator_strength = row["Oscillator Strength"]
            spectrum += self._gaussian(self.x_grid, oscillator_strength, energy)
        return spectrum

    def generate_spectrum(self) -> NDArray[np.float64]:
        """
        Generate full spectrum across all states.

        Returns
        -------
        NDArray[np.float64]
            Complete convolved spectrum
        """
        spectrum = np.zeros(len(self.x_grid))
        for state in pd.unique(self.energy_intensities_df["State"]):
            filtered_df = self.energy_intensities_df[
                self.energy_intensities_df["State"] == state
            ]
            spectrum += self._convolute_spectrum(filtered_df)
        self.full_spectrum = spectrum
        return spectrum

    def return_semiclassical_spectrum(self):
        spec = np.copy(self.generate_spectrum())
        spec /= max(spec)
        return [np.copy(self.x_grid), spec]

    def save_to_file(self, filename: str) -> None:
        """
        Save generated spectrum to file.

        Parameters
        ----------
        filename : str
            Output file path
        """
        with open(filename, "w") as f:
            for i, x in enumerate(self.x_grid):
                y = self.full_spectrum[i]
                f.write(f"{x} {y}\n")
