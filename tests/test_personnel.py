# -*- coding: utf-8 -*-
import pytest
import numpy as np
from barista.personnel.alberto import Alberto


def test_alberto_parsing():
    a = Alberto(
        "tests/alberto/alberto_raw.in.out"
    )  # Replace with your actual file

    expected_matrix = np.loadtxt("tests/alberto/alberto_absolute_test.dat")

    assert np.all(a.bulk_energy_array == expected_matrix)
