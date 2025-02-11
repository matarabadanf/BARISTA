# -*- coding: utf-8 -*-
import pytest
import numpy as np
from barista.personnel.alberto import Alberto

a = Alberto("tests/alberto/alberto_raw.in.out") 

def test_alberto_parsing():

    expected_matrix = np.loadtxt("tests/alberto/alberto_absolute_test.dat")

    assert np.all(a.bulk_energy_array == expected_matrix)

def test_alberto_units():
    a = Alberto("tests/alberto/alberto_raw.in.out", units='eV')

    expected_matrix = np.loadtxt('tests/alberto/absolute_relative_eV.dat')

    assert np.all(a.energy_array == expected_matrix)

def test_alberto_reference_energy():
    a = Alberto("tests/alberto/alberto_raw.in.out", reference_energy=-600.0, units='Hartree')

    expected_matrix = np.loadtxt("tests/alberto/reference_energy_hartree.dat")

    assert np.all(np.isclose(a.energy_array, expected_matrix))
    
