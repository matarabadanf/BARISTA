import pytest
import numpy as np
from barista.personnel.javi import Javi

def test_javi_ethylene_a():
    # Initialize Javi with the data files
    j = Javi(
        'ethylene_a/engrad0.dat',
        'ethylene_a/engrad1.dat',
        'ethylene_a/nacme2.dat'
    )
    
    # Load data directly from files
    ga = np.loadtxt('ethylene_a/engrad0.dat').flatten()
    gb = np.loadtxt('ethylene_a/engrad1.dat').flatten()
    hab = np.loadtxt('ethylene_a/nacme2.dat').flatten()
    
    # Test gradients match
    assert np.all(np.isclose(j._ga, ga)), "g_a and j.ga should be equal"
    assert np.all(np.isclose(j._gb, gb)), "g_b and j.gb should be equal"
    assert np.all(np.isclose(j.h_ab, hab)), "nac and j.nac should be equal"
    
    # Test g_ab calculation
    g_ab = (gb - ga) / 2
    assert np.all(np.isclose(j.g_ab, g_ab)), "g_ab and j.g_ab should be equal"
    
    # Test beta calculation
    tan_2_beta = 2 * np.dot(g_ab, hab) / (np.dot(g_ab, g_ab) - np.dot(hab, hab))
    beta_2 = np.atan(tan_2_beta)
    beta = beta_2 / 2
    assert np.all(np.isclose(j.beta, beta)), "beta and j.beta should be equal"
    
    # Reference values
    ref_pitch = 0.0949
    ref_asymm = 0.5320
    ref_sigma = 0.9550
    ref_theta = 0
    
    # Test against reference values with some tolerance
    tolerance = 1e-4
    assert abs(j.pitch - ref_pitch) < tolerance, f"Pitch difference exceeds tolerance: {abs(j.pitch - ref_pitch)}"
    assert abs(j.asymmetry - ref_asymm) < tolerance, f"Asymmetry difference exceeds tolerance: {abs(j.asymmetry - ref_asymm)}"
    assert abs(j.sigma - ref_sigma) < tolerance, f"Sigma difference exceeds tolerance: {abs(j.sigma - ref_sigma)}"
    # assert abs(j.theta_s % np.pi - ref_theta % np.pi) < tolerance, f"Theta_s difference exceeds tolerance: {abs(j.theta_s % np.pi - ref_theta % np.pi)}"
    
    # Test conical characterization values
    assert abs(j.p[0] - 0.60) < tolerance, f"p[0] should be approximately 0.60, got {j.p[0]}"
    assert abs(j.b[0] - 1.07) < tolerance, f"b[0] should be approximately 1.07, got {j.b[0]}"
