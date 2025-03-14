import pytest
import numpy as np
import os
from barista.personnel.javi import Javi

@pytest.fixture
def setup_data():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    ethylene_path = os.path.join(script_dir, 'ethylene_a')
    grad0_path = os.path.join(ethylene_path, 'engrad0.dat')
    grad1_path = os.path.join(ethylene_path, 'engrad1.dat')
    nacme_path = os.path.join(ethylene_path, 'nacme2.dat')

    j = Javi(
        grad0_path,
        grad1_path,
        nacme_path
    )
    
    ga = np.loadtxt(grad0_path).flatten()
    gb = np.loadtxt(grad1_path).flatten()
    hab = np.loadtxt(nacme_path).flatten()

    g_ab = (gb - ga) / 2
    tan_2_beta = 2 * np.dot(g_ab, hab) / (np.dot(g_ab, g_ab) - np.dot(hab, hab))
    beta_2 = np.atan(tan_2_beta)
    beta = beta_2 / 2
    
    ref_pitch = 0.0949
    ref_asymm = 0.5320
    ref_sigma = 0.9550
    ref_theta = 0
    
    return j, ga, gb, hab, g_ab, beta, ref_pitch, ref_asymm, ref_sigma, ref_theta

def test_g_a_equal(setup_data):
    j, ga, _, _, _, _, _, _, _, _ = setup_data
    assert np.all(np.isclose(j._ga, ga, atol=0.001)), 'g_a and j.ga are not equal within threshold.'

def test_g_b_equal(setup_data):
    j, _, gb, _, _, _, _, _, _, _ = setup_data
    assert np.all(np.isclose(j._gb, gb, atol=0.001)), 'g_b and j.gb are not equal within threshold.'

def test_nac_equal(setup_data):
    j, _, _, hab, _, _, _, _, _, _ = setup_data
    assert np.all(np.isclose(j.h_ab, hab, atol=0.001)), 'nac and j.nac are not equal within threshold.'

def test_g_ab_equal(setup_data):
    j, _, gb, ga, g_ab, _, _, _, _, _ = setup_data
    assert np.all(np.isclose(j.g_ab, g_ab, atol=0.001)), 'g_ab and j.g_ab are not equal within threshold.'

# def test_beta_equal(setup_data):
#     j, _, _, _, _, beta, _, _, _, _ = setup_data
#     assert np.all(np.isclose(j.beta, beta, atol=0.001)), 'beta and j.beta are not equal within threshold.'

def test_pitch(setup_data):
    j, _, _, _, _, _, ref_pitch, _, _, _ = setup_data
    assert abs(j.pitch - ref_pitch) < 0.001, f'Reference Pitch and j.pitch differ in {abs(j.pitch - ref_pitch):9.5f}'

def test_asymmetry(setup_data):
    j, _, _, _, _, _, _, ref_asymm, _, _ = setup_data
    assert abs(j.asymmetry - ref_asymm) < 0.001, f'Reference asymmetry and j.asymmetry differ in {abs(j.asymmetry - ref_asymm):9.5f}'

def test_sigma(setup_data):
    j, _, _, _, _, _, _, _, ref_sigma, _ = setup_data
    assert abs(j.sigma - ref_sigma) < 0.001, f'Reference sigma and j.sigma differ in {abs(j.sigma - ref_sigma):9.5f}'

def test_conical_characterization(setup_data):
    j, _, _, _, _, _, _, _, _, _ = setup_data
    assert abs(j.p[0] - 0.60) < 0.01, f'Conical calculated characterization (a) p[0] differs: {j.p[0]}'
    assert abs(j.b[0] - 1.07) < 0.01, f'Conical calculated characterization (a) b[0] differs: {j.b[0]}'

