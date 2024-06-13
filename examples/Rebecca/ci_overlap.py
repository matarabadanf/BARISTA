'''
Overlap of two CISD wave functions (they can be obtained from different
geometries).
'''

from functools import reduce
import numpy
from pyscf import gto, scf, ci

#
# RCISD wavefunction overlap
#
myhf1 = gto.M(atom='H 0 0 0; F 0 0 1.1', basis='sto-3g', verbose=0).apply(scf.RHF).run()
ci1 = ci.CISD(myhf1).run()
print('CISD energy of mol1', ci1.e_tot)

myhf2 = gto.M(atom='H 0 0 0; F 0 0 1.2', basis='sto-3g', verbose=0).apply(scf.RHF).run()
ci2 = ci.CISD(myhf2).run()
print('CISD energy of mol2', ci2.e_tot)

s12 = gto.intor_cross('cint1e_ovlp_sph', myhf1.mol, myhf2.mol)
print(s12)
s12 = reduce(numpy.dot, (myhf1.mo_coeff.T, s12, myhf2.mo_coeff))
print(s12)
nmo = myhf2.mo_energy.size
nocc = myhf2.mol.nelectron // 2
print('<CISD-mol1|CISD-mol2> = ', ci.cisd.overlap(ci1.ci, ci2.ci, nmo, nocc, s12))
