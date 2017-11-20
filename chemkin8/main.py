import numpy as np
import chemkin8
from chemkin8 import chemkin
import os

# Input variables
x1 = [2., 1., .5, 1., 1., 1.]
x2 = [2., 1., .5, 1., 1., 1., .5, 1.]
T = 1500

test_data_dir = os.path.join(os.path.dirname(chemkin8.__file__), 'tests/')
fname1 = os.path.join(test_data_dir, 'rxns.xml')
fname2 = os.path.join(test_data_dir, 'rxnset_long.xml')
fname3 = os.path.join(test_data_dir, 'rxns_reversible.xml')

c = chemkin.chemkin(fname1)
print(c.reaction_rates(x1, T))

c = chemkin.chemkin(fname2)
print(c.reaction_rates(x2, T))

c = chemkin.chemkin(fname3)
print(c.reaction_rates(x2, T))

print("NASA Coeffs:", c.NASAcoeffs())

