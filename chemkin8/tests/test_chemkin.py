import os
from chemkin8 import chemkin
import chemkin8
import numpy as np

test_data_dir = os.path.join(os.path.dirname(chemkin8.__file__), 'tests/')
fname = os.path.join(test_data_dir, 'rxns.xml')
fname_r = os.path.join(test_data_dir, 'rxns_reversible.xml')


def test_k_const():
    c = chemkin.chemkin(fname)
    try:
        c.k_constant(['a'])
    except TypeError as err:
        assert(type(err) == TypeError)

def test_k_arr_t():
    c = chemkin.chemkin(fname)
    try:
        c.k_arrhenius(1., 2., -1500)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_k_arr_a():
    c = chemkin.chemkin(fname)
    try:
        c.k_arrhenius(-1., 2., 1500)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_k_arr_float():
    c = chemkin.chemkin(fname)
    try:
        c.k_arrhenius(1., 'e', 1500)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_k_mod_t():
    c = chemkin.chemkin(fname)
    try:
        c.k_modified(1., 1., 2., -1500)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_k_mod_complex():
    c = chemkin.chemkin(fname)
    try:
        c.k_modified(1., 0+8j, 2., 1500)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_k_mod_a():
    c = chemkin.chemkin(fname)
    try:
        c.k_modified(-1., 0., 2., 1500)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_k_mod_float():
    c = chemkin.chemkin(fname)
    try:
        c.k_modified(1., 1., 'e', 1500)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_progress_reac_k0():
    c = chemkin.chemkin(fname)
    try:
        c.progress_reaction([1,2,3],0,[2,1,0])
    except ValueError as err:
        assert(type(err) == ValueError)

def test_progress_reac_result():
    c = chemkin.chemkin(fname)
    assert c.progress_reaction([1,2,3],10,[2,1,0]) == 20

def test_progress_reac_k_length():
    c = chemkin.chemkin(fname)
    try:
        c.progress_reaction([1,2,3],10,[2,1])
    except ValueError as err:
        assert(type(err) == ValueError)

def test_progress_reac_k_zero():
    c = chemkin.chemkin(fname)
    try:
        c.progress_reaction([1,2,3],0,[2,1])
    except ValueError as err:
        assert(type(err) == ValueError)

def test_progress_reac_reac_conc():
    c = chemkin.chemkin(fname)
    try:
        c.progress_reaction([0,0,0],10,[2,1,1])
    except ValueError as err:
        assert(type(err) == ValueError)

def test_progress_reac_reac_coef():
    c = chemkin.chemkin(fname)
    try:
        c.progress_reaction([2,1,0],10,[0,0,0])
    except ValueError as err:
        assert(type(err) == ValueError)

def test_progress_system_k_exist():
    c = chemkin.chemkin(fname)
    try:
        c.progress_system([2., 1., .5, 1., 1., 1.], 1500)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_progress_system_k_dim():
    c = chemkin.chemkin(fname)
    c.k_system(1500)
    try:
        c.kf = [1., 2.]
        c.progress_system([2., 1., .5, 1., 1., 1.], 1500)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_progress_system_v2():
    c = chemkin.chemkin(fname)
    c.k_system(1500)
    try:
        c.v1 = [[0., 1., 1., 0., 0., 0+8j],
                [1., 0., 1., 0., 0., 0.],
                [1., 0., 0., 0., 1., 0.]]
        c.progress_system([2., 1., .5, 1., 1., 1.])
    except ValueError as err:
        assert(type(err) == ValueError)

def test_reaction_result():
    c = chemkin.chemkin(fname)
    assert c.reaction_rates([2., 1., .5, 1., 1., 1.], 1500) == [-227364086.53073898, 227364586.53073898, 231985198.37073097, -2311055.9199959813, 500.0, -229675142.45073497]

def test_reaction_v1_v2_dim_1():
    c = chemkin.chemkin(fname)
    try:
        c.v1 = [[1,2]]
        c.v2 = [[1,2], [3,4]]
        c.reaction_rates([2., 1., .5, 1., 1., 1.], 1500)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_reaction_v1_v2_dim_2():
    c = chemkin.chemkin(fname)
    try:
        c.v1 = [[1,2,3]]
        c.v2 = [[1,2], [3,4]]
        c.reaction_rates([2., 1., .5, 1., 1., 1.], 1500)
    except ValueError as err:
        assert(type(err) == ValueError)


def test_reaction_prod_coef_complex():
    c = chemkin.chemkin(fname)
    try:
        c.v2 = [[0., 1., 1., 0., 0., 0+8j],
                [1., 0., 1., 0., 0., 0.],
                [1., 0., 0., 0., 1., 0.]]
        c.reaction_rates([2., 1., .5, 1., 1., 1.], 1500)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_reaction_concentration_dim():
    c = chemkin.chemkin(fname)
    try:
        c.reaction_rates([2., .5, 1., 1., 1.], 1500)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_reaction_reversible():
    c = chemkin.chemkin(fname_r)
    assert(c.reaction_rates([2., 1., .5, 1., 1., 1., .5, 1.], 1500) == [132279654839465.03, -312877268431120.12, -137621783394319.91, 60704918199955.906, 64529905746598.25, 336486898686318.19, -41877131955541.195, -101625193691356.08])

def test_reaction_low_temp():
    c = chemkin.chemkin(fname_r)
    assert(c.reaction_rates([2., 1., .5, 1., 1., 1., .5, 1.], 300) == [2.1355696871822284e+22, -2.124491832491173e+22, -2.1577254121208575e+22, 15944750005197.543, 1.1077867058624066e+20, 2.1355696997475826e+22, -95741906595346.828, -13966889975866.303])

def test_nasa_feed():
    c = chemkin.chemkin(fname_r)
    c.parseNASA(1500, feed=[1,2,3])
    assert(c.nasa == [1,2,3])

def test_nasa_missing():
    c = chemkin.chemkin(fname_r)
    c.nasa = None
    try:
        c.reaction_rates([2., 1., .5, 1., 1., 1., .5, 1.], 300)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_repr():
    c = chemkin.chemkin(fname)
    assert(repr(c) == 'chemkin()')

'''
def test_nuclear():
    1. Test against known reaction steps
    2. Check against known half-life values
'''
