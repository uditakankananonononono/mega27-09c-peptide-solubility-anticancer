import numpy as np

from pepx.encoders import aac, dpc, descriptor_array, one_hot


def test_one_hot_shape_and_content():
    m = one_hot("ACD", 5)
    assert m.shape == (5, 21)
    assert m[0].argmax() == 0 and m[1].argmax() == 1 and m[2].argmax() == 2
    assert m[3].sum() == 0 and m[4].sum() == 0  # zero padding


def test_aac_sums_to_one():
    v = aac("AACCDE")
    assert v.sum() == np.float32(1.0)
    assert v[0] == np.float32(2 / 6)


def test_dpc_dim_and_norm():
    v = dpc("KLAKLAK")
    assert v.shape == (400,)
    assert abs(v.sum() - 1.0) < 1e-5


def test_descriptor_array_finite():
    v = descriptor_array("KLAKLAKKLAKLAK")
    assert v.dtype == np.float32 and np.all(np.isfinite(v))
