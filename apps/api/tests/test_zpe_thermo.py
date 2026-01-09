from __future__ import annotations

from services.zpe.thermo import calc_zpe_and_s_vib


def test_calc_zpe_and_s_vib_units():
    zpe_ev, s_vib_jmol_k = calc_zpe_and_s_vib(
        [1000.0],
        low_cut_cm=50.0,
        temperature=298.15,
    )
    assert abs(zpe_ev - 0.0619920992) < 1.0e-6
    assert s_vib_jmol_k > 0.1
