import unittest
from dataclasses import replace
import numpy as np
from scipy.integrate import quad
from model import Parameters, simulate

class ExposureTests(unittest.TestCase):
    def test_decaying_clearance_against_total_sflt1_solution(self):
        p=Parameters(degrader_destruction_probability=0)
        ks=np.log(2)*24/p.sflt1_half_hours
        for half in [1,3,7]:
            df,_=simulate(p,degrader_half_days=half)
            lam=np.log(2)/half
            def integral_rate(t):
                return ks*t+(p.degrader_initial_pm*p.degrader_encounter_pm_day)/lam*(-np.expm1(-lam*t))
            for i in [50,350,700]:
                t=df.day.iloc[i]
                expected=np.exp(-integral_rate(t))*(10+ks*10*quad(lambda u:np.exp(integral_rate(u)),0,t)[0])
                self.assertAlmostEqual(df.total_sflt1_ng_mL.iloc[i],expected,places=5)

    def test_very_long_half_life_matches_constant_exposure(self):
        a,_=simulate(Parameters(degrader_destruction_probability=0),degrader_half_days=None)
        b,_=simulate(Parameters(degrader_destruction_probability=0),degrader_half_days=1e10)
        for col in ['free_vegf_pg_mL','free_plgf_pg_mL','total_sflt1_ng_mL']:
            np.testing.assert_allclose(a[col],b[col],rtol=1e-6)

    def test_zero_clearance_and_validation(self):
        p=Parameters(degrader_encounter_pm_day=0)
        a,_=simulate(p,degrader_half_days=1)
        np.testing.assert_allclose(a.free_vegf_pg_mL,a.free_vegf_pg_mL.iloc[0],rtol=1e-6)
        for half in [0,-1,np.nan,np.inf]:
            with self.assertRaises(ValueError):simulate(degrader_half_days=half)

    def test_default_is_single_dose_with_one_day_half_life(self):
        p=Parameters(degrader_destruction_probability=0)
        self.assertEqual(p.degrader_half_days,1)
        default,meta=simulate(p)
        explicit,_=simulate(p,degrader_half_days=1)
        np.testing.assert_allclose(default,explicit)
        self.assertEqual(meta['degrader_half_days'],1)
        constant,_=simulate(p,degrader_half_days=None)
        self.assertGreater(default.total_sflt1_ng_mL.iloc[-1],constant.total_sflt1_ng_mL.iloc[-1])
