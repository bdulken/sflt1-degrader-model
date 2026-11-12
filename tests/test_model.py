import unittest
from dataclasses import replace
import numpy as np
from model import Parameters as DefaultParameters, setup, simulate

def Parameters(**kwargs):
    return DefaultParameters(**dict(dict(degrader_half_days=1e12,vegf_pg_ml=100,plgf_pg_ml=100,platelet_vegf_pg_ml=0,other_capacity_pm=100,platelet_kon_pm_day=0,plgf_other_capacity_pm=0,platelet_plgf_pg_ml=0,plgf_platelet_kon_pm_day=0),**kwargs))

class ModelTests(unittest.TestCase):
    def test_untreated_is_steady_with_independent_complex_clearance(self):
        p=Parameters(sv_half_hours=2,sp_half_hours=80)
        df,_=simulate(p,treated=False)
        for name in ['free_vegf_pg_mL','free_plgf_pg_mL','total_sflt1_ng_mL']:
            np.testing.assert_allclose(df[name],df[name].iloc[0],rtol=1e-6)

    def test_sflt1_total_matches_analytic_turnover(self):
        p=Parameters()
        df,_=simulate(p)
        ks=np.log(2)*24/p.sflt1_half_hours
        expected=p.sflt1_ng_ml*(ks/(ks+p.degrader_per_day)+p.degrader_per_day/(ks+p.degrader_per_day)*np.exp(-(ks+p.degrader_per_day)*df.day))
        np.testing.assert_allclose(df.total_sflt1_ng_mL,expected,rtol=2e-6)

    def test_release_conserves_ligand_when_natural_losses_vanish(self):
        p=Parameters(other_complex_half_hours=1e12,vegf_half_hours=1e12,plgf_half_hours=1e12,sv_half_hours=1e12,sp_half_hours=1e12,co_v=0,co_p=0)
        df,_=simulate(p,production=[0,0,0])
        np.testing.assert_allclose(df.total_vegf_pg_mL,100,rtol=1e-7)
        np.testing.assert_allclose(df.total_plgf_pg_mL,100,rtol=1e-7)

    def test_treatment_window_and_zero_input(self):
        p=Parameters(start_day=3,stop_day=5)
        df,_=simulate(p)
        np.testing.assert_allclose(df.loc[df.day<=3,'total_sflt1_ng_mL'],10,rtol=1e-6)
        self.assertGreater(df.total_sflt1_ng_mL.iloc[-1],df.loc[df.day==5,'total_sflt1_ng_mL'].iloc[0])
        empty,_=simulate(Parameters(sflt1_ng_ml=0,vegf_pg_ml=0,plgf_pg_ml=0))
        self.assertEqual(empty.drop(columns='day').to_numpy().max(),0)

    def test_full_codegradation_mass_balance(self):
        p=Parameters(other_complex_half_hours=1e12,vegf_half_hours=1e12,plgf_half_hours=1e12,sv_half_hours=1e12,sp_half_hours=1e12)
        df,_=simulate(p,production=[0,0,0])
        np.testing.assert_allclose(df.total_vegf_pg_mL+df.cumulative_vegf_codegraded_pg_mL,100,rtol=1e-7)
        np.testing.assert_allclose(df.total_plgf_pg_mL+df.cumulative_plgf_codegraded_pg_mL,100,rtol=1e-7)

class ExtraPoolTests(unittest.TestCase):
    def test_competing_pool_and_platelets_untreated_steady(self):
        p=Parameters(other_capacity_pm=100,platelet_vegf_pg_ml=200,platelet_release_per_day=0.01)
        df,_=simulate(p,treated=False)
        for name in ['free_vegf_pg_mL','other_bound_vegf_pg_mL','platelet_vegf_pg_mL','total_vegf_pg_mL','total_sflt1_ng_mL']:
            np.testing.assert_allclose(df[name],df[name].iloc[0],rtol=1e-6)
        self.assertAlmostEqual(df.total_vegf_pg_mL.iloc[0],100)
        self.assertAlmostEqual(df.all_pools_vegf_pg_mL.iloc[0],300)
        original,_=simulate(Parameters(other_capacity_pm=0),treated=False)
        self.assertLess(df.free_vegf_pg_mL.iloc[0],original.free_vegf_pg_mL.iloc[0])

    def test_platelet_storage_is_inaccessible(self):
        df,_=simulate(Parameters(platelet_vegf_pg_ml=500))
        np.testing.assert_allclose(df.platelet_vegf_pg_mL,500,rtol=1e-7)
        original,_=simulate(Parameters())
        np.testing.assert_allclose(df.free_vegf_pg_mL,original.free_vegf_pg_mL,rtol=1e-6)

    def test_all_pool_mass_balance_with_platelet_release(self):
        p=Parameters(other_capacity_pm=80,platelet_vegf_pg_ml=200,platelet_release_per_day=0.4,platelet_input_multiplier=0,platelet_half_hours=1e12,other_complex_half_hours=1e12,vegf_half_hours=1e12,sv_half_hours=1e12,sp_half_hours=1e12)
        df,_=simulate(p,production=[0,0,0])
        np.testing.assert_allclose(df.all_pools_vegf_pg_mL+df.cumulative_vegf_codegraded_pg_mL,300,rtol=1e-7)
        self.assertLess(df.platelet_vegf_pg_mL.iloc[-1],1)
        self.assertTrue((df.other_bound_vegf_pg_mL/p.vegf_mw_kda <= p.other_capacity_pm+1e-7).all())

    def test_invalid_calibration_and_manual_override(self):
        p=Parameters(platelet_vegf_pg_ml=10000,platelet_release_per_day=10)
        with self.assertRaisesRegex(ValueError,'negative VEGF production'):
            simulate(p)
        df,_=simulate(p,production=[0,0,0])
        self.assertTrue(np.isfinite(df.to_numpy()).all())

if __name__=='__main__': unittest.main()

class DefaultAllocationTests(unittest.TestCase):
    def test_equal_nonfree_pools_and_free_remainder(self):
        p=DefaultParameters()
        df,_=simulate(p,treated=False)
        first=df.iloc[0]
        np.testing.assert_allclose([first.bound_vegf_pg_mL,first.other_bound_vegf_pg_mL,first.platelet_vegf_pg_mL],31.17339717903927,rtol=1e-7)
        self.assertAlmostEqual(first.all_pools_vegf_pg_mL,100)
        self.assertAlmostEqual(first.free_vegf_pg_mL,6.479808462882173)
        for col in ['free_vegf_pg_mL','platelet_vegf_pg_mL','other_bound_vegf_pg_mL','free_plgf_pg_mL','platelet_plgf_pg_mL','other_bound_plgf_pg_mL']:
            np.testing.assert_allclose(df[col],first[col],rtol=1e-6)

    def test_exchange_conserves_mass(self):
        p=DefaultParameters(vegf_half_hours=1e12,sv_half_hours=1e12,sp_half_hours=1e12,other_complex_half_hours=1e12,platelet_half_hours=1e12,platelet_input_multiplier=0)
        df,_=simulate(p,production=[0,0,0])
        expected=df.all_pools_vegf_pg_mL.iloc[0]
        np.testing.assert_allclose(df.all_pools_vegf_pg_mL+df.cumulative_vegf_codegraded_pg_mL,expected,rtol=1e-7)

class PlgfPoolTests(unittest.TestCase):
    def test_default_allocation(self):
        df,_=simulate(DefaultParameters(),treated=False)
        first=df.iloc[0]
        np.testing.assert_allclose([first.bound_plgf_pg_mL,first.other_bound_plgf_pg_mL,first.platelet_plgf_pg_mL],22.091093341564427,rtol=1e-7)
        self.assertAlmostEqual(first.all_pools_plgf_pg_mL,100)
        self.assertAlmostEqual(first.free_plgf_pg_mL,33.726719975306715)

    def test_both_ligands_mass_balance_with_release_and_exchange(self):
        p=DefaultParameters(vegf_half_hours=1e12,plgf_half_hours=1e12,sv_half_hours=1e12,sp_half_hours=1e12,other_complex_half_hours=1e12,plgf_other_complex_half_hours=1e12,platelet_half_hours=1e12,plgf_platelet_half_hours=1e12,platelet_input_multiplier=0,plgf_platelet_input_multiplier=0,platelet_release_per_day=0.2,plgf_platelet_release_per_day=0.3)
        df,_=simulate(p,production=[0,0,0])
        for ligand in ['vegf','plgf']:
            total=df[f'all_pools_{ligand}_pg_mL']
            np.testing.assert_allclose(total+df[f'cumulative_{ligand}_codegraded_pg_mL'],total.iloc[0],rtol=1e-7)
            self.assertTrue((df[f'other_bound_{ligand}_pg_mL']>=0).all())

    def test_plgf_calibration_rejects_negative_input(self):
        with self.assertRaisesRegex(ValueError,'negative PlGF production'):
            simulate(DefaultParameters(platelet_plgf_pg_ml=10000,plgf_platelet_release_per_day=10))
