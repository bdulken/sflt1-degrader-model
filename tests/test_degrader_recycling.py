import unittest
from dataclasses import replace
import numpy as np
from model import Parameters, simulate

class RecyclingTests(unittest.TestCase):
    def test_degrader_mass_balance_and_event_fates(self):
        for probability in [0,.05,.25,.5,1]:
            p=Parameters(degrader_destruction_probability=probability)
            df,_=simulate(p)
            np.testing.assert_allclose(df.degrader_pM+df.cumulative_degrader_destroyed_pM+df.cumulative_degrader_background_loss_pM,p.degrader_initial_pm,rtol=1e-7)
            np.testing.assert_allclose(df.cumulative_degrader_destroyed_pM,probability*df.cumulative_target_degradation_pM,atol=1e-8)
            np.testing.assert_allclose(df.cumulative_degrader_recycling_events_pM,(1-probability)*df.cumulative_target_degradation_pM,atol=1e-8)

    def test_full_consumption_limits_target_removal_to_dose(self):
        p=Parameters(degrader_destruction_probability=1)
        df,_=simulate(p,degrader_half_days=None)
        self.assertLessEqual(df.cumulative_target_degradation_pM.iloc[-1],p.degrader_initial_pm+1e-7)
        np.testing.assert_allclose(df.degrader_pM+df.cumulative_target_degradation_pM,p.degrader_initial_pm,rtol=1e-7)

    def test_perfect_recycling_can_exceed_one_target_per_degrader(self):
        p=Parameters(degrader_destruction_probability=0)
        df,_=simulate(p,degrader_half_days=None)
        np.testing.assert_allclose(df.degrader_pM,p.degrader_initial_pm,rtol=1e-7)
        self.assertGreater(df.cumulative_target_degradation_pM.iloc[-1],p.degrader_initial_pm)

    def test_zero_dose_is_untreated(self):
        df,_=simulate(Parameters(degrader_initial_pm=0))
        for col in ['total_sflt1_ng_mL','free_vegf_pg_mL','free_plgf_pg_mL']:
            np.testing.assert_allclose(df[col],df[col].iloc[0],rtol=1e-6)

    def test_delayed_dose_and_stopped_activity(self):
        p=Parameters(start_day=3,stop_day=5)
        df,_=simulate(p)
        self.assertTrue((df.loc[df.day<3,'degrader_pM']==0).all())
        self.assertEqual(df.loc[df.day==3,'degrader_pM'].iloc[0],p.degrader_initial_pm)
        np.testing.assert_allclose(df.loc[df.day>=3,'degrader_pM']+df.loc[df.day>=3,'cumulative_degrader_destroyed_pM']+df.loc[df.day>=3,'cumulative_degrader_background_loss_pM'],100,rtol=1e-7)
        np.testing.assert_allclose(df.loc[df.day>=5,'cumulative_target_degradation_pM'],df.loc[df.day==5,'cumulative_target_degradation_pM'].iloc[0],rtol=1e-7)
        self.assertTrue((df.loc[df.day>=5,'degrader_clearance_per_day']==0).all())

    def test_probability_validation(self):
        for value in [-.1,1.1,np.nan]:
            with self.assertRaises(ValueError):simulate(Parameters(degrader_destruction_probability=value))
