import numpy as np
from turbine_costsse.NEWturbine_costsse_2015 import *
import unittest
from openmdao.api import pyOptSparseDriver, ExecComp, IndepVarComp, Problem


class test_HubSystemCostAdder2015(unittest.TestCase):

    def setUp(self):

        self.rtol = 1E-5
        self.atol = 1E-5

        hub_cost = np.random.rand(1)*150000.+100000.
        pitch_system_cost = np.random.rand(1)*150000.+100000.
        spinner_cost = np.random.rand(1)*150000.+100000.

        prob = Problem()
        root = prob.root = Group()

        root.add('hub_cost', IndepVarComp('hub_cost', float(hub_cost)), promotes=['*'])
        root.add('pitch_system_cost', IndepVarComp('pitch_system_cost', float(pitch_system_cost)), promotes=['*'])
        root.add('spinner_cost', IndepVarComp('spinner_cost', float(spinner_cost)), promotes=['*'])
        root.add('HubSystemCostAdder2015', HubSystemCostAdder2015(), promotes=['*'])

        # set up optimizer
        prob.driver = pyOptSparseDriver()
        prob.driver.options['optimizer'] = 'SNOPT'
        prob.driver.add_objective('hub_system_cost')

        # select design variables
        prob.driver.add_desvar('hub_cost')
        prob.driver.add_desvar('pitch_system_cost')
        prob.driver.add_desvar('spinner_cost')

        # initialize problem
        prob.setup()

        # run problem
        prob.run_once()

        # pass results to self for use with unit test
        self.J = prob.check_total_derivatives(out_stream=None)

    def test(self):
        np.testing.assert_allclose(self.J[('hub_system_cost', 'hub_cost')]['J_fwd'], self.J[('hub_system_cost', 'hub_cost')]['J_fd'], self.rtol, self.atol)
        np.testing.assert_allclose(self.J[('hub_system_cost', 'pitch_system_cost')]['J_fwd'], self.J[('hub_system_cost', 'pitch_system_cost')]['J_fd'], self.rtol, self.atol)
        np.testing.assert_allclose(self.J[('hub_system_cost', 'spinner_cost')]['J_fwd'], self.J[('hub_system_cost', 'spinner_cost')]['J_fd'], self.rtol, self.atol)



if __name__ == "__main__":
    unittest.main()
