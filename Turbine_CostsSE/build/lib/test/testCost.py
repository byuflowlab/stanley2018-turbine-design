#!/usr/bin/env python
# encoding: utf-8
"""
test_environment_gradients.py
Created by Andrew Ning on 2013-12-18.
Copyright (c) NREL. All rights reserved.
"""


import unittest
import numpy as np
from openmdao.api import pyOptSparseDriver, Problem, Group
from openmdao.api import ScipyOptimizer, IndepVarComp
import time
from turbine_costsse.NEWnrel_csm_tcc_2015 import nrel_csm_tcc_2015



class TestRotorStuff(unittest.TestCase):

    def setUp(self):

        prob = Problem()
        root = prob.root = Group()

        root.add('nrel_csm_tcc_2015', nrel_csm_tcc_2015(), promotes=['*'])
        root.add('machine_rating', IndepVarComp('machine_rating', 5000.), promotes=['*'])

        prob.driver = pyOptSparseDriver()
        prob.driver.options['optimizer'] = 'SNOPT'

        obj = 'turbine_cost'

        #ROTOR MASS and ROTOR COST

        prob.driver.add_objective('%s'%obj)


        prob.driver.add_desvar('machine_rating')

        prob.setup()

        prob['rotor_diameter'] = 100.0
        prob['turbine_class'] = 'II/III'
        prob['blade_has_carbon'] = True
        prob['blade_number'] = 3
        # prob['machine_rating'] = 5000.0
        prob['hub_height'] = 80.0
        prob['bearing_number'] = 2
        prob['crane'] = True
        # prob['offshore'] = False

        # Rotor force calculations for nacelle inputs
        maxTipSpd = 80.0
        maxEfficiency = 0.9

        ratedHubPower  = prob['machine_rating']*1000. / maxEfficiency
        rotorSpeed     = (maxTipSpd/(0.5*prob['rotor_diameter'])) * (60.0 / (2.*np.pi))
        prob['rotor_torque'] = ratedHubPower/(rotorSpeed*(np.pi/30))

        prob.run_once()

        self.J = prob.check_total_derivatives(out_stream=None)
        self.obj = obj

        print 'analytic: '
        print self.J[('%s'%obj, 'machine_rating')]['J_fwd']
        print 'fd: '
        print self.J[('%s'%obj, 'machine_rating')]['J_fd']

    def test_rotor(self):
        obj = self.obj
        np.testing.assert_allclose(self.J[('%s'%obj, 'machine_rating')]['J_fwd'], self.J[('%s'%obj, 'machine_rating')]['J_fd'], 1e-6, 1e-6)


if __name__ == '__main__':
    unittest.main()
