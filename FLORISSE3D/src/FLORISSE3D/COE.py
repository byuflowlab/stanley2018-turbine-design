import numpy as np
from math import pi, log
import time
from openmdao.api import Group, Component, Problem, ScipyGMRES, IndepVarComp
#, pyOptSparseDriver
from FLORISSE3D.floris import AEPGroup
from FLORISSE3D.BOS import BOSgroup
from turbine_costsse.NEWnrel_csm_tcc_2015 import nrel_csm_tcc_2015
from turbine_costsse.NEWnrel_csm_tcc_2015 import cost_group_variableRotorStudy
from FLORISSE3D.simpleTower import calcMass, TowerDiscretization

class rotorCostComponent(Component):
    """
    Component to calculate the cost of the wind farm
    """

    def __init__(self, nTurbines, nGroups):

        super(rotorCostComponent, self).__init__()

        self.nTurbines = nTurbines
        self.nGroups = nGroups
        self.add_param('hGroup', np.zeros(nTurbines), desc='An array indicating which turbines are of each height')
        for i in range(nGroups):
            self.add_param('rotorCost%s'%i, 0.0, desc='costs of the rotors')

        self.add_output('rotorCost', np.zeros(nTurbines), desc='Array cost of rotors')

    def solve_nonlinear(self, params, unknowns, resids):
        nGroups = self.nGroups
        nTurbs = self.nTurbines
        hGroup = params['hGroup']

        rotorCost = np.zeros(self.nTurbines)
        for j in range(nGroups):
            for k in range(nTurbs):
                if j == hGroup[k]:
                    rotorCost[k] = params['rotorCost%s'%j]

        unknowns['rotorCost'] = rotorCost


    def linearize(self, params, unknowns, resids):

        J = {}

        for j in range(self.nGroups):
            J['rotorCost','rotorCost%s'%j] = np.zeros(self.nTurbines)
            for k in range(self.nTurbines):
                if j == params['hGroup'][k]:
                    J['rotorCost','rotorCost%s'%j][k] = 1.
        return J


class nacelleCostComponent(Component):
    """
    Component to calculate the cost of the wind farm
    """

    def __init__(self, nTurbines, nGroups):

        super(nacelleCostComponent, self).__init__()

        self.nTurbines = nTurbines
        self.nGroups = nGroups
        self.add_param('hGroup', np.zeros(nTurbines), desc='An array indicating which turbines are of each height')
        for i in range(nTurbines):
            self.add_param('nacelleCost%s'%i, 0.0, desc='costs of the nacelles')

        self.add_output('nacelleCost', np.zeros(nTurbines), desc='Array cost of nacelles')

    def solve_nonlinear(self, params, unknowns, resids):
        nGroups = self.nGroups
        nTurbs = self.nTurbines
        hGroup = params['hGroup']

        nacelleCost = np.zeros(self.nTurbines)
        for j in range(nGroups):
            for k in range(nTurbs):
                if j == hGroup[k]:
                    nacelleCost[k] = params['nacelleCost%s'%j]

        unknowns['nacelleCost'] = nacelleCost


    def linearize(self, params, unknowns, resids):

        J = {}

        for j in range(self.nGroups):
            J['nacelleCost','nacelleCost%s'%j] = np.zeros(self.nTurbines)
            for k in range(self.nTurbines):
                if j == params['hGroup'][k]:
                    J['nacelleCost','nacelleCost%s'%j][k] = 1.
        return J


class farmCost(Component):
    """
    Component to calculate the cost of the wind farm
    #TODO check gradients
    """

    def __init__(self, nTurbines, nGroups):

        super(farmCost, self).__init__()

        self.nTurbines = nTurbines
        self.nGroups = nGroups
        self.add_param('hGroup', np.zeros(nTurbines), desc='An array indicating which turbines are of each height')
        for i in range(nGroups):
            self.add_param('mass%s'%i, 0.0, units='kg',
                        desc='mass of each tower')
        self.add_param('rotorCost', np.zeros(nTurbines), desc='costs of the rotors')
        self.add_param('nacelleCost', np.zeros(nTurbines), desc='costs of the nacelles')

        self.add_output('cost', 0.0, desc='Cost of the wind farm')

    def solve_nonlinear(self, params, unknowns, resids):

        nTurbines = self.nTurbines
        nGroups = self.nGroups
        hGroup = params['hGroup']

        """these are for the NREL 5 MW reference turbine"""
        # nacelle_cost = 2446465.19*3/4.#turbine_costsse_2015.py run
        # rotor_cost = 1658752.71*3/4. #turbine_costsse_2015.py run
        # nacelle_cost = 1715919.90 #nrel_csm_tcc_2015.py run
        # rotor_cost = 1206984.20 #nrel_csm_tcc_2015.py run
        rotor_cost = 0.0
        nacelle_cost = 0.0
        for i in range(nTurbines):
            rotor_cost += params['rotorCost'][i]
            nacelle_cost += params['nacelleCost'][i]

        tower_mass_cost_coefficient = 3.08 #$/kg
        self.tower_mass_cost_coefficient = tower_mass_cost_coefficient

        tower_mass = np.zeros(nTurbines)
        for j in range(nGroups):
            for k in range(nTurbines):
                if j == hGroup[k]:
                    tower_mass[k] = params['mass%s'%j]

        tower_cost = np.zeros(nTurbines)
        for i in range(nTurbines):
            tower_cost[i] = tower_mass_cost_coefficient*tower_mass[i] #new cost from Katherine


        # print 'ROTOR COST: ', rotor_cost
        # print 'NACELLE COST: ', nacelle_cost
        # print 'TOWER COST: ', tower_cost

        parts_cost_farm = nacelle_cost + rotor_cost + np.sum(tower_cost) #parts cost for the entire wind farm
        turbine_multiplier = 4./3.
        turbine_cost = turbine_multiplier * parts_cost_farm
        unknowns['cost'] = turbine_cost

    def linearize(self, params, unknowns, resids):
        nTurbs = self.nTurbines
        tower_mass_cost_coefficient = self.tower_mass_cost_coefficient

        J = {}
        group_count = np.zeros(self.nGroups)
        for i in range(self.nGroups):
            for j in range(nTurbs):
                if params['hGroup'][j]==i:
                    group_count[i] += 1

        for i in range(self.nGroups):
            J['cost','mass%s'%i] = group_count[i]*tower_mass_cost_coefficient*4./3.

        J['cost','rotorCost'] = np.zeros((1,self.nTurbines))
        J['cost','nacelleCost'] = np.zeros((1,self.nTurbines))
        for i in range(nTurbs):
            J['cost','rotorCost'][0][i] = 4./3.
            J['cost','nacelleCost'][0][i] = 4./3.

        return J


class COEComponent(Component):
    """
    # Componenet to calculate the cost of energy (COE)
    """

    def __init__(self, nTurbines):

        super(COEComponent, self).__init__()

        # self.deriv_options['form'] = 'forward'
        # self.deriv_options['step_size'] = 1.E-6
        # self.deriv_options['step_calc'] = 'relative'

        self.nTurbines = nTurbines
        self.add_param('cost', 0.0, desc='Cost of the wind farm')
        self.add_param('AEP', 0.0, desc='AEP of the wind farm')
        self.add_param('BOS', 0.0, desc='BOS cost')

        self.add_output('COE', 0.0, desc='Cost of Energy for the wind farm')


    def solve_nonlinear(self, params, unknowns, resids):

        cost = params['cost']
        AEP = params['AEP']
        BOS = params['BOS']

        fixed_charge_rate = 0.102 #http://www.nrel.gov/docs/fy15osti/63267.pdf pg 58: 2013 COE
        tax_rate = 0.389 #http://www.nrel.gov/docs/fy15osti/63267.pdf pg 54 tax credit calculation
        O_M_coeff = 0.01 #operating and maintainence cost per kWh

        # unknowns['COE'] = 1000.*(fixed_charge_rate*(cost+bos)+ 0.0122*19566000*(1-tax_rate))/AEP # $/MWh
        unknowns['COE'] = 1000.*(fixed_charge_rate*(cost+BOS)+ O_M_coeff*AEP*(1.-tax_rate))/AEP # $/MWh

    def linearize(self, params, unknowns, resids):

        cost = params['cost']
        AEP = params['AEP']
        BOS = params['BOS']

        fixed_charge_rate = 0.102
        tax_rate = 0.389
        O_M_coeff = 0.01

        J = {}
        J['COE', 'cost'] = 1000.*fixed_charge_rate/AEP
        J['COE', 'AEP'] = -1000.*((BOS+cost)*fixed_charge_rate+AEP*O_M_coeff*(1.-tax_rate))/(AEP**2) + \
                            1000.*O_M_coeff*(1.-tax_rate)/AEP
        J['COE', 'BOS'] = 1000.*fixed_charge_rate/AEP


        return J


class COEGroup(Group):
    """
    Group containing components ot calculate COEGroup
    """
    def __init__(self, nTurbines, nGroups, nDirections, nPoints, nFull):

        super(COEGroup, self).__init__()

        self.add('farmCost', farmCost(nTurbines, nGroups), promotes=['*'])
        self.add('COEComponent', COEComponent(nTurbines), promotes=['*'])
        for i in range(nGroups):
            self.add('rotor_nacelle_costs%s'%i, cost_group_variableRotorStudy(), promotes=['bearing_number'])#, promotes=['turbine_class','blade_has_carbon','bearing_number'])
            self.add('calcMass%s'%i, calcMass(nFull))
            self.add('TowerDiscretization%s'%i, TowerDiscretization(nPoints, nFull))
        self.add('rotorCostComponent', rotorCostComponent(nTurbines, nGroups), promotes=['*'])
        self.add('nacelleCostComponent', nacelleCostComponent(nTurbines, nGroups), promotes=['*'])
        self.add('BOSgroup', BOSgroup(nTurbines), promotes=['*'])
        self.add('AEPGroup', AEPGroup(nTurbines, nDirections=nDirections,
                    use_rotor_components=False, differentiable=True,
                    optimizingLayout=False, nSamples=0), promotes=['*'])

        for i in range(nGroups):
            self.connect('rotor_nacelle_costs%s.rotor_cost'%i,'rotorCost%s'%i)
            self.connect('rotor_nacelle_costs%s.nacelle_cost'%i,'nacelleCost%s'%i)
            self.connect('TowerDiscretization%s.z_full'%i,'calcMass%s.z_full'%i)
            self.connect('TowerDiscretization%s.d_full'%i,'calcMass%s.d_full'%i)
            self.connect('TowerDiscretization%s.t_full'%i,'calcMass%s.t_full'%i)
            self.connect('calcMass%s.mass'%i,'mass%s'%i)


if __name__=="__main__":
    """
    This is just to test during development
    """

    # H1 = 150.
    # H2 = 75.
    #
    # H1_H2 = np.array([0,0,1,0,1,0,1,1,1])
    #
    # nTurbines = len(H1_H2)
    #
    # prob = Problem()
    # root = prob.root = Group()
    #
    # root.add('H1', IndepVarComp('turbineH1', H1), promotes=['*'])
    # root.add('H2', IndepVarComp('turbineH2', H2), promotes=['*'])
    # root.add('getTurbineZ', getTurbineZ(nTurbines), promotes=['*'])
    # root.add('farmCost', farmCost(nTurbines), promotes=['*'])
    #
    # prob.driver = pyOptSparseDriver()
    # prob.driver.options['optimizer'] = 'SNOPT'
    # prob.driver.opt_settings['Major iterations limit'] = 1000
    #
    # prob.driver.add_objective('cost')
    #
    # # --- Design Variables ---
    # prob.driver.add_desvar('turbineH1', lower=60., upper=None)
    # prob.driver.add_desvar('turbineH2', lower=60., upper=None)
    #
    # prob.setup()
    #
    # #prob['turbineH1'] = H1
    # #prob['turbineH2'] = H2
    # prob['H1_H2'] = H1_H2
    #
    # prob.run()
    #
    # print "Cost: ", prob['cost']
    # print 'H1: ', prob['turbineH1']
    # print 'H2: ', prob['turbineH2']


    prob = Problem()
    root = prob.root = Group()

    mass = np.array([100000.])
    nTurbines = len(mass)
    root.add('COEGroup', farmCost(nTurbines), promotes=['*'])

    prob.setup()

    #prob['turbineH1'] = H1
    #prob['turbineH2'] = H2
    for i in range(nTurbines):
        prob['mass%s'%i] = mass[i]

    prob.run()

    print "Cost: ", prob['cost']
    for i in range(nTurbines):
        print 'Mass ', i, ': ', prob['mass%s'%i]
