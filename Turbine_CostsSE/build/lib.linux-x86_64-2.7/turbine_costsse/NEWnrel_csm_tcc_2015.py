"""
tcc_csm_component.py

Created by NWTC Systems Engineering Sub-Task on 2012-08-01.
Copyright (c) NREL. All rights reserved.
"""

import numpy as np
from openmdao.api import Group, Component, Problem, ScipyGMRES, IndepVarComp
from NEWturbine_costsse_2015 import Turbine_CostsSE_2015

# --------------------------------------------------------------------
class BladeMass(Component):

    def __init__(self):

        super(BladeMass, self).__init__()


        # Variables
        self.add_param('rotor_diameter', 0.0, desc= 'rotor diameter of the machine')
        self.add_param('turbine_class', 'I', desc= 'turbine class')
        self.add_param('blade_has_carbon', False, desc= 'does the blade have carbon?')
        self.add_param('blade_mass_coefficient', 0.5, desc= 'A in the blade mass equation: A*(rotor_diameter/B)^exp')
        self.add_param('rotor_diameter_denominator', 2.0, desc= 'B in the blade mass equation: A*(rotor_diameter/B)^exp')
        self.add_param('blade_user_exponent', 2.5, desc= 'optional user-entered exponent for the blade mass equation')

        # Outputs
        self.add_output('blade_mass', 0.0, desc= 'component mass [kg]')

    def solve_nonlinear(self, params, unknowns, resids):

        # select the exponent for the blade mass equation
        exponent = 0.0
        if params['turbine_class'] == 'I':
            if params['blade_has_carbon']:
              exponent = 2.47
            else:
              exponent = 2.54
        elif params['turbine_class'] == 'II/III':
            if params['blade_has_carbon']:
              exponent = 2.44
            else:
              exponent = 2.50
        else:
            exponent = params['blade_user_exponent']

        # calculate the blade mass
        unknowns['blade_mass'] = params['blade_mass_coefficient']*(params['rotor_diameter']/params['rotor_diameter_denominator'])**exponent

    def linearize(self, params, unknowns, resids):

        exponent = 0.0
        if params['turbine_class'] == 'I':
            if params['blade_has_carbon']:
              exponent = 2.47
            else:
              exponent = 2.54
        elif params['turbine_class'] == 'II/III':
            if params['blade_has_carbon']:
              exponent = 2.44
            else:
              exponent = 2.50
        else:
            exponent = params['blade_user_exponent']

        J = {}
        J['blade_mass','rotor_diameter'] = params['blade_mass_coefficient']*exponent*params['rotor_diameter']**(exponent-1)/(params['rotor_diameter_denominator']**exponent)

        return J

  # --------------------------------------------------------------------
class HubMass(Component):

    def __init__(self):

        super(HubMass, self).__init__()

        # Variables
        self.add_param('blade_mass', 0.0, desc= 'component mass [kg]')
        self.add_param('hub_mass_coeff', 2.3, desc= 'A in the hub mass equation: A*blade_mass + B')
        self.add_param('hub_mass_intercept', 1320., desc= 'B in the hub mass equation: A*blade_mass + B')

        # Outputs
        self.add_output('hub_mass', 0.0, desc= 'component mass [kg]')

    def solve_nonlinear(self, params, unknowns, resids):

        # calculate the hub mass
        unknowns['hub_mass'] = params['hub_mass_coeff']*params['blade_mass']+params['hub_mass_intercept']

    def linearize(self, params, unknowns, resids):

        J = {}
        J['hub_mass','blade_mass'] = params['hub_mass_coeff']

        return J


# --------------------------------------------------------------------
class PitchSystemMass(Component):

    def __init__(self):

        super(PitchSystemMass, self).__init__()

        # Variables
        self.add_param('blade_mass', 0.0, desc= 'component mass [kg]')
        self.add_param('blade_number', 3, desc= 'number of rotor blades')
        self.add_param('pitch_bearing_mass_coeff', 0.1295, desc= 'A in the pitch bearing mass equation: A*blade_mass*blade_number + B')
        self.add_param('pitch_bearing_mass_intercept', 491.31, desc= 'B in the pitch bearing mass equation: A*blade_mass*blade_number + B')
        self.add_param('bearing_housing_percent', 0.3280, desc= 'bearing housing percentage (in decimal form: ex 10% is 0.10)')
        self.add_param('mass_sys_offset', 555.0, desc= 'mass system offset')

        # Outputs
        self.add_output('pitch_system_mass', 0.0, desc= 'component mass [kg]')

    def solve_nonlinear(self, params, unknowns, resids):

        # calculate the hub mass
        pitchBearingMass = params['pitch_bearing_mass_coeff']* params['blade_mass'] * params['blade_number'] + params['pitch_bearing_mass_intercept']
        unknowns['pitch_system_mass'] = pitchBearingMass * (1 + params['bearing_housing_percent']) + params['mass_sys_offset']

    def linearize(self, params, unknowns, resids):

        J = {}
        J['pitch_system_mass','blade_mass'] = params['pitch_bearing_mass_coeff']*params['blade_number']*(1 + params['bearing_housing_percent'])


        return J

# --------------------------------------------------------------------
class SpinnerMass(Component):

    def __init__(self):

        super(SpinnerMass, self).__init__()

        # Variables
        self.add_param('rotor_diameter', 0.0, desc= 'rotor diameter of the machine')
        self.add_param('spinner_mass_coeff', 15.5, desc= 'A in the spinner mass equation: A*rotor_diameter + B')
        self.add_param('spinner_mass_intercept', -980.0, desc= 'B in the spinner mass equation: A*rotor_diameter + B')

        # Outputs
        self.add_output('spinner_mass', 0.0, desc='component mass [kg]')

    def solve_nonlinear(self, params, unknowns, resids):

        # calculate the spinner mass
        unknowns['spinner_mass'] = params['spinner_mass_coeff'] * params['rotor_diameter'] + params['spinner_mass_intercept']

    def linearize(self, params, unknowns, resids):

        J = {}
        J['spinner_mass','rotor_diameter'] = params['spinner_mass_coeff']

        return J

# --------------------------------------------------------------------
class LowSpeedShaftMass(Component):

    def __init__(self):

        super(LowSpeedShaftMass, self).__init__()

        # self.deriv_options['step_size'] = 1.E-4
        # self.deriv_options['step_calc'] = 'relative'
        # self.deriv_options['form'] = 'central'

        # Variables
        self.add_param('blade_mass', 0.0, desc='mass for a single wind turbine blade')
        self.add_param('machine_rating', 0.0, desc='machine rating') #kW
        self.add_param('lss_mass_coeff', 13., desc='A in the lss mass equation: A*(blade_mass*rated_power)^exp + B')
        self.add_param('lss_mass_exp', 0.65, desc='exp in the lss mass equation: A*(blade_mass*rated_power)^exp + B')
        self.add_param('lss_mass_intercept', 775., desc='B in the lss mass equation: A*(blade_mass*rated_power)^exp + B')

        # Outputs
        self.add_output('low_speed_shaft_mass', 0.0, desc='component mass [kg]')

    def solve_nonlinear(self, params, unknowns, resids):

        # calculate the lss mass
        # print 'Turbine Costs: LowSpeedShaftMass: ', params['machine_rating']
        # print 'LSS: ', params['machine_rating']
        unknowns['low_speed_shaft_mass'] = params['lss_mass_coeff'] * (params['blade_mass'] * params['machine_rating']/1000.)**params['lss_mass_exp'] + params['lss_mass_intercept']

    def linearize(self, params, unknowns, resids):

        J = {}
        # print 'blade mass: ', params['blade_mass']
        # print 'wrt mass: ', params['lss_mass_coeff']*params['lss_mass_exp'] * (params['blade_mass'] * params['machine_rating']/1000.)**params['lss_mass_exp']/params['blade_mass']
        # print 'wrt rating: ', params['lss_mass_coeff']*params['lss_mass_exp'] * (params['blade_mass'] * params['machine_rating']/1000.)**params['lss_mass_exp']/params['machine_rating']

        J['low_speed_shaft_mass','blade_mass'] = params['lss_mass_coeff']*params['lss_mass_exp'] * (params['blade_mass'] * params['machine_rating']/1000.)**params['lss_mass_exp']/params['blade_mass']
        J['low_speed_shaft_mass','machine_rating'] = params['lss_mass_coeff']*params['lss_mass_exp'] * (params['blade_mass'] * params['machine_rating']/1000.)**params['lss_mass_exp']/params['machine_rating']

        return J
# --------------------------------------------------------------------
class BearingMass(Component):

    def __init__(self):

        super(BearingMass, self).__init__()

        # self.deriv_options['step_size'] = 1.E-4
        # self.deriv_options['step_calc'] = 'relative'
        # self.deriv_options['form'] = 'central'


        # Variables
        self.add_param('rotor_diameter', 0.0, desc= 'rotor diameter of the machine')
        self.add_param('bearing_mass_coeff', 0.0001, desc= 'A in the bearing mass equation: A*rotor_diameter^exp') #default from ppt
        self.add_param('bearing_mass_exp', 3.5, desc= 'exp in the bearing mass equation: A*rotor_diameter^exp') #default from ppt

        # Outputs
        self.add_output('main_bearing_mass', 0.0, desc='component mass [kg]')

    def solve_nonlinear(self, params, unknowns, resids):

        # calculates the mass of a SINGLE bearing
        unknowns['main_bearing_mass'] = params['bearing_mass_coeff'] * params['rotor_diameter']**params['bearing_mass_exp']

    def linearize(self, params, unknowns, resids):

        J = {}
        # print 'Comp stuff'
        # print params['bearing_mass_exp']*unknowns['main_bearing_mass']/params['rotor_diameter']
        # print 'rotor_diameter: ', params['rotor_diameter']
        J['main_bearing_mass','rotor_diameter'] = params['bearing_mass_exp']*unknowns['main_bearing_mass']/params['rotor_diameter']
        # print 'J: ', J

        return J
# --------------------------------------------------------------------
class GearboxMass(Component):

    def __init__(self):

        super(GearboxMass, self).__init__()

        # Variables
        self.add_param('rotor_torque', 0.0, desc = 'torque from rotor at rated power') #JMF do we want this default? N*m
        self.add_param('gearbox_mass_coeff', 113., desc= 'A in the gearbox mass equation: A*rotor_torque^exp')
        self.add_param('gearbox_mass_exp', 0.71, desc= 'exp in the gearbox mass equation: A*rotor_torque^exp')

        # Outputs
        self.add_output('gearbox_mass', 0.0, desc='component mass [kg]')

    def solve_nonlinear(self, params, unknowns, resids):

        # calculate the gearbox mass
        unknowns['gearbox_mass'] = params['gearbox_mass_coeff'] * (params['rotor_torque']/1000.0)**params['gearbox_mass_exp']

    def linearize(self, params, unknowns, resids):

        J = {}
        J['gearbox_mass','rotor_torque'] = params['gearbox_mass_coeff']/(1000.**params['gearbox_mass_exp']) * params['gearbox_mass_exp'] * params['rotor_torque']**(params['gearbox_mass_exp']-1.)

        return J
# --------------------------------------------------------------------
class HighSpeedSideMass(Component):

    def __init__(self):

        super(HighSpeedSideMass, self).__init__()

        # Variables
        self.add_param('machine_rating', 0.0, desc='machine rating') #kW
        self.add_param('hss_mass_coeff', 0.19894, desc= 'NREL CSM hss equation; removing intercept since it is negligible')

        # Outputs
        self.add_output('high_speed_side_mass', 0.0, desc='component mass [kg]')

    def solve_nonlinear(self, params, unknowns, resids):

        # TODO: this is in DriveSE; replace this with code in DriveSE and have DriveSE use this code??
        unknowns['high_speed_side_mass'] = params['hss_mass_coeff'] * params['machine_rating']
        # print 'Turbine Costs: HighSpeedSideMass: ', params['machine_rating']

    def linearize(self, params, unknowns, resids):

        J = {}
        J['high_speed_side_mass','machine_rating'] = params['hss_mass_coeff']

        return J
# --------------------------------------------------------------------
class GeneratorMass(Component):

    def __init__(self):

        super(GeneratorMass, self).__init__()

        # Variables
        self.add_param('machine_rating', 0.0, desc='machine rating')#kW
        self.add_param('generator_mass_coeff', 2300., desc= 'A in the generator mass equation: A*rated_power + B') #default from ppt
        self.add_param('generator_mass_intercept', 3400., desc= 'B in the generator mass equation: A*rated_power + B') #default from ppt

        # Outputs
        self.add_output('generator_mass', 0.0, desc='component mass [kg]')

    def solve_nonlinear(self, params, unknowns, resids):

        # calculate the generator mass
        unknowns['generator_mass'] = params['generator_mass_coeff'] * params['machine_rating']/1000. + params['generator_mass_intercept']
        # print 'Turbine Costs: GeneratorMass', params['machine_rating']

    def linearize(self, params, unknowns, resids):

        J = {}
        J['generator_mass','machine_rating'] = params['generator_mass_coeff']/1000.
        return J
# --------------------------------------------------------------------
class BedplateMass(Component):

    def __init__(self):

        super(BedplateMass, self).__init__()

        # Variables
        self.add_param('rotor_diameter', 0.0, desc= 'rotor diameter of the machine')
        self.add_param('bedplate_mass_exp', 2.2, desc= 'exp in the bedplate mass equation: rotor_diameter^exp')

        # Outputs
        self.add_output('bedplate_mass', 0.0, desc='component mass [kg]')

    def solve_nonlinear(self, params, unknowns, resids):

        # calculate the bedplate mass
        unknowns['bedplate_mass'] = params['rotor_diameter']**params['bedplate_mass_exp']

    def linearize(self, params, unknowns, resids):

        J = {}
        J['bedplate_mass','rotor_diameter'] = params['bedplate_mass_exp']*params['rotor_diameter']**(params['bedplate_mass_exp']-1.)
        return J
# --------------------------------------------------------------------
class YawSystemMass(Component):

    def __init__(self):

        super(YawSystemMass, self).__init__()

        # Variables
        self.add_param('rotor_diameter', 0.0, desc= 'rotor diameter of the machine')
        self.add_param('yaw_mass_coeff', 0.0009, desc= 'A in the yaw mass equation: A*rotor_diameter^exp') #NREL CSM
        self.add_param('yaw_mass_exp', 3.314, desc= 'exp in the yaw mass equation: A*rotor_diameter^exp') #NREL CSM

        # Outputs
        self.add_output('yaw_system_mass', 0.0, desc='component mass [kg]')

    def solve_nonlinear(self, params, unknowns, resids):

        # calculate yaw system mass #TODO - 50% adder for non-bearing mass
        unknowns['yaw_system_mass'] = 1.5 * (params['yaw_mass_coeff'] * params['rotor_diameter'] ** params['yaw_mass_exp']) #JMF do we really want to expose all these?

#TODO: no variable speed mass; ignore for now

    def linearize(self, params, unknowns, resids):

        J = {}
        J['yaw_system_mass','rotor_diameter'] = 1.5 * params['yaw_mass_coeff'] * params['yaw_mass_exp'] * params['rotor_diameter'] ** (params['yaw_mass_exp']-1.)
        return J
# --------------------------------------------------------------------
class HydraulicCoolingMass(Component):

    def __init__(self):

        super(HydraulicCoolingMass, self).__init__()

        # Variables
        self.add_param('machine_rating', 0.0, desc='machine rating')
        self.add_param('hvac_mass_coeff', 0.08, desc= 'hvac linear coefficient') #NREL CSM

        # Outputs
        self.add_output('hydraulic_cooling_mass', 0.0, desc='component mass [kg]')

    def solve_nonlinear(self, params, unknowns, resids):

        # calculate hvac system mass
        unknowns['hydraulic_cooling_mass'] = params['hvac_mass_coeff'] * params['machine_rating']
        # print 'Turbine Costs: HydraulicCoolingMass: ', params['machine_rating']

    def linearize(self, params, unknowns, resids):

        J = {}
        J['hydraulic_cooling_mass','machine_rating'] = params['hvac_mass_coeff']
        return J

# --------------------------------------------------------------------
class NacelleCoverMass(Component):

    def __init__(self):

        super(NacelleCoverMass, self).__init__()

        # Variables
        self.add_param('machine_rating', 0.0, desc='machine rating')
        self.add_param('cover_mass_coeff', 1.2817, desc= 'A in the spinner mass equation: A*rotor_diameter + B')
        self.add_param('cover_mass_intercept', 428.19, desc= 'B in the spinner mass equation: A*rotor_diameter + B')

        # Outputs
        self.add_output('nacelle_cover_mass', 0.0, desc='component mass [kg]')

    def solve_nonlinear(self, params, unknowns, resids):

        # calculate nacelle cover mass
        unknowns['nacelle_cover_mass'] = params['cover_mass_coeff'] * params['machine_rating'] + params['cover_mass_intercept']
        # print 'Turbine Costs: NacelleCoverMass: ', params['machine_rating']

    def linearize(self, params, unknowns, resids):

        J = {}
        J['nacelle_cover_mass','machine_rating'] = params['cover_mass_coeff']
        return J

# TODO: ignoring controls and electronics mass for now

# --------------------------------------------------------------------
class OtherMainframeMass(Component):
    # nacelle platforms, service crane, base hardware

    def __init__(self):

        super(OtherMainframeMass, self).__init__()

        # Variables
        self.add_param('bedplate_mass', 0.0, desc='component mass [kg]')
        self.add_param('nacelle_platforms_mass_coeff', 0.125, desc='nacelle platforms mass coefficient as a function of bedplate mass [kg/kg]') #default from old CSM
        self.add_param('crane', True, desc='flag for presence of onboard crane')
        self.add_param('crane_weight', 3000., desc='weight of onboard crane')
        #TODO: there is no base hardware mass model in the old model. Cost is not dependent on mass.

        # Outputs
        self.add_output('other_mass', 0.0, desc='component mass [kg]')

    def solve_nonlinear(self, params, unknowns, resids):

        # calculate nacelle cover mass
        nacelle_platforms_mass = params['nacelle_platforms_mass_coeff'] * params['bedplate_mass']

        # --- crane ---
        if (params['crane']):
            crane_mass =  params['crane_weight']
        else:
            crane_mass = 0.

        unknowns['other_mass'] = nacelle_platforms_mass + crane_mass

    def linearize(self, params, unknowns, resids):

        # --- crane ---
        if (params['crane']):
            crane_mass =  params['crane_weight'] #this is confusing: mass or weight?
        else:
            crane_mass = 0.

        J = {}
        J['other_mass','bedplate_mass'] = params['nacelle_platforms_mass_coeff']
        return J

# --------------------------------------------------------------------
class TransformerMass(Component):

    def __init__(self):

        super(TransformerMass, self).__init__()

        # Variables
        self.add_param('machine_rating', 0.0, desc='machine rating')
        self.add_param('transformer_mass_coeff', 1915., desc= 'A in the transformer mass equation: A*rated_power + B') #default from ppt
        self.add_param('transformer_mass_intercept', 1910., desc= 'B in the transformer mass equation: A*rated_power + B') #default from ppt

        # Outputs
        self.add_output('transformer_mass', 0.0, desc='component mass [kg]')

    def solve_nonlinear(self, params, unknowns, resids):

        # calculate the transformer mass
        unknowns['transformer_mass'] = params['transformer_mass_coeff'] * params['machine_rating']/1000. + params['transformer_mass_intercept']
        # print 'Turbine Costs: TransformerMass: ', params['machine_rating']

    def linearize(self, params, unknowns, resids):

        J = {}
        J['transformer_mass','machine_rating'] = params['transformer_mass_coeff']/1000.
        return J

# --------------------------------------------------------------------

#Won't use this
class TowerMass(Component):

    def __init__(self):

        super(TowerMass, self).__init__()

        # Variables
        self.add_param('hub_height', 0.0, desc= 'hub height of wind turbine above ground / sea level')
        self.add_param('tower_mass_coeff', 19.828, desc= 'A in the tower mass equation: A*hub_height + B') #default from ppt
        self.add_param('tower_mass_exp', 2.0282, desc= 'B in the tower mass equation: A*hub_height + B') #default from ppt

        # Outputs
        self.add_output('tower_mass', 0.0, desc='component mass [kg]')

    def solve_nonlinear(self, params, unknowns, resids):

        # calculate the tower mass
        unknowns['tower_mass'] = params['tower_mass_coeff'] * params['hub_height'] ** params['tower_mass_exp']

    def linearize(self, params, unknowns, resids):

        J = {}
        J['tower_mass','hub_height'] = params['tower_mass_coeff'] * params['tower_mass_exp'] * params['hub_height'] ** (params['tower_mass_exp']-1.)
        return J


# Turbine mass adder
class turbine_mass_adder(Component):

    def __init__(self):

        super(turbine_mass_adder, self).__init__()

        # Inputs
        # rotor
        self.add_param('blade_mass', 0.0, desc= 'component mass [kg]')
        self.add_param('hub_mass', 0.0, desc='component mass [kg]')
        self.add_param('pitch_system_mass', 0.0, desc='component mass [kg]')
        self.add_param('spinner_mass', 0.0, desc='component mass [kg]')
        # nacelle
        self.add_param('low_speed_shaft_mass', 0.0, desc='component mass [kg]')
        self.add_param('main_bearing_mass', 0.0, desc='component mass [kg]')
        self.add_param('gearbox_mass', 0.0, desc='component mass [kg]')
        self.add_param('high_speed_side_mass', 0.0, desc='component mass [kg]')
        self.add_param('generator_mass', 0.0, desc='component mass [kg]')
        self.add_param('bedplate_mass', 0.0, desc='component mass [kg]')
        self.add_param('yaw_system_mass', 0.0, desc='component mass [kg]')
        self.add_param('hydraulic_cooling_mass', 0.0, desc='component mass [kg]')
        self.add_param('nacelle_cover_mass', 0.0, desc='component mass [kg]')
        self.add_param('other_mass', 0.0, desc='component mass [kg]')
        self.add_param('transformer_mass', 0.0, desc='component mass [kg]')
        # tower
        self.add_param('tower_mass', 0.0, desc='component mass [kg]')

        # Parameters
        self.add_param('blade_number', 3, desc = 'number of rotor blades')
        self.add_param('bearing_number', 2, desc = 'number of main bearings')

        # Outputs
        self.add_output('hub_system_mass', 0.0, desc='hub system mass')
        self.add_output('rotor_mass', 0.0, desc='hub system mass')
        self.add_output('nacelle_mass', 0.0, desc='nacelle mass')
        self.add_output('turbine_mass', 0.0, desc='turbine mass')

    def solve_nonlinear(self, params, unknowns, resids):

        unknowns['hub_system_mass'] = params['hub_mass'] + params['pitch_system_mass'] + params['spinner_mass']
        unknowns['rotor_mass'] = params['blade_mass'] * params['blade_number'] + unknowns['hub_system_mass']
        unknowns['nacelle_mass'] = params['low_speed_shaft_mass'] + params['bearing_number'] * params['main_bearing_mass'] + \
                            params['gearbox_mass'] + params['high_speed_side_mass'] + params['generator_mass'] + \
                            params['bedplate_mass'] + params['yaw_system_mass'] + params['hydraulic_cooling_mass'] + \
                            params['nacelle_cover_mass'] + params['other_mass'] + params['transformer_mass']
        unknowns['turbine_mass'] = unknowns['rotor_mass'] + unknowns['nacelle_mass'] + params['tower_mass']

    def linearize(self, params, unknowns, resids):

        J = {}

        J['hub_system_mass','blade_mass'] = 0.
        J['hub_system_mass','hub_mass'] = 1.
        J['hub_system_mass','pitch_system_mass'] = 1.
        J['hub_system_mass','spinner_mass'] = 1.
        J['hub_system_mass','low_speed_shaft_mass'] = 0.
        J['hub_system_mass','main_bearing_mass'] = 0.
        J['hub_system_mass','gearbox_mass'] = 0.
        J['hub_system_mass','high_speed_side_mass'] = 0.
        J['hub_system_mass','generator_mass'] = 0.
        J['hub_system_mass','bedplate_mass'] = 0.
        J['hub_system_mass','yaw_system_mass'] = 0.
        J['hub_system_mass','hydraulic_cooling_mass'] = 0.
        J['hub_system_mass','nacelle_cover_mass'] = 0.
        J['hub_system_mass','other_mass'] = 0.
        J['hub_system_mass','transformer_mass'] = 0.
        J['hub_system_mass','tower_mass'] = 0.

        J['rotor_mass','blade_mass'] = params['blade_number']
        J['rotor_mass','hub_mass'] = 1.
        J['rotor_mass','pitch_system_mass'] = 1.
        J['rotor_mass','spinner_mass'] = 1.
        J['rotor_mass','low_speed_shaft_mass'] = 0.
        J['rotor_mass','main_bearing_mass'] = 0.
        J['rotor_mass','gearbox_mass'] = 0.
        J['rotor_mass','high_speed_side_mass'] = 0.
        J['rotor_mass','generator_mass'] = 0.
        J['rotor_mass','bedplate_mass'] = 0.
        J['rotor_mass','yaw_system_mass'] = 0.
        J['rotor_mass','hydraulic_cooling_mass'] = 0.
        J['rotor_mass','nacelle_cover_mass'] = 0.
        J['rotor_mass','other_mass'] = 0.
        J['rotor_mass','transformer_mass'] = 0.
        J['rotor_mass','tower_mass'] = 0.

        J['nacelle_mass','blade_mass'] = 0.
        J['nacelle_mass','hub_mass'] = 0.
        J['nacelle_mass','pitch_system_mass'] = 0.
        J['nacelle_mass','spinner_mass'] = 0.
        J['nacelle_mass','low_speed_shaft_mass'] = 1.
        J['nacelle_mass','main_bearing_mass'] = params['bearing_number']
        J['nacelle_mass','gearbox_mass'] = 1.
        J['nacelle_mass','high_speed_side_mass'] = 1.
        J['nacelle_mass','generator_mass'] = 1.
        J['nacelle_mass','bedplate_mass'] = 1.
        J['nacelle_mass','yaw_system_mass'] = 1.
        J['nacelle_mass','hydraulic_cooling_mass'] = 1.
        J['nacelle_mass','nacelle_cover_mass'] = 1.
        J['nacelle_mass','other_mass'] = 1.
        J['nacelle_mass','transformer_mass'] = 1.
        J['nacelle_mass','tower_mass'] = 0.

        J['turbine_mass','blade_mass'] = params['blade_number']
        J['turbine_mass','hub_mass'] = 1.
        J['turbine_mass','pitch_system_mass'] = 1.
        J['turbine_mass','spinner_mass'] = 1.
        J['turbine_mass','low_speed_shaft_mass'] = 1.
        J['turbine_mass','main_bearing_mass'] = params['bearing_number']
        J['turbine_mass','gearbox_mass'] = 1.
        J['turbine_mass','high_speed_side_mass'] = 1.
        J['turbine_mass','generator_mass'] = 1.
        J['turbine_mass','bedplate_mass'] = 1.
        J['turbine_mass','yaw_system_mass'] = 1.
        J['turbine_mass','hydraulic_cooling_mass'] = 1.
        J['turbine_mass','nacelle_cover_mass'] = 1.
        J['turbine_mass','other_mass'] = 1.
        J['turbine_mass','transformer_mass'] = 1.
        J['turbine_mass','tower_mass'] = 1.


        return J

# -------------------------------------------------------------------


class nrel_csm_mass_2015(Group):

    def __init__(self):

        super(nrel_csm_mass_2015, self).__init__()

        self.add('BladeMass', BladeMass(), promotes=['*'])
        self.add('HubMass', HubMass(), promotes=['*'])
        self.add('PitchSystemMass', PitchSystemMass(), promotes=['*'])
        self.add('SpinnerMass', SpinnerMass(), promotes=['*'])
        self.add('LowSpeedShaftMass', LowSpeedShaftMass(), promotes=['*'])
        self.add('BearingMass', BearingMass(), promotes=['*'])
        self.add('GearboxMass', GearboxMass(), promotes=['*'])
        self.add('HighSpeedSideMass', HighSpeedSideMass(), promotes=['*'])
        self.add('GeneratorMass', GeneratorMass(), promotes=['*'])
        self.add('BedplateMass', BedplateMass(), promotes=['*'])
        self.add('YawSystemMass', YawSystemMass(), promotes=['*'])
        self.add('HydraulicCoolingMass', HydraulicCoolingMass(), promotes=['*'])
        self.add('NacelleCoverMass', NacelleCoverMass(), promotes=['*'])
        self.add('OtherMainframeMass', OtherMainframeMass(), promotes=['*'])
        self.add('TransformerMass', TransformerMass(), promotes=['*'])
        self.add('TowerMass', TowerMass(), promotes=['*'])
        self.add('turbine_mass_adder', turbine_mass_adder(), promotes=['*'])




class nrel_csm_tcc_2015(Group):

    def __init__(self):

        super(nrel_csm_tcc_2015, self).__init__()

        self.add('trb_mass',nrel_csm_mass_2015(), promotes=['*'])
        self.add('tcc', Turbine_CostsSE_2015(), promotes=['*'])


class mass_group_variableRotorStudy(Group):

    def __init__(self):

        super(mass_group_variableRotorStudy, self).__init__()

        # self.add('BladeMass', BladeMass(), promotes=['*'])
        self.add('HubMass', HubMass(), promotes=['*'])
        self.add('PitchSystemMass', PitchSystemMass(), promotes=['*'])
        self.add('SpinnerMass', SpinnerMass(), promotes=['*'])
        self.add('LowSpeedShaftMass', LowSpeedShaftMass(), promotes=['*'])
        self.add('BearingMass', BearingMass(), promotes=['*'])
        self.add('GearboxMass', GearboxMass(), promotes=['*'])
        self.add('HighSpeedSideMass', HighSpeedSideMass(), promotes=['*'])
        self.add('GeneratorMass', GeneratorMass(), promotes=['*'])
        self.add('BedplateMass', BedplateMass(), promotes=['*'])
        self.add('YawSystemMass', YawSystemMass(), promotes=['*'])
        self.add('HydraulicCoolingMass', HydraulicCoolingMass(), promotes=['*'])
        self.add('NacelleCoverMass', NacelleCoverMass(), promotes=['*'])
        self.add('OtherMainframeMass', OtherMainframeMass(), promotes=['*'])
        self.add('TransformerMass', TransformerMass(), promotes=['*'])
        # self.add('TowerMass', TowerMass(), promotes=['*'])
        self.add('turbine_mass_adder', turbine_mass_adder(), promotes=['*'])




class cost_group_variableRotorStudy(Group):

    def __init__(self):

        super(cost_group_variableRotorStudy, self).__init__()

        self.add('mass', mass_group_variableRotorStudy(), promotes=['*'])
        self.add('tcc', Turbine_CostsSE_2015(), promotes=['*'])






# def mass_example():
#
#     # simple test of module
#     trb = nrel_csm_mass_2015()
#     trb.rotor_diameter = 126.0
#     trb.turbine_class = 'I'
#     trb.blade_has_carbon = False
#     trb.blade_number = 3
#     trb.machine_rating = 5000.0
#     trb.hub_height = 90.0
#     trb.bearing_number = 2
#     trb.crane = True
#
#     # Rotor force calculations for nacelle inputs
#     maxTipSpd = 80.0
#     maxEfficiency = 0.90
#
#     ratedHubPower  = trb.machine_rating*1000. / maxEfficiency
#     rotorSpeed     = (maxTipSpd/(0.5*trb.rotor_diameter)) * (60.0 / (2*np.pi))
#     trb.rotor_torque = ratedHubPower/(rotorSpeed*(np.pi/30))
#
#     trb.run()
#
#     print "The results for the NREL 5 MW Reference Turbine in an offshore 20 m water depth location are:"
#     print "Overall turbine mass is {0:.2f} kg".format(trb.turbine_mass)
#     for io in trb.list_outputs():
#         val = getattr(trb, io)
#         print io + ' ' + str(val)
#
# def cost_example():
#
#     # simple test of module
#     trb = nrel_csm_tcc_2015()
#     trb.rotor_diameter = 100.0
#     trb.turbine_class = 'II/III'
#     trb.blade_has_carbon = True
#     trb.blade_number = 3
#     trb.machine_rating = 5000.0
#     trb.hub_height = 80.0
#     trb.bearing_number = 2
#     trb.crane = True
#     trb.offshore = False
#
#     # Rotor force calculations for nacelle inputs
#     maxTipSpd = 80.0
#     maxEfficiency = 0.9
#
#     ratedHubPower  = trb.machine_rating*1000. / maxEfficiency
#     rotorSpeed     = (maxTipSpd/(0.5*trb.rotor_diameter)) * (60.0 / (2*np.pi))
#     trb.rotor_torque = ratedHubPower/(rotorSpeed*(np.pi/30))
#
#     trb.run()
#
#     print "The results for the NREL 5 MW Reference Turbine in an offshore 20 m water depth location are:"
#     print "Overall turbine mass is {0:.2f} kg".format(trb.turbine_mass)
#     print "Overall turbine cost is ${0:.2f} USD".format(trb.turbine_cost)

    # for io in trb.list_inputs():
    #     val = getattr(trb, io)
    #     print io + ' ' + str(val)
    # for io in trb.list_outputs():
    #     val = getattr(trb, io)
    #     print io + ' ' + str(val)

if __name__ == "__main__":

    #mass_example()

    # cost_example()

    prob = Problem()
    root = prob.root = Group()

    root.add('nrel_csm_tcc_2015', nrel_csm_tcc_2015(), promotes=['*'])
    prob.setup()

    prob['rotor_diameter'] = 126.4
    prob['turbine_class'] = 'II/III'
    prob['blade_has_carbon'] = True
    prob['blade_number'] = 3
    prob['machine_rating'] = 5000.
    prob['hub_height'] = 80.0
    prob['bearing_number'] = 2
    prob['crane'] = True
    # prob['offshore'] = True

    # Rotor force calculations for nacelle inputs
    maxTipSpd = 80.0
    maxEfficiency = 0.9

    ratedHubPower  = prob['machine_rating']*1000. / maxEfficiency
    rotorSpeed     = (maxTipSpd/(0.5*prob['rotor_diameter'])) * (60.0 / (2.*np.pi))
    prob['rotor_torque'] = ratedHubPower/(rotorSpeed*(np.pi/30))


    prob.run()

    # print "The results for the NREL 5 MW Reference Turbine in an offshore 20 m water depth location are:"
    # print "Overall turbine mass is {0:.2f} kg".format(prob['turbine_mass'])
    # print "Overall turbine cost is ${0:.2f} USD".format(prob['turbine_cost'])
    #
    # print "Blade Cost: ", prob['blade_cost']
    # print "Hub Cost: ", prob['hub_cost']
    # print 'Pitch System Cost: ', prob['pitch_system_cost']
    # print 'spinner cost: ', prob['spinner_cost']
    # print 'hub_system_cost: ', prob['hub_system_cost']
    # print 'rotor cost: ', prob['rotor_cost']
    # print 'lss cost: ', prob['lss_cost']
    # print 'bearings cost: ', prob['bearings_cost']
    # print 'gearbox cost: ', prob['gearbox_cost']
    # print 'hss cost: ', prob['hss_cost']
    # print 'generator cost: ', prob['generator_cost']
    # print 'bedplate cost: ', prob['bedplate_cost']
    # print 'yaw system cost: ', prob['yaw_system_cost']
    # print 'variable speed elctrics cost: ', prob['variable_speed_elec_cost']
    # print 'hydraulics cost: ', prob['hydraulic_cooling_cost']
    # print 'nacelle cover cost: ', prob['nacelle_cover_cost']
    # print 'electric connections cost: ', prob['elec_connec_cost']
    # print 'controls cost: ', prob['controls_cost']
    # print 'other costs: ', prob['other_mainframe_cost']
    # print 'Transformer costs: ', prob['transformer_cost']
    # print 'nacelle system: ', prob['nacelle_cost']
    # print 'tower cost: ', prob['tower_cost']
    # print 'towerCost: ', prob['TowerCost']

    print 'rotor cost: ', prob['rotor_cost']
    print 'nacelle cost: ', prob['nacelle_cost']
    print
    print 'rotor_diameter: ', prob['rotor_diameter']


    # for io in trb.list_inputs():
    #     val = getattr(trb, io)
    #     print io + ' ' + str(val)
    # for io in trb.list_outputs():
    #     val = getattr(trb, io)
    #     print io + ' ' + str(val)
