"""
turbine_costsse_2015.py

Created by Janine Freeman 2015 based on turbine_costsse.py 2012.
Copyright (c) NREL. All rights reserved.
"""

import numpy as np
from openmdao.api import Group, Component, Problem, ScipyGMRES, IndepVarComp

###### Rotor

#-------------------------------------------------------------------------------
class BladeCost2015(Component):

    def __init__(self):

        super(BladeCost2015, self).__init__()

        # variables
        self.add_param('blade_mass', 0.0, desc='component mass [kg]')
        self.add_param('blade_mass_cost_coeff', 14.6, desc='blade mass-cost coefficient [$/kg]') #mass-cost coefficient with default from ppt

        # Outputs
        self.add_output('blade_cost', 0.0, desc='Overall wind turbine component capital costs excluding transportation costs')

    def solve_nonlinear(self, params, unknowns, resids):

        # calculate component cost
        BladeCost2015 = params['blade_mass_cost_coeff'] * params['blade_mass']
        unknowns['blade_cost'] = BladeCost2015

    def linearize(self, params, unknowns, resids):

        J = {}
        J['blade_cost','blade_mass'] = params['blade_mass_cost_coeff']

        return J


# -----------------------------------------------------------------------------------------------
class HubCost2015(Component):

    def __init__(self):

        super(HubCost2015, self).__init__()

        # variables
        self.add_param('hub_mass', 0.0, desc='component mass [kg]')
        self.add_param('hub_mass_cost_coeff', 3.90, desc='hub mass-cost coefficient [$/kg]')

        # Outputs
        self.add_output('hub_cost', 0.0, desc='Overall wind turbine component capital costs excluding transportation costs')

    def solve_nonlinear(self, params, unknowns, resids):

        # calculate component cost
        HubCost2015 = params['hub_mass_cost_coeff'] * params['hub_mass']
        unknowns['hub_cost'] = HubCost2015

    def linearize(self, params, unknowns, resids):

        J = {}
        J['hub_cost','hub_mass'] = params['hub_mass_cost_coeff']

        return J


#-------------------------------------------------------------------------------
class PitchSystemCost2015(Component):

    def __init__(self):

        super(PitchSystemCost2015, self).__init__()

        # variables
        self.add_param('pitch_system_mass', 0.0, desc='component mass [kg]')
        self.add_param('pitch_system_mass_cost_coeff', 22.1, desc='pitch system mass-cost coefficient [$/kg]') #mass-cost coefficient with default from list

        # Outputs
        self.add_output('pitch_system_cost', 0.0, desc='Overall wind turbine component capital costs excluding transportation costs')

    def solve_nonlinear(self, params, unknowns, resids):

        #calculate system costs
        PitchSystemCost2015 = params['pitch_system_mass_cost_coeff'] * params['pitch_system_mass']
        unknowns['pitch_system_cost'] = PitchSystemCost2015

    def linearize(self, params, unknowns, resids):

        J = {}
        J['pitch_system_cost','pitch_system_mass'] = params['pitch_system_mass_cost_coeff']

        return J


#-------------------------------------------------------------------------------
class SpinnerCost2015(Component):

    def __init__(self):

        super(SpinnerCost2015, self).__init__()

        # variables
        self.add_param('spinner_mass', 0.0, desc='component mass [kg]')
        self.add_param('spinner_mass_cost_coeff', 11.1, desc='spinner/nose cone mass-cost coefficient [$/kg]') #mass-cost coefficient with default from ppt

        # Outputs
        self.add_output('spinner_cost', 0.0, desc='Overall wind turbine component capital costs excluding transportation costs')

    def solve_nonlinear(self, params, unknowns, resids):

        #calculate system costs
        SpinnerCost2015 = params['spinner_mass_cost_coeff'] * params['spinner_mass']
        unknowns['spinner_cost'] = SpinnerCost2015

    def linearize(self, params, unknowns, resids):

        J = {}
        J['spinner_cost','spinner_mass'] = params['spinner_mass_cost_coeff']

        return J


#-------------------------------------------------------------------------------
class HubSystemCostAdder2015(Component):

    def __init__(self):

        super(HubSystemCostAdder2015, self).__init__()

        # variables
        self.add_param('hub_cost', 0.0, desc='hub component cost')
        self.add_param('pitch_system_cost', 0.0, desc='pitch system cost')
        self.add_param('spinner_cost', 0.0, desc='spinner component cost')

        # multipliers
        self.add_param('hub_assemblyCostMultiplier', 0.0, desc='rotor assembly cost multiplier')
        self.add_param('hub_overheadCostMultiplier', 0.0, desc='rotor overhead cost multiplier')
        self.add_param('hub_profitMultiplier', 0.0, desc='rotor profit multiplier')
        self.add_param('hub_transportMultiplier', 0.0, desc='rotor transport multiplier')

        # Outputs
        self.add_output('hub_system_cost', 0.0, desc='Overall wind sub-assembly capial costs including transportation costs')

    def solve_nonlinear(self, params, unknowns, resids):

        partsCost = params['hub_cost']+ params['pitch_system_cost'] + params['spinner_cost']

        # updated calculations below to account for assembly, transport, overhead and profit
        unknowns['hub_system_cost'] = (1 + params['hub_transportMultiplier'] + params['hub_profitMultiplier']) * ((1 + params['hub_overheadCostMultiplier'] + params['hub_assemblyCostMultiplier']) * partsCost)

    def linearize(self, params, unknowns, resids):
        A = (1 + params['hub_transportMultiplier'] + params['hub_profitMultiplier'])
        B = (1 + params['hub_overheadCostMultiplier'] + params['hub_assemblyCostMultiplier'])
        J = {}
        J['hub_system_cost','hub_cost'] = A * B
        J['hub_system_cost','pitch_system_cost'] = A * B
        J['hub_system_cost','spinner_cost'] = A * B

        return J

#-------------------------------------------------------------------------------
class RotorCostAdder2015(Component):
    """
    RotorCostAdder adds up individual rotor system and component costs to get overall rotor cost.
    """

    def __init__(self):

        super(RotorCostAdder2015, self).__init__()

        # variables
        self.add_param('blade_cost', 0.0, desc='individual blade cost')
        self.add_param('hub_system_cost', 0.0, desc='cost for hub system')

        # parameters
        self.add_param('blade_number', 3, desc='number of rotor blades')

        # Outputs
        self.add_output('rotor_cost', 0.0, desc='Overall wind sub-assembly capial costs including transportation costs')

    def solve_nonlinear(self, params, unknowns, resids):

        unknowns['rotor_cost'] = params['blade_cost'] * params['blade_number'] + params['hub_system_cost']

    def linearize(self, params, unknowns, resids):
        J = {}
        J['rotor_cost','blade_cost'] = params['blade_number']
        J['rotor_cost','hub_system_cost'] = 1.

        return J


###### Nacelle

# -------------------------------------------------
class LowSpeedShaftCost2015(Component):

    def __init__(self):

        super(LowSpeedShaftCost2015, self).__init__()

        # variables
        self.add_param('low_speed_shaft_mass', 0.0, desc='component mass [kg]') #mass input
        self.add_param('lss_mass_cost_coeff', 11.9, desc='low speed shaft mass-cost coefficient [$/kg]') #mass-cost coefficient with default from ppt

        # Outputs
        self.add_output('lss_cost', 0.0, desc='Overall wind turbine component capital costs excluding transportation costs') #initialize cost output

    def solve_nonlinear(self, params, unknowns, resids):

        # calculate component cost
        LowSpeedShaftCost2015 = params['lss_mass_cost_coeff'] * params['low_speed_shaft_mass']
        unknowns['lss_cost'] = LowSpeedShaftCost2015 #assign the cost to this object so don't have to return anything

    def linearize(self, params, unknowns, resids):
        J = {}
        J['lss_cost','low_speed_shaft_mass'] = params['lss_mass_cost_coeff']

        return J


#-------------------------------------------------------------------------------
class BearingsCost2015(Component):

    def __init__(self):

        super(BearingsCost2015, self).__init__()

        # variables
        self.add_param('main_bearing_mass', 0.0, desc='component mass [kg]') #mass input
        self.add_param('bearing_number', 2, desc='number of main bearings []') #number of main bearings- defaults to 2
        self.add_param('bearings_mass_cost_coeff', 4.5, desc='main bearings mass-cost coefficient [$/kg]') #mass-cost coefficient- HALF of the 12.70 in powerpoint because it was based on TWO bearings

        # Outputs
        self.add_output('bearings_cost', 0.0, desc='Overall wind turbine component capital costs excluding transportation costs')

    def solve_nonlinear(self, params, unknowns, resids):

        #calculate component cost
        BearingsCost2015 = params['bearings_mass_cost_coeff'] * params['main_bearing_mass'] * params['bearing_number']
        unknowns['bearings_cost'] = BearingsCost2015

    def linearize(self, params, unknowns, resids):
        J = {}
        J['bearings_cost','main_bearing_mass'] = params['bearings_mass_cost_coeff'] * params['bearing_number']

        return J


#-------------------------------------------------------------------------------
class GearboxCost2015(Component):

    def __init__(self):

        super(GearboxCost2015, self).__init__()

        # variables
        self.add_param('gearbox_mass', 0.0, desc='component mass')
        self.add_param('gearbox_mass_cost_coeff', 12.9, desc='gearbox mass-cost coefficient [$/kg]') #mass-cost coefficient with default from ppt

        # Outputs
        self.add_output('gearbox_cost', 0.0, desc='Overall wind turbine component capital costs excluding transportation costs')

    def solve_nonlinear(self, params, unknowns, resids):

        # calculate component cost
        GearboxCost2015 = params['gearbox_mass_cost_coeff'] * params['gearbox_mass']
        unknowns['gearbox_cost'] = GearboxCost2015

    def linearize(self, params, unknowns, resids):
        J = {}
        J['gearbox_cost','gearbox_mass'] = params['gearbox_mass_cost_coeff']

        return J


#-------------------------------------------------------------------------------
class HighSpeedSideCost2015(Component):

    def __init__(self):

        super(HighSpeedSideCost2015, self).__init__()

        # variables
        self.add_param('high_speed_side_mass', 0.0, desc='component mass [kg]')
        self.add_param('high_speed_side_mass_cost_coeff', 6.8, desc='high speed side mass-cost coefficient [$/kg]') #mass-cost coefficient with default from list

        # Outputs
        self.add_output('hss_cost', 0.0, desc='Overall wind turbine component capital costs excluding transportation costs')

    def solve_nonlinear(self, params, unknowns, resids):

        # calculate component cost
        MechBrakeCost2015 = params['high_speed_side_mass_cost_coeff'] * params['high_speed_side_mass']
        unknowns['hss_cost'] = MechBrakeCost2015

    def linearize(self, params, unknowns, resids):
        J = {}
        J['hss_cost','high_speed_side_mass'] = params['high_speed_side_mass_cost_coeff']

        return J


#-------------------------------------------------------------------------------
class GeneratorCost2015(Component):

    def __init__(self):

        super(GeneratorCost2015, self).__init__()

        # variables
        self.add_param('generator_mass', 0.0, desc='component mass [kg]')
        self.add_param('generator_mass_cost_coeff', 12.4, desc='generator mass cost coefficient [$/kg]') #mass-cost coefficient with default from ppt

        # Outputs
        self.add_output('generator_cost', 0.0, desc='Overall wind turbine component capital costs excluding transportation costs')

    def solve_nonlinear(self, params, unknowns, resids):

        #calculate component cost
        GeneratorCost2015 = params['generator_mass_cost_coeff'] * params['generator_mass']
        unknowns['generator_cost'] = GeneratorCost2015

    def linearize(self, params, unknowns, resids):
        J = {}
        J['generator_cost','generator_mass'] = params['generator_mass_cost_coeff']

        return J


#-------------------------------------------------------------------------------
class BedplateCost2015(Component):

    def __init__(self):

        super(BedplateCost2015, self).__init__()

        # variables
        self.add_param('bedplate_mass', 0.0, desc='component mass [kg]')
        self.add_param('bedplate_mass_cost_coeff', 2.9, desc='bedplate mass-cost coefficient [$/kg]') #mass-cost coefficient with default from ppt

        # Outputs
        self.add_output('bedplate_cost', 0.0, desc='Overall wind turbine component capital costs excluding transportation costs')

    def solve_nonlinear(self, params, unknowns, resids):

        # calculate component cost
        BedplateCost2015 = params['bedplate_mass_cost_coeff'] * params['bedplate_mass']
        unknowns['bedplate_cost'] = BedplateCost2015

    def linearize(self, params, unknowns, resids):
        J = {}
        J['bedplate_cost','bedplate_mass'] = params['bedplate_mass_cost_coeff']

        return J


#---------------------------------------------------------------------------------
class YawSystemCost2015(Component):

    def __init__(self):

        super(YawSystemCost2015, self).__init__()

        # variables
        self.add_param('yaw_system_mass', 0.0, desc='component mass [kg]')
        self.add_param('yaw_system_mass_cost_coeff', 8.3, desc='yaw system mass cost coefficient [$/kg]') #mass-cost coefficient with default from list

        # Outputs
        self.add_output('yaw_system_cost', 0.0, desc='Overall wind turbine component capital costs excluding transportation costs')

    def solve_nonlinear(self, params, unknowns, resids):

        # calculate cost
        YawSystemCost2015 = params['yaw_system_mass_cost_coeff'] * params['yaw_system_mass']
        unknowns['yaw_system_cost'] = YawSystemCost2015

    def linearize(self, params, unknowns, resids):
        J = {}
        J['yaw_system_cost','yaw_system_mass'] = params['yaw_system_mass_cost_coeff']

        return J


#---------------------------------------------------------------------------------
class VariableSpeedElecCost2015(Component):

    def __init__(self):

        super(VariableSpeedElecCost2015, self).__init__()

        # variables
        self.add_param('variable_speed_elec_mass', 0.0, desc='component mass [kg]')
        self.add_param('variable_speed_elec_mass_cost_coeff', 18.8, desc='variable speed electronics mass cost coefficient [$/kg]') #mass-cost coefficient with default from list

        # Outputs
        self.add_output('variable_speed_elec_cost', 0.0, desc='Overall wind turbine component capital costs excluding transportation costs')

    def solve_nonlinear(self, params, unknowns, resids):

        # calculate cost
        VariableSpeedElecCost2015 = params['variable_speed_elec_mass_cost_coeff'] * params['variable_speed_elec_mass']
        unknowns['variable_speed_elec_cost'] = VariableSpeedElecCost2015

    def linearize(self, params, unknowns, resids):
        J = {}
        J['variable_speed_elec_cost','variable_speed_elec_mass'] = params['variable_speed_elec_mass_cost_coeff']

        return J


#---------------------------------------------------------------------------------
class HydraulicCoolingCost2015(Component):

    def __init__(self):

        super(HydraulicCoolingCost2015, self).__init__()

        # variables
        self.add_param('hydraulic_cooling_mass', 0.0, desc='component mass [kg]')
        self.add_param('hydraulic_cooling_mass_cost_coeff', 124., desc='hydraulic and cooling system mass cost coefficient [$/kg]') #mass-cost coefficient with default from list

        # Outputs
        self.add_output('hydraulic_cooling_cost', 0.0, desc='Overall wind turbine component capital costs excluding transportation costs')

    def solve_nonlinear(self, params, unknowns, resids):

        # calculate cost
        HydraulicCoolingCost2015 = params['hydraulic_cooling_mass_cost_coeff'] * params['hydraulic_cooling_mass']
        unknowns['hydraulic_cooling_cost'] = HydraulicCoolingCost2015

    def linearize(self, params, unknowns, resids):
        J = {}
        J['hydraulic_cooling_cost','hydraulic_cooling_mass'] = params['hydraulic_cooling_mass_cost_coeff']

        return J


#---------------------------------------------------------------------------------
class NacelleCoverCost2015(Component):

    def __init__(self):

        super(NacelleCoverCost2015, self).__init__()

        # variables
        self.add_param('nacelle_cover_mass', 0.0, desc='component mass [kg]')
        self.add_param('nacelle_cover_mass_cost_coeff', 5.7, desc='nacelle cover mass cost coefficient [$/kg]') #mass-cost coefficient with default from list

        # Outputs
        self.add_output('nacelle_cover_cost', 0.0, desc='Overall wind turbine component capital costs excluding transportation costs')

    def solve_nonlinear(self, params, unknowns, resids):

        # calculate cost
        NacelleCoverCost2015 = params['nacelle_cover_mass_cost_coeff'] * params['nacelle_cover_mass']
        unknowns['nacelle_cover_cost'] = NacelleCoverCost2015

    def linearize(self, params, unknowns, resids):
        J = {}
        J['nacelle_cover_cost','nacelle_cover_mass'] = params['nacelle_cover_mass_cost_coeff']

        return J


#---------------------------------------------------------------------------------
class ElecConnecCost2015(Component):

    def __init__(self):

        super(ElecConnecCost2015, self).__init__()

        # variables
        self.add_param('machine_rating', 0.0, desc='machine rating')
        self.add_param('elec_connec_machine_rating_cost_coeff', 41.85, desc='electrical connections cost coefficient per kW') #default from old CSM

        # Outputs
        self.add_output('elec_connec_cost', 0.0, desc='Overall wind turbine component capital costs excluding transportation costs')

    def solve_nonlinear(self, params, unknowns, resids):

        # electronic systems, hydraulics and controls
        unknowns['elec_connec_cost'] = params['elec_connec_machine_rating_cost_coeff']* params['machine_rating']

    def linearize(self, params, unknowns, resids):
        J = {}
        J['elec_connec_cost','machine_rating'] = params['elec_connec_machine_rating_cost_coeff']

        return J


#---------------------------------------------------------------------------------
class ControlsCost2015(Component):

    def __init__(self):

        super(ControlsCost2015, self).__init__()

        # variables
        self.add_param('machine_rating', 0.0, desc='machine rating')
        self.add_param('controls_machine_rating_cost_coeff', 21.15, desc='controls cost coefficient per kW') #default from old CSM

        # Outputs
        self.add_output('controls_cost', 0.0, desc='Overall wind turbine component capital costs excluding transportation costs')

    def solve_nonlinear(self, params, unknowns, resids):

        unknowns['controls_cost'] = params['machine_rating'] * params['controls_machine_rating_cost_coeff']

    def linearize(self, params, unknowns, resids):
        J = {}
        J['controls_cost','machine_rating'] = params['controls_machine_rating_cost_coeff']

        return J


#---------------------------------------------------------------------------------
class OtherMainframeCost2015(Component):

    def __init__(self):

        super(OtherMainframeCost2015, self).__init__()

        #model all three (nacelle platform, service crane, and base hardware) from old model

        # variables
        self.add_param('other_mass', 0.0, desc='component mass [kg]')
        self.add_param('nacelle_platforms_mass_cost_coeff', 17.1, desc='nacelle platforms mass cost coefficient [$/kg]') #default from old CSM
        self.add_param('crane', True, desc='flag for presence of onboard crane')
        self.add_param('crane_cost', 12000., desc='crane cost if present [$]') #default from old CSM
        self.add_param('bedplate_cost', 0.0, desc='component cost [USD]')
        self.add_param('base_hardware_cost_coeff', 0.7, desc='base hardware cost coefficient based on bedplate cost') #default from old CSM

        # Outputs
        self.add_output('other_mainframe_cost', 0.0, desc='Overall wind turbine component capital costs excluding transportation costs')

    def solve_nonlinear(self, params, unknowns, resids):

        # nacelle platform cost
        # print 'platform mass: ', params['other_mass']
        NacellePlatformsCost = params['nacelle_platforms_mass_cost_coeff'] * params['other_mass']

        # crane cost
        if (params['crane']):
            NacellePlatformsCost = params['nacelle_platforms_mass_cost_coeff'] * (params['other_mass'] - 3000.0)
            craneCost  = params['crane_cost']
        else:
            NacellePlatformsCost = params['nacelle_platforms_mass_cost_coeff'] * params['other_mass']
            craneCost  = 0.0

        # base hardware cost
        #BaseHardwareCost = self.bedplate_cost * self.base_hardware_cost_coeff

        #aggregate all three mainframe costs
        MainFrameCost = (NacellePlatformsCost + craneCost) # + BaseHardwareCost)
        unknowns['other_mainframe_cost']  = MainFrameCost

    def linearize(self, params, unknowns, resids):

        J = {}
        J['other_mainframe_cost','other_mass'] = params['nacelle_platforms_mass_cost_coeff']

        return J


#-------------------------------------------------------------------------------
class TransformerCost2015(Component):

    def __init__(self):

        super(TransformerCost2015, self).__init__()

        # variables
        self.add_param('transformer_mass', 0.0, desc='component mass [kg]')
        self.add_param('transformer_mass_cost_coeff', 18.8, desc='transformer mass cost coefficient [$/kg]') #mass-cost coefficient with default from ppt

        # Outputs
        self.add_output('transformer_cost', 0.0, desc='Overall wind turbine component capital costs excluding transportation costs')

    def solve_nonlinear(self, params, unknowns, resids):

        #calculate component cost
        TransformerCost2015 = params['transformer_mass_cost_coeff'] * params['transformer_mass']
        unknowns['transformer_cost'] = TransformerCost2015

    def linearize(self, params, unknowns, resids):

        J = {}
        J['transformer_cost','transformer_mass'] = params['transformer_mass_cost_coeff']

        return J

#-------------------------------------------------------------------------------
class NacelleSystemCostAdder2015(Component):

    def __init__(self):

        super(NacelleSystemCostAdder2015, self).__init__()

        # variables
        self.add_param('lss_cost', 0.0, desc='component cost')
        self.add_param('bearings_cost', 0.0, desc='component cost')
        self.add_param('gearbox_cost', 0.0, desc='component cost')
        self.add_param('hss_cost', 0.0, desc='component cost')
        self.add_param('generator_cost', 0.0, desc='component cost')
        self.add_param('bedplate_cost', 0.0, desc='component cost')
        self.add_param('yaw_system_cost', 0.0, desc='component cost')
        self.add_param('variable_speed_elec_cost', 0.0, desc='component cost')
        self.add_param('hydraulic_cooling_cost', 0.0, desc='component cost')
        self.add_param('nacelle_cover_cost', 0.0, desc='component cost')
        self.add_param('elec_connec_cost', 0.0, desc='component cost')
        self.add_param('controls_cost', 0.0, desc='component cost')
        self.add_param('other_mainframe_cost', 0.0, desc='component cost')
        self.add_param('transformer_cost', 0.0, desc='component cost')

        #multipliers
        self.add_param('nacelle_assemblyCostMultiplier', 0.0, desc='nacelle assembly cost multiplier')
        self.add_param('nacelle_overheadCostMultiplier', 0.0, desc='nacelle overhead cost multiplier')
        self.add_param('nacelle_profitMultiplier', 0.0, desc='nacelle profit multiplier')
        self.add_param('nacelle_transportMultiplier', 0.0, desc='nacelle transport multiplier')

        # returns
        self.add_output('nacelle_cost', 0.0, desc='component cost')

    def solve_nonlinear(self, params, unknowns, resids):

        # aggregation of nacelle costs
        partsCost = params['lss_cost'] + \
                    params['bearings_cost'] + \
                    params['gearbox_cost'] + \
                    params['hss_cost'] + \
                    params['generator_cost'] + \
                    params['bedplate_cost'] + \
                    params['yaw_system_cost'] + \
                    params['variable_speed_elec_cost'] + \
                    params['hydraulic_cooling_cost'] + \
                    params['nacelle_cover_cost'] + \
                    params['elec_connec_cost'] + \
                    params['controls_cost'] + \
                    params['other_mainframe_cost'] + \
                    params['transformer_cost']

        #apply multipliers for assembly, transport, overhead, and profits
        unknowns['nacelle_cost'] = (1. + params['nacelle_transportMultiplier'] + params['nacelle_profitMultiplier']) * ((1 + params['nacelle_overheadCostMultiplier'] + params['nacelle_assemblyCostMultiplier']) * partsCost)

    def linearize(self, params, unknowns, resids):
        A = (1. + params['nacelle_transportMultiplier'] + params['nacelle_profitMultiplier'])
        B = (1 + params['nacelle_overheadCostMultiplier'] + params['nacelle_assemblyCostMultiplier'])
        J = {}
        J['nacelle_cost','lss_cost'] = A * B
        J['nacelle_cost','bearings_cost'] = A * B
        J['nacelle_cost','gearbox_cost'] = A * B
        J['nacelle_cost','hss_cost'] = A * B
        J['nacelle_cost','generator_cost'] = A * B
        J['nacelle_cost','bedplate_cost'] = A * B
        J['nacelle_cost','yaw_system_cost'] = A * B
        J['nacelle_cost','variable_speed_elec_cost'] = A * B
        J['nacelle_cost','hydraulic_cooling_cost'] = A * B
        J['nacelle_cost','nacelle_cover_cost'] = A * B
        J['nacelle_cost','elec_connec_cost'] = A * B
        J['nacelle_cost','controls_cost'] = A * B
        J['nacelle_cost','other_mainframe_cost'] = A * B
        J['nacelle_cost','transformer_cost'] = A * B

        return J
#---------------------------------------------------------------------------------------------
# class Nacelle_CostsSE_2015(FullNacelleCostModel):
#
#     def __init__(self):
#
#         super(Nacelle_CostsSE_2015, self).__init__()
#
#         '''
#            Nacelle_CostsSE class
#               The Rotor_costsSE class is used to represent the rotor costs of a wind turbine.
#         '''
#
#         # variables
#         low_speed_shaft_mass = Float(iotype='in', units='kg', desc='component mass [kg]')
#         main_bearing_mass = Float(iotype='in', units='kg', desc='component mass [kg]') #mass input
#         gearbox_mass = Float(iotype='in', units='kg', desc='component mass')
#         high_speed_side_mass = Float(iotype='in', units='kg', desc='component mass [kg]')
#         generator_mass = Float(iotype='in', units='kg', desc='component mass [kg]')
#         bedplate_mass = Float(iotype='in', units='kg', desc='component mass [kg]')
#         yaw_system_mass = Float(iotype='in', units='kg', desc='component mass [kg]')
#         variable_speed_elec_mass = Float(iotype='in', units='kg', desc='component mass [kg]')
#         hydraulic_cooling_mass = Float(iotype='in', units='kg', desc='component mass [kg]')
#         nacelle_cover_mass = Float(iotype='in', units='kg', desc='component mass [kg]')
#         other_mass = Float(iotype='in', units='kg', desc='component mass [kg]')
#         transformer_mass = Float(iotype='in', units='kg', desc='component mass [kg]')
#
#         # parameters / high level inputs
#         machine_rating = Float(iotype='in', units='kW', desc='machine rating')
#         offshore = Bool(iotype='in', desc='flag for offshore project')
#         crane = Bool(iotype='in', desc='flag for presence of onboard crane')
#         bearing_number = Float(2, iotype='in', desc='number of main bearings []') #number of main bearings- defaults to 2
#
#         # coefficients
#         lss_mass_cost_coeff = Float(11.9, iotype='in', units='USD/kg', desc='low speed shaft mass-cost coefficient [$/kg]')
#         bearings_mass_cost_coeff = Float(4.5, iotype='in', units='USD/kg', desc='main bearings mass-cost coefficient [$/kg]') #mass-cost coefficient- HALF of the 12.70 in powerpoint because it was based on TWO bearings
#         gearbox_mass_cost_coeff = Float(12.9, iotype='in', units='USD/kg', desc='gearbox mass-cost coefficient [$/kg]')
#         high_speed_side_mass_cost_coeff = Float(6.8, iotype='in', units='USD/kg', desc='high speed side mass-cost coefficient [$/kg]') #mass-cost coefficient with default from list
#         generator_mass_cost_coeff = Float(12.4, iotype='in', units= 'USD/kg', desc='generator mass cost coefficient [$/kg]')
#         bedplate_mass_cost_coeff = Float(2.9, iotype='in', units='USD/kg', desc='bedplate mass-cost coefficient [$/kg]') #mass-cost coefficient with default from ppt
#         yaw_system_mass_cost_coeff = Float(8.3, iotype='in', units='USD/kg', desc='yaw system mass cost coefficient [$/kg]') #mass-cost coefficient with default from list
#         variable_speed_elec_mass_cost_coeff = Float(18.8, iotype='in', units='USD/kg', desc='variable speed electronics mass cost coefficient [$/kg]') #mass-cost coefficient with default from list
#         hydraulic_cooling_mass_cost_coeff = Float(124., iotype='in', units='USD/kg', desc='hydraulic and cooling system mass cost coefficient [$/kg]') #mass-cost coefficient with default from list
#         nacelle_cover_mass_cost_coeff = Float(5.7, iotype='in', units='USD/kg', desc='nacelle cover mass cost coefficient [$/kg]') #mass-cost coefficient with default from list
#         elec_connec_machine_rating_cost_coeff = Float(41.85, iotype='in', units='USD/kW', desc='2002 electrical connections cost coefficient per kW')
#         controls_machine_rating_cost_coeff = Float(21.15, iotype='in', units='USD/kW', desc='controls cost coefficient per kW') #default from old CSM
#         nacelle_platforms_mass_cost_coeff = Float(17.1, iotype='in', units='USD/kg', desc='nacelle platforms mass cost coefficient [$/kg]') #default from old CSM
#         crane_cost = Float(12000.0, iotype='in', units='USD', desc='crane cost if present [$]') #default from old CSM
#         base_hardware_cost_coeff = Float(0.7, iotype='in', desc='base hardware cost coefficient based on bedplate cost') #default from old CSM
#         transformer_mass_cost_coeff = Float(18.8, iotype='in', units= 'USD/kg', desc='transformer mass cost coefficient [$/kg]') #mass-cost coefficient with default from ppt
#
#         #multipliers
#         nacelle_assemblyCostMultiplier = Float(0.0, iotype='in', desc='nacelle assembly cost multiplier')
#         nacelle_overheadCostMultiplier = Float(0.0, iotype='in', desc='nacelle overhead cost multiplier')
#         nacelle_profitMultiplier = Float(0.0, iotype='in', desc='nacelle profit multiplier')
#         nacelle_transportMultiplier = Float(0.0, iotype='in', desc='nacelle transport multiplier')
#
#         # outputs
#         cost = Float(iotype='out', units='USD', desc='component cost')
#
#     def configure(self):
#
#         configure_full_ncc(self)
#
#         # select components #this creates an instance of each of these components with the name in ''
#         self.replace('lssCC', LowSpeedShaftCost2015())
#         self.replace('bearingsCC', BearingsCost2015())
#         self.replace('gearboxCC', GearboxCost2015())
#         self.replace('hssCC', HighSpeedSideCost2015())
#         self.replace('generatorCC', GeneratorCost2015())
#         self.replace('bedplateCC', BedplateCost2015())
#         self.replace('yawSysCC', YawSystemCost2015())
#         self.add('vsCC', VariableSpeedElecCost2015())
#         self.add('hydraulicCC', HydraulicCoolingCost2015())
#         self.add('nacelleCC', NacelleCoverCost2015())
#         self.add('elecCC', ElecConnecCost2015())
#         self.add('controlsCC', ControlsCost2015())
#         self.add('mainframeCC', OtherMainframeCost2015())
#         self.add('transformerCC', TransformerCost2015())
#         self.replace('ncc', NacelleSystemCostAdder2015())
#
#         self.driver.workflow.add(['vsCC', 'hydraulicCC', 'nacelleCC', 'elecCC', 'controlsCC', 'mainframeCC', 'transformerCC'])
#
#         # connect inputs
#         self.connect('low_speed_shaft_mass', 'lssCC.low_speed_shaft_mass')
#         self.connect('lss_mass_cost_coeff', 'lssCC.lss_mass_cost_coeff')
#         self.connect('main_bearing_mass', 'bearingsCC.main_bearing_mass')
#         self.connect('bearing_number', 'bearingsCC.bearing_number')
#         self.connect('bearings_mass_cost_coeff', 'bearingsCC.bearings_mass_cost_coeff')
#         self.connect('gearbox_mass', 'gearboxCC.gearbox_mass')
#         self.connect('gearbox_mass_cost_coeff', 'gearboxCC.gearbox_mass_cost_coeff')
#         self.connect('high_speed_side_mass', 'hssCC.high_speed_side_mass')
#         self.connect('high_speed_side_mass_cost_coeff', 'hssCC.high_speed_side_mass_cost_coeff')
#         self.connect('generator_mass', 'generatorCC.generator_mass')
#         self.connect('generator_mass_cost_coeff', 'generatorCC.generator_mass_cost_coeff')
#         self.connect('bedplate_mass', ['bedplateCC.bedplate_mass'])
#         self.connect('bedplate_mass_cost_coeff', 'bedplateCC.bedplate_mass_cost_coeff')
#         self.connect('yaw_system_mass', 'yawSysCC.yaw_system_mass')
#         self.connect('yaw_system_mass_cost_coeff', 'yawSysCC.yaw_system_mass_cost_coeff')
#         self.connect('variable_speed_elec_mass', 'vsCC.variable_speed_elec_mass')
#         self.connect('variable_speed_elec_mass_cost_coeff', 'vsCC.variable_speed_elec_mass_cost_coeff')
#         self.connect('hydraulic_cooling_mass', 'hydraulicCC.hydraulic_cooling_mass')
#         self.connect('hydraulic_cooling_mass_cost_coeff', 'hydraulicCC.hydraulic_cooling_mass_cost_coeff')
#         self.connect('nacelle_cover_mass', 'nacelleCC.nacelle_cover_mass')
#         self.connect('nacelle_cover_mass_cost_coeff', 'nacelleCC.nacelle_cover_mass_cost_coeff')
#         self.connect('machine_rating','elecCC.machine_rating')
#         self.connect('elec_connec_machine_rating_cost_coeff', 'elecCC.elec_connec_machine_rating_cost_coeff')
#         self.connect('controls_machine_rating_cost_coeff', 'controlsCC.controls_machine_rating_cost_coeff')
#         self.connect('machine_rating','controlsCC.machine_rating')
#         self.connect('other_mass', 'mainframeCC.nacelle_platforms_mass')
#         self.connect('nacelle_platforms_mass_cost_coeff', 'mainframeCC.nacelle_platforms_mass_cost_coeff')
#         self.connect('crane', 'mainframeCC.crane')
#         self.connect('crane_cost', 'mainframeCC.crane_cost')
#         self.connect('base_hardware_cost_coeff', 'mainframeCC.base_hardware_cost_coeff')
#         self.connect('transformer_mass', 'transformerCC.transformer_mass')
#         self.connect('transformer_mass_cost_coeff', 'transformerCC.transformer_mass_cost_coeff')
#
#         # internal connections
#         self.connect('bedplateCC.cost', 'mainframeCC.bedplate_cost')
#         self.connect('vsCC.cost','ncc.variable_speed_elec_cost')
#         self.connect('hydraulicCC.cost','ncc.hydraulic_cooling_cost')
#         self.connect('nacelleCC.cost','ncc.nacelle_cover_cost')
#         self.connect('elecCC.cost','ncc.elec_connec_cost')
#         self.connect('controlsCC.cost','ncc.controls_cost')
#         self.connect('mainframeCC.cost','ncc.other_mainframe_cost')
#         self.connect('transformerCC.cost','ncc.transformer_cost')
#
#         # connect multipliers
#         self.connect('nacelle_assemblyCostMultiplier', 'ncc.nacelle_assemblyCostMultiplier')
#         self.connect('nacelle_overheadCostMultiplier', 'ncc.nacelle_overheadCostMultiplier')
#         self.connect('nacelle_profitMultiplier', 'ncc.nacelle_profitMultiplier')
#         self.connect('nacelle_transportMultiplier', 'ncc.nacelle_transportMultiplier')
#
#
# ###### Tower

#-------------------------------------------------------------------------------
class TowerCost2015(Component):

    def __init__(self):

        super(TowerCost2015, self).__init__()

        # variables
        self.add_param('tower_mass', 0.0, desc='tower mass [kg]')
        self.add_param('tower_mass_cost_coefficient', 2.9, desc='tower mass-cost coefficient [$/kg]') #mass-cost coefficient with default from ppt

        # Outputs
        self.add_output('tower_cost', 0.0, desc='Overall wind turbine component capital costs excluding transportation costs')

    def solve_nonlinear(self, params, unknowns, resids):

        # calculate component cost
        TowerCost2015 = params['tower_mass_cost_coefficient'] * params['tower_mass']
        unknowns['tower_cost'] = TowerCost2015

    def linearize(self, params, unknowns, resids):

        J = {}
        J['tower_cost','tower_mass'] = params['tower_mass_cost_coefficient']

        return J


#-------------------------------------------------------------------------------
class TowerCostAdder2015(Component):

    def __init__(self):

        super(TowerCostAdder2015, self).__init__()

        # variables
        self.add_param('tower_cost', 0.0, desc='component cost')

        # multipliers
        self.add_param('tower_assemblyCostMultiplier', 0.0, desc='tower assembly cost multiplier')
        self.add_param('tower_overheadCostMultiplier', 0.0, desc='tower overhead cost multiplier')
        self.add_param('tower_profitMultiplier', 0.0, desc='tower profit cost multiplier')
        self.add_param('tower_transportMultiplier', 0.0, desc='tower transport cost multiplier')

        # returns
        self.add_output('TowerCost', 0.0, desc='component cost')

    def solve_nonlinear(self, params, unknowns, resids):

        partsCost = params['tower_cost']
        unknowns['TowerCost'] = (1 + params['tower_transportMultiplier'] + params['tower_profitMultiplier']) * ((1 + params['tower_overheadCostMultiplier'] + params['tower_assemblyCostMultiplier']) * partsCost)

    def linearize(self, params, unknowns, resids):

        J = {}
        J['TowerCost','tower_cost'] = (1 + params['tower_transportMultiplier'] + params['tower_profitMultiplier']) * ((1 + params['tower_overheadCostMultiplier'] + params['tower_assemblyCostMultiplier']))

        return J

#-------------------------------------------------------------------------------
class Turbine_CostsSE_2015(Group):

    def __init__(self):

        super(Turbine_CostsSE_2015, self).__init__()

        self.add('BladeCost2015', BladeCost2015(), promotes=['*'])
        self.add('HubCost2015', HubCost2015(), promotes=['*'])
        self.add('PitchSystemCost2015', PitchSystemCost2015(), promotes=['*'])
        self.add('SpinnerCost2015', SpinnerCost2015(), promotes=['*'])
        self.add('HubSystemCostAdder2015', HubSystemCostAdder2015(), promotes=['*'])
        self.add('RotorCostAdder2015', RotorCostAdder2015(), promotes=['*'])
        self.add('LowSpeedShaftCost2015', LowSpeedShaftCost2015(), promotes=['*'])
        self.add('BearingsCost2015', BearingsCost2015(), promotes=['*'])
        self.add('GearboxCost2015', GearboxCost2015(), promotes=['*'])
        self.add('HighSpeedSideCost2015', HighSpeedSideCost2015(), promotes=['*'])
        self.add('GeneratorCost2015', GeneratorCost2015(), promotes=['*'])
        self.add('BedplateCost2015', BedplateCost2015(), promotes=['*'])
        self.add('YawSystemCost2015', YawSystemCost2015(), promotes=['*'])
        self.add('VariableSpeedElecCost2015', VariableSpeedElecCost2015(), promotes=['*'])
        self.add('HydraulicCoolingCost2015', HydraulicCoolingCost2015(), promotes=['*'])
        self.add('NacelleCoverCost2015', NacelleCoverCost2015(), promotes=['*'])
        self.add('ElecConnecCost2015', ElecConnecCost2015(), promotes=['*'])
        self.add('ControlsCost2015', ControlsCost2015(), promotes=['*'])
        self.add('OtherMainframeCost2015', OtherMainframeCost2015(), promotes=['*'])
        self.add('TransformerCost2015', TransformerCost2015(), promotes=['*'])
        self.add('NacelleSystemCostAdder2015', NacelleSystemCostAdder2015(), promotes=['*'])
        self.add('TowerCost2015', TowerCost2015(), promotes=['*'])
        self.add('TowerCostAdder2015', TowerCostAdder2015(), promotes=['*'])
        # self.add('Tower_CostsSE_2015', Tower_CostsSE_2015(), promotes=['*'])
        self.add('TurbineCostAdder2015', TurbineCostAdder2015(), promotes=['*'])



#-------------------------------------------------------------------------------
class TurbineCostAdder2015(Component):

    def __init__(self):

        super(TurbineCostAdder2015, self).__init__()

        # Variables
        self.add_param('rotor_cost', 0.0, desc='rotor cost')
        self.add_param('nacelle_cost', 0.0, desc='nacelle cost')
        self.add_param('TowerCost', 0.0, desc='tower cost')

        # parameters
        self.add_param('turbine_assemblyCostMultiplier', 0.0, desc='turbine multiplier for assembly cost in manufacturing')
        self.add_param('turbine_overheadCostMultiplier', 0.0, desc='turbine multiplier for overhead')
        self.add_param('turbine_profitMultiplier', 0.0, desc='turbine multiplier for profit markup')
        self.add_param('turbine_transportMultiplier', 0.0, desc='turbine multiplier for transport costs')

        # Outputs
        self.add_output('turbine_cost', 0.0, desc='Overall wind turbine capial costs including transportation costs')

    def solve_nonlinear(self, params, unknowns, resids):

        partsCost = params['rotor_cost'] + params['nacelle_cost'] + params['TowerCost']

        unknowns['turbine_cost'] = (1. + params['turbine_transportMultiplier'] + params['turbine_profitMultiplier']) * ((1. + params['turbine_overheadCostMultiplier'] + params['turbine_assemblyCostMultiplier']) * partsCost)

    def linearize(self, params, unknowns, resids):
        A = (1. + params['turbine_transportMultiplier'] + params['turbine_profitMultiplier'])
        B = (1. + params['turbine_overheadCostMultiplier'] + params['turbine_assemblyCostMultiplier'])

        J = {}
        J['turbine_cost','rotor_cost'] = A * B
        J['turbine_cost','nacelle_cost'] = A * B
        J['turbine_cost','TowerCost'] = A * B

        return J

#-------------------------------------------------------------------------------

def example():

    # simple test of module

    turbine = Turbine_CostsSE_2015()

    turbine.blade_mass = 17650.67  # inline with the windpact estimates
    turbine.hub_mass = 31644.5
    turbine.pitch_system_mass = 17004.0
    turbine.spinner_mass = 1810.5
    turbine.low_speed_shaft_mass = 31257.3
    #bearingsMass = 9731.41
    turbine.main_bearing_mass = 9731.41 / 2
    turbine.second_bearing_mass = 9731.41 / 2 #KLD - revisit this in new model
    turbine.gearbox_mass = 30237.60
    turbine.high_speed_side_mass = 1492.45
    turbine.generator_mass = 16699.85
    turbine.bedplate_mass = 93090.6
    turbine.yaw_system_mass = 11878.24
    turbine.tower_mass = 434559.0
    turbine.variable_speed_elec_mass = 1000. #Float(iotype='in', units='kg', desc='component mass [kg]')
    turbine.hydraulic_cooling_mass = 1000. #Float(iotype='in', units='kg', desc='component mass [kg]')
    turbine.nacelle_cover_mass = 1000. #Float(iotype='in', units='kg', desc='component mass [kg]')
    turbine.nacelle_platforms_mass = 1000. #Float(iotype='in', units='kg', desc='component mass [kg]')
    turbine.transformer_mass = 1000. #Float(iotype='in', units='kg', desc='component mass [kg]')

    # other inputs
    turbine.machine_rating = 5000.0
    turbine.blade_number = 3
    turbine.crane = True
    turbine.offshore = True
    turbine.bearing_number = 2

    turbine.run()

    print "The results for the NREL 5 MW Reference Turbine in an offshore 20 m water depth location are:"
    print
    print "Overall rotor cost with 3 advanced blades is ${0:.2f} USD".format(turbine.rotorCC.cost)
    print "Blade cost is ${0:.2f} USD".format(turbine.rotorCC.bladeCC.cost)
    print "Hub cost is ${0:.2f} USD".format(turbine.rotorCC.hubCC.cost)
    print "Pitch system cost is ${0:.2f} USD".format(turbine.rotorCC.pitchSysCC.cost)
    print "Spinner cost is ${0:.2f} USD".format(turbine.rotorCC.spinnerCC.cost)
    print
    print "Overall nacelle cost is ${0:.2f} USD".format(turbine.nacelleCC.cost)
    print "LSS cost is ${0:.2f} USD".format(turbine.nacelleCC.lssCC.cost)
    print "Main bearings cost is ${0:.2f} USD".format(turbine.nacelleCC.bearingsCC.cost)
    print "Gearbox cost is ${0:.2f} USD".format(turbine.nacelleCC.gearboxCC.cost)
    print "High speed side cost is ${0:.2f} USD".format(turbine.nacelleCC.hssCC.cost)
    print "Generator cost is ${0:.2f} USD".format(turbine.nacelleCC.generatorCC.cost)
    print "Bedplate cost is ${0:.2f} USD".format(turbine.nacelleCC.bedplateCC.cost)
    print "Yaw system cost is ${0:.2f} USD".format(turbine.nacelleCC.yawSysCC.cost)
    print "Variable speed electronics cost is ${0:.2f} USD".format(turbine.nacelleCC.vsCC.cost)
    print "HVAC cost is ${0:.2f} USD".format(turbine.nacelleCC.hydraulicCC.cost)
    print "Electrical connections cost is ${0:.2f} USD".format(turbine.nacelleCC.elecCC.cost)
    print "Controls cost is ${0:.2f} USD".format(turbine.nacelleCC.controlsCC.cost)
    print "Mainframe cost is ${0:.2f} USD".format(turbine.nacelleCC.mainframeCC.cost)
    print "Transformer cost is ${0:.2f} USD".format(turbine.nacelleCC.transformerCC.cost)
    print
    print "Tower cost is ${0:.2f} USD".format(turbine.towerCC.cost)
    print
    print "The overall turbine cost is ${0:.2f} USD".format(turbine.turbine_cost)
    print


if __name__ == "__main__":

    example()
