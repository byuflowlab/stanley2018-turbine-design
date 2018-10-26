import numpy as np
from math import pi, log
from openmdao.api import Group, Component, Problem, ScipyGMRES, IndepVarComp, ExecComp
from FLORISSE3D.simpleTower import Tower
from FLORISSE3D.GeneralWindFarmComponents import calculate_boundary, SpacingComp,\
            BoundaryComp, get_z, getTurbineZ, AEPobj, DeMUX, hGroups, randomStart,\
            getRotorDiameter, getRatedPower, DeMUX, Myy_estimate, bladeLengthComp,\
	    minHeight,SpacingConstraint
from FLORISSE3D.SimpleRotorSE import SimpleRotorSE, create_rotor_functions
from FLORISSE3D.COE import COEGroup

class getRating(Component):
    """
    find turbine rating from rotor diameter
    """

    def __init__(self, nTurbines):

        super(getRating, self).__init__()

        # self.deriv_options['form'] = 'forward'
        # self.deriv_options['step_size'] = 1.E-6
        # self.deriv_options['step_calc'] = 'relative'

        self.nTurbines = nTurbines
        self.add_param('rotorDiameter', np.zeros(nTurbines), desc='array of rotor diameters')

        self.add_output('ratedPower', np.zeros(nTurbines), units='kW',  desc='rated power array')
        for i in range(nTurbines):
            self.add_output('rated_powers%s'%i, 0.0, units='kW', desc='rated power of each turbine')


    def solve_nonlinear(self, params, unknowns, resids):

        rotorDiameter = params['rotorDiameter']
        ratedPower = np.zeros(self.nTurbines)

        for i in range(self.nTurbines):
            ratedPower[i] = 5000.*rotorDiameter[i]**2/126.4**2
            unknowns['rated_powers%s'%i] = 5000.*rotorDiameter[i]**2/126.4**2

        unknowns['ratedPower'] = ratedPower

    def linearize(self, params, unknowns, resids):

        J = {}
        J['ratedPower', 'rotorDiameter'] = np.zeros((self.nTurbines, self.nTurbines))
        for i in range(nTurbines):
            J['ratedPower', 'rotorDiameter'][i][i] = 2.*5000.*params['rotorDiameter'][i]/126.4**2

        return J


class getMinFreq(Component):
    """
    linear relation for the frequency constraint/rotation rate of the turbine wrt rotor diameter
    """

    def __init__(self):

        super(getMinFreq, self).__init__()
        self.add_param('diameter', 0.0, desc='rotor diameter')
        self.add_output('minFreq', 0.0, desc='frequency constraint')


    def solve_nonlinear(self, params, unknowns, resids):
        unknowns['minFreq'] = -0.002305*params['diameter']+0.4873

    def linearize(self, params, unknowns, resids):
        J = {}
        J['minFreq', 'diameter'] = -0.002305
        return J


class freqConstraint(Component):
    """
    linear relation for the frequency constraint/rotation rate of the turbine wrt rotor diameter
    """

    def __init__(self):

        super(freqConstraint, self).__init__()
        self.add_param('freq', 0.0, desc='frequency')
        self.add_param('minFreq', 0.0, desc='upper bound for freq')
        self.add_output('freqConstraint', 0.0, desc='frequency constraint')


    def solve_nonlinear(self, params, unknowns, resids):
        unknowns['freqConstraint'] = params['freq']-1.1*params['minFreq']

    def linearize(self, params, unknowns, resids):

        J = {}
        J['freqConstraint', 'freq'] = 1.
        J['freqConstraint', 'minFreq'] = -1.1
        return J

class freqConstraintGroup(Group):
    """
    linear relation for the frequency constraint/rotation rate of the turbine wrt rotor diameter
    """

    def __init__(self):

        super(freqConstraintGroup, self).__init__()

        self.add('getMinFreq', getMinFreq(), promotes=['*'])
        self.add('freqConstraint', freqConstraint(), promotes=['*'])


class optCOE(Group):

    def __init__(self,nGroups, nPoints, nFull, nTurbs, nDirections, nVertices, minSpacing):

        super(optCOE, self).__init__()

        interp_spline_ratedQ, interp_spline_blade_mass, interp_spline_Vrated, interp_spline_I1, interp_spline_I2, interp_spline_I3, interp_spline_ratedT, interp_spline_extremeT = create_rotor_functions()

        for i in range(nGroups):
            self.add('get_z_param%s'%i, get_z(nPoints)) #have derivatives
            self.add('get_z_full%s'%i, get_z(nFull)) #have derivatives
            self.add('Tower%s_max_thrust'%i, Tower(nPoints, nFull), promotes=['L_reinforced','mrhox','E','sigma_y','gamma_f','gamma_b','rhoAir','z0','zref','shearExp','rho'])
            self.add('Tower%s_max_speed'%i, Tower(nPoints, nFull), promotes=['L_reinforced','mrhox','E','sigma_y','gamma_f','gamma_b','rhoAir','z0','zref','shearExp','rho'])
            self.add('bladeLengthComp%s'%i, bladeLengthComp()) #have derivatives
            self.add('minHeight%s'%i, minHeight()) #have derivatives
            self.add('freqConstraintGroup%s'%i, freqConstraintGroup())


            self.add('Rotor%s'%i, SimpleRotorSE(interp_spline_ratedQ, interp_spline_blade_mass, interp_spline_Vrated, interp_spline_I1, interp_spline_I2, interp_spline_I3, interp_spline_ratedT, interp_spline_extremeT))
            self.add('split_I%s'%i, DeMUX(6)) #have derivatives
            self.add('Myy_estimate%s'%i, Myy_estimate()) #have derivatives

        self.add('Zs', DeMUX(nTurbs)) #have derivatives
        self.add('hGroups', hGroups(nTurbs, nGroups), promotes=['*']) #have derivatives
        self.add('getRotorDiameter', getRotorDiameter(nTurbs, nGroups), promotes=['*']) #have derivatives
        self.add('getRatedPower', getRatedPower(nTurbs, nGroups), promotes=['*'])    #have derivatives

        self.add('COEGroup', COEGroup(nTurbs, nGroups, nDirections, nPoints, nFull), promotes=['*']) #TODO check derivatives?




        self.connect('turbineZ', 'Zs.Array')

        for i in range(nGroups):
            self.connect('Rotor%s.ratedQ'%i, 'rotor_nacelle_costs%s.rotor_torque'%i)

            self.connect('Rotor%s.blade_mass'%i, 'rotor_nacelle_costs%s.blade_mass'%i)
            self.connect('Rotor%s.Vrated'%i,'Tower%s_max_thrust.Vel'%i)
            self.connect('Rotor%s.I'%i, 'split_I%s.Array'%i)
            self.connect('split_I%s.output%s'%(i,2),'Tower%s_max_speed.It'%i)
            self.connect('Tower%s_max_speed.It'%i,'Tower%s_max_thrust.It'%i)
            self.connect('Rotor%s.ratedT'%i,'Tower%s_max_thrust.Fx'%i)
            self.connect('Rotor%s.extremeT'%i,'Tower%s_max_speed.Fx'%i)

            self.connect('Myy_estimate%s.Myy'%i,'Tower%s_max_thrust.Myy'%i)
            self.connect('Myy_estimate%s.Myy'%i,'Tower%s_max_speed.Myy'%i)

            self.connect('Tower%s_max_thrust.freq'%i,'freqConstraintGroup%s.freq'%i)

        for i in range(nGroups):
            self.connect('rotor_nacelle_costs%s.rotor_mass'%i, 'Tower%s_max_speed.rotor_mass'%i)
            self.connect('rotor_nacelle_costs%s.nacelle_mass'%i, 'Tower%s_max_speed.nacelle_mass'%i)

            self.connect('Tower%s_max_speed.rotor_mass'%i, 'Tower%s_max_thrust.rotor_mass'%i)
            self.connect('Tower%s_max_speed.nacelle_mass'%i, 'Tower%s_max_thrust.nacelle_mass'%i)

        # for j in range(nGroups):
        #     self.connect('rotor_diameters%s'%j,'rotor_nacelle_costs%s.rotor_diameter'%j)
        #     self.connect('rated_powers%s'%j,'rotor_nacelle_costs%s.machine_rating'%j)

        for j in range(nGroups):
            self.connect('rotorDiameter%s'%j,'rotor_nacelle_costs%s.rotor_diameter'%j)
            self.connect('ratedPower%s'%j,'rotor_nacelle_costs%s.machine_rating'%j)

        for i in range(nGroups):
            self.connect('get_z_param%s.z_param'%i, 'Tower%s_max_thrust.z_param'%i)
            self.connect('get_z_full%s.z_param'%i, 'Tower%s_max_thrust.z_full'%i)
            self.connect('get_z_param%s.z_param'%i, 'Tower%s_max_speed.z_param'%i)
            self.connect('get_z_full%s.z_param'%i, 'Tower%s_max_speed.z_full'%i)

            self.connect('Zs.output%s'%i, 'get_z_param%s.turbineZ'%i)
            self.connect('Zs.output%s'%i, 'get_z_full%s.turbineZ'%i)
            self.connect('Zs.output%s'%i, 'Tower%s_max_thrust.L'%i)
            self.connect('Zs.output%s'%i, 'Tower%s_max_speed.L'%i)

            self.connect('get_z_param%s.z_param'%i, 'Tower%s_max_thrust.z_param'%i)
            self.connect('get_z_full%s.z_param'%i, 'Tower%s_max_thrust.z_full'%i)

            self.connect('get_z_param%s.z_param'%i, 'TowerDiscretization%s.z_param'%i)
            self.connect('get_z_full%s.z_param'%i, 'TowerDiscretization%s.z_full'%i)
            self.connect('rho', 'calcMass%s.rho'%i)


	self.add('spacing_comp', SpacingComp(nTurbines=nTurbs), promotes=['*'])

        # add constraint definitions
	self.add('spacing_con', SpacingConstraint(nTurbs), promotes=['*'])

        if nVertices > 0:
            # add component that enforces a convex hull wind farm boundary
            self.add('boundary_con', BoundaryComp(nVertices=nVertices, nTurbines=nTurbs), promotes=['*'])

if __name__=="__main__":
    """
    This is just to test during development
    """
