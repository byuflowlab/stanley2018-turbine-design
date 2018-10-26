import numpy as np
from openmdao.api import Group, Problem, IndepVarComp, pyOptSparseDriver
from FLORISSE3D.setupOptimization import *
from FLORISSE3D.simpleTower import Tower
from FLORISSE3D.GeneralWindFarmComponents import calculate_boundary, SpacingComp,\
            BoundaryComp, get_z, getTurbineZ, AEPobj, DeMUX, hGroups, randomStart,\
            getRotorDiameter, getRatedPower, DeMUX, Myy_estimate, bladeLengthComp, minHeight
from FLORISSE3D.COE import COEGroup
from FLORISSE3D.floris import AEPGroup
from FLORISSE3D.rotorComponents import getRating, freqConstraintGroup, optCOE
from FLORISSE3D.SimpleRotorSE import SimpleRotorSE, create_rotor_functions
import os
from sys import argv

if __name__ == '__main__':

    """0.075"""
    # rotor diameter 0:  159.99
    # turbineH 0:  89.995
    # rating 0:  9999.99
    # diam 0:  [ 5.16509528  5.08902842  3.87      ]
    # t 0:  [ 0.03077212  0.01961062  0.01222393]
    """"""
    #0.5 spacing
    # rotor diameter 0:  144.0
    # turbineH 0:  123.453852558
    # rating 0:  9657.61991137
    # diam 0:  [ 6.3         6.19176708  3.87      ]
    # t 0:  [ 0.02851286  0.01933791  0.01140267]
    """0.175"""
    # rotor diameter 0:  159.99
    # turbineH 0:  149.423006514
    # rating 0:  9999.99
    # diam 0:  [ 6.3         6.3         4.90224649]
    # t 0:  [ 0.03569724  0.02243389  0.01169662]
    """"""
    #0.5
    # rotor diameter 0:  144.0
    # turbineH 0:  124.419405527
    # rating 0:  9999.99
    # diam 0:  [ 6.3   6.3   3.87]
    # t 0:  [ 0.02925928  0.01957163  0.01152075]

    """0.275"""
    # rotor diameter 0:  159.99
    # turbineH 0:  154.540469583
    # rating 0:  9999.99
    # diam 0:  [ 6.3  6.3  6.3]
    # t 0:  [ 0.03688503  0.02285054  0.01151209]
    """"""
    #0.5
    # rotor diameter 0:  144.0
    # turbineH 0:  131.400125401
    # rating 0:  9999.99
    # diam 0:  [ 6.3  6.3  6.3]
    # t 0:  [ 0.03064379  0.01974818  0.01094128]

    shearExp = 0.275

    """setup the turbine locations"""
    nGroups = 1

    turbX = np.array([-9000000.,9000000.])
    turbY = np.array([9000000.,-9000000.])
    nTurbs = len(turbX)

    minSpacing = 2.0

    nDirs = 23
    nSpeeds = 5
    dirs, freqs, speeds = alturasRose(nDirs)
    windDirections, windFrequencies, windSpeeds = setup_weibull(dirs,freqs,speeds,nSpeeds)
    nDirections = len(windDirections)

    """initial yaw values"""
    yaw = np.zeros((nDirections, nTurbs))

    nPoints = 3
    nFull = 15

    nVertices = 0

    d_param = np.array([ 5.12457725,  5.04837357,  3.87      ])
    t_param = np.array([ 0.03051178,  0.01943812,  0.01211713])

    rotorDiameter = np.array([159.99])
    turbineZ = np.array([147.019339595])
    ratedPower = np.array([9999.99])


    """OpenMDAO"""
    prob = Problem()
    root = prob.root = Group()

    # Design Variables
    for i in range(nGroups):
        root.add('ratedPower%s'%i, IndepVarComp('ratedPower%s'%i, float(ratedPower[i]), units='kW'), promotes=['*'])
        root.add('d_param%s'%i, IndepVarComp('d_param%s'%i, d_param), promotes=['*'])
        root.add('t_param%s'%i, IndepVarComp('t_param%s'%i, t_param), promotes=['*'])
        root.add('turbineH%s'%i, IndepVarComp('turbineH%s'%i, float(turbineZ[i])), promotes=['*'])
        root.add('rotorDiameter%s'%i, IndepVarComp('rotorDiameter%s'%i, float(rotorDiameter[i])), promotes=['*'])

    root.add('optCOE', optCOE(nGroups,nPoints,nFull,nTurbs,nDirections,nVertices,minSpacing),promotes=['*'])

    print 'added'

    for i in range(nGroups):
        root.connect('rotorDiameter%s'%i, 'Rotor%s.rotorDiameter'%i)
        root.connect('ratedPower%s'%i, 'Rotor%s.turbineRating'%i)

        root.connect('rotorDiameter%s'%i, 'Myy_estimate%s.rotor_diameter'%i)
        root.connect('rotorDiameter%s'%i,'bladeLengthComp%s.rotor_diameter'%i)
        root.connect('rotorDiameter%s'%i,'freqConstraintGroup%s.diameter'%i)

        root.connect('turbineH%s'%i, 'minHeight%s.height'%i)
        root.connect('rotorDiameter%s'%i, 'minHeight%s.diameter'%i)

    for i in range(nGroups):
        root.connect('d_param%s'%i, 'Tower%s_max_thrust.d_param'%i)
        root.connect('t_param%s'%i, 'Tower%s_max_thrust.t_param'%i)
        root.connect('d_param%s'%i, 'Tower%s_max_speed.d_param'%i)
        root.connect('t_param%s'%i, 'Tower%s_max_speed.t_param'%i)
        root.connect('d_param%s'%i, 'TowerDiscretization%s.d_param'%i)
        root.connect('t_param%s'%i, 'TowerDiscretization%s.t_param'%i)
    print 'connected'


    prob.driver = pyOptSparseDriver()
    prob.driver.options['optimizer'] = 'SNOPT'
    prob.driver.opt_settings['Major iterations limit'] = 1000
    prob.driver.opt_settings['Major optimality tolerance'] = 1.E-4
    prob.driver.opt_settings['Major feasibility tolerance'] = 1.E-4
    prob.driver.opt_settings['Verify level'] = 0
    prob.driver.opt_settings['Scale option'] = 1
    prob.driver.opt_settings['Scale tolerance'] = .95

    prob.driver.opt_settings['Summary file'] = 'testSum.out'

    prob.driver.add_objective('COE', scaler=100.)

    for i in range(nGroups):
        prob.driver.add_desvar('d_param%s'%i, lower=3.87, upper=6.3, scaler=0.1)
        prob.driver.add_desvar('t_param%s'%i, lower=0.001, upper=None, scaler=1.)
        prob.driver.add_desvar('turbineH%s'%i, lower=10., scaler=0.01)
        prob.driver.add_desvar('rotorDiameter%s'%i, lower=10., upper=159.99, scaler=1.)
        prob.driver.add_desvar('ratedPower%s'%i, lower=500., upper=9999.99, scaler=0.00001)

    for i in range(nGroups):
        prob.driver.add_constraint('Tower%s_max_thrust.shell_buckling'%i, upper=np.ones(nFull), scaler=1E2)
        prob.driver.add_constraint('Tower%s_max_speed.shell_buckling'%i, upper=np.ones(nFull), scaler=1E2)
        prob.driver.add_constraint('freqConstraintGroup%s.freqConstraint'%i, lower=0.0)
        prob.driver.add_constraint('minHeight%s.minHeight'%i, lower=0.0)

    prob.root.ln_solver.options['single_voi_relevance_reduction'] = True
    prob.root.ln_solver.options['mode'] = 'rev'

    prob.setup(check=False)

    setupTower(nFull, prob)
    simpleSetup(nTurbs, prob)
    prob['Uref'] = windSpeeds
    prob['windDirections'] = windDirections
    prob['windFrequencies'] = windFrequencies

    for i in range(nDirections):
        prob['yaw%s'%i] = yaw[i]
    prob['turbineX'] = turbX
    prob['turbineY'] = turbY
    prob['shearExp'] = shearExp

    for i in range(nGroups):
        prob['Tower%s_max_speed.Vel'%i] = 60.

    prob.run()
    print 'AEP: ', prob['AEP']
    print 'COE: ', prob['COE']

    print 'rotor diameter 0: ', prob['rotorDiameter0']
    print 'turbineH 0: ', prob['turbineH0']
    print 'rating 0: ', prob['ratedPower0']
    print 'diam 0: ', prob['d_param0']
    print 't 0: ', prob['t_param0']

    print 'hGroup: ', prob['hGroup']

    print 'turbineX: ', prob['turbineX']
    print 'turbineY: ', prob['turbineY']
