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

    # AEP:  57192467.9891
    # COE:  80.5873879983
    # rotor diameter 0:  159.99
    # turbineH 0:  89.995
    # rating 0:  9999.99
    # diam 0:  [ 5.16511715  5.08905243  3.87      ]
    # t 0:  [ 0.03077193  0.01961053  0.01222393]
    # hGroup:  [ 0.  0.]
    # turbineX:  [-9000000.  9000000.]
    # turbineY:  [ 9000000. -9000000.]

    """"""
    #0.5 spacing
    # AEP:  35408290.961
    # COE:  114.788690603
    # rotor diameter 0:  113.0
    # turbineH 0:  107.429567248
    # rating 0:  9999.99
    # diam 0:  [ 6.3  6.3  6.3]
    # t 0:  [ 0.02355099  0.01612411  0.00984985]
    # hGroup:  [ 0.  0.]
    # turbineX:  [-9000000.  9000000.]
    # turbineY:  [ 9000000. -9000000.]




    """0.175"""
    # AEP:  67530852.3369
    # COE:  73.2488805509
    # rotor diameter 0:  159.99
    # turbineH 0:  147.054032843
    # rating 0:  9999.99
    # diam 0:  [ 6.3         6.23996357  4.14459639]
    # t 0:  [ 0.03520666  0.02290331  0.01202712]
    # hGroup:  [ 0.  0.]
    # turbineX:  [-9000000.  9000000.]
    # turbineY:  [ 9000000. -9000000.]

    """"""
    #0.5
    # AEP:  41447119.2375
    # COE:  95.7255865188
    # rotor diameter 0:  113.0
    # turbineH 0:  114.497605088
    # rating 0:  7533.37099361
    # diam 0:  [ 6.3  6.3  6.3]
    # t 0:  [ 0.03455852  0.01509497  0.00919154]
    # hGroup:  [ 0.  0.]
    # turbineX:  [-9000000.  9000000.]
    # turbineY:  [ 9000000. -9000000.]

    # AEP:  43287664.4449
    # COE:  96.6792278452
    # rotor diameter 0:  113.0
    # turbineH 0:  113.031085579
    # rating 0:  9999.99
    # diam 0:  [ 6.3  6.3  6.3]
    # t 0:  [ 0.03727617  0.0165213   0.00985223]
    # hGroup:  [ 0.  0.]
    # turbineX:  [-9000000.  9000000.]
    # turbineY:  [ 9000000. -9000000.]




    """0.275"""
    # AEP:  79761138.7223
    # COE:  63.7390614962
    # rotor diameter 0:  159.99
    # turbineH 0:  154.540591754
    # rating 0:  9999.99
    # diam 0:  [ 6.3  6.3  6.3]
    # t 0:  [ 0.03688522  0.0228506   0.01151211]
    # hGroup:  [ 0.  0.]
    # turbineX:  [-9000000.  9000000.]
    # turbineY:  [ 9000000. -9000000.]

    """"""
    #0.5
    # AEP:  55740095.5698
    # COE:  78.6878351692
    # rotor diameter 0:  113.0
    # turbineH 0:  120.146324962
    # rating 0:  9999.99
    # diam 0:  [ 6.3  6.3  6.3]
    # t 0:  [ 0.06060468  0.01703565  0.00985564]
    # hGroup:  [ 0.  0.]
    # turbineX:  [-9000000.  9000000.]
    # turbineY:  [ 9000000. -9000000.]


    shearExp = 0.075
    num = 3

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

    # d_param = np.array([ 5.12457725,  5.04837357,  3.87      ])
    # t_param = np.array([ 0.03051178,  0.01943812,  0.01211713])
    #
    # # rotorDiameter = np.array([159.99])
    # rotorDiameter = np.array([100.])
    # turbineZ = np.array([200.])
    # # ratedPower = np.array([9999.99])
    # ratedPower = np.array([6000.99])
    #
    # rotorDiameter = np.array([113.0])
    # turbineZ = np.array([114.496593904])
    # # turbineZ = np.array([107.496593904])
    # ratedPower = np.array([7533.32849056])
    # # ratedPower = np.array([9999.99])
    # d_param = np.array([ 6.3,  6.3,  6.3])
    # t_param = np.array([ 0.03455555,  0.01509493,  0.00919154])

    rotorDiameter = np.random.rand(1)*50.+63.
    turbineZ = np.random.rand(1)*50.+63.
    ratedPower = np.random.rand(1)*5000.+2000.
    d_param = np.random.rand(3)*4.+2.3
    t_param = np.random.rand(3)*0.02+0.009

    # rotor diameter 0:  113.0
    # turbineH 0:  107.429567248
    # rating 0:  9999.99
    # diam 0:  [ 6.3  6.3  6.3]
    # t 0:  [ 0.02355099  0.01612411  0.00984985]


    # rotorDiameter = np.array([113.0])
    # turbineZ = np.array([107.429567248])
    # # ratedPower = np.array([7533.37099361])
    # ratedPower = np.array([7000.99])
    # d_param = np.array([ 6.3 , 6.3 , 6.3])
    # t_param = np.array([ 0.02355099 , 0.01612411 , 0.00984985])




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

    if shearExp == 0.075:
        prob.driver.opt_settings['Summary file'] = 'IsoSum075_05_%s.out'%num

    if shearExp == 0.175:
        prob.driver.opt_settings['Summary file'] = 'IsoSum175_05_%s.out'%num

    if shearExp == 0.275:
        prob.driver.opt_settings['Summary file'] = 'IsoSum275_05_%s.out'%num

    prob.driver.add_objective('COE', scaler=1.)

    for i in range(nGroups):
        prob.driver.add_desvar('d_param%s'%i, lower=3.87, upper=6.3, scaler=0.1)
        prob.driver.add_desvar('t_param%s'%i, lower=0.001, upper=None, scaler=1.)
        prob.driver.add_desvar('turbineH%s'%i, lower=10., scaler=0.01)
        # prob.driver.add_desvar('rotorDiameter%s'%i, lower=10., upper=159.99, scaler=1.)
        # prob.driver.add_desvar('rotorDiameter%s'%i, lower=10., upper=113., scaler=1.)
        prob.driver.add_desvar('rotorDiameter%s'%i, lower=10., upper=115., scaler=1.)
        prob.driver.add_desvar('ratedPower%s'%i, lower=500., upper=9999.99, scaler=0.00001)
        # prob.driver.add_desvar('ratedPower%s'%i, lower=500., upper=7600., scaler=0.00001)

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

    print 'Constraints: '
    print 'Tower0_max_thrust.shell_buckling: ', prob['Tower0_max_thrust.shell_buckling']
    print 'Tower0_max_speed.shell_buckling: ', prob['Tower0_max_speed.shell_buckling']
    print 'freqConstraintGroup0.freqConstraint: ', prob['freqConstraintGroup0.freqConstraint']
    print 'minHeight0.minHeight: ', prob['minHeight0.minHeight']

    print
    print
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
