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

    shearExp = float(argv[1])/1000.
    spacing_multiplier = float(argv[2])/10.
    jobID = int(argv[3])
    taskID = int(argv[4])

    """setup the turbine locations"""
    nGroups = 1

    rotor_diameter = 80.
    w = 8.*rotor_diameter
    h = w/2.
    turbineX = np.array([-h,h,-w,w,-w-h,w+h,-2.*w,2.*w,w,-w,h,-h,0.,0.,-h,h,-w,w,-w-h,
            w+h,2.*w,-2.*w,w+h,-w-h,w,-w,h,-h,0.,0.,-h,h])
    turbineY = np.array([-2.*w,2.*w,-w-h,w+h,-w,w,-h,h,-w-h,w+h,-w,w,-h,h,0.,0.,h,-h,w,
            -w,-h,h,0.,0.,h,-h,w,-w,w+h,-w-h,2.*w,-2.*w])
    turbineX = turbineX*spacing_multiplier
    turbineY = turbineY*spacing_multiplier

    nTurbs = len(turbineX)

    if spacing_multiplier == 0.5:
        turbineXpert = np.random.rand(nTurbs)*rotor_diameter*1.0 - 0.5*rotor_diameter
        turbineYpert = np.random.rand(nTurbs)*rotor_diameter*1.0 - 0.5*rotor_diameter
    else:
        turbineXpert = np.random.rand(nTurbs)*rotor_diameter*2.8 - 1.4*rotor_diameter
        turbineYpert = np.random.rand(nTurbs)*rotor_diameter*2.8 - 1.4*rotor_diameter

    turbineXstart = turbineX + turbineXpert*spacing_multiplier
    turbineYstart = turbineY + turbineYpert*spacing_multiplier

    if spacing_multiplier == 0.5:
        turbineXstart = turbineX
        turbineYstart = turbineY

    R = np.sqrt((2.*w*spacing_multiplier)**2+(h*spacing_multiplier)**2)
    print 'R: ', R

    nDirs = 23
    nSpeeds = 5
    dirs, freqs, speeds = alturasRose(nDirs)
    windDirections, windFrequencies, windSpeeds = setup_weibull(dirs,freqs,speeds,nSpeeds)
    nDirections = len(windDirections)

    """initial yaw values"""
    yaw = np.zeros((nDirections, nTurbs))

    nPoints = 3
    nFull = 15

    minSpacing = 2.0

    if shearExp == 0.075:
        if spacing_multiplier > 0.5:
            d_param = np.array([ 5.16511715,  5.08905243,  3.87      ])
            t_param = np.array([ 0.03077193,  0.01961053,  0.01222393])
            rotorDiameter = np.array([159.99])
            turbineZ = np.array([89.995])
            ratedPower = np.array([9999.99])
            print 'shearExp: ', shearExp
        elif spacing_multiplier == 0.5:
            d_param = np.array([ 6.3 , 6.3,  6.3    ])
            t_param = np.array([ 0.02355099,  0.01612411,  0.00984985])
            rotorDiameter = np.array([113.])
            turbineZ = np.array([107.429567248])
            ratedPower = np.array([9999.99])

    elif shearExp == 0.175:
        if spacing_multiplier > 0.5:
            d_param = np.array([ 6.3       ,  6.23996357,  4.14459639])
            t_param = np.array([ 0.03520666,  0.02290331,  0.01202712])
            rotorDiameter = np.array([159.99])
            turbineZ = np.array([147.054032843])
            ratedPower = np.array([9999.99])
            print 'shearExp: ', shearExp
        elif spacing_multiplier == 0.5:
            d_param = np.array([ 6.3 , 6.3,  6.3    ])
            t_param = np.array([0.03727617,  0.0165213 ,  0.00985223])
            rotorDiameter = np.array([113.])
            turbineZ = np.array([113.031085579])
            ratedPower = np.array([9999.99])

    elif shearExp == 0.275:
        if spacing_multiplier > 0.5:
            d_param = np.array([ 6.3,  6.3,  6.3])
            t_param = np.array([ 0.03688522,  0.0228506,   0.01151211])
            rotorDiameter = np.array([159.99])
            turbineZ = np.array([154.540591754])
            ratedPower = np.array([9999.99])
            print 'shearExp: ', shearExp
        elif spacing_multiplier == 0.5:
            d_param = np.array([ 6.3 , 6.3,  6.3    ])
            t_param = np.array([0.06060468,  0.01703565,  0.00985564])
            rotorDiameter = np.array([113.])
            turbineZ = np.array([120.146324962])
            ratedPower = np.array([9999.99])

    nVertices = 1

    try:
        os.makedirs('sequential')
    except OSError:
        if not os.path.isdir('sequential'):
            raise

    try:
        os.makedirs('sequential/%s'%shearExp)
    except OSError:
        if not os.path.isdir('sequential/%s'%shearExp):
            raise

    try:
        os.makedirs('sequential/%s/%s'%(shearExp,spacing_multiplier))
    except OSError:
        if not os.path.isdir('sequential/%s/%s'%(shearExp,spacing_multiplier)):
            raise

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

    prob.driver.opt_settings['Summary file'] = 'sequential/%s/%s/summary_%s_%s.out'%(shearExp,spacing_multiplier,jobID,taskID)
    prob.driver.opt_settings['Print file'] = 'sequential/%s/%s/verify_%s_%s.out'%(shearExp,spacing_multiplier,jobID,taskID)

    prob.driver.add_objective('COE', scaler=10.)

    if spacing_multiplier > 0.5:
        prob.driver.add_desvar('turbineX', scaler=0.1)
        prob.driver.add_desvar('turbineY', scaler=0.1)

    for i in range(nGroups):
        prob.driver.add_desvar('d_param%s'%i, lower=3.87, upper=6.3, scaler=0.1)
        prob.driver.add_desvar('t_param%s'%i, lower=0.001, upper=None, scaler=1.)
    #     prob.driver.add_desvar('turbineH%s'%i, lower=10., scaler=0.01)
    #     prob.driver.add_desvar('rotorDiameter%s'%i, lower=10., upper=159.99, scaler=1.)
    #     prob.driver.add_desvar('ratedPower%s'%i, lower=500., upper=9999.99, scaler=0.00001)

    for i in range(nGroups):
        prob.driver.add_constraint('Tower%s_max_thrust.shell_buckling'%i, upper=np.ones(nFull))
        prob.driver.add_constraint('Tower%s_max_speed.shell_buckling'%i, upper=np.ones(nFull))
        prob.driver.add_constraint('freqConstraintGroup%s.freqConstraint'%i, lower=0.0)
    #     prob.driver.add_constraint('minHeight%s.minHeight'%i, lower=0.0)

    if spacing_multiplier > 0.6:
        prob.driver.add_constraint('spacing_con', lower=np.zeros(int(((nTurbs - 1.) * nTurbs / 2.))), scaler=1E-1)

        prob.driver.add_constraint('boundaryDistances', lower=(np.zeros(1*turbineX.size)), scaler=1E-1)
                                      #     active_tol=2. * rotor_diameter)
    prob.root.ln_solver.options['single_voi_relevance_reduction'] = True
    prob.root.ln_solver.options['mode'] = 'rev'

    prob.setup(check=True)

    setupTower(nFull, prob)
    simpleSetup(nTurbs, prob)
    prob['Uref'] = windSpeeds
    prob['windDirections'] = windDirections
    prob['windFrequencies'] = windFrequencies

    for i in range(nDirections):
        prob['yaw%s'%i] = yaw[i]
    prob['turbineX'] = turbineXstart
    prob['turbineY'] = turbineYstart
    prob['shearExp'] = shearExp

    for i in range(nGroups):
        prob['Tower%s_max_speed.Vel'%i] = 60.

    # provide values for hull constraint
    prob['boundary_radius'] = R
    prob['boundary_center'] = np.array([0.,0.])

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

    diameters = np.zeros((nGroups, 3))
    thicknesses = np.zeros((nGroups, 3))
    ratedPowers = np.zeros(nGroups)
    rotorDiameters = np.zeros(nGroups)
    for i in range(nGroups):
        diameters[i] = prob['d_param%s'%i]
        thicknesses[i] = prob['t_param%s'%i]
        ratedPowers[i] = prob['ratedPower%s'%i]
        rotorDiameters[i] = prob['rotorDiameter%s'%i]

    np.savetxt('sequential/%s/%s/XYZ_%s_%s.txt'%(shearExp,spacing_multiplier,jobID, taskID), np.c_[prob['turbineX'], prob['turbineY'], prob['turbineZ']], header="turbineX, turbineY, turbineZ")
    np.savetxt('sequential/%s/%s/XYstart_%s_%s.txt'%(shearExp,spacing_multiplier,jobID, taskID), np.c_[turbineXstart, turbineYstart], header="turbineXstart, turbineYstart")
    np.savetxt('sequential/%s/%s/diameter_%s_%s.txt'%(shearExp,spacing_multiplier,jobID, taskID), np.c_[diameters], header="diameters")
    np.savetxt('sequential/%s/%s/thickness_%s_%s.txt'%(shearExp,spacing_multiplier,jobID, taskID), np.c_[thicknesses], header="thicknesses")
    np.savetxt('sequential/%s/%s/COE_%s_%s.txt'%(shearExp,spacing_multiplier,jobID, taskID), np.c_[prob['COE'], prob['AEP']])
    np.savetxt('sequential/%s/%s/power_diam_%s_%s.txt'%(shearExp,spacing_multiplier,jobID, taskID), np.c_[ratedPowers,rotorDiameters])
