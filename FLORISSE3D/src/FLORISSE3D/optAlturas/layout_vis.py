import numpy as np
import matplotlib.pyplot as plt
from math import sin,cos,radians

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
import matplotlib.gridspec as gridspec

nTurbs = 32
shearExp = 0.08
spacing_multiplier = 0.5
nGroups = 2
windDirections = np.array([0.])
windFrequencies = np.array([1.])
windSpeeds = np.array([10.])
nDirections = len(windDirections)
"""initial yaw values"""
yaw = np.zeros((nDirections, nTurbs))
nPoints = 3
nFull = 15
nVertices = 1
minSpacing = 2.
Rad = 659.6969000988257

resolution = 250
x = np.linspace(-750.,750.,resolution)
y = np.linspace(-750.,750.,resolution)
XX, YY = np.meshgrid(x, y)
xx = XX.flatten()
yy = YY.flatten()
n = len(xx)
plot_x = np.zeros((resolution,resolution))
plot_y = np.zeros((resolution,resolution))

for i in range(resolution):
    for j in range(resolution):
        plot_x[i][j] = xx[resolution*i+j]
        plot_y[i][j] = yy[resolution*i+j]

"""OpenMDAO"""
prob = Problem()
root = prob.root = Group()

# Design Variables
for i in range(nGroups):
    root.add('ratedPower%s'%i, IndepVarComp('ratedPower%s'%i, 0., units='kW'), promotes=['*'])
    root.add('d_param%s'%i, IndepVarComp('d_param%s'%i, np.zeros(3)), promotes=['*'])
    root.add('t_param%s'%i, IndepVarComp('t_param%s'%i, np.zeros(3)), promotes=['*'])
    root.add('turbineH%s'%i, IndepVarComp('turbineH%s'%i, 0.), promotes=['*'])
    root.add('rotorDiameter%s'%i, IndepVarComp('rotorDiameter%s'%i, 0.), promotes=['*'])

root.add('optCOE', optCOE(nGroups,nPoints,nFull,nTurbs,nDirections,nVertices,minSpacing,nSamples=n),promotes=['*'])
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

prob.setup(check=False)

setupTower(nFull, prob)
simpleSetup(nTurbs, prob)
prob['Uref'] = windSpeeds
prob['windDirections'] = windDirections
prob['windFrequencies'] = windFrequencies

for i in range(nDirections):
    prob['yaw%s'%i] = yaw[i]
prob['shearExp'] = shearExp

for i in range(nGroups):
    prob['Tower%s_max_speed.Vel'%i] = 60.

# provide values for hull constraint
prob['boundary_radius'] = Rad
prob['boundary_center'] = np.array([0.,0.])

prob['wsPositionX'] = xx
prob['wsPositionY'] = yy


salmon = (244./255.,159./255.,163./255.)
blue_purple = (149./255.,156./255.,254./255.)


nSteps = 50

file = 'z_diameter4.txt'
opt = open(file)
optimized = np.loadtxt(opt)
data = optimized
D = np.zeros((nSteps,nTurbs))
for i in range(nSteps):
    D[i][:] = data[i*nTurbs:(i+1)*nTurbs]

file = 'z_z4.txt'
opt = open(file)
optimized = np.loadtxt(opt)
data = optimized
Z = np.zeros((nSteps,nTurbs))
for i in range(nSteps):
    Z[i][:] = data[i*nTurbs:(i+1)*nTurbs]

file = 'z_x4.txt'
opt = open(file)
optimized = np.loadtxt(opt)
data = optimized
X = np.zeros((nSteps,nTurbs))
for i in range(nSteps):
    X[i][:] = data[i*nTurbs:(i+1)*nTurbs]

file = 'z_y4.txt'
opt = open(file)
optimized = np.loadtxt(opt)
data = optimized
Y = np.zeros((nSteps,nTurbs))
for i in range(nSteps):
    Y[i][:] = data[i*nTurbs:(i+1)*nTurbs]

file = 'z_rated4.txt'
opt = open(file)
optimized = np.loadtxt(opt)
data = optimized
R = np.zeros((nSteps,nTurbs))
for i in range(nSteps):
    R[i][:] = data[i*nTurbs:(i+1)*nTurbs]

file = 'z_tparam4.txt'
opt = open(file)
optimized = np.loadtxt(opt)
data = optimized
T_param = np.zeros((nSteps,18))
for i in range(nSteps):
    # print i
    T_param[i][:] = data[i*18:(i+1)*18]

T_param1 = np.zeros((nSteps,3))
T_param2 = np.zeros((nSteps,3))

for j in range(nSteps):
    t1 = T_param[j][0:3]
    t2 = T_param[j][3:6]
    t3 = T_param[j][6:9]
    t4 = T_param[j][9:12]
    if t1[0] < t2[0]:
        T_param1[j][:] = t1[:]
        T_param2[j][:] = t2[:]
    elif t1[0] < t3[0]:
        T_param1[j][:] = t1[:]
        T_param2[j][:] = t3[:]
    elif t1[0] < t4[0]:
        T_param1[j][:] = t1[:]
        T_param2[j][:] = t4[:]

    elif t2[0] < t1[0]:
        T_param1[j][:] = t2[:]
        T_param2[j][:] = t1[:]
    elif t2[0] < t3[0]:
        T_param1[j][:] = t2[:]
        T_param2[j][:] = t3[:]
    elif t2[0] < t4[0]:
        T_param1[j][:] = t2[:]
        T_param2[j][:] = t4[:]

    elif t3[0] < t1[0]:
        T_param1[j][:] = t3[:]
        T_param2[j][:] = t1[:]
    elif t3[0] < t2[0]:
        T_param1[j][:] = t3[:]
        T_param2[j][:] = t2[:]
    elif t3[0] < t4[0]:
        T_param1[j][:] = t3[:]
        T_param2[j][:] = t4[:]

    elif t4[0] < t1[0]:
        T_param1[j][:] = t4[:]
        T_param2[j][:] = t1[:]
    elif t4[0] < t2[0]:
        T_param1[j][:] = t4[:]
        T_param2[j][:] = t2[:]
    elif t4[0] < t3[0]:
        T_param1[j][:] = t4[:]
        T_param2[j][:] = t3[:]



file = 'z_dparam4.txt'
opt = open(file)
optimized = np.loadtxt(opt)
data = optimized
D_param = np.zeros((nSteps,18))
for i in range(nSteps):
    D_param[i][:] = data[i*18:(i+1)*18]
    # print data[i*18:(i+1)*18]

D_param1 = np.zeros((nSteps,3))
D_param2 = np.zeros((nSteps,3))

for j in range(nSteps):
    d1 = D_param[j][0:3]
    d2 = D_param[j][3:6]
    d3 = D_param[j][6:9]
    d4 = D_param[j][9:12]
    if d1[0] < d2[0]:
        D_param1[j][:] = d1[:]
        D_param2[j][:] = d2[:]
    elif d1[0] < d3[0]:
        D_param1[j][:] = d1[:]
        D_param2[j][:] = d3[:]
    elif d1[0] < d4[0]:
        D_param1[j][:] = d1[:]
        D_param2[j][:] = d4[:]

    elif d2[0] < d1[0]:
        D_param1[j][:] = d2[:]
        D_param2[j][:] = d1[:]
    elif d2[0] < d3[0]:
        D_param1[j][:] = d2[:]
        D_param2[j][:] = d3[:]
    elif d2[0] < d4[0]:
        D_param1[j][:] = d2[:]
        D_param2[j][:] = d4[:]

    elif d3[0] < d1[0]:
        D_param1[j][:] = d3[:]
        D_param2[j][:] = d1[:]
    elif d3[0] < d2[0]:
        D_param1[j][:] = d3[:]
        D_param2[j][:] = d2[:]
    elif d3[0] < d4[0]:
        D_param1[j][:] = d3[:]
        D_param2[j][:] = d4[:]

    elif d4[0] < d1[0]:
        D_param1[j][:] = d4[:]
        D_param2[j][:] = d1[:]
    elif d4[0] < d2[0]:
        D_param1[j][:] = d4[:]
        D_param2[j][:] = d2[:]
    elif d4[0] < d3[0]:
        D_param1[j][:] = d4[:]
        D_param2[j][:] = d3[:]
    D_param1[0] = np.array([6.,6.,6.])
    D_param1[1] = np.array([6.,6.,6.])
    D_param2[0] = np.array([6.,6.,6.])
    D_param2[1] = np.array([6.,6.,6.])


n = 15
x = np.zeros((nSteps+n*(nSteps-1),nTurbs))
y = np.zeros((nSteps+n*(nSteps-1),nTurbs))
d = np.zeros((nSteps+n*(nSteps-1),nTurbs))
z = np.zeros((nSteps+n*(nSteps-1),nTurbs))
r = np.zeros((nSteps+n*(nSteps-1),nTurbs))
d_param1 = np.zeros((nSteps+n*(nSteps-1),3))
d_param2 = np.zeros((nSteps+n*(nSteps-1),3))
t_param1 = np.zeros((nSteps+n*(nSteps-1),3))
t_param2 = np.zeros((nSteps+n*(nSteps-1),3))

x[0] = X[0]
y[0] = Y[0]
d[0] = D[0]
z[0] = Z[0]
r[0] = R[0]
d_param1[0] = D_param1[0]
d_param2[0] = D_param2[0]
t_param1[0] = T_param1[0]
t_param2[0] = T_param2[0]

for k in range(nSteps-1):
    x[(k+1)*(n+1)-1] = X[k+1]
    y[(k+1)*(n+1)-1] = Y[k+1]
    d[(k+1)*(n+1)-1] = D[k+1]
    z[(k+1)*(n+1)-1] = Z[k+1]
    r[(k+1)*(n+1)-1] = R[k+1]
    d_param1[(k+1)*(n+1)-1] = D_param1[k+1]
    d_param2[(k+1)*(n+1)-1] = D_param2[k+1]
    t_param1[(k+1)*(n+1)-1] = T_param1[k+1]
    t_param2[(k+1)*(n+1)-1] = T_param2[k+1]
    for i in range(n):
        for j in range(nTurbs):
            x[k*(n+1)+i][j] = (X[k+1][j]-X[k][j])/float(n+1.)*float(i+1)+X[k][j]
            y[k*(n+1)+i][j] = (Y[k+1][j]-Y[k][j])/float(n+1.)*float(i+1)+Y[k][j]
            d[k*(n+1)+i][j] = (D[k+1][j]-D[k][j])/float(n+1.)*float(i+1)+D[k][j]
            z[k*(n+1)+i][j] = (Z[k+1][j]-Z[k][j])/float(n+1.)*float(i+1)+Z[k][j]
            r[k*(n+1)+i][j] = (R[k+1][j]-R[k][j])/float(n+1.)*float(i+1)+R[k][j]
        for j in range(3):
            d_param1[k*(n+1)+i][j] = (D_param1[k+1][j]-D_param1[k][j])/float(n+1.)*float(i+1)+D_param1[k][j]
            d_param2[k*(n+1)+i][j] = (D_param2[k+1][j]-D_param2[k][j])/float(n+1.)*float(i+1)+D_param2[k][j]
            t_param1[k*(n+1)+i][j] = (T_param1[k+1][j]-T_param1[k][j])/float(n+1.)*float(i+1)+T_param1[k][j]
            t_param2[k*(n+1)+i][j] = (T_param2[k+1][j]-T_param2[k][j])/float(n+1.)*float(i+1)+T_param2[k][j]










f = plt.figure(figsize=(14.5,7.5))
# gridspec.GridSpec(14,12)
# f, ((a1, a2, a3, a4), (a5,a6,a7,a8)) = plt.subplots(2,4,figsize=(14.5,12.), gridspec_kw = {'width_ratios':[4,5,2.5,2.5], 'height_ratios':[1,2]})

a1 = plt.subplot2grid((16, 28), (0, 0),rowspan=8,colspan=8)
a2 = plt.subplot2grid((16, 28), (0, 8),rowspan=8,colspan=10)
a3 = plt.subplot2grid((16, 28), (0, 18),rowspan=8,colspan=5)
a4 = plt.subplot2grid((16, 28), (0, 23),rowspan=8,colspan=5)
a5 = plt.subplot2grid((16, 28), (8, 0),rowspan=8,colspan=8)
a6 = plt.subplot2grid((16, 28), (8, 8),rowspan=8,colspan=8)

iteration = np.array([])
convergence = np.array([])
COE_arr = np.array([])
AEP_arr = np.array([])

deg = 0.
counter = 1
# for k in range(nSteps+n*(nSteps-1)):
for k in range(501):
    print k
    a1.cla()
    a2.cla()
    a3.cla()
    a4.cla()
    a5.cla()
    a6.cla()

    #layout, top left
    a1.axis('off')
    circ = plt.Circle((0.,0.), Rad,facecolor="None",edgecolor="black",alpha=0.8)
    a1.add_patch(circ)
    for i in range(nTurbs):
        if i%2 == 0:
            color=blue_purple
        else:
            color=salmon
        circ = plt.Circle((x[k][i],y[k][i]), d[k][i]/2.,facecolor=color,edgecolor=color)
        a1.add_patch(circ)
    a1.set_xlim(-800.,800.)
    a1.set_ylim(-750.,750.)
    # plt.axis('equal')




    #turbines, top middle
    a2.axis('off')

    R1 = d[k][0]/2.
    R2 = d[k][1]/2.
    H1 = z[k][0]
    H2 = z[k][1]

    x1 = 70.
    x2 = 200.

    a2.set_xlim(0.,275.)
    a2.set_ylim(0.,200.)
    # plt.axis('equal')

    bladeX = np.array([3.,7.,10.,15.,20.,25.,30.,35.,30.,25.,20.,15.,10.,5.,3.,3.])
    bladeY = np.array([0.,0.,0.8,1.5,1.7,1.9,2.1,2.3,2.4,2.4,2.4,2.4,2.4,2.4,2.4,0.])-1.5

    c1 = R1/35.
    d1 = d_param1[k]
    px1 = np.array([x1-d1[0]/2,x1-d1[1]/2,x1-d1[2]/2,x1+d1[2]/2,x1+d1[1]/2,x1+d1[0]/2,x1-d1[0]/2])
    py1 = np.array([0,H1/2,H1-3.*c1,H1-3.*c1,H1/2,0,0])
    a2.plot(px1,py1,color=blue_purple, linewidth=2)
    #add blades
    hub1 = plt.Circle((x1,H1), 3*c1, color=blue_purple, fill=False, linewidth=2)
    a2.add_artist(hub1)

    angle1 = -65.+deg
    blade1X = bladeX*cos(radians(angle1))-bladeY*sin(radians(angle1))
    blade1Y = bladeX*sin(radians(angle1))+bladeY*cos(radians(angle1))

    blade2X = bladeX*cos(radians(angle1+120.))-bladeY*sin(radians(angle1+120.))
    blade2Y = bladeX*sin(radians(angle1+120.))+bladeY*cos(radians(angle1+120.))

    blade3X = bladeX*cos(radians(angle1+240.))-bladeY*sin(radians(angle1+240.))
    blade3Y = bladeX*sin(radians(angle1+240.))+bladeY*cos(radians(angle1+240.))

    a2.plot(blade1X*c1+x1, blade1Y*c1+H1, linewidth=2, color=blue_purple)
    a2.plot(blade2X*c1+x1, blade2Y*c1+H1, linewidth=2, color=blue_purple)
    a2.plot(blade3X*c1+x1, blade3Y*c1+H1, linewidth=2, color=blue_purple)


    c2 = R2/35.
    d2 = d_param2[k]
    px2 = np.array([x2-d2[0]/2,x2-d2[1]/2,x2-d2[2]/2,x2+d2[2]/2,x2+d2[1]/2,x2+d2[0]/2,x2-d2[0]/2])
    py2 = np.array([0,H2/2,H2-3.*c2,H2-3.*c2,H2/2,0,0])
    a2.plot(px2,py2,color=salmon, linewidth=2)
    #add blades
    hub2 = plt.Circle((x2,H2), 3*c2, color=salmon, fill=False, linewidth=2)
    a2.add_artist(hub2)

    angle2 = -20.+deg*.6
    blade1X = bladeX*cos(radians(angle2))-bladeY*sin(radians(angle2))
    blade1Y = bladeX*sin(radians(angle2))+bladeY*cos(radians(angle2))

    blade2X = bladeX*cos(radians(angle2+120.))-bladeY*sin(radians(angle2+120.))
    blade2Y = bladeX*sin(radians(angle2+120.))+bladeY*cos(radians(angle2+120.))

    blade3X = bladeX*cos(radians(angle2+240.))-bladeY*sin(radians(angle2+240.))
    blade3Y = bladeX*sin(radians(angle2+240.))+bladeY*cos(radians(angle2+240.))

    a2.plot(blade1X*c2+x2, blade1Y*c2+H2, linewidth=2, color=salmon)
    a2.plot(blade2X*c2+x2, blade2Y*c2+H2, linewidth=2, color=salmon)
    a2.plot(blade3X*c2+x2, blade3Y*c2+H2, linewidth=2, color=salmon)
    a2.axis('off')

    plt.draw()

    deg += 5.

    prob['turbineX'] = x[k]
    prob['turbineY'] = y[k]
    prob['ratedPower0'] = r[k][0]
    prob['ratedPower1'] = r[k][1]
    prob['rotorDiameter0'] = d[k][0]
    prob['rotorDiameter1'] = d[k][1]
    prob['turbineH0'] = z[k][0]
    prob['turbineH1'] = z[k][1]
    prob['d_param0'] = d_param1[k]
    prob['d_param1'] = d_param2[k]
    prob['t_param0'] = t_param1[k]
    prob['t_param1'] = t_param2[k]


    prob['wsPositionZ'] = np.ones_like(xx)*H2
    prob.run()
    vmin = 0.
    vmax = 10.
    wind_speed_array1 = prob['wsArray0']
    plot_ws1 = np.zeros((resolution,resolution))
    for i in range(resolution):
        for j in range(resolution):
            plot_ws1[i][j] = wind_speed_array1[resolution*i+j]
    im = a6.pcolormesh(plot_x, plot_y, plot_ws1, cmap='Greys_r', vmin=vmin, vmax=vmax)

    # print prob['wsPositionZ']


    prob['wsPositionZ'] = np.ones_like(xx)*H1
    prob.run()
    wind_speed_array2 = prob['wsArray0']
    plot_ws2 = np.zeros((resolution,resolution))
    for i in range(resolution):
        for j in range(resolution):
            plot_ws1[i][j] = wind_speed_array1[resolution*i+j]
            plot_ws2[i][j] = wind_speed_array2[resolution*i+j]
    im = a5.pcolormesh(plot_x, plot_y, plot_ws2, cmap='Greys_r', vmin=vmin, vmax=vmax,shading='gouraud')


    # print prob['wsPositionZ']

    COE = prob['COE']
    AEP = prob['AEP']/1E6

    COE_arr = np.append(COE_arr,COE)
    AEP_arr = np.append(AEP_arr,AEP)
    iteration = np.append(iteration,float(k)/16.)

    # print 'COE: ', repr(COE_arr)
    # print 'AEP: ', repr(AEP_arr)
    # print 'iteration: ', repr(iteration)

    a3.plot(iteration,COE_arr,'k',alpha=0.5)
    a3.set_xlim(0.,32.)
    a3.set_ylim(0.,100.)
    a3.set_xticks((10,20,30))
    a3.set_ylabel('cost of energy ($/MWh)')
    a3.set_xlabel('iteration')

    a4.plot(iteration,AEP_arr,'k',alpha=0.75)
    a4.set_xlim(0.,32.)
    a4.set_ylim(0.,600.)
    a4.set_xticks((10,20,30))
    a4.set_ylabel('annual energy production (GWh)')
    a4.set_xlabel('iteration')

    # plt.tight_layout()

    circ = plt.Circle((0.,0.), Rad,facecolor="None",edgecolor="black",alpha=0.8)
    a5.add_patch(circ)

    circ = plt.Circle((0.,0.), Rad,facecolor="None",edgecolor="black",alpha=0.8)
    a6.add_patch(circ)


    for i in range(nTurbs):
        if i%2 == 0:
            color=blue_purple
        else:
            color = salmon
        R = d[k][i]/2.
        tx = x[k][i]
        ty = y[k][i]
        theta = (np.deg2rad(90.))
        num = 1.
        num1 = 0.
        xt = tx-R*np.sin(theta)
        xb = tx+R*np.sin(theta)
        yt = ty+R*np.cos(theta)
        yb = ty-R*np.cos(theta)
        a5.plot([xt,xb],[yt*num-num1,yb*num-num1],color=color,linewidth=2)
        a6.plot([xt,xb],[yt*num-num1,yb*num-num1],color=color,linewidth=2)


    a5.set_ylim(-Rad-50.,Rad+50.)
    a5.set_xlim(-Rad-100.,Rad+100.)
    a6.set_ylim(-Rad-50.,Rad+50.)
    a6.set_xlim(-Rad-100.,Rad+100.)

    a5.axis('off')
    a6.axis('off')

    plt.tight_layout()
    plt.pause(0.0001)



    # if k < 10:
    #     plt.savefig('vid/lay_00%s.png'%k,transparent=True)
    # elif k >= 10 and k < 100:
    #     plt.savefig('vid/lay_0%s.png'%k,transparent=True)
    # else:
    #     plt.savefig('vid/lay_%s.png'%k,transparent=True)
