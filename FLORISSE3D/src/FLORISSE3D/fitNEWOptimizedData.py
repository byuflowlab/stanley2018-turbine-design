import numpy as np
from scipy.interpolate import LinearNDInterpolator,SmoothBivariateSpline
from scipy import ndimage
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
import math
np.set_printoptions(threshold='nan')

#loading data
filename = 'src/florisse3D/optRotor/BEST_DATA.txt'
opened = open(filename)
data = np.loadtxt(opened)
"ratedPower, rotorDiameter, ratedQ, blade_mass, Vrated, I1, I2, I3, ratedT, extremeT"
ratedPower = data[:,0]
rotorDiameter = data[:,1]
ratedQ = data[:,2]
blade_mass = data[:,3]
Vrated = data[:,4]
I1 = data[:,5]
I2 = data[:,6]
I3 = data[:,7]
ratedT = data[:,8]
extremeT = data[:,9]

print max(ratedPower)
print max(rotorDiameter)
print max(ratedQ)
print max(blade_mass)
print max(Vrated)
print max(I1)
print max(I2)
print max(I3)
print max(ratedT)
print max(extremeT)

newRatedQ = np.zeros(30)
newRatedPower = np.zeros(30)
newRotorDiameter = np.zeros(30)

for i in range(30):
    newRatedQ[i] = ratedQ[i*10]/max(ratedQ)
    newRatedPower[i] = ratedPower[i*10]/max(ratedPower)
    newRotorDiameter[i] = rotorDiameter[i*10]/max(rotorDiameter)

results_rated_power = np.zeros((20,20))
results_rotor_diameter = np.zeros((20,20))
results_ratedQ = np.zeros((20,20))
results_blade_mass = np.zeros((20,20))
results_Vrated = np.zeros((20,20))
results_I1 = np.zeros((20,20))
results_I2 = np.zeros((20,20))
results_I3 = np.zeros((20,20))
results_ratedT = np.zeros((20,20))
results_extremeT = np.zeros((20,20))

n = 0
# print ratedPower
for i in range(20):
    for j in range(20):
        for k in range(len(ratedPower)):
            results_rated_power[i][j] = float(i*500.+500.)
            results_rotor_diameter[i][j] = float(j*6.+46.)
            if ratedPower[k]==i*500.+500. and rotorDiameter[k]==float(j*6.+46.):
                results_ratedQ[i][j] = ratedQ[n]
                results_blade_mass[i][j] = blade_mass[n]
                results_Vrated[i][j] = Vrated[n]
                results_I1[i][j] = I1[n]
                results_I2[i][j] = I2[n]
                results_I3[i][j] = I3[n]
                results_ratedT[i][j] = ratedT[n]
                results_extremeT[i][j] = extremeT[n]
                n += 1


# for i in range(20):
#     for j in range(20):
#         if results_blade_mass[i][j] != 0:
#             plt.plot(i*500.+500.,j*6.+46.,'ob')
#         else:
#             plt.plot(i*500.+500.,j*6.+46.,'or')
# plt.show()

ratedPower = ratedPower/max(ratedPower)
rotorDiameter = rotorDiameter/max(rotorDiameter)
ratedQ = ratedQ/max(ratedQ)
blade_mass = blade_mass/max(blade_mass)
Vrated = Vrated/max(Vrated)
I1 = I1/max(I1)
I2 = I2/max(I2)
I3 = I3/max(I3)
ratedT = ratedT/max(ratedT)
extremeT = extremeT/max(extremeT)


cartcoord = list(zip(ratedPower,rotorDiameter))
w = np.ones(len(ratedPower))*2.

order = 2

interp_spline_ratedQ = SmoothBivariateSpline(ratedPower,rotorDiameter,ratedQ,w,kx=order,ky=order)
interp_spline_blade_mass = SmoothBivariateSpline(ratedPower,rotorDiameter,blade_mass,w,kx=order,ky=order)#, s=1000.)
interp_spline_Vrated = SmoothBivariateSpline(ratedPower,rotorDiameter,Vrated,w,kx=order,ky=order)
interp_spline_I1 = SmoothBivariateSpline(ratedPower,rotorDiameter,I1,w,kx=order,ky=order)
interp_spline_I2 = SmoothBivariateSpline(ratedPower,rotorDiameter,I2,w,kx=order,ky=order)
interp_spline_I3 = SmoothBivariateSpline(ratedPower,rotorDiameter,I3,w,kx=order,ky=order)
interp_spline_ratedT = SmoothBivariateSpline(ratedPower,rotorDiameter,ratedT,w,kx=order,ky=order)
interp_spline_extremeT = SmoothBivariateSpline(ratedPower,rotorDiameter,extremeT,w,kx=order,ky=order)

num = 100
x = np.linspace(500.,10000.,num)/10000.
y = np.linspace(40.,160.,num)/160.
# x = np.linspace(0.,1.,num)
# y = np.linspace(0.,1.,num)
X,Y = np.meshgrid(x,y)
#
Z_ratedQ = np.zeros((num,num))
Z_blade_mass = np.zeros((num,num))
Z_Vrated = np.zeros((num,num))
Z_I1 = np.zeros((num,num))
Z_I2 = np.zeros((num,num))
Z_I3 = np.zeros((num,num))
Z_ratedT = np.zeros((num,num))
Z_extremeT = np.zeros((num,num))

for i in range(num):
    for j in range(num):
        Z_ratedQ[j][i] = interp_spline_ratedQ(x[i],y[j])
        Z_blade_mass[j][i] = interp_spline_blade_mass(x[i],y[j])
        Z_Vrated[j][i] = interp_spline_Vrated(x[i],y[j])
        Z_I1[j][i] = interp_spline_I1(x[i],y[j])
        Z_I2[j][i] = interp_spline_I2(x[i],y[j])
        Z_I3[j][i] = interp_spline_I3(x[i],y[j])
        Z_ratedT[j][i] = interp_spline_ratedT(x[i],y[j])
        Z_extremeT[j][i] = interp_spline_extremeT(x[i],y[j])

# fig, ax = plt.subplots(nrows=2, ncols=4, subplot_kw={'projection': '3d'})
fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# ax[0][0].plot_surface(X, Y, Z_blade_mass)
# ax[0][0].plot(ratedPower, rotorDiameter, blade_mass, 'or')
#
# ax[0][1].plot_surface(X, Y, Z_Vrated)
# ax[0][1].plot(ratedPower, rotorDiameter, Vrated, 'or')
#
# ax[0][2].plot_surface(X, Y, Z_ratedQ)
# ax[0][2].plot(ratedPower, rotorDiameter, ratedQ, 'or')
#
# ax[0][3].plot_surface(X, Y, Z_ratedT)
# ax[0][3].plot(ratedPower, rotorDiameter, ratedT, 'or')


# ax[1][0].plot_surface(X, Y, Z_I1)
# ax[1][0].plot(ratedPower, rotorDiameter, I1, 'or')

ax.plot_surface(X, Y, Z_ratedQ,cmap=plt.cm.gray)
# ax.plot(ratedPower, rotorDiameter, ratedQ, 'ow')
ax.plot(newRatedPower, newRotorDiameter, newRatedQ, 'ow')

# ax[1][1].plot_surface(X, Y, Z_I2)
# ax[1][1].plot(ratedPower, rotorDiameter, I2, 'or')
#
# ax[1][2].plot_surface(X, Y, Z_I3)
# ax[1][2].plot(ratedPower, rotorDiameter, I3, 'or')
#
# ax[1][3].plot_surface(X, Y, Z_extremeT)
# ax[1][3].plot(ratedPower, rotorDiameter, extremeT, 'or')
#
# print 'max, max: ', interp_spline_blade_mass(10000.,160.), results_blade_mass[19][19]
#
# #
# # ax[1].plot_wireframe(X, Y, Z_smooth_blade_mass, color='k')
# #
# # # for axes in ax:
# # #     # axes.set_zlim(-0.2,1)
# # #     axes.set_axis_off()
# #
plt.axis('off')
fig.tight_layout()
plt.savefig('fit.pdf', transparent=True)


# plt.xlabel('Turbine Rating')
# plt.ylabel('Rotor Diameter')
# plt.title('Blade Mass')

#
# #
# # fig, ax = plt.subplots(nrows=1, ncols=2, subplot_kw={'projection': '3d'})
# # ax[0].plot_wireframe(X, Y, Z_Vrated, color='k')
# #
# # ax[1].plot_wireframe(X, Y, Z_smooth_Vrated, color='k')
# #
# # # for axes in ax:
# # #     # axes.set_zlim(-0.2,1)
# # #     axes.set_axis_off()
# #
# # fig.tight_layout()
# # plt.xlabel('Turbine Rating')
# # plt.ylabel('Rotor Diameter')
print len(w)
plt.show()

print interp_spline_blade_mass.get_coeffs()
print interp_spline_blade_mass.ev(0.5,0.5,dx=0,dy=1)
print interp_spline_blade_mass.ev(0.5,0.5,dx=1,dy=0)
print interp_spline_blade_mass.ev(0.5,0.5,dx=1,dy=1)
