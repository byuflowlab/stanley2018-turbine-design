from openmdao.api import Component, Group, Problem, IndepVarComp
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.interpolate import SmoothBivariateSpline


class SimpleRotorSE(Component):
    def __init__(self, ratedQfunc, blade_massfunc, Vratedfunc, I1func, I2func, I3func, ratedTfunc, extremeTfunc):

        super(SimpleRotorSE, self).__init__()

        self.deriv_options['step_size'] = 1.E-4
        self.deriv_options['step_calc'] = 'relative'
        self.deriv_options['form'] = 'central'

        self.ratedTfunc = ratedTfunc
        self.ratedQfunc = ratedQfunc
        self.blade_massfunc = blade_massfunc
        self.Vratedfunc = Vratedfunc
        self.extremeTfunc = extremeTfunc
        self.I1func = I1func
        self.I2func = I2func
        self.I3func = I3func

        self.add_param('turbineRating', 0.0, units='kW', desc='turbine rating (kW)')
        self.add_param('rotorDiameter', 0.0, units='m', desc='rotor diameter (m)')

        self.add_output('ratedT', 0., desc='rated thrust')
        self.add_output('ratedQ', 0., desc='rated torque')
        self.add_output('blade_mass', 0., desc='blade mass')
        self.add_output('Vrated', 0., desc='rated wind speed')
        self.add_output('extremeT', 0., desc='extreme thrust')
        self.add_output('I', np.zeros(6), desc='I all blades')

    def solve_nonlinear(self, params, unknowns, resids):

        ratedTfunc = self.ratedTfunc
        ratedQfunc = self.ratedQfunc
        blade_massfunc = self.blade_massfunc
        Vratedfunc = self.Vratedfunc
        extremeTfunc = self.extremeTfunc
        I1func = self.I1func
        I2func = self.I2func
        I3func = self.I3func

        rating = params['turbineRating']
        diam = params['rotorDiameter']

        ratedPower_normal = 10000.0
        rotorDiameter_normal = 160.0
        ratedQ_normal = 6111001.99443
        blade_mass_normal = 26426.7546363
        Vrated_normal = 28.6866746136
        I1_normal = 70849909.3796
        I2_normal = 34069906.949
        I3_normal = 28001484.5761
        ratedT_normal = 1475530.07363
        extremeT_normal = 240326.697284

        # print 'IN COMP: RATING: ', rating
        # print 'IN COMP: DIAMETER: ', diam/rotorDiameter_normal
        # print 'IN COMP: blade_mass: ', blade_massfunc(rating/ratedPower_normal,diam/rotorDiameter_normal) * blade_mass_normal
        # print 'IN COMP: ratedT: ', ratedTfunc(rating/ratedPower_normal,diam/rotorDiameter_normal) * ratedT_normal
        # print 'IN COMP: Vrated: ', Vratedfunc(rating/ratedPower_normal,diam/rotorDiameter_normal) * Vrated_normal
        # print blade_massfunc(0.5,0.5)
        # print blade_massfunc(0.5,0.6)
        # print blade_massfunc(0.5,0.7)
        # print blade_massfunc(0.5,0.8)


        unknowns['ratedT'] = ratedTfunc(rating/ratedPower_normal,diam/rotorDiameter_normal) * ratedT_normal
        unknowns['ratedQ'] = ratedQfunc(rating/ratedPower_normal,diam/rotorDiameter_normal) * ratedQ_normal
        unknowns['blade_mass'] = blade_massfunc(rating/ratedPower_normal,diam/rotorDiameter_normal) * blade_mass_normal
        unknowns['Vrated'] = Vratedfunc(rating/ratedPower_normal,diam/rotorDiameter_normal) * Vrated_normal
        unknowns['extremeT'] = extremeTfunc(rating/ratedPower_normal,diam/rotorDiameter_normal) * extremeT_normal
        unknowns['I'] = np.array([I1func(rating/ratedPower_normal,diam/rotorDiameter_normal) * I1_normal,\
                        I2func(rating/ratedPower_normal,diam/rotorDiameter_normal) * I2_normal,\
                        I3func(rating/ratedPower_normal,diam/rotorDiameter_normal) * I3_normal,0.,0.,0.])


    def linearize(self, params, unknowns, resids):

        ratedTfunc = self.ratedTfunc
        ratedQfunc = self.ratedQfunc
        blade_massfunc = self.blade_massfunc
        Vratedfunc = self.Vratedfunc
        extremeTfunc = self.extremeTfunc
        I1func = self.I1func
        I2func = self.I2func
        I3func = self.I3func

        rating = params['turbineRating']
        diam = params['rotorDiameter']

        ratedPower_normal = 10000.0
        rotorDiameter_normal = 160.0
        ratedQ_normal = 6111001.99443
        blade_mass_normal = 26426.7546363
        Vrated_normal = 28.6866746136
        I1_normal = 70849909.3796
        I2_normal = 34069906.949
        I3_normal = 28001484.5761
        ratedT_normal = 1475530.07363
        extremeT_normal = 240326.697284

        dratedT_dratedPower = ratedTfunc(rating/ratedPower_normal,diam/rotorDiameter_normal,dx=1,dy=0) * ratedT_normal/ratedPower_normal
        dratedT_drotorDiameter = ratedTfunc(rating/ratedPower_normal,diam/rotorDiameter_normal,dx=0,dy=1) * ratedT_normal/rotorDiameter_normal

        dratedQ_dratedPower = ratedQfunc(rating/ratedPower_normal,diam/rotorDiameter_normal,dx=1,dy=0) * ratedQ_normal/ratedPower_normal
        dratedQ_drotorDiameter = ratedQfunc(rating/ratedPower_normal,diam/rotorDiameter_normal,dx=0,dy=1) * ratedQ_normal/rotorDiameter_normal

        dblade_mass_dratedPower = blade_massfunc(rating/ratedPower_normal,diam/rotorDiameter_normal,dx=1,dy=0) * blade_mass_normal/ratedPower_normal
        dblade_mass_drotorDiameter = blade_massfunc(rating/ratedPower_normal,diam/rotorDiameter_normal,dx=0,dy=1) * blade_mass_normal/rotorDiameter_normal

        dVrated_dratedPower = Vratedfunc(rating/ratedPower_normal,diam/rotorDiameter_normal,dx=1,dy=0) * Vrated_normal/ratedPower_normal
        dVrated_drotorDiameter = Vratedfunc(rating/ratedPower_normal,diam/rotorDiameter_normal,dx=0,dy=1) * Vrated_normal/rotorDiameter_normal

        dextremeT_dratedPower = extremeTfunc(rating/ratedPower_normal,diam/rotorDiameter_normal,dx=1,dy=0) * extremeT_normal/ratedPower_normal
        dextremeT_drotorDiameter = extremeTfunc(rating/ratedPower_normal,diam/rotorDiameter_normal,dx=0,dy=1) * extremeT_normal/rotorDiameter_normal

        dI1_dratedPower = I1func(rating/ratedPower_normal,diam/rotorDiameter_normal,dx=1,dy=0) * I1_normal/ratedPower_normal
        dI1_drotorDiameter = I1func(rating/ratedPower_normal,diam/rotorDiameter_normal,dx=0,dy=1) * I1_normal/rotorDiameter_normal

        dI2_dratedPower = I2func(rating/ratedPower_normal,diam/rotorDiameter_normal,dx=1,dy=0) * I2_normal/ratedPower_normal
        dI2_drotorDiameter = I2func(rating/ratedPower_normal,diam/rotorDiameter_normal,dx=0,dy=1) * I2_normal/rotorDiameter_normal

        dI3_dratedPower = I3func(rating/ratedPower_normal,diam/rotorDiameter_normal,dx=1,dy=0) * I3_normal/ratedPower_normal
        dI3_drotorDiameter = I3func(rating/ratedPower_normal,diam/rotorDiameter_normal,dx=0,dy=1) * I3_normal/rotorDiameter_normal

        J = {}
        J['ratedT', 'turbineRating'] = dratedT_dratedPower
        J['ratedT', 'rotorDiameter'] = dratedT_drotorDiameter

        J['ratedQ', 'turbineRating'] = dratedQ_dratedPower
        J['ratedQ', 'rotorDiameter'] = dratedQ_drotorDiameter

        J['blade_mass', 'turbineRating'] = dblade_mass_dratedPower
        J['blade_mass', 'rotorDiameter'] = dblade_mass_drotorDiameter

        J['Vrated', 'turbineRating'] = dVrated_dratedPower
        J['Vrated', 'rotorDiameter'] = dVrated_drotorDiameter

        J['extremeT', 'turbineRating'] = dextremeT_dratedPower
        J['extremeT', 'rotorDiameter'] = dextremeT_drotorDiameter

        J['I', 'turbineRating'] = np.zeros(6)
        J['I', 'turbineRating'][0] = dI1_dratedPower
        J['I', 'turbineRating'][1] = dI2_dratedPower
        J['I', 'turbineRating'][2] = dI3_dratedPower

        J['I', 'rotorDiameter'] = np.zeros(6)
        J['I', 'rotorDiameter'][0] = dI1_drotorDiameter
        J['I', 'rotorDiameter'][1] = dI2_drotorDiameter
        J['I', 'rotorDiameter'][2] = dI3_drotorDiameter


        return J



def create_rotor_functions():

    #loading data
    filename = '/home/flowlab/PJ/FLORISSE3D/doc/BEST_DATA.txt'
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

    w = np.ones(len(ratedPower))*2.
    order = 2

    interp_spline_ratedQ = SmoothBivariateSpline(ratedPower,rotorDiameter,ratedQ,w,kx=order,ky=order)
    interp_spline_blade_mass = SmoothBivariateSpline(ratedPower,rotorDiameter,blade_mass,w,kx=order,ky=order)
    interp_spline_Vrated = SmoothBivariateSpline(ratedPower,rotorDiameter,Vrated,w,kx=order,ky=order)
    interp_spline_I1 = SmoothBivariateSpline(ratedPower,rotorDiameter,I1,w,kx=order,ky=order)
    interp_spline_I2 = SmoothBivariateSpline(ratedPower,rotorDiameter,I2,w,kx=order,ky=order)
    interp_spline_I3 = SmoothBivariateSpline(ratedPower,rotorDiameter,I3,w,kx=order,ky=order)
    interp_spline_ratedT = SmoothBivariateSpline(ratedPower,rotorDiameter,ratedT,w,kx=order,ky=order)
    interp_spline_extremeT = SmoothBivariateSpline(ratedPower,rotorDiameter,extremeT,w,kx=order,ky=order)

    return interp_spline_ratedQ, interp_spline_blade_mass, interp_spline_Vrated, interp_spline_I1, interp_spline_I2, interp_spline_I3, interp_spline_ratedT, interp_spline_extremeT


if __name__=="__main__":
    interp_spline_ratedQ, interp_spline_blade_mass, interp_spline_Vrated, interp_spline_I1, interp_spline_I2, interp_spline_I3, interp_spline_ratedT, interp_spline_extremeT = create_rotor_functions()

    num = 100
    x = np.linspace(500.,10000.,num)/10000.
    y = np.linspace(50.,160.,num)/160.
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

    fig, ax = plt.subplots(nrows=2, ncols=4, subplot_kw={'projection': '3d'})
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')

    ax[0][0].plot_surface(X, Y, Z_blade_mass)
    ax[0][1].plot_surface(X, Y, Z_Vrated)
    ax[0][2].plot_surface(X, Y, Z_ratedQ)
    ax[0][3].plot_surface(X, Y, Z_ratedT)
    ax[1][0].plot_surface(X, Y, Z_I1)
    ax[1][1].plot_surface(X, Y, Z_I2)
    ax[1][2].plot_surface(X, Y, Z_I3)
    ax[1][3].plot_surface(X, Y, Z_extremeT)

    rp = 0.8
    diam = np.linspace(0.,1.,1000)
    mass = np.zeros(1000)
    for i in range(1000):
        mass[i] = interp_spline_blade_mass(rp,diam[i])
    plt.figure(2)
    plt.plot(diam,mass)

    print interp_spline_blade_mass(0.5,0.5)
    print interp_spline_blade_mass(0.5,0.6)
    print interp_spline_blade_mass(0.5,0.7)
    print interp_spline_blade_mass(0.5,0.8)

    # fig.tight_layout()
    # plt.xlabel('Turbine Rating')
    # plt.ylabel('Rotor Diameter')

    plt.show()
