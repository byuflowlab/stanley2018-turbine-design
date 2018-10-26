from openmdao.api import Component, Group, Problem, IndepVarComp
import numpy as np
from scipy.optimize import fsolve
from math import sin, cos, sinh, cosh, sqrt, pi
from scipy.optimize import brentq, root
import _shellbuckling
import _axialShear
from akima import Akima
from commonse.WindWaveDrag import TowerWindDrag
from commonse.utilities import interp_with_deriv
from math import pi


class topMassAdder(Component):
    def __init__(self):

        super(topMassAdder, self).__init__()

        self.add_param('rotor_mass', 0.0, units='m', desc='rotor mass')
        self.add_param('nacelle_mass', 0.0, units='m', desc='nacelle mass')

        self.add_output('m', np.array([0.]), desc='top mass')

    def solve_nonlinear(self, params, unknowns, resids):
        m = params['rotor_mass']+params['nacelle_mass']
        unknowns['m'] = np.array([m])

    def linearize(self, params, unknowns, resids):

        J = {}
        J['m','rotor_mass'] = np.array([1.])
        J['m','nacelle_mass'] = np.array([1.])

        return J


class Tower(Group):
    """group to calculate mass, frequency, stresses, and buckling for a tower"""

    def __init__(self, nPoints, nFull):

        super(Tower, self).__init__()

        # set finite difference options (fd used for testing only)
        # self.deriv_options['form'] = 'central'
        # self.deriv_options['step_size'] = 1.E-8
        # self.deriv_options['step_calc'] = 'relative'

        self.add('TowerDiscretization', TowerDiscretization(nPoints, nFull), promotes=['*'])
        self.add('calcMass', calcMass(nFull), promotes=['*'])
        self.add('dynamicQgroup', dynamicQgroup(nFull), promotes=['*'])
        self.add('shellBucklingGroup', shellBucklingGroup(nPoints, nFull), promotes=['*'])
        self.add('averageI', averageI(nFull), promotes=['*'])
        self.add('m_L', m_L(), promotes=['*'])
        self.add('freq', freq(nFull), promotes=['*'])
        self.add('windLoads', TowerWindDrag(nFull))
        self.add('topMassAdder', topMassAdder(), promotes=['*'])

        self.connect('towerSpeeds', 'windLoads.U')
        self.connect('z_full', 'windLoads.z')
        self.connect('d_full', 'windLoads.d')
        self.connect('rhoAir', 'windLoads.rho')
        self.connect('windLoads.windLoads:Px', 'qx')
        self.connect('windLoads.windLoads:Py', 'qy')
        self.connect('m','Mt')


class calcMass(Component):
    """
    Calculate the mass of the cylinder tower
    """

    def __init__(self, nFull):

        super(calcMass, self).__init__()

        self.nFull = nFull
        self.add_param('z_full', np.zeros(nFull), units='m', desc='Height of the Tower')
        self.add_param('d_full', np.zeros(nFull), units='m', desc='parameterized diameter')
        self.add_param('t_full', np.zeros(nFull), units='m', desc='parameterized thickness')
        self.add_param('rho', np.zeros(nFull), desc='density of the material')

        self.add_output('mass', 0.0, desc='tower mass')


    def solve_nonlinear(self, params, unknowns, resids):
        z_full = params['z_full']
        d = params['d_full']
        t = params['t_full']
        rho = params['rho']
        nFull = self.nFull

        sectionMass = np.zeros(nFull-1)
        for i in range(nFull-1):
            outer = 1./3.*pi*((d[i]/2.)**2+(d[i]/2.)*(d[i+1]/2.)+ \
                                (d[i+1]/2.)**2)*(z_full[i+1]-z_full[i])
            inner = 1./3.*pi*(((d[i]/2.)-t[i])**2+((d[i]/2.)-t[i])* \
                                ((d[i+1]/2.)-t[i+1])+((d[i+1]/2.)-t[i+1])**2)*(z_full[i+1]-z_full[i])
            sectionMass[i] = (outer-inner)*rho[i]

        # bottom_outer = 1./3.*3.141592653589793*(r[0]**2+r[0]*r[1]+r[1]**2)*H/2.
        # bottom_inner = 1./3.*3.141592653589793*((r[0]-t[0])**2+(r[0]-t[0])*(r[1]-t[1])+(r[1]-t[1])**2)*H/2.
        # top_outer = 1./3.*3.141592653589793*(r[1]**2+r[1]*r[2]+r[2]**2)*H/2.
        # top_inner = 1./3.*3.141592653589793*((r[1]-t[1])**2+(r[1]-t[1])*(r[2]-t[2])+(r[2]-t[2])**2)*H/2.

        # unknowns['mass'] = (bottom_outer + top_outer - bottom_inner - top_inner)*rho
        unknowns['mass'] = np.sum(sectionMass)

    def linearize(self, params, unknowns, resids):
        z_full = params['z_full']
        d = params['d_full']
        t = params['t_full']
        rho = params['rho'][0] #assuming density remains constant
        nFull = self.nFull
        mass = unknowns['mass']
        #
        # J = {}
        # J['mass', 'turbineH'] = mass/H
        #
        dmass_dd = np.zeros((1,nFull))
        dmass_dt = np.zeros((1,nFull))
        dmass_dz = np.zeros((1,nFull))


        dmass_dd[0][0] = (z_full[1]-z_full[0])*pi/12.*((2.*d[0]+d[1])-(2.*(d[0]-2.*t[0])+(d[1]-2.*t[1])))*rho
        dmass_dt[0][0] = -1.*(z_full[1]-z_full[0])*pi/12.*(-4.*(d[0]-2.*t[0])-2.*(d[1]-2.*t[1]))*rho
        dmass_dz[0][0] = (-rho/12.*pi*(d[0]**2+d[0]*d[1]+d[1]**2))-(-rho/12.*pi*((d[0]-2.*t[0])**2+ \
                                (d[0]-2.*t[0])*(d[1]-2.*t[1])+(d[1]-2.*t[1])**2))

        """For the dmass_dd and dmass_dt, the SECOND TO LAST value is failing the test. This
        loop is working for every value except that one"""
        for i in range(1, nFull-1):
            dmass_dd[0][i] = (z_full[i+1]-z_full[i])*pi/12.*((d[i-1]+2.*d[i])+(2.*d[i]+d[i+1])-((d[i-1]-2.*t[i-1])+2.* \
                                (d[i]-2.*t[i]))-(2.*(d[i]-2.*t[i])+(d[i+1]-2.*t[i+1])))*rho
            dmass_dt[0][i] = -1.*(z_full[i+1]-z_full[i])*pi/12.*(-4.*(d[i]-2.*t[i])-2.*(d[i-1]-2.*t[i-1])+-4.*(d[i]- \
                                2.*t[i])-2.*(d[i+1]-2*t[i+1]))*rho
            z1 = (rho/12.*pi*(d[i-1]**2+d[i-1]*d[i]+d[i]**2))-(rho/12.*pi*((d[i-1]-2.*t[i-1])**2+ \
                                    (d[i-1]-2.*t[i-1])*(d[i]-2.*t[i])+(d[i]-2.*t[i])**2))
            z2 = (rho/12.*pi*(d[i]**2+d[i]*d[i+1]+d[i+1]**2))-(rho/12.*pi*((d[i]-2.*t[i])**2+ \
                                    (d[i]-2.*t[i])*(d[i+1]-2.*t[i+1])+(d[i+1]-2.*t[i+1])**2))
            dmass_dz[0][i] = z1-z2
            # dmass_dz[0][i+1] =

        # dmass_dd[0][nFull-2] = (z_full[nFull-1]-z_full[nFull-2])*3.141592653589793/12.*((d[nFull-3]+2.*d[nFull-2])+(2.*d[nFull-2]+d[nFull-1])-((d[nFull-3]-2.*t[nFull-3])+2.* \
        #                     (d[nFull-2]-2.*t[nFull-2]))-(2.*(d[nFull-2]-2.*t[nFull-2])+(d[nFull-1]-2.*t[nFull-1])))*rho
        # dmass_dt[0][nFull-2] = -1.*(z_full[nFull-1]-z_full[nFull-2])*3.141592653589793/12.*(-4.*(d[nFull-2]-2.*t[nFull-2])-2.*(d[nFull-3]-2.*t[nFull-3])+-4.*(d[nFull-2]- \
        #                     2.*t[nFull-2])-2.*(d[nFull-1]-2*t[nFull-1]))*rho

        dmass_dd[0][-1] = (z_full[-1]-z_full[nFull-2])*pi/12.*((2.*d[-1]+d[nFull-2])-(2.*(d[-1]-2.*t[-1])+(d[nFull-2]-2.*t[nFull-2])))*rho
        dmass_dt[0][-1] = -1.*(z_full[-1]-z_full[nFull-2])*pi/12.*(-4.*(d[-1]-2.*t[-1])-2.*(d[nFull-2]-2.*t[nFull-2]))*rho
        dmass_dz[0][-1] = (rho/12.*pi*(d[nFull-2]**2+d[nFull-2]*d[-1]+d[-1]**2))-(rho/12.*pi*((d[nFull-2]-2.*t[nFull-2])**2+ \
                                (d[nFull-2]-2.*t[nFull-2])*(d[-1]-2.*t[-1])+(d[-1]-2.*t[-1])**2))

        J = {}
        J['mass', 'd_full'] = dmass_dd
        J['mass', 't_full'] = dmass_dt
        J['mass', 'z_full'] = dmass_dz

        return J


class TowerDiscretization(Component):
    """discretize geometry into finite element nodes"""

    #inputs

    def __init__(self, nPoints, nFull):

        super(TowerDiscretization, self).__init__()

        self.nFull = nFull
        self.nPoints = nPoints

         # variables
        self.add_param('z_param', np.zeros(nPoints), units='m', desc='parameterized locations along tower, linear lofting between')
        self.add_param('d_param', np.zeros(nPoints), units='m', desc='tower diameter at corresponding locations')
        self.add_param('t_param', np.zeros(nPoints), units='m', desc='shell thickness at corresponding locations')
        self.add_param('z_full', np.zeros(nFull), units='m', desc='locations along tower')

        #out
        self.add_output('d_full', np.zeros(nFull), units='m', desc='tower diameter at corresponding locations')
        self.add_output('t_full', np.zeros(nFull), units='m', desc='shell thickness at corresponding locations')


    def solve_nonlinear(self, params, unknowns, resids):
        z_param = params['z_param']
        d_param = params['d_param']
        t_param = params['t_param']
        z_full = params['z_full']

        d_full, self.ddfull_dzfull, self.ddfull_dzparam, self.ddfull_ddparam = interp_with_deriv(z_full, z_param, d_param)
        t_full, self.dtfull_dzfull, self.dtfull_dzparam, self.dtfull_dtparam = interp_with_deriv(z_full, z_param, t_param)

        unknowns['d_full'] = d_full
        unknowns['t_full'] = t_full

    def linearize(self, params, unknowns, resids):

        J = {}

        J['d_full', 'd_param'] = self.ddfull_ddparam
        J['d_full', 't_param'] = np.zeros((self.nFull, self.nPoints))
        J['d_full', 'z_param'] = self.ddfull_dzparam
        J['d_full', 'z_full'] = self.ddfull_dzfull

        J['t_full', 'd_param'] = np.zeros((self.nFull, self.nPoints))
        J['t_full', 't_param'] = self.dtfull_dtparam
        J['t_full', 'z_param'] = self.dtfull_dzparam
        J['t_full', 'z_full'] = self.dtfull_dzfull

        return J


class dynamicQgroup(Group):
    """group to calculate dynamic pressure at each point of the tower"""

    def __init__(self, nFull):

        super(dynamicQgroup, self).__init__()

        self.add('speed', powWindTower(nFull), promotes=['*'])
        self.add('q_dyn', dynamic_q(nFull), promotes=['*'])


class powWindTower(Component):
    """calculate the wind speed up the height of the tower"""
    def __init__(self, nFull):

        super(powWindTower, self).__init__()

        self.nFull = nFull
        self.add_param('Vel', 0.0, desc='reference wind speed')
        self.add_param('zref', 90.0, units='m', desc='reference height')
        self.add_param('z0', 0.0, units='m', desc='reference zero height')
        self.add_param('shearExp', 0.15, desc='wind shear exponent')
        self.add_param('z_full', np.zeros(nFull), units='m', desc='heights of interest')

        self.add_output('towerSpeeds', np.zeros(nFull), desc='speeds at different heights')

    def solve_nonlinear(self, params, unknowns, resids):

        Vel = params['Vel']
        zref = params['zref']
        z0 = params['z0']
        shearExp = params['shearExp']
        z = params['z_full']

        speeds = Vel*((z-z0)/(zref-z0))**shearExp

        unknowns['towerSpeeds'] = speeds

    def linearize(self, params, unknowns, resids):

        nFull = self.nFull
        Vel = params['Vel']
        zref = params['zref']
        z0 = params['z0']
        shearExp = params['shearExp']
        z = params['z_full']

        speeds = Vel*((z-z0)/(zref-z0))**shearExp

        J = {}
        dspeeds_dV = ((z-z0)/(zref-z0))**shearExp
        dspeeds_dz = np.zeros((nFull, nFull))
        for i in range(nFull):
            if z[i] == 0:
                z[i] = .000000001 #at z = 0, there is an infinite gradient
            dspeeds_dz[i][i] = Vel/((zref-z0)**shearExp)*shearExp*(z[i])**(shearExp-1.)

        J['towerSpeeds', 'z_full'] = dspeeds_dz
        J['towerSpeeds', 'Vel'] = dspeeds_dV

        return J


class dynamic_q(Component):
    """calculate the dynamic pressure"""
    def __init__(self, nFull):

        super(dynamic_q, self).__init__()

        self.nFull = nFull

        self.add_param('rhoAir', 1.225, units='kg/m**3', desc='density of air')
        self.add_param('towerSpeeds', np.zeros(nFull), desc='speeds of interest')

        self.add_output('q_dyn', np.zeros(nFull), desc='dynamic pressure at each point')

    def solve_nonlinear(self, params, unknowns, resids):

        rho = params['rhoAir']
        speeds = params['towerSpeeds']

        q_dyn = 0.5*rho*speeds**2

        unknowns['q_dyn'] = q_dyn

    def linearize(self, params, unknowns, resids):

        rho = params['rhoAir']
        speeds = params['towerSpeeds']
        q_dyn = unknowns['q_dyn']
        nFull = self.nFull

        J = {}
        dq_dspeeds = np.zeros((nFull,nFull))
        for i in range(nFull):
            dq_dspeeds[i][i] = rho*speeds[i] #The garient will be wrong at z=0, but this shouldn't matter as z0 remains fixed

        J['q_dyn', 'towerSpeeds'] = dq_dspeeds
        return J


class hoopStressEurocode(Component):
    """Hoop stress at each point"""

    #inputs

    def __init__(self, nFull):

        super(hoopStressEurocode, self).__init__()

        self.nFull = nFull

        self.add_param('d_full', np.zeros(nFull), desc='diameter at each point')
        self.add_param('t_full', np.zeros(nFull), desc='thickness at each point')
        self.add_param('L_reinforced', np.zeros(nFull), units='m')
        self.add_param('q_dyn', np.zeros(nFull), desc='dynamic pressure at each point')

        self.add_output('hoop_stress', np.zeros(nFull), desc='hoop stress at each point')

    def solve_nonlinear(self, params, unknowns, resids):

        d = params['d_full']
        t = params['t_full']
        L_reinforced = params['L_reinforced']

        q_dyn = params['q_dyn']
        r = d/2.0-t/2.0  # radius of cylinder middle surface
        omega = L_reinforced/np.sqrt(r*t)

        C_theta = 1.5  # clamped-clamped
        k_w = 0.46*(1.0 + 0.1*np.sqrt(C_theta/omega*r/t))
        Peq = k_w*q_dyn
        hoop_stress = -Peq*r/t
        dhoop_dq = -k_w*r/t
        self.dhoop_dq = dhoop_dq

        unknowns['hoop_stress'] = hoop_stress

    def linearize(self, params, unknowns, resids):
        #TODO fix this

        nFull = self.nFull

        d = params['d_full']
        t = params['t_full']
        L_reinforced = params['L_reinforced']
        q_dyn = params['q_dyn']

        C_theta = 1.5
        r = d/2.0-t/2.0
        omega = L_reinforced/np.sqrt(r*t)

        hoop_stress = unknowns['hoop_stress']

        J = {}
        dhoop_dD_full = np.zeros((nFull, nFull))
        dhoop_dT_full = np.zeros((nFull, nFull))

        d1 = 0.023*q_dyn*r
        d2 = (0.25*C_theta*r/(L_reinforced*np.sqrt(r*t)))+(0.5*C_theta*np.sqrt(r*t)/(L_reinforced*t))
        d3 = t*np.sqrt(C_theta*r*np.sqrt(r*t)/(L_reinforced*t))
        d4 = 0.23*q_dyn*(1.+0.1*np.sqrt(C_theta*r*np.sqrt(r*t)/(L_reinforced*t)))/t
        der_d = -d1*d2/d3-d4

        t1 = 0.023*q_dyn*r
        t2 = (C_theta*(0.5*d-t)*r/(2*L_reinforced*t*np.sqrt(r*t)))-\
            (C_theta*r*np.sqrt(r*t)/(L_reinforced*t**2))-\
            (0.5*C_theta*np.sqrt(r*t)/(L_reinforced*t))
        t3 = t*np.sqrt(C_theta*r*np.sqrt(r*t)/(L_reinforced*t))
        t4 = 0.46*q_dyn*r*(1.+0.1*np.sqrt(C_theta*r*np.sqrt(r*t)/(L_reinforced*t)))/(t**2)
        t5 = 0.23*q_dyn*(1.+0.1*np.sqrt(C_theta*r*np.sqrt(r*t)/(L_reinforced*t)))/t
        der_t = -t1*t2/t3+t4+t5

        dhoop_dD_full = np.zeros((nFull,nFull))
        dhoop_dT_full = np.zeros((nFull,nFull))

        for i in range(nFull):
            dhoop_dD_full[i][i] = der_d[i]
            dhoop_dT_full[i][i] = der_t[i]

        dhoop_dq = np.zeros((nFull,nFull))
        for i in range(nFull):
            dhoop_dq[i][i] = self.dhoop_dq[i]
        J['hoop_stress', 'd_full'] = dhoop_dD_full
        J['hoop_stress', 't_full'] = dhoop_dT_full
        J['hoop_stress', 'q_dyn'] = dhoop_dq

        return J


class Fz_comp(Component):
    """Fz from gravity"""

    def __init__(self):

        super(Fz_comp, self).__init__()
        self.add_param('m', np.array([0.]), units='kg', desc='mass at top')
        self.add_output('Fz', 0.0, units='kg', desc='z force')

    def solve_nonlinear(self,params,unknowns,resids):
        unknowns['Fz'] = -9.81*float(params['m'])

    def linearize(self,params,unknowns,resids):
        J = {}
        J['Fz','m'] = -9.81

        return J

class axial_and_shear(Component):
    """axial stress at each point"""

    def __init__(self, nFull):

        super(axial_and_shear, self).__init__()

        self.nFull = nFull
        self.add_param('m', np.array([0.]), units='kg', desc='mass at top')
        self.add_param('d_full', np.zeros(nFull), desc='diameter at each point')
        self.add_param('t_full', np.zeros(nFull), desc='thickness at each point')
        self.add_param('z_full', np.zeros(nFull), units='m', desc='location on tower')
        self.add_param('Fx', 0.0, desc='fx force at top of the tower')
        self.add_param('Fy', 0.0, desc='fy at the top of the tower')
        self.add_param('Fz', 0.0, desc='z force at top of tower')
        self.add_param('qx', np.zeros(nFull), desc='wind load in x')
        self.add_param('qy', np.zeros(nFull), desc='wind load in y')
        self.add_param('Mxx', 0.0, desc='moments at the top of the tower, xx')
        self.add_param('Myy', 0.0, desc='moments at the top of the tower, yy')
        self.add_param('rho', np.zeros(nFull), desc='density of the tower')
        self.add_param('mrhox', np.zeros([1]), units='m', desc='center of mass displacement in x') #need to change the length of this array as necessary

        self.add_output('axial_stress', np.zeros(nFull), desc='axial stress at each point')
        self.add_output('shear_stress', np.zeros(nFull), desc='shear stress at each point')

    def solve_nonlinear(self, params, unknowns, resids):

        m = params['m']
        d_full = params['d_full']
        t_full = params['t_full']
        z_full = params['z_full']
        Fz = params['Fz']
        Fx = params['Fx']
        Fy = params['Fy']
        Mxx = params['Mxx']
        Myy = params['Myy']
        qx = params['qx']
        qy = params['qy']
        rho = params['rho']
        mrhox = params['mrhox']

        # Fz = m * -9.81

        # print 'Fz: ', Fz
        # print 'm: ', m

        axial_stress, shear_stress = _axialShear.axial_and_shear(d_full, t_full, z_full, Fx, Fy, Fz, qx, qy, \
                                    Mxx, Myy, rho, mrhox, m)

        unknowns['axial_stress'] = axial_stress
        unknowns['shear_stress'] = shear_stress

    #TODO need gradients wrt other inputs (m, Fx, Mxx, Myy, qx)
    def linearize(self, params, unknowns, resids):

        m = params['m']
        d_full = params['d_full']
        t_full = params['t_full']
        z_full = params['z_full']
        Fz = m * -9.81
        Fx = params['Fx']
        Fy = params['Fy']
        Mxx = params['Mxx']
        Myy = params['Myy']
        qx = params['qx']
        qy = params['qy']
        rho = params['rho']
        mrhox = params['mrhox']
        nFull = self.nFull

        #wrt diameter
        dd = np.eye(nFull)
        td = np.zeros((nFull, nFull))
        zd = np.zeros((nFull, nFull))
        Fxd = np.zeros(nFull)
        Fyd = np.zeros(nFull)
        Fzd = np.zeros(nFull)
        qxd = np.zeros((nFull, nFull))
        qyd = np.zeros((nFull, nFull))
        Mxxd = np.zeros(nFull)
        Myyd = np.zeros(nFull)
        md = np.zeros(nFull)
        axial_stress,da_dd,shear_stress,ds_dd = _axialShear.axial_and_shear_dv(d_full, \
                            dd,t_full,td,z_full,zd,Fx,Fxd,Fy,Fyd,Fz,Fzd,qx,qxd,qy,qyd,Mxx,   \
                            Mxxd,Myy,Myyd,rho,mrhox,m,md)

        #wrt thickness
        dd = np.zeros((nFull, nFull))
        td = np.eye(nFull)
        zd = np.zeros((nFull, nFull))
        Fxd = np.zeros(nFull)
        Fyd = np.zeros(nFull)
        Fzd = np.zeros(nFull)
        qxd = np.zeros((nFull, nFull))
        qyd = np.zeros((nFull, nFull))
        Mxxd = np.zeros(nFull)
        Myyd = np.zeros(nFull)
        md = np.zeros(nFull)
        axial_stress,da_dt,shear_stress,ds_dt = _axialShear.axial_and_shear_dv(d_full, \
                            dd,t_full,td,z_full,zd,Fx,Fxd,Fy,Fyd,Fz,Fzd,qx,qxd,qy,qyd,Mxx,   \
                            Mxxd,Myy,Myyd,rho,mrhox,m,md)

        #wrt z
        dd = np.zeros((nFull, nFull))
        td = np.zeros((nFull, nFull))
        zd = np.eye(nFull)
        Fxd = np.zeros(nFull)
        Fyd = np.zeros(nFull)
        Fzd = np.zeros(nFull)
        qxd = np.zeros((nFull, nFull))
        qyd = np.zeros((nFull, nFull))
        Mxxd = np.zeros(nFull)
        Myyd = np.zeros(nFull)
        md = np.zeros(nFull)
        axial_stress,da_dz,shear_stress,ds_dz = _axialShear.axial_and_shear_dv(d_full, \
                            dd,t_full,td,z_full,zd,Fx,Fxd,Fy,Fyd,Fz,Fzd,qx,qxd,qy,qyd,Mxx,   \
                            Mxxd,Myy,Myyd,rho,mrhox,m,md)

        #wrt Fx
        dd = np.zeros((nFull, nFull))
        td = np.zeros((nFull, nFull))
        zd = np.zeros((nFull, nFull))
        Fxd = np.ones(nFull)
        Fyd = np.zeros(nFull)
        Fzd = np.zeros(nFull)
        qxd = np.zeros((nFull, nFull))
        qyd = np.zeros((nFull, nFull))
        Mxxd = np.zeros(nFull)
        Myyd = np.zeros(nFull)
        md = np.zeros(nFull)
        axial_stress,da_dFx,shear_stress,ds_dFx = _axialShear.axial_and_shear_dv(d_full, \
                            dd,t_full,td,z_full,zd,Fx,Fxd,Fy,Fyd,Fz,Fzd,qx,qxd,qy,qyd,Mxx,   \
                            Mxxd,Myy,Myyd,rho,mrhox,m,md)

        #wrt Fy
        dd = np.zeros((nFull, nFull))
        td = np.zeros((nFull, nFull))
        zd = np.zeros((nFull, nFull))
        Fxd = np.zeros(nFull)
        Fyd = np.ones(nFull)
        Fzd = np.zeros(nFull)
        qxd = np.zeros((nFull, nFull))
        qyd = np.zeros((nFull, nFull))
        Mxxd = np.zeros(nFull)
        Myyd = np.zeros(nFull)
        md = np.zeros(nFull)
        axial_stress,da_dFy,shear_stress,ds_dFy = _axialShear.axial_and_shear_dv(d_full, \
                            dd,t_full,td,z_full,zd,Fx,Fxd,Fy,Fyd,Fz,Fzd,qx,qxd,qy,qyd,Mxx,   \
                            Mxxd,Myy,Myyd,rho,mrhox,m,md)

        #wrt Fz
        dd = np.zeros((nFull, nFull))
        td = np.zeros((nFull, nFull))
        zd = np.zeros((nFull, nFull))
        Fxd = np.zeros(nFull)
        Fyd = np.zeros(nFull)
        Fzd = np.ones(nFull)
        qxd = np.zeros((nFull, nFull))
        qyd = np.zeros((nFull, nFull))
        Mxxd = np.zeros(nFull)
        Myyd = np.zeros(nFull)
        md = np.zeros(nFull)
        axial_stress,da_dFz,shear_stress,ds_dFz = _axialShear.axial_and_shear_dv(d_full, \
                            dd,t_full,td,z_full,zd,Fx,Fxd,Fy,Fyd,Fz,Fzd,qx,qxd,qy,qyd,Mxx,   \
                            Mxxd,Myy,Myyd,rho,mrhox,m,md)

        #wrt qx (wind loads in x)
        dd = np.zeros((nFull, nFull))
        td = np.zeros((nFull, nFull))
        zd = np.zeros((nFull, nFull))
        Fxd = np.zeros(nFull)
        Fyd = np.zeros(nFull)
        Fzd = np.zeros(nFull)
        qxd = np.eye(nFull)
        qyd = np.zeros((nFull, nFull))
        Mxxd = np.zeros(nFull)
        Myyd = np.zeros(nFull)
        md = np.zeros(nFull)
        axial_stress,da_dqx,shear_stress,ds_dqx = _axialShear.axial_and_shear_dv(d_full, \
                            dd,t_full,td,z_full,zd,Fx,Fxd,Fy,Fyd,Fz,Fzd,qx,qxd,qy,qyd,Mxx,   \
                            Mxxd,Myy,Myyd,rho,mrhox,m,md)

        #wrt qy (wind loads in y)
        dd = np.zeros((nFull, nFull))
        td = np.zeros((nFull, nFull))
        zd = np.zeros((nFull, nFull))
        Fxd = np.zeros(nFull)
        Fyd = np.zeros(nFull)
        Fzd = np.zeros(nFull)
        qxd = np.zeros((nFull, nFull))
        qyd = np.eye(nFull)
        Mxxd = np.zeros(nFull)
        Myyd = np.zeros(nFull)
        md = np.zeros(nFull)
        axial_stress,da_dqy,shear_stress,ds_dqy = _axialShear.axial_and_shear_dv(d_full, \
                            dd,t_full,td,z_full,zd,Fx,Fxd,Fy,Fyd,Fz,Fzd,qx,qxd,qy,qyd,Mxx,   \
                            Mxxd,Myy,Myyd,rho,mrhox,m,md)

        #wrt Mxx
        dd = np.zeros((nFull, nFull))
        td = np.zeros((nFull, nFull))
        zd = np.zeros((nFull, nFull))
        Fxd = np.zeros(nFull)
        Fyd = np.zeros(nFull)
        Fzd = np.zeros(nFull)
        qxd = np.zeros((nFull, nFull))
        qyd =  np.zeros((nFull, nFull))
        Mxxd = np.ones(nFull)
        Myyd = np.zeros(nFull)
        md = np.zeros(nFull)
        axial_stress,da_dMxx,shear_stress,ds_dMxx = _axialShear.axial_and_shear_dv(d_full, \
                            dd,t_full,td,z_full,zd,Fx,Fxd,Fy,Fyd,Fz,Fzd,qx,qxd,qy,qyd,Mxx,   \
                            Mxxd,Myy,Myyd,rho,mrhox,m,md)

        #wrt Myy
        dd = np.zeros((nFull, nFull))
        td = np.zeros((nFull, nFull))
        zd = np.zeros((nFull, nFull))
        Fxd = np.zeros(nFull)
        Fyd = np.zeros(nFull)
        Fzd = np.zeros(nFull)
        qxd = np.zeros((nFull, nFull))
        qyd =  np.zeros((nFull, nFull))
        Mxxd = np.zeros(nFull)
        Myyd = np.ones(nFull)
        md = np.zeros(nFull)
        axial_stress,da_dMyy,shear_stress,ds_dMyy = _axialShear.axial_and_shear_dv(d_full, \
                            dd,t_full,td,z_full,zd,Fx,Fxd,Fy,Fyd,Fz,Fzd,qx,qxd,qy,qyd,Mxx,   \
                            Mxxd,Myy,Myyd,rho,mrhox,m,md)

        #wrt m
        dd = np.zeros((nFull, nFull))
        td = np.zeros((nFull, nFull))
        zd = np.zeros((nFull, nFull))
        Fxd = np.zeros(nFull)
        Fyd = np.zeros(nFull)
        Fzd = np.zeros(nFull)
        qxd = np.zeros((nFull, nFull))
        qyd =  np.zeros((nFull, nFull))
        Mxxd = np.zeros(nFull)
        Myyd = np.zeros(nFull)
        md = np.ones(nFull)
        axial_stress,da_dm,shear_stress,ds_dm = _axialShear.axial_and_shear_dv(d_full, \
                            dd,t_full,td,z_full,zd,Fx,Fxd,Fy,Fyd,Fz,Fzd,qx,qxd,qy,qyd,Mxx,   \
                            Mxxd,Myy,Myyd,rho,mrhox,m,md)

        J = {}

        J['axial_stress', 'd_full'] = da_dd.T
        J['axial_stress', 't_full'] = da_dt.T
        J['axial_stress', 'z_full'] = da_dz.T
        J['axial_stress', 'Fx'] = np.diagonal(da_dFx)
        J['axial_stress', 'Fy'] = np.diagonal(da_dFy)
        J['axial_stress', 'Fz'] = np.diagonal(da_dFz)
        J['axial_stress', 'qx'] = da_dqx.T
        J['axial_stress', 'qy'] = da_dqy.T
        J['axial_stress', 'Mxx'] = np.diagonal(da_dMxx)
        J['axial_stress', 'Myy'] = np.diagonal(da_dMyy)
        J['axial_stress', 'm'] = np.diagonal(da_dm)

        J['shear_stress', 'd_full'] = ds_dd.T
        J['shear_stress', 't_full'] = ds_dt.T
        J['shear_stress', 'z_full'] = ds_dz.T
        J['shear_stress', 'Fx'] = np.diagonal(ds_dFx)
        J['shear_stress', 'Fy'] = np.diagonal(ds_dFy)
        J['shear_stress', 'Fz'] = np.diagonal(ds_dFz)
        J['shear_stress', 'qx'] = ds_dqx.T
        J['shear_stress', 'qy'] = ds_dqy.T
        J['shear_stress', 'Mxx'] = np.diagonal(ds_dMxx)
        J['shear_stress', 'Myy'] = np.diagonal(ds_dMyy)
        J['shear_stress', 'm'] = np.diagonal(ds_dm)

        return J


class shellBuckling(Component):
    def __init__(self, nFull):

        super(shellBuckling, self).__init__()

        self.deriv_options['form'] = 'central'
        self.deriv_options['step_size'] = 1.E-6
        self.deriv_options['step_type'] = 'relative'

        self.nFull = nFull
        self.add_param('d_full', np.zeros(nFull), desc='diameter at specified locations')
        self.add_param('t_full', np.zeros(nFull), desc='thickness at specified locations')
        self.add_param('axial_stress', np.zeros(nFull), desc='axial stress at specified locations')
        self.add_param('hoop_stress', np.zeros(nFull), desc='shoop stress at specified locations')
        self.add_param('shear_stress', np.zeros(nFull), desc='shear stress at specified locations')
        self.add_param('L_reinforced', np.zeros(nFull), units='m')
        self.add_param('E', np.zeros(nFull), units='N/m**2', desc='modulus of elasticity')
        self.add_param('sigma_y', np.zeros(nFull), units='N/m**2', desc='yield stress')
        self.add_param('gamma_f', 1.35)
        self.add_param('gamma_b', 1.1)

        self.add_output('shell_buckling', np.zeros(nFull), desc='shell buckling at each point')


    def solve_nonlinear(self, params, unknowns, resids):
        d_full = params['d_full']
        t_full = params['t_full']
        axial_stress = params['axial_stress']
        hoop_stress = params['hoop_stress']
        shear_stress = params['shear_stress']
        L_reinforced = params['L_reinforced']
        E = params['E']
        sigma_y = params['sigma_y']
        gamma_f = params['gamma_f']
        gamma_b = params['gamma_b']

        # axial_stress = np.array([ -1.79803856e+08, -1.18724763e+08, -3.99967970e+07])
        # # hoop_stress = np.array([     -0.    ,     -103157.96160052 ,-140123.42531843])
        # shear_stress = np.array([ 5051257.78080602,  5593468.41665133 , 6363968.7620388 ])
        unknowns['shell_buckling'] = _shellbuckling.shellbucklingeurocode(d_full,
                                    t_full, axial_stress, hoop_stress, shear_stress,
                                    L_reinforced, E, sigma_y, gamma_f, gamma_b)


    def linearize(self, params, unknowns, resids):
        d_full = params['d_full']
        t_full = params['t_full']
        axial_stress = params['axial_stress']
        hoop_stress = params['hoop_stress']
        shear_stress = params['shear_stress']
        L_reinforced = params['L_reinforced']
        E = params['E']
        sigma_y = params['sigma_y']
        gamma_f = params['gamma_f']
        gamma_b = params['gamma_b']
        nFull = self.nFull

        #wrt diameter
        dd = np.eye(nFull)
        td = np.zeros((nFull, nFull))
        sigma_zd = np.zeros((nFull, nFull))
        sigma_td = np.zeros((nFull, nFull))
        tau_ztd = np.zeros((nFull, nFull))
        eu_utilization, dbuckling_dd = _shellbuckling.shellbucklingeurocode_dv(d_full,dd, \
                t_full,td,axial_stress,sigma_zd,hoop_stress,sigma_td,shear_stress,tau_ztd, \
                L_reinforced,E,sigma_y,gamma_f,gamma_b)

        #wrt thickness
        dd = np.zeros((nFull, nFull))
        td = np.eye(nFull)
        sigma_zd = np.zeros((nFull, nFull))
        sigma_td = np.zeros((nFull, nFull))
        tau_ztd = np.zeros((nFull, nFull))
        eu_utilization, dbuckling_dt = _shellbuckling.shellbucklingeurocode_dv(d_full,dd, \
                t_full,td,axial_stress,sigma_zd,hoop_stress,sigma_td,shear_stress,tau_ztd, \
                L_reinforced,E,sigma_y,gamma_f,gamma_b)

        #wrt axial stress
        dd = np.zeros((nFull, nFull))
        td = np.zeros((nFull, nFull))
        sigma_zd = np.eye(nFull)
        sigma_td = np.zeros((nFull, nFull))
        tau_ztd = np.zeros((nFull, nFull))
        eu_utilization, dbuckling_dsigmaz = _shellbuckling.shellbucklingeurocode_dv(d_full,dd, \
                t_full,td,axial_stress,sigma_zd,hoop_stress,sigma_td,shear_stress,tau_ztd, \
                L_reinforced,E,sigma_y,gamma_f,gamma_b)

        #wrt hoop stress
        dd = np.zeros((nFull, nFull))
        td = np.zeros((nFull, nFull))
        sigma_zd = np.zeros((nFull, nFull))
        sigma_td = np.eye(nFull)
        tau_ztd = np.zeros((nFull, nFull))
        eu_utilization, dbuckling_dsigmat = _shellbuckling.shellbucklingeurocode_dv(d_full,dd, \
                t_full,td,axial_stress,sigma_zd,hoop_stress,sigma_td,shear_stress,tau_ztd, \
                L_reinforced,E,sigma_y,gamma_f,gamma_b)

        #wrt shear stress
        dd = np.zeros((nFull, nFull))
        td = np.zeros((nFull, nFull))
        sigma_zd = np.zeros((nFull, nFull))
        sigma_td = np.zeros((nFull, nFull))
        tau_ztd = np.eye(nFull)
        eu_utilization, dbuckling_dtauzt = _shellbuckling.shellbucklingeurocode_dv(d_full,dd, \
                t_full,td,axial_stress,sigma_zd,hoop_stress,sigma_td,shear_stress,tau_ztd, \
                L_reinforced,E,sigma_y,gamma_f,gamma_b)

        J = {}
        J['shell_buckling', 'd_full'] = dbuckling_dd
        J['shell_buckling', 't_full'] = dbuckling_dt
        J['shell_buckling', 'axial_stress'] = dbuckling_dsigmaz
        J['shell_buckling', 'hoop_stress'] = dbuckling_dsigmat
        J['shell_buckling', 'shear_stress'] = dbuckling_dtauzt

        return J


class shellBucklingGroup(Group):
    """Group to calculate the shell buckling"""

    def __init__(self, nPoints, nFull):

        super(shellBucklingGroup, self).__init__()

        self.add('Fz_comp', Fz_comp(), promotes=['*'])
        self.add('hoopStressEurocode', hoopStressEurocode(nFull), promotes=['*'])
        self.add('axial_and_shear', axial_and_shear(nFull), promotes=['*'])
        self.add('shellBuckling', shellBuckling(nFull), promotes=['*'])


class averageI(Component):

    def __init__(self, nFull):

        super(averageI, self).__init__()

        self.nFull = nFull
        self.add_param('d_full', np.zeros(nFull), units='m', desc='diameters of the tower')
        self.add_param('t_full', np.zeros(nFull), units='m', desc='thickness of the tower')

        self.add_output('I', 0.0, desc='average moment of inertia of the tower')

    def solve_nonlinear(self, params, unknowns, resids):

        d_full = params['d_full']
        t_full = params['t_full']
        nFull = self.nFull

        """average d and t then find I, thin wall cylinder"""
        # d_avg = np.sum(d_full)/nFull
        # t_avg = np.sum(t_full)/nFull
        #
        # I = pi*(d_avg/2.)**3*t_avg
        # print 'Average D and T: ', I
        # unknowns['I'] = I

        """average the moments of inertia"""
        I = pi*(d_full/2.)**3*t_full
        # print 'Average Moments of Inertia: ', np.sum(I)/nFull
        unknowns['I'] = np.sum(I)/nFull

    def linearize(self, params, unknowns, resids):

        d_full = params['d_full']
        t_full = params['t_full']
        nFull = self.nFull

        J = {}
        dd = pi*t_full/120.*3.*d_full**2
        dt = pi/120.*(d_full)**3
        dI_dd = np.zeros((1,nFull))
        dI_dt = np.zeros((1,nFull))
        for i in range(nFull):
            dI_dd[0][i] = dd[i]
            dI_dt[0][i] = dt[i]
        J['I','d_full'] = dI_dd
        J['I','t_full'] = dI_dt

        return J


class m_L(Component):

    def __init__(self):

        super(m_L, self).__init__()

        self.add_param('L', 0.0, units='m', desc='Height of the Tower')
        self.add_param('mass', 0.0, units='kg', desc='mass of the tower') #needs to be generalized for multiple m's

        self.add_output('m_L', 0.0, desc='first natural frequency')

    def solve_nonlinear(self, params, unknowns, resids):

        L = params['L']
        m = params['mass']

        unknowns['m_L'] = m/L

    def linearize(self, params, unknowns, resids):

        L = params['L']
        m = params['mass']

        J = {}
        J['m_L', 'mass'] = 1./L
        J['m_L', 'L'] = -1.*m/L**2

        return J


class freq(Component):

    def __init__(self, nFull):

        super(freq, self).__init__()

        # self.deriv_options['type'] = 'fd'
        # self.deriv_options['form'] = 'central'
        # self.deriv_options['step_size'] = 1.E-4
        # self.deriv_options['step_type'] = 'relative'

        self.add_param('L', 0.0, units='m', desc='Height of the Tower')
        self.add_param('m_L', 0.0, units='kg', desc='mass of the tower per unit length') #needs to be generalized for multiple m's
        self.add_param('I', 0.0, desc='averaged moment of inertia of the tower')
        self.add_param('E', np.zeros(nFull), units = 'N/m**2', desc='Modulus of Elasticity')
        self.add_param('Mt', np.zeros(1), units='kg', desc='mass at the top of the tower')
        self.add_param('It', 0.0, desc='moment of inertia at the top of the tower')

        self.add_output('freq', 0.0, desc='first natural frequency')


    def solve_nonlinear(self, params, unknowns, resids):

        L = params['L']
        m = params['m_L']
        I = params['I']
        E = params['E'][0]
        Mt = params['Mt']
        It = params['It']


        # # constant
        # Mt = 3.4
        # It = 6.4
        # E = 1.2
        #
        # # variable
        # m = 2.6
        # L = 5.3
        # I = 3.6

        def R(lam, m, L, I):
            return 1 + cos(lam)*cosh(lam) + lam*Mt/(m*L)*(cos(lam)*sinh(lam) - sin(lam)*cosh(lam)) \
            - lam**3*It/(m*L**3)*(cosh(lam)*sin(lam) + sinh(lam)*cos(lam)) \
            + lam**4*Mt*It/(m**2*L**4)*(1 - cos(lam)*cosh(lam))

        def freq(x):
            m = x[0]
            L = x[1]
            I = x[2]

            # l = brentq(R, 0.0, 100., args=(m, L, I))
            l = root(R, 1.0, args=(m, L, I))['x'][0]
            # print 'm: ', m
            omega = l**2*sqrt(E*I/(m*L**4))
            f = omega/(2*pi)  # divided by 2pi to give Hz

            sl = sin(l)
            cl = cos(l)
            shl = sinh(l)
            chl = cosh(l)

            pfpm = -omega/(4*pi*m)
            pfpl = omega/(l*pi)
            prpm = -l*Mt/(m**2*L)*(cl*shl - sl*chl) \
                + l**3*It/(m**2*L**3)*(chl*sl + shl*cl) \
                - 2*l**4*Mt*It/(m**3*L**4)*(1 - cl*chl)
            prpl = -2*It*l**3*cos(l)*cosh(l)/(L**3*m) - 3*It*l**2*(sin(l)*cosh(l) + cos(l)*sinh(l))/(L**3*m) + It*Mt*l**4*(sin(l)*cosh(l) - cos(l)*sinh(l))/(L**4*m**2) + 4*It*Mt*l**3*(-cos(l)*cosh(l) + 1)/(L**4*m**2) - sin(l)*cosh(l) + cos(l)*sinh(l) - 2*Mt*l*sin(l)*sinh(l)/(L*m) + Mt*(-sin(l)*cosh(l) + cos(l)*sinh(l))/(L*m)

            """PJ adding this stuff"""
            prpMt = It*l**4*(1.-cos(l)*cosh(l))/(L**4*m**2) \
                + l*(-cosh(l)*sin(l)+cos(l)*sinh(l))/(L*m)
            prpIt = l**4*Mt*(1.-cos(l)*cosh(l))/(L**4*m**2) \
                - l**3*(cosh(l)*sin(l)+cos(l)*sinh(l))/(L**3*m)
            """"""


            pfpL = -omega/(pi*L)
            prpL = -l*Mt/(m*L**2)*(cl*shl - sl*chl) \
                + 3*l**3*It/(m*L**4)*(chl*sl + shl*cl) \
                - 4*l**4*Mt*It/(m**2*L**5)*(1 - cl*chl)

            pfpI = omega/(4*pi*I)


            dfdm = pfpm - pfpl*prpm/prpl
            dfdL = pfpL - pfpl*prpL/prpl
            dfdI = pfpI

            """PJ adding this stuff"""
            dfdMt = -pfpl*prpMt/prpl
            dfdIt = -pfpl*prpIt/prpl
            """"""

            dfdx = np.array([dfdm, dfdL, dfdI, dfdMt, dfdIt])

            return f, dfdx

        x = np.array([m, L, I])
        f, self.dfdx = freq(x)

        unknowns['freq'] = f

    def linearize(self, params, unknowns, resids):

        dfdx = self.dfdx

        J = {}
        J['freq', 'L'] = dfdx[1]
        J['freq', 'm_L'] = dfdx[0]
        J['freq', 'I'] = dfdx[2]
        J['freq', 'Mt'] = dfdx[3]
        J['freq', 'It'] = dfdx[4]

        return J
