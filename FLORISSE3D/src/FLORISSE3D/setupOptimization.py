import numpy as np
from openmdao.api import Problem, Group, IndepVarComp, pyOptSparseDriver, ExecComp, ScipyOptimizer
from scipy.spatial import ConvexHull
from FLORISSE3D.GeneralWindFarmComponents import calculate_boundary
import os
from rotorse.precomp import Profile, Orthotropic2DMaterial, CompositeSection, _precomp
from scipy.interpolate import RectBivariateSpline
from akima import Akima, akima_interp

def setupTower(n, prob):
    prob['L_reinforced'] = 30.0*np.ones(n)  # [m] buckling length
    # prob['Toweryaw'] = 0.0

    # --- material props ---
    prob['E'] = 210.e9*np.ones(n)
    # prob['G'] = 80.8e9*np.ones(n)
    prob['rho'] = 8500.0*np.ones(n)
    prob['sigma_y'] = 450.0e6*np.ones(n)

    # --- extra mass ----
    # prob['midx'] = np.array([n-1], dtype=int)  # RNA mass at top
    # prob['m'] = np.array([285598.8])
    # prob['mIxx'] = np.array([1.14930678e+08])
    # prob['mIyy'] = np.array([2.20354030e+07])
    # prob['mIzz'] = np.array([1.87597425e+07])
    # prob['mIxy'] = np.array([0.00000000e+00])
    # prob['mIxz'] = np.array([5.03710467e+05])
    # prob['mIyz'] = np.array([0.00000000e+00])
    prob['mrhox'] = np.array([-1.13197635]) # Does not change with rotor_diameter
    # prob['mrhoy'] = np.array([0.])
    # prob['mrhoz'] = np.array([0.50875268])
    # prob['nMass'] = len(prob['midx'])
    # prob['addGravityLoadForExtraMass'] = True
    # -----------

    # --- wind ---
    prob['zref'] = 50.0
    prob['z0'] = 0.0
    # ---------------

    # if addGravityLoadForExtraMass=True be sure not to double count by adding those force here also
    # # --- loading case 1: max Thrust ---
    # wind_Uref1 = 11.73732
    # prob['plidx1'] = np.array([n-1], dtype=int)  # at  top
    # prob['Fx1'] = np.array([1284744.19620519])
    # prob['Fy1'] = np.array([0.])
    # prob['Fz1'] = np.array([-2914124.84400512])
    # prob['Mxx1'] = np.array([3963732.76208099])
    # prob['Myy1'] = np.array([-2275104.79420872])
    # prob['Mzz1'] = np.array([-346781.68192839])
    # prob['nPL'] = len(prob['plidx1'])
    # # ---------------

    # # --- loading case 2: max wind speed ---
    # wind_Uref2 = 70.0
    # prob['plidx2'] = np.array([n-1], dtype=int)  # at  top
    # prob['Fx2'] = np.array([930198.60063279])
    # prob['Fy2'] = np.array([0.])
    # prob['Fz2'] = np.array([-2883106.12368949])
    # prob['Mxx2'] = np.array([-1683669.22411597])
    # prob['Myy2'] = np.array([-2522475.34625363])
    # prob['Mzz2'] = np.array([147301.97023764])
    # # ---------------

    # --- safety factors ---
    prob['gamma_f'] = 1.35
    # prob['gamma_m'] = 1.3
    # prob['gamma_n'] = 1.0
    prob['gamma_b'] = 1.1
    # ---------------

    # --- constraints ---
    # prob['min_d_to_t'] = 120.0
    # prob['min_taper'] = 0.4
    # ---------------


def amaliaWind(prob):
    windSpeeds = np.array([6.53163342, 6.11908394, 6.13415514, 6.0614625,  6.21344602,
                                5.87000793, 5.62161519, 5.96779107, 6.33589422, 6.4668016,
                                7.9854581,  7.6894432,  7.5089221,  7.48638098, 7.65764618,
                                6.82414044, 6.36728201, 5.95982999, 6.05942132, 6.1176321,
                                5.50987893, 4.18461796, 4.82863115, 0.,         0.,         0.,
                                5.94115843, 5.94914252, 5.59386528, 6.42332524, 7.67904937,
                                7.89618066, 8.84560463, 8.51601497, 8.40826823, 7.89479475,
                                7.86194762, 7.9242645,  8.56269962, 8.94563889, 9.82636368,
                               10.11153102, 9.71402212, 9.95233636,  10.35446959, 9.67156182,
                                9.62462527, 8.83545158, 8.18011771, 7.9372492,  7.68726143,
                                7.88134508, 7.31394723, 7.01839896, 6.82858346, 7.06213432,
                                7.01949894, 7.00575122, 7.78735165, 7.52836352, 7.21392201,
                                7.4356621,  7.54099962, 7.61335262, 7.90293531, 7.16021596,
                                7.19617087, 7.5593657,  7.03278586, 6.76105501, 6.48004694,
                                6.94716392])

    windFrequencies = np.array([1.17812570e-02, 1.09958570e-02, 9.60626600e-03, 1.21236860e-02,
                               1.04722450e-02, 1.00695140e-02, 9.68687400e-03, 1.00090550e-02,
                               1.03715390e-02, 1.12172280e-02, 1.52249700e-02, 1.56279300e-02,
                               1.57488780e-02, 1.70577560e-02, 1.93535770e-02, 1.41980570e-02,
                               1.20632100e-02, 1.20229000e-02, 1.32111160e-02, 1.74605400e-02,
                               1.72994400e-02, 1.43993790e-02, 7.87436000e-03, 0.00000000e+00,
                               2.01390000e-05, 0.00000000e+00, 3.42360000e-04, 3.56458900e-03,
                               7.18957000e-03, 8.80068000e-03, 1.13583200e-02, 1.41576700e-02,
                               1.66951900e-02, 1.63125500e-02, 1.31709000e-02, 1.09153300e-02,
                               9.48553000e-03, 1.01097900e-02, 1.18819700e-02, 1.26069900e-02,
                               1.58895900e-02, 1.77021600e-02, 2.04208100e-02, 2.27972500e-02,
                               2.95438600e-02, 3.02891700e-02, 2.69861000e-02, 2.21527500e-02,
                               2.12465500e-02, 1.82861400e-02, 1.66147400e-02, 1.90111800e-02,
                               1.90514500e-02, 1.63932050e-02, 1.76215200e-02, 1.65341460e-02,
                               1.44597600e-02, 1.40370300e-02, 1.65745000e-02, 1.56278200e-02,
                               1.53459200e-02, 1.75210100e-02, 1.59702700e-02, 1.51041500e-02,
                               1.45201100e-02, 1.34527800e-02, 1.47819600e-02, 1.33923300e-02,
                               1.10562900e-02, 1.04521380e-02, 1.16201970e-02, 1.10562700e-02])

    # windSpeeds = windSpeeds* 1.714285714 #to get the average speed higher, close to 12 m/s
    nDirections = len(windSpeeds)
    windDirections = np.linspace(0,360-360/nDirections, nDirections)

    index = np.where(windSpeeds==0.0)
    windSpeeds = np.delete(windSpeeds, index[0])
    windFrequencies = np.delete(windFrequencies, index[0])
    windDirections = np.delete(windDirections, index[0])
    nDirections = len(windSpeeds)

    prob['Uref'] = windSpeeds
    prob['windDirections'] = windDirections
    prob['windFrequencies'] = windFrequencies

    return nDirections

    # return windSpeeds, windFrequencies, windDirections, nDirections


def simpleSetup(nTurbs, prob):
    axialInduction = np.zeros(nTurbs)
    Ct = np.zeros(nTurbs)
    Cp = np.zeros(nTurbs)
    generatorEfficiency = np.zeros(nTurbs)
    # define initial values
    for turbI in range(0, nTurbs):
        axialInduction[turbI] = 1.0/3.0
        Ct[turbI] = 4.0*axialInduction[turbI]*(1.0-axialInduction[turbI])
        # Cp[turbI] = 0.7737/0.944 * 4.0 * 1.0/3.0 * np.power((1 - 1.0/3.0), 2)
        # Cp[turbI] = 0.7737 * 4.0 * 1.0/3.0 * np.power((1. - 1.0/3.0), 2)
        Cp[turbI] = 0.42
        generatorEfficiency[turbI] = 1.0#0.944

        # print Cp[turbI]

    prob['axialInduction'] = axialInduction # kg/m^3
    prob['generatorEfficiency'] = generatorEfficiency
    prob['air_density'] = 1.1716
    prob['Ct_in'] = Ct
    prob['Cp_in'] = Cp
    prob['floris_params:cos_spread'] = 1E12 # turns off cosine spread (just needs to be very large)

    # prob['turbine_class'] = 'II/III'
    # prob['blade_has_carbon'] = True
    # prob['bearing_number'] = 2


def setupGrid(nRows, rotor_diameter, spacing):

    points = np.linspace(start=spacing*rotor_diameter, stop=nRows*spacing*rotor_diameter, num=nRows)
    xpoints, ypoints = np.meshgrid(points, points)
    turbineX = np.ndarray.flatten(xpoints)
    turbineY = np.ndarray.flatten(ypoints)

    return turbineX, turbineY


def setupBoundaryConstraints(turbineX_bounds, turbineY_bounds):
    # generate boundary constraint
    locations = np.zeros((len(turbineX_bounds),2))
    for i in range(len(turbineX_bounds)):
        locations[i][0] = turbineX_bounds[i]
        locations[i][1] = turbineY_bounds[i]

    boundaryVertices, boundaryNormals = calculate_boundary(locations)
    nVertices = boundaryVertices.shape[0]

    return nVertices, boundaryVertices, boundaryNormals


def setupRotor(nGroups, prob):
    # ROTOR stuff
        # === blade grid ===
    for i in range(nGroups):

        """CHANGE HERE"""
        prob['Rotor%s.initial_aero_grid'%i] = np.array([0.02222276, 0.06666667, 0.11111057, 0.16666667, 0.23333333, 0.3, 0.36666667,
            0.43333333, 0.5, 0.56666667, 0.63333333, 0.7, 0.76666667, 0.83333333, 0.88888943, 0.93333333,
            0.97777724])  # (Array): initial aerodynamic grid on unit radius
        # prob['Rotor%s.initial_aero_grid'%i] = np.array([0.0453591772, 0.886075949, 0.131850127, 0.1859177245, 0.2507911392, 0.315664557, 0.3805379747,
        #     0.4454113924, 0.5102848101, 0.5751582278, 0.6400316456, 0.7049050633, 0.769778481, 0.8346518987, 0.8887136076, 0.9319620253,
        #     0.975210443])
        prob['Rotor%s.initial_str_grid'%i] = np.array([0.0, 0.00492790457512, 0.00652942887106, 0.00813095316699, 0.00983257273154,
            0.0114340970275, 0.0130356213234, 0.02222276, 0.024446481932, 0.026048006228, 0.06666667, 0.089508406455,
            0.11111057, 0.146462614229, 0.16666667, 0.195309105255, 0.23333333, 0.276686558545, 0.3, 0.333640766319,
            0.36666667, 0.400404310407, 0.43333333, 0.5, 0.520818918408, 0.56666667, 0.602196371696, 0.63333333,
            0.667358391486, 0.683573824984, 0.7, 0.73242031601, 0.76666667, 0.83333333, 0.88888943, 0.93333333, 0.97777724,
            1.0])  # (Array): initial structural grid on unit radius
        prob['Rotor%s.idx_cylinder_aero'%i] = 3  # (Int): first idx in r_aero_unit of non-cylindrical section, constant twist inboard of here
        prob['Rotor%s.idx_cylinder_str'%i] = 14  # (Int): first idx in r_str_unit of non-cylindrical section
        prob['Rotor%s.hubFraction'%i] = 0.025  # (Float): hub location as fraction of radius
        # ------------------

        # === blade geometry ===
        """CHANGE HERE"""
        prob['Rotor%s.r_aero'%i] = np.array([0.02222276, 0.06666667, 0.11111057, 0.2, 0.23333333, 0.3, 0.36666667, 0.43333333,
            0.5, 0.56666667, 0.63333333, 0.64, 0.7, 0.83333333, 0.88888943, 0.93333333,
            0.97777724])  # (Array): new aerodynamic grid on unit radius
        # prob['Rotor%s.r_aero'%i] = np.array([0.0453591772, 0.886075949, 0.131850127, 0.1859177245, 0.2507911392, 0.315664557, 0.3805379747,
        #     0.4454113924, 0.5102848101, 0.5751582278, 0.6400316456, 0.7049050633, 0.769778481, 0.8346518987, 0.8887136076, 0.9319620253,
        #     0.975210443])
        """CHANGE HERE"""
        prob['Rotor%s.r_max_chord'%i] = 0.23577  # (Float): location of max chord on unit radius
        # prob['Rotor%s.r_max_chord'%i] = 0.2507911392
        """CHANGE HERE"""
        # prob['Rotor%s.chord_sub'%i] = np.array([3.542, 4.652, 3.256, 1.419])
        prob['Rotor%s.chord_sub'%i] = np.array([3.2612, 4.5709, 3.3178, 1.4621])  # (Array, m): chord at control points. defined at hub, then at linearly spaced locations from r_max_chord to tip
        # prob['chord_sub'] = np.array([2.2612, 4.5709, 3.3178, 1.4621])
        prob['Rotor%s.theta_sub'%i] = np.array([13.2783, 7.46036, 2.89317, -0.0878099])  # (Array, deg): twist at control points.  defined at linearly spaced locations from r[idx_cylinder] to tip
        prob['Rotor%s.precurve_sub'%i] = np.array([0.0, 0.0, 0.0])  # (Array, m): precurve at control points.  defined at same locations at chord, starting at 2nd control point (root must be zero precurve)
        prob['Rotor%s.delta_precurve_sub'%i] = np.array([0.0, 0.0, 0.0])  # (Array, m): adjustment to precurve to account for curvature from loading
        prob['Rotor%s.sparT'%i] = np.array([0.05, 0.047754, 0.045376, 0.031085, 0.0061398])  # (Array, m): spar cap thickness parameters
        prob['Rotor%s.teT'%i] = np.array([0.1, 0.09569, 0.06569, 0.02569, 0.00569])  # (Array, m): trailing-edge thickness parameters
        # prob['Rotor%s.bladeLength'%i] = 61.5  # (Float, m): blade length (if not precurved or swept) otherwise length of blade before curvature
        # prob['Rotor%s.bladeLength'%i] = 31.5  # (Float, m): blade length (if not precurved or swept) otherwise length of blade before curvature

        prob['Rotor%s.delta_bladeLength'%i] = 0.0  # (Float, m): adjustment to blade length to account for curvature from loading
        prob['Rotor%s.precone'%i] = 2.5  # (Float, deg): precone angle
        prob['Rotor%s.tilt'%i] = 5.0  # (Float, deg): shaft tilt
        prob['Rotor%s.yaw'%i] = 0.0  # (Float, deg): yaw error
        prob['Rotor%s.nBlades'%i] = 3  # (Int): number of blades
        # ------------------

        # === airfoil files ===
        basepath = '/Users/ningrsrch/Dropbox/Programs/RotorSE/src/rotorse/5MW_AFFiles'

        # load all airfoils
        airfoil_types = [0]*8
        airfoil_types[0] = os.path.join(basepath, 'Cylinder1.dat')
        airfoil_types[1] = os.path.join(basepath, 'Cylinder2.dat')
        airfoil_types[2] = os.path.join(basepath, 'DU40_A17.dat')
        airfoil_types[3] = os.path.join(basepath, 'DU35_A17.dat')
        airfoil_types[4] = os.path.join(basepath, 'DU30_A17.dat')
        airfoil_types[5] = os.path.join(basepath, 'DU25_A17.dat')
        airfoil_types[6] = os.path.join(basepath, 'DU21_A17.dat')
        airfoil_types[7] = os.path.join(basepath, 'NACA64_A17.dat')

        # place at appropriate radial stations
        prob['Rotor%s.af_idx'%i] = np.array([0, 0, 1, 2, 3, 3, 4, 5, 5, 6, 6, 7, 7, 7, 7, 7, 7])
        prob['Rotor%s.af_str_idx'%i] = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4, \
                5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7])

        n = len(prob['Rotor%s.af_idx'%i])
        af = [0]*n
        for j in range(n):
            af[j] = airfoil_types[int(prob['Rotor%s.af_idx'%i][j])]

        print 'af: ', af
        print 'airfoil_types: ', airfoil_types
        prob['Rotor%s.airfoil_types'%i] = airfoil_types  # (List): names of airfoil file
        # ----------------------

        # === atmosphere ===
        prob['Rotor%s.rho'%i] = 1.225  # (Float, kg/m**3): density of air
        prob['Rotor%s.mu'%i] = 1.81206e-5  # (Float, kg/m/s): dynamic viscosity of air
        prob['Rotor%s.shearExp'%i] = 0.25  # (Float): shear exponent
        prob['Rotor%s.hubHt'%i] = np.array([90.0])  # (Float, m): hub height
        prob['Rotor%s.turbine_class'%i] = 'I'  # (Enum): IEC turbine class
        prob['Rotor%s.turbulence_class'%i] = 'B'  # (Enum): IEC turbulence class class
        prob['Rotor%s.cdf_reference_height_wind_speed'%i] = 50.0  # (Float): reference hub height for IEC wind speed (used in CDF calculation)
        prob['Rotor%s.g'%i] = 9.81  # (Float, m/s**2): acceleration of gravity
        # ----------------------

        # === control ===
        prob['Rotor%s.control:Vin'%i] = 3.0  # (Float, m/s): cut-in wind speed
        prob['Rotor%s.control:Vout'%i] = 25.0  # (Float, m/s): cut-out wind speed
        # prob['Rotor%s.control:ratedPower'%i] = 5e6  # (Float, W): rated power
        prob['Rotor%s.control:minOmega'%i] = 0.0  # (Float, rpm): minimum allowed prob rotation speed
        prob['Rotor%s.control:maxOmega'%i] = 12.0  # (Float, rpm): maximum allowed prob rotation speed
        prob['Rotor%s.control:tsr'%i] = 7.55  # (Float): tip-speed ratio in Region 2 (should be optimized externally)
        prob['Rotor%s.control:pitch'%i] = 0.0  # (Float, deg): pitch angle in region 2 (and region 3 for fixed pitch machines)
        prob['Rotor%s.pitch_extreme'%i] = 0.0  # (Float, deg): worst-case pitch at survival wind condition
        prob['Rotor%s.azimuth_extreme'%i] = 0.0  # (Float, deg): worst-case azimuth at survival wind condition
        prob['Rotor%s.VfactorPC'%i] = 0.7  # (Float): fraction of rated speed at which the deflection is assumed to representative throughout the power curve calculation
        # ----------------------

        # === aero and structural analysis options ===
        prob['Rotor%s.nSector'%i] = 4  # (Int): number of sectors to divide prob face into in computing thrust and power
        prob['Rotor%s.npts_coarse_power_curve'%i] = 20  # (Int): number of points to evaluate aero analysis at
        prob['Rotor%s.npts_spline_power_curve'%i] = 200  # (Int): number of points to use in fitting spline to power curve
        prob['Rotor%s.AEP_loss_factor'%i] = 1.0  # (Float): availability and other losses (soiling, array, etc.)
        prob['Rotor%s.drivetrainType'%i] = 'geared'  # (Enum)
        prob['Rotor%s.nF'%i] = 5  # (Int): number of natural frequencies to compute
        prob['Rotor%s.dynamic_amplication_tip_deflection'%i] = 1.35  # (Float): a dynamic amplification factor to adjust the static deflection calculation
        # ----------------------

        # === materials and composite layup  ===
        basepath = '/Users/ningrsrch/Dropbox/Programs/RotorSE/src/rotorse/5MW_PrecompFiles'

        materials = Orthotropic2DMaterial.listFromPreCompFile(os.path.join(basepath, 'materials.inp'))

        ncomp = len(prob['Rotor%s.initial_str_grid'%i])
        upper = [0]*ncomp
        lower = [0]*ncomp
        webs = [0]*ncomp
        profile = [0]*ncomp

        prob['Rotor%s.leLoc'%i] = np.array([0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.498, 0.497, 0.465, 0.447, 0.43, 0.411,
            0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
            0.4, 0.4, 0.4, 0.4])    # (Array): array of leading-edge positions from a reference blade axis (usually blade pitch axis). locations are normalized by the local chord length. e.g. leLoc[i] = 0.2 means leading edge is 0.2*chord[i] from reference axis.  positive in -x direction for airfoil-aligned coordinate system
        prob['Rotor%s.sector_idx_strain_spar'%i] = [2]*ncomp  # (Array): index of sector for spar (PreComp definition of sector)
        prob['Rotor%s.sector_idx_strain_te'%i] = [3]*ncomp  # (Array): index of sector for trailing-edge (PreComp definition of sector)
        web1 = np.array([-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 0.4114, 0.4102, 0.4094, 0.3876, 0.3755, 0.3639, 0.345, 0.3342, 0.3313, 0.3274, 0.323, 0.3206, 0.3172, 0.3138, 0.3104, 0.307, 0.3003, 0.2982, 0.2935, 0.2899,\
                0.2867, 0.2833, 0.2817, 0.2799, 0.2767, 0.2731, 0.2664, 0.2607, 0.2562, 0.1886, -1.0])
        web2 = np.array([-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 0.5886, 0.5868, 0.5854, 0.5508, 0.5315, 0.5131, 0.4831, 0.4658, 0.4687, 0.4726, 0.477, 0.4794, 0.4828, 0.4862, 0.4896, 0.493, 0.4997, 0.5018, 0.5065, 0.5101, \
                0.5133, 0.5167, 0.5183, 0.5201, 0.5233, 0.5269, 0.5336, 0.5393, 0.5438, 0.6114, -1.0])
        web3 = np.array([-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, \
                1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0])
        # prob['Rotor%s.chord_str_ref'%i] = np.array([3.2612, 3.3100915356, 3.32587052924, 3.34159388653, 3.35823798667, 3.37384375335,
        #     3.38939112914, 3.4774055542, 3.49839685, 3.51343645709, 3.87017220335, 4.04645623801, 4.19408216643,
        #      4.47641008477, 4.55844487985, 4.57383098262, 4.57285771934, 4.51914315648, 4.47677655262, 4.40075650022,
        #      4.31069949379, 4.20483735936, 4.08985563932, 3.82931757126, 3.74220276467, 3.54415796922, 3.38732428502,
        #      3.24931446473, 3.23421422609, 3.22701537997, 3.21972125648, 3.08979310611, 2.95152261813, 2.330753331,
        #      2.05553464181, 1.82577817774, 1.5860853279, 1.4621])  # (Array, m): chord distribution for reference section, thickness of structural layup scaled with reference thickness
        prob['Rotor%s.chord_str_ref'%i] = np.array([3.2612, 3.3100915356, 3.32587052924, 3.34159388653, 3.35823798667, 3.37384375335,
            3.38939112914, 3.4774055542, 3.49839685, 3.51343645709, 3.87017220335, 4.04645623801, 4.19408216643,
             4.47641008477, 4.55844487985, 4.57383098262, 4.57285771934, 4.51914315648, 4.47677655262, 4.40075650022,
             4.31069949379, 4.20483735936, 4.08985563932, 3.82931757126, 3.74220276467, 3.54415796922, 3.38732428502,
             3.24931446473, 3.23421422609, 3.22701537997, 3.21972125648, 3.08979310611, 2.95152261813, 2.330753331,
             2.05553464181, 1.82577817774, 1.5860853279, 1.4621])  # (Array, m): chord distribution for reference section, thickness of structural layup scaled with reference thickness
        prob['Rotor%s.thick_str_ref'%i] = np.array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 0.404457084248, 0.404457084248,
                                       0.349012780126, 0.349012780126, 0.349012780126, 0.349012780126, 0.29892003076, 0.29892003076, 0.25110545018, 0.25110545018, 0.25110545018, 0.25110545018,
                                       0.211298863564, 0.211298863564, 0.211298863564, 0.211298863564, 0.17933792591, 0.17933792591, 0.17933792591, 0.17933792591, 0.17933792591, 0.17933792591,
                                       0.17933792591, 0.17933792591])  # (Array, m): airfoil thickness distribution for reference section, thickness of structural layup scaled with reference thickness

        prob['Rotor%s.capTriaxThk'%i] = np.array([0.30, 0.29, 0.28, 0.275, 0.27])
        prob['Rotor%s.capCarbThk'%i] = np.array([4.2, 2.5, 1.0, 0.90, 0.658])
        prob['Rotor%s.tePanelTriaxThk'%i] = np.array([0.30, 0.29, 0.28, 0.275, 0.27])
        prob['Rotor%s.tePanelFoamThk'%i] = np.array([9.00, 7.00, 5.00, 3.00, 2.00])

        for j in range(ncomp):

            webLoc = []
            if web1[j] != -1:
                webLoc.append(web1[j])
            if web2[j] != -1:
                webLoc.append(web2[j])
            if web3[j] != -1:
                webLoc.append(web3[j])

            upper[j], lower[j], webs[j] = CompositeSection.initFromPreCompLayupFile(os.path.join(basepath, 'layup_' + str(j+1) + '.inp'), webLoc, materials)
            profile[j] = Profile.initFromPreCompFile(os.path.join(basepath, 'shape_' + str(j+1) + '.inp'))

        prob['Rotor%s.materials'%i] = materials  # (List): list of all Orthotropic2DMaterial objects used in defining the geometry
        prob['Rotor%s.upperCS'%i] = upper  # (List): list of CompositeSection objections defining the properties for upper surface
        prob['Rotor%s.lowerCS'%i] = lower  # (List): list of CompositeSection objections defining the properties for lower surface
        prob['Rotor%s.websCS'%i] = webs  # (List): list of CompositeSection objections defining the properties for shear webs
        prob['Rotor%s.profile'%i] = profile  # (List): airfoil shape at each radial position
        # --------------------------------------


        # === fatigue ===
        prob['Rotor%s.rstar_damage'%i] = np.array([0.000, 0.022, 0.067, 0.111, 0.167, 0.233, 0.300, 0.367, 0.433, 0.500,
            0.567, 0.633, 0.700, 0.767, 0.833, 0.889, 0.933, 0.978])  # (Array): nondimensional radial locations of damage equivalent moments
        prob['Rotor%s.Mxb_damage'%i] = 1e3*np.array([2.3743E+003, 2.0834E+003, 1.8108E+003, 1.5705E+003, 1.3104E+003,
            1.0488E+003, 8.2367E+002, 6.3407E+002, 4.7727E+002, 3.4804E+002, 2.4458E+002, 1.6339E+002,
            1.0252E+002, 5.7842E+001, 2.7349E+001, 1.1262E+001, 3.8549E+000, 4.4738E-001])  # (Array, N*m): damage equivalent moments about blade c.s. x-direction
        prob['Rotor%s.Myb_damage'%i] = 1e3*np.array([2.7732E+003, 2.8155E+003, 2.6004E+003, 2.3933E+003, 2.1371E+003,
            1.8459E+003, 1.5582E+003, 1.2896E+003, 1.0427E+003, 8.2015E+002, 6.2449E+002, 4.5229E+002,
            3.0658E+002, 1.8746E+002, 9.6475E+001, 4.2677E+001, 1.5409E+001, 1.8426E+000])  # (Array, N*m): damage equivalent moments about blade c.s. y-direction
        prob['Rotor%s.strain_ult_spar'%i] = 1.0e-2  # (Float): ultimate strain in spar cap
        prob['Rotor%s.strain_ult_te'%i] = 2500*1e-6 * 2   # (Float): uptimate strain in trailing-edge panels, note that I am putting a factor of two for the damage part only.
        prob['Rotor%s.eta_damage'%i] = 1.35*1.3*1.0  # (Float): safety factor for fatigue
        prob['Rotor%s.m_damage'%i] = 10.0  # (Float): slope of S-N curve for fatigue analysis
        prob['Rotor%s.N_damage'%i] = 365*24*3600*20.0  # (Float): number of cycles used in fatigue analysis  TODO: make function of rotation speed


def setupSimpleRotorSE():
    num = 100
    x = np.linspace(500.,10000.,num)
    y = np.linspace(25.,160.,num)

    X,Y = np.meshgrid(x,y)

    filename = 'src/florisse3D/optRotor/ratedT100.txt'
    openedFile = open(filename)
    loadedData = np.loadtxt(openedFile)
    ratedT = loadedData[:]

    filename = 'src/florisse3D/optRotor/ratedQ100.txt'
    openedFile = open(filename)
    loadedData = np.loadtxt(openedFile)
    ratedQ = loadedData[:]

    filename = 'src/florisse3D/optRotor/blade_mass100.txt'
    openedFile = open(filename)
    loadedData = np.loadtxt(openedFile)
    blade_mass = loadedData[:]

    filename = 'src/florisse3D/optRotor/Vrated100.txt'
    openedFile = open(filename)
    loadedData = np.loadtxt(openedFile)
    Vrated = loadedData[:]

    filename = 'src/florisse3D/optRotor/extremeT100.txt'
    openedFile = open(filename)
    loadedData = np.loadtxt(openedFile)
    extremeT = loadedData[:]

    global ratedT_func
    ratedT_func = RectBivariateSpline(x, y, ratedT)
    global ratedQ_func
    ratedQ_func = RectBivariateSpline(x, y, ratedQ)
    global blade_mass_func
    blade_mass_func = RectBivariateSpline(x, y, blade_mass)
    global Vrated_func
    Vrated_func = RectBivariateSpline(x, y, Vrated)
    global extremeT_func
    extremeT_func = RectBivariateSpline(x, y, extremeT)


def amaliaWind_23(prob):
    windSpeeds = np.array([6.53163342, 6.11908394, 6.13415514, 6.0614625,  6.21344602,
                                5.87000793, 5.62161519, 5.96779107, 6.33589422, 6.4668016,
                                7.9854581,  7.6894432,  7.5089221,  7.48638098, 7.65764618,
                                6.82414044, 6.36728201, 5.95982999, 6.05942132, 6.1176321,
                                5.50987893, 4.18461796, 4.82863115, 0.,         0.,         0.,
                                5.94115843, 5.94914252, 5.59386528, 6.42332524, 7.67904937,
                                7.89618066, 8.84560463, 8.51601497, 8.40826823, 7.89479475,
                                7.86194762, 7.9242645,  8.56269962, 8.94563889, 9.82636368,
                               10.11153102, 9.71402212, 9.95233636,  10.35446959, 9.67156182,
                                9.62462527, 8.83545158, 8.18011771, 7.9372492,  7.68726143,
                                7.88134508, 7.31394723, 7.01839896, 6.82858346, 7.06213432,
                                7.01949894, 7.00575122, 7.78735165, 7.52836352, 7.21392201,
                                7.4356621,  7.54099962, 7.61335262, 7.90293531, 7.16021596,
                                7.19617087, 7.5593657,  7.03278586, 6.76105501, 6.48004694,
                                6.94716392])

    windFrequencies = np.array([1.17812570e-02, 1.09958570e-02, 9.60626600e-03, 1.21236860e-02,
                               1.04722450e-02, 1.00695140e-02, 9.68687400e-03, 1.00090550e-02,
                               1.03715390e-02, 1.12172280e-02, 1.52249700e-02, 1.56279300e-02,
                               1.57488780e-02, 1.70577560e-02, 1.93535770e-02, 1.41980570e-02,
                               1.20632100e-02, 1.20229000e-02, 1.32111160e-02, 1.74605400e-02,
                               1.72994400e-02, 1.43993790e-02, 7.87436000e-03, 0.00000000e+00,
                               2.01390000e-05, 0.00000000e+00, 3.42360000e-04, 3.56458900e-03,
                               7.18957000e-03, 8.80068000e-03, 1.13583200e-02, 1.41576700e-02,
                               1.66951900e-02, 1.63125500e-02, 1.31709000e-02, 1.09153300e-02,
                               9.48553000e-03, 1.01097900e-02, 1.18819700e-02, 1.26069900e-02,
                               1.58895900e-02, 1.77021600e-02, 2.04208100e-02, 2.27972500e-02,
                               2.95438600e-02, 3.02891700e-02, 2.69861000e-02, 2.21527500e-02,
                               2.12465500e-02, 1.82861400e-02, 1.66147400e-02, 1.90111800e-02,
                               1.90514500e-02, 1.63932050e-02, 1.76215200e-02, 1.65341460e-02,
                               1.44597600e-02, 1.40370300e-02, 1.65745000e-02, 1.56278200e-02,
                               1.53459200e-02, 1.75210100e-02, 1.59702700e-02, 1.51041500e-02,
                               1.45201100e-02, 1.34527800e-02, 1.47819600e-02, 1.33923300e-02,
                               1.10562900e-02, 1.04521380e-02, 1.16201970e-02, 1.10562700e-02])

    # windSpeeds = windSpeeds* 1.714285714 #to get the average speed higher, close to 12 m/s
    nDirections = len(windSpeeds)
    windDirections = np.linspace(0,360-360/nDirections, nDirections)

    index = np.where(windSpeeds==0.0)
    windSpeeds = np.delete(windSpeeds, index[0])
    windFrequencies = np.delete(windFrequencies, index[0])
    windDirections = np.delete(windDirections, index[0])

    wd = np.zeros(23)
    ws = np.zeros(23)
    wf = np.zeros(23)
    for i in range(len(wd)):
        first = 3*i
        second = 3*i+1
        third = 3*i+2
        wd[i] = (windDirections[first]+windDirections[second]+windDirections[third])/3.
        wf[i] = windFrequencies[first]+windFrequencies[second]+windFrequencies[third]
        ws[i] = (windSpeeds[first]*windFrequencies[first]+windSpeeds[second]*windFrequencies[second]+windSpeeds[third]*windFrequencies[third])/wf[i]

    nDirections = len(ws)


    prob['Uref'] = ws
    prob['windDirections'] = wd
    prob['windFrequencies'] = wf

    return nDirections


def amaliaLayout():
    locations = np.loadtxt('/home/flowlab/PJ/FLORISSE3D/doc/layout_amalia.txt')
    turbineX = locations[:, 0]
    turbineY = locations[:, 1]

    return turbineX, turbineY


def alturasRose(nDirections):
    # wf = np.array([247.,280.,220.,153.,40.,92.,90.,91.,128.,181.,173.,116.,78.,62.,53.,49.,
    #    44.,44.,49.,52.,66.,73.,82.,78.,87.,93.,99.,78.,61.,70.,90.,109.,
    #    112.,112.,124.,170.,])
    wf = np.array([78.,99.,93.,87.,78.,82.,73.,66.,52.,49.,44.,44.,49.,53.,62.,78.,
       116.,173.,181.,128.,91.,90.,92.,40.,153.,220.,280.,247.,170.,124.,112.,112.,
       109.,90.,70.,61.])

    wf = wf/sum(wf)

    nDirections = len(wf)
    wd = np.linspace(0.,360.-360/nDirections, nDirections)
    ws = np.ones(nDirections)*8.

    windDirections = wd
    windSpeeds = ws
    windFrequencies = wf

    spline_freq = Akima(windDirections, windFrequencies)
    spline_speed = Akima(windDirections, windSpeeds)
    num = nDirections
    dirs = np.linspace(0.,360.-360./float(num), num)
    ddir = dirs[1]-dirs[0]

    frequencies = np.zeros(num)
    speeds = np.zeros(num)

    num_int = 1000

    dir_int1 = np.linspace(dirs[0],dirs[0]+ddir/2.,num_int/2)
    freq_int1 = np.zeros(num_int/2)
    speed_freq_int1 = np.zeros(num_int/2)
    for j in range(num_int/2):
        freq_int1[j],_,_,_ = spline_freq.interp(dir_int1[j])
        ws,_,_,_ = spline_speed.interp(dir_int1[j])
        speed_freq_int1[j] = freq_int1[j]*ws

    dir_int2 = np.linspace(dirs[0],dirs[0]+ddir/2.,num_int/2)
    freq_int2 = np.zeros(num_int/2)
    speed_freq_int2 = np.zeros(num_int/2)
    for j in range(num_int/2):
        freq_int2[j],_,_,_ = spline_freq.interp(dir_int2[j])
        ws,_,_,_ = spline_speed.interp(dir_int2[j])
        speed_freq_int2[j] = freq_int2[j]*ws

    frequencies[0] = np.trapz(freq_int1,dir_int1)+np.trapz(freq_int2,dir_int2)
    speeds[0] = (np.trapz(speed_freq_int1,dir_int1)+np.trapz(speed_freq_int2,dir_int2))/\
        (np.trapz(freq_int1,dir_int1)+np.trapz(freq_int2,dir_int2))

    for i in range(1,num):
        dir_int = np.linspace(dirs[i]-ddir/2.,dirs[i]+ddir/2.,num_int)
        freq_int = np.zeros(num_int)
        speed_freq_int = np.zeros(num_int)
        for j in range(num_int):
            freq_int[j],_,_,_ = spline_freq.interp(dir_int[j])
            ws,_,_,_ = spline_speed.interp(dir_int[j])
            speed_freq_int[j] = freq_int[j]*ws
        frequencies[i] = np.trapz(freq_int,dir_int)
        speeds[i] = np.trapz(speed_freq_int,dir_int)/np.trapz(freq_int,dir_int)

    frequencies = frequencies/sum(frequencies)
    for i in range(len(frequencies)):
        if speeds[i] < 0.:
            speeds[i] = 0.
        if frequencies[i] < 0.:
            frequencies[i] = 0.

    return dirs, frequencies, speeds



def amaliaRose(nDirections):
    windSpeeds = np.array([6.53163342, 6.11908394, 6.13415514, 6.0614625,  6.21344602,
                                5.87000793, 5.62161519, 5.96779107, 6.33589422, 6.4668016,
                                7.9854581,  7.6894432,  7.5089221,  7.48638098, 7.65764618,
                                6.82414044, 6.36728201, 5.95982999, 6.05942132, 6.1176321,
                                5.50987893, 4.18461796, 4.82863115, 0.,         0.,         0.,
                                5.94115843, 5.94914252, 5.59386528, 6.42332524, 7.67904937,
                                7.89618066, 8.84560463, 8.51601497, 8.40826823, 7.89479475,
                                7.86194762, 7.9242645,  8.56269962, 8.94563889, 9.82636368,
                               10.11153102, 9.71402212, 9.95233636,  10.35446959, 9.67156182,
                                9.62462527, 8.83545158, 8.18011771, 7.9372492,  7.68726143,
                                7.88134508, 7.31394723, 7.01839896, 6.82858346, 7.06213432,
                                7.01949894, 7.00575122, 7.78735165, 7.52836352, 7.21392201,
                                7.4356621,  7.54099962, 7.61335262, 7.90293531, 7.16021596,
                                7.19617087, 7.5593657,  7.03278586, 6.76105501, 6.48004694,
                                6.94716392])

    windFrequencies = np.array([1.17812570e-02, 1.09958570e-02, 9.60626600e-03, 1.21236860e-02,
                               1.04722450e-02, 1.00695140e-02, 9.68687400e-03, 1.00090550e-02,
                               1.03715390e-02, 1.12172280e-02, 1.52249700e-02, 1.56279300e-02,
                               1.57488780e-02, 1.70577560e-02, 1.93535770e-02, 1.41980570e-02,
                               1.20632100e-02, 1.20229000e-02, 1.32111160e-02, 1.74605400e-02,
                               1.72994400e-02, 1.43993790e-02, 7.87436000e-03, 0.00000000e+00,
                               2.01390000e-05, 0.00000000e+00, 3.42360000e-04, 3.56458900e-03,
                               7.18957000e-03, 8.80068000e-03, 1.13583200e-02, 1.41576700e-02,
                               1.66951900e-02, 1.63125500e-02, 1.31709000e-02, 1.09153300e-02,
                               9.48553000e-03, 1.01097900e-02, 1.18819700e-02, 1.26069900e-02,
                               1.58895900e-02, 1.77021600e-02, 2.04208100e-02, 2.27972500e-02,
                               2.95438600e-02, 3.02891700e-02, 2.69861000e-02, 2.21527500e-02,
                               2.12465500e-02, 1.82861400e-02, 1.66147400e-02, 1.90111800e-02,
                               1.90514500e-02, 1.63932050e-02, 1.76215200e-02, 1.65341460e-02,
                               1.44597600e-02, 1.40370300e-02, 1.65745000e-02, 1.56278200e-02,
                               1.53459200e-02, 1.75210100e-02, 1.59702700e-02, 1.51041500e-02,
                               1.45201100e-02, 1.34527800e-02, 1.47819600e-02, 1.33923300e-02,
                               1.10562900e-02, 1.04521380e-02, 1.16201970e-02, 1.10562700e-02])

    windDirections = np.linspace(0.,360.-360./float(len(windSpeeds)), len(windSpeeds))
    spline_freq = Akima(windDirections, windFrequencies)
    spline_speed = Akima(windDirections, windSpeeds)
    num = nDirections
    dirs = np.linspace(0.,360.-360./float(num), num)
    ddir = dirs[1]-dirs[0]

    frequencies = np.zeros(num)
    speeds = np.zeros(num)

    num_int = 1000

    dir_int1 = np.linspace(dirs[0],dirs[0]+ddir/2.,num_int/2)
    freq_int1 = np.zeros(num_int/2)
    speed_freq_int1 = np.zeros(num_int/2)
    for j in range(num_int/2):
        freq_int1[j],_,_,_ = spline_freq.interp(dir_int1[j])
        ws,_,_,_ = spline_speed.interp(dir_int1[j])
        speed_freq_int1[j] = freq_int1[j]*ws

    dir_int2 = np.linspace(dirs[0],dirs[0]+ddir/2.,num_int/2)
    freq_int2 = np.zeros(num_int/2)
    speed_freq_int2 = np.zeros(num_int/2)
    for j in range(num_int/2):
        freq_int2[j],_,_,_ = spline_freq.interp(dir_int2[j])
        ws,_,_,_ = spline_speed.interp(dir_int2[j])
        speed_freq_int2[j] = freq_int2[j]*ws

    frequencies[0] = np.trapz(freq_int1,dir_int1)+np.trapz(freq_int2,dir_int2)
    speeds[0] = (np.trapz(speed_freq_int1,dir_int1)+np.trapz(speed_freq_int2,dir_int2))/\
        (np.trapz(freq_int1,dir_int1)+np.trapz(freq_int2,dir_int2))

    for i in range(1,num):
        dir_int = np.linspace(dirs[i]-ddir/2.,dirs[i]+ddir/2.,num_int)
        freq_int = np.zeros(num_int)
        speed_freq_int = np.zeros(num_int)
        for j in range(num_int):
            freq_int[j],_,_,_ = spline_freq.interp(dir_int[j])
            ws,_,_,_ = spline_speed.interp(dir_int[j])
            speed_freq_int[j] = freq_int[j]*ws
        frequencies[i] = np.trapz(freq_int,dir_int)
        speeds[i] = np.trapz(speed_freq_int,dir_int)/np.trapz(freq_int,dir_int)

    frequencies = frequencies/sum(frequencies)
    for i in range(len(frequencies)):
        if speeds[i] < 0.:
            speeds[i] = 0.
        if frequencies[i] < 0.:
            frequencies[i] = 0.

    return dirs, frequencies, speeds


def setup_weibull(windDirections,windFrequencies,windSpeeds,nSpeeds):

    def Weibull(x,L):
        k = 1.76
        if L < 0.0001:
            L = 0.0001
        return (k/L)*(x/L)**(k-1)*np.exp(-(x/L)**k)

    nDirections = len(windDirections)
    dirs = np.zeros(nDirections*nSpeeds)
    freqs = np.zeros(nDirections*nSpeeds)
    speeds = np.zeros(nDirections*nSpeeds)

    #direction loops
    for i in range(nDirections):
        for j in range(nSpeeds):
            dirs[i*nSpeeds+j] = windDirections[i]

    #speed and frequency loops
    for i in range(nDirections):
        avg_speed = windSpeeds[i]
        speed_dist = np.linspace((25.)/(2.*float(nSpeeds))+0.001,25.-(25.)/(2.*float(nSpeeds)),nSpeeds)
        # print 'speed_dist: ', speed_dist
        dspeed = speed_dist[1]-speed_dist[0]
        num_int = 1000
        for j in range(nSpeeds):
            speed_int = np.linspace(speed_dist[j]-dspeed/2.,speed_dist[j]+dspeed/2.,num_int)
            freq_int = Weibull(speed_int,avg_speed)
            speed_freq = np.trapz(freq_int,speed_int)
            speeds[i*nSpeeds+j] = speed_dist[j]
            freqs[i*nSpeeds+j] = speed_freq*windFrequencies[i]

    return dirs, freqs, speeds

if __name__=="__main__":
    windDirections = np.array([0.,90.,180.,270.])
    windSpeeds = np.array([10.,12.,5.,5.8])
    windFrequencies = np.array([0.2,0.5,0.1,0.2])
    nSpeeds = 5

    dirs, freqs, speeds = setup_weibull(windDirections,windFrequencies,windSpeeds,nSpeeds)
    print 'dirs: ', dirs
    print 'freqs: ', freqs
    print 'speeds: ', speeds
    print sum(freqs)
