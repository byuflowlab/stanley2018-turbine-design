import numpy as np
import os


if __name__=="__main__":
    ratedPower = np.linspace(500.,10000.,20)
    rotorDiameter = np.linspace(46.,160.,20)
    number = np.linspace(0.,100.,101)

    BEST_DATA = np.zeros((400,10))

    rated_power = 0.
    rotor_diameter = 0.
    ratedQ = 0.
    blade_mass = 0.
    Vrated = 0.
    I1 = 0.
    I2 = 0.
    I3 = 0.
    ratedT = 0.
    extremeT = 0.

    # for power in range(20):
    #     for diameter in range(20):
    #         for num in range(101):
    for power in range(20):
        for diameter in range(20):

            OPT_ratedQ = 0.
            OPT_blade_mass = 1.E10
            OPT_Vrated = 0.
            OPT_I1 = 0.
            OPT_I2 = 0.
            OPT_I3 = 0.
            OPT_ratedT = 0.
            OPT_extremeT = 0.

            for num in range(100):
                filename = 'src/florisse3D/optRotor/passed/results_%s_%s_%s.txt'%(ratedPower[power],rotorDiameter[diameter],number[num])
                if os.path.isfile(filename):
                    print filename
                    opened = open(filename)
                    data = np.loadtxt(opened)
                    rated_power = ratedPower[power]
                    rotor_diameter = rotorDiameter[diameter]
                    ratedQ = data[0]
                    blade_mass = data[1]
                    Vrated = data[2]
                    I1 = data[3]
                    I2 = data[4]
                    I3 = data[5]
                    ratedT = data[6]
                    extremeT = data[7]

                    if ratedQ > 0. and ratedT > 0. and blade_mass < OPT_blade_mass:
                        OPT_ratedQ = ratedQ
                        OPT_blade_mass = blade_mass
                        OPT_Vrated = Vrated
                        OPT_I1 = I1
                        OPT_I2 = I2
                        OPT_I3 = I3
                        OPT_ratedT = ratedT
                        OPT_extremeT = extremeT
                        # print rated_power, rotor_diameter, OPT_ratedQ, OPT_blade_mass, OPT_Vrated, OPT_I1, OPT_I2, OPT_I3, OPT_ratedT, OPT_extremeT

            BEST_DATA[20*power+diameter][0] = rated_power
            BEST_DATA[20*power+diameter][1] = rotor_diameter
            BEST_DATA[20*power+diameter][2] = OPT_ratedQ
            BEST_DATA[20*power+diameter][3] = OPT_blade_mass
            BEST_DATA[20*power+diameter][4] = OPT_Vrated
            BEST_DATA[20*power+diameter][5] = OPT_I1
            BEST_DATA[20*power+diameter][6] = OPT_I2
            BEST_DATA[20*power+diameter][7] = OPT_I3
            BEST_DATA[20*power+diameter][8] = OPT_ratedT
            BEST_DATA[20*power+diameter][9] = OPT_extremeT

    np.savetxt('src/florisse3D/optRotor/BEST_DATA.txt', BEST_DATA)
