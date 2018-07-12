import numpy as np
from openmdao.api import Group, Component

class transportationCost(Component):
    def __init__(self, nTurbines):

        super(transportationCost, self).__init__()

        self.nTurbines = nTurbines
        self.add_param('transportDist', 0.0, desc='transportation distance')
        self.add_param('cost', 0.0, desc='TCC (whole farm)')

        self.add_output('transportation_cost', 0.0, desc='transportation cost')

    def solve_nonlinear(self, params, unknowns, resids):
        nTurbs = float(self.nTurbines)
        turbineCost = params['cost']

        transportation_cost = turbineCost

        transportation_cost += 1867.*params['transportDist']**(0.726)*nTurbs

        unknowns['transportation_cost'] = transportation_cost

    def linearize(self, params, unknowns, resids):
        nTurbs = float(self.nTurbines)

        J = {}
        J['transportation_cost', 'cost'] = 1.
        return J

"""
"""
class powerPerformanceCost(Component):
    #TODO what to do for gradients on this?
    def __init__(self, nTurbines):

        super(powerPerformanceCost, self).__init__()

        self.nTurbines = nTurbines
        self.add_param('turbineZ', np.zeros(nTurbines), units='m', desc='the hub heights of each turbine')

        self.add_output('power_performance_cost', 0.0, desc='power performance cost')

    def solve_nonlinear(self, params, unknowns, resids):

        nTurbs = float(self.nTurbines)
        avgHeight = np.sum(params['turbineZ'])/nTurbs
        # avgHeight = 110.

        cost = np.zeros(int(nTurbs))

        for i in range(int(nTurbs)):

            hL = 85.0
            hU = 95.0

            c3 = -114.8
            c2 = 30996.0
            c1 = -2781030.0
            c0 = 83175600.0

            mL1 = 232600.0
            mU1 = 290000.0

            if avgHeight <= hL:
                multiplier1 = mL1
            elif avgHeight >= hU:
                multiplier1 = mU1
            else:
                multiplier1 = c3*avgHeight**3 + c2*avgHeight**2 + c1*avgHeight + c0

            c3 = -48.4
            c2 = 13068.0
            c1 = -1172490.0
            c0 = 35061600.0

            mL2 = 92600.
            mU2 = 116800.

            if avgHeight <= hL:
                multiplier2 = mL2
            elif avgHeight >= hU:
                multiplier2 = mU2
            else:
                multiplier2 = c3*avgHeight**3 + c2*avgHeight**2 + c1*avgHeight + c0

            permanent = 2. #number of permanent met towers
            temporary = 2.

            cost[i] = 200000. + permanent*multiplier1 + temporary*multiplier2

        power_perf_cost = np.sum(cost)/nTurbs

        unknowns['power_performance_cost'] = power_perf_cost

    def linearize(self, params, unknowns, resids):
        nTurbs = self.nTurbines

        J = {}
        J['power_performance_cost', 'turbineZ'] = np.zeros([1,nTurbs])
        return J


class accessRoadCost(Component):
    def __init__(self, nTurbines):

        super(accessRoadCost, self).__init__()

        self.nTurbines = nTurbines
        self.add_param('rotorDiameter', np.ones(nTurbines)*126.4, desc='rotor diameter')

        self.add_output('access_road_cost', 0.0, desc='access road cost')

    def solve_nonlinear(self, params, unknowns, resids):
        nTurbs = float(self.nTurbines)

        factor1 = 62653.6 #complex layout, flat terrain
        factor2 = 30.9
        rotorDiameter = params['rotorDiameter']
        constructionTime = 5.
        accessRoadEntrances = 1.

        roads_cost = 0.0
        for i in range(int(nTurbs)):
            roads_cost += factor1 + rotorDiameter[i]*factor2

        roads_cost = roads_cost*1.05
        roads_cost += 1.05*(constructionTime*55500. + accessRoadEntrances*3800.)

        unknowns['access_road_cost'] = roads_cost

    def linearize(self, params, unknowns, resids):
        factor2 = 30.9
        J = {}
        J['access_road_cost','rotorDiameter'] = np.ones((1,self.nTurbines))*factor2*1.05

        return J


class foundationCost(Component):
    def __init__(self, nTurbines):

        super(foundationCost, self).__init__()

        self.nTurbines = nTurbines
        self.add_param('rotorDiameter', np.ones(nTurbines)*126.4, desc='rotor diameter')
        self.add_param('turbineZ', np.zeros(nTurbines), units='m', desc='the hub heights of each turbine')
        self.add_param('ratedPower', np.ones(nTurbines)*5000., desc='rated power array')

        self.add_output('foundation_cost', 0.0, desc='foundation cost')

    def solve_nonlinear(self, params, unknowns, resids):
        nTurbs = float(self.nTurbines)
        ratedPower = params['ratedPower']
        turbineZ = params['turbineZ']
        rotorDiameter = params['rotorDiameter']

        # topMass = 24.05066 #??
        # #topMass = 88. #??

        topMass = np.zeros(int(nTurbs))
        for i in range(int(nTurbs)):
            topMass[i] = 50. #TODO don't really know what this is: need to fix for sure
        self.topMass = topMass

        cost = 0.
        for i in range(int(nTurbs)):
            cost += ratedPower[i]*rotorDiameter[i]*topMass[i]/1000. + 163421.5*nTurbs**(-0.1458) + (turbineZ[i]-80.)*500.

        foundation_cost = np.sum(cost)

        unknowns['foundation_cost'] = foundation_cost

    def linearize(self, params, unknowns, resids):
        nTurbs = self.nTurbines
        topMass = self.topMass

        J = {}
        J['foundation_cost', 'turbineZ'] = np.ones((1,nTurbs))*500.
        J['foundation_cost', 'rotorDiameter'] = np.zeros((1,nTurbs))
        J['foundation_cost', 'ratedPower'] = np.zeros((1,nTurbs))


        for i in range(nTurbs):
            J['foundation_cost', 'rotorDiameter'][0][i] = params['ratedPower'][i]*topMass[i]/1000.
            J['foundation_cost', 'ratedPower'][0][i] = params['rotorDiameter'][i]*topMass[i]/1000.
        return J


class erectionCost(Component):
    def __init__(self, nTurbines):

        super(erectionCost, self).__init__()

        self.nTurbines = nTurbines
        self.add_param('turbineZ', np.zeros(nTurbines), units='m', desc='the hub heights of each turbine')
        self.add_param('ratedPower', np.ones(nTurbines)*5000., desc='rated power array')

        self.add_output('erection_cost', 0.0, desc='erection cost')

    def solve_nonlinear(self, params, unknowns, resids):
        nTurbs = float(self.nTurbines)
        turbineZ = params['turbineZ']
        ratedPower = params['ratedPower']

        weatherDelayDays = 5.
        craneBreakdowns = 1.

        cost = np.zeros(int(nTurbs))
        for i in range(int(nTurbs)):
            cost[i] = (37.*ratedPower[i] + 27000.*nTurbs**(-0.42145) + (turbineZ[i]-80.)*500.)

        erection_cost = np.sum(cost)+ 20000.*weatherDelayDays + 35000.*craneBreakdowns + 181.*nTurbs + 1834.

        unknowns['erection_cost'] = erection_cost

    def linearize(self, params, unknowns, resids):
        nTurbs = self.nTurbines

        J = {}
        J['erection_cost', 'turbineZ'] = np.ones([1, nTurbs])*500.
        J['erection_cost', 'ratedPower'] = np.ones([1, nTurbs])*37.

        return J


class electircalMaterialsCost(Component):
    def __init__(self, nTurbines):

        super(electircalMaterialsCost, self).__init__()

        self.nTurbines = nTurbines
        self.add_param('rotorDiameter', np.ones(nTurbines)*126.4, desc='rotor diameter')

        self.add_output('electrical_materials_cost', 0.0, desc='electrical materials cost')

    def solve_nonlinear(self, params, unknowns, resids):
        nTurbs = float(self.nTurbines)

        factor1 = 67519.4 #complex layout, flat terrain
        factor2 = 27874.4
        factor3 = 681.7
        thermalBackfill = 2. #TODO (what is this)

        elec_mat_cost = nTurbs*factor2

        rotor_diameter = np.sum(params['rotorDiameter'])/nTurbs #avg rotor diameter

        #TODO (area)
        elec_mat_cost += 5.*35375. + 1.*50000. + rotor_diameter*nTurbs*factor3 + \
                    thermalBackfill*5. + 41945.

        unknowns['electrical_materials_cost'] = elec_mat_cost

    def linearize(self, params, unknowns, resids):
        factor3 = 681.7
        J = {}
        J['electrical_materials_cost','rotorDiameter'] = np.ones((1,self.nTurbines))*factor3

        return J


class electircalInstallationCost(Component):
    def __init__(self, nTurbines):

        super(electircalInstallationCost, self).__init__()

        self.nTurbines = nTurbines
        self.add_param('rotorDiameter', np.ones(nTurbines)*126.4, desc='rotor diameter')

        self.add_output('electrical_installation_cost', 0.0, desc='electrical installation cost')

    def solve_nonlinear(self, params, unknowns, resids):
        nTurbs = float(self.nTurbines)

        factor1 = 7683.5
        factor2 = 564.9
        factor3 = 446.0
        rockTrenchingLength = 10.

        elec_instal_cost = 5.*14985.+155000.
        rotor_diameter = np.sum(params['rotorDiameter'])/nTurbs #avg rotor diameter

        #TODO area
        elec_instal_cost += nTurbs*(factor1 + rotor_diameter*(factor2 + \
                factor3*rockTrenchingLength/100.0)) + 10000.

        unknowns['electrical_installation_cost'] = elec_instal_cost

    def linearize(self, params, unknowns, resids):
        factor2 = 564.9
        factor3 = 446.0
        rockTrenchingLength = 10.
        J = {}
        J['electrical_installation_cost','rotorDiameter'] = np.ones((1,self.nTurbines))*\
                factor2 + factor3*rockTrenchingLength/100.

        return J


class insuranceCost(Component):
    def __init__(self, nTurbines):

        super(insuranceCost, self).__init__()

        self.nTurbines = nTurbines
        self.add_param('cost', 0.0, desc='TCC (whole farm)')
        self.add_param('foundation_cost', 0.0, desc='foundation costs')
        self.add_param('ratedPower', np.ones(nTurbines)*5000., desc='rated power array')

        self.add_output('insurance_cost', 0.0, desc='insurance cost')
        self.add_output('alpha_insurance', 0.0, desc='alpha insurance')

    def solve_nonlinear(self, params, unknowns, resids):
        nTurbs = float(self.nTurbines)

        turbineCost = params['cost']/(np.sum(params['ratedPower']))

        alpha_insurance = 3.5 + 0.7 + 0.4 + 1.0
        insurance_cost = (0.7 + 0.4 + 1.0) * turbineCost * 37.5

        alpha_insurance /= 1000.0
        insurance_cost += 0.02*params['foundation_cost'] + 20000.

        unknowns['insurance_cost'] = insurance_cost
        unknowns['alpha_insurance'] = alpha_insurance

    def linearize(self, params, unknowns, resids):
        nTurbs = float(self.nTurbines)

        J = {}
        J['insurance_cost', 'cost'] = (0.7+0.4+1.0)*37.5/(np.sum(params['ratedPower']))
        J['insurance_cost', 'foundation_cost'] = 0.02
        J['insurance_cost', 'ratedPower'] = np.ones((1,int(nTurbs)))*-1.*(0.7+0.4+1.0)*37.58*params['cost']/(np.sum(params['ratedPower'])**2)

        return J


class markupCost(Component):
    def __init__(self, nTurbines):

        super(markupCost, self).__init__()

        self.nTurbines = nTurbines
        self.add_param('transportation_cost', 0.0, desc='transportation costs')

        self.add_output('markup_cost', 0.0, desc='markup cost')
        self.add_output('alpha_markup', 0.0, desc='alpha markup')

    def solve_nonlinear(self, params, unknowns, resids):
        nTurbs = float(self.nTurbines)

        contingency = 3.0
        warranty = 0.02
        useTax = 0.0
        overhead = 5.0
        profitMargin = 5.0

        alpha_markup = (contingency + warranty + useTax + overhead + profitMargin)/100.0
        markup_cost = -alpha_markup * params['transportation_cost']

        unknowns['markup_cost'] = markup_cost
        unknowns['alpha_markup'] = alpha_markup

    def linearize(self, params, unknowns, resids):
        alpha_markup = unknowns['alpha_markup']

        J = {}
        J['markup_cost', 'transportation_cost'] = -alpha_markup

        return J


class BOScalc(Component):
    def __init__(self, nTurbines):

        super(BOScalc, self).__init__()

        self.nTurbines = nTurbines
        self.add_param('transportation_cost', 0.0, desc='transportation cost')
        self.add_param('power_performance_cost', 0.0, desc='power performance cost')
        self.add_param('access_road_cost', 0.0, desc='access road cost')
        self.add_param('foundation_cost', 0.0, desc='foundations cost')
        self.add_param('erection_cost', 0.0, desc='erection cost')
        self.add_param('electrical_materials_cost', 0.0, desc='electrical materials cost')
        self.add_param('electrical_installation_cost', 0.0, desc='electrical installation cost')
        self.add_param('insurance_cost', 0.0, desc='insurance cost')
        self.add_param('markup_cost', 0.0, desc='markup cost')
        self.add_param('alpha_insurance', 0.0, desc='alpha insurance')
        self.add_param('alpha_markup', 0.0, desc='alpha markup')
        self.add_param('cost', 0.0, desc='TCC (whole farm)')

        self.add_output('BOS', 0.0, desc='BOS costs')

    def solve_nonlinear(self, params, unknowns, resids):
        nTurbs = float(self.nTurbines)
        #constant = 16867866.8985 #BOS costs that remain constant
        constant = 18507985.8237
        total_cost = constant+params['transportation_cost']+params['power_performance_cost']+\
                    params['access_road_cost']+params['foundation_cost']+params['erection_cost']+\
                    params['electrical_materials_cost']+params['electrical_installation_cost']+\
                    params['insurance_cost']+params['markup_cost']
        self.total_cost = total_cost

        # print 'total: ', total_cost
        # print 'access_road_cost: ', params['access_road_cost']
        # print 'electrical_materials_cost: ', params['electrical_materials_cost']
        # print 'electrical_installation_cost: ', params['electrical_installation_cost']
	# print 'foundation: ', params['foundation_cost']

        alpha = params['alpha_markup'] + params['alpha_insurance']
        self.alpha = alpha

        #multiplier
        total_cost /= (1.0-alpha)

        #remove TCC
        total_cost -= params['cost']

        unknowns['BOS'] = total_cost

    def linearize(self, params, unknowns, resids):
        nTurbs = float(self.nTurbines)
        alpha = self.alpha

        J = {}
        J['BOS', 'transportation_cost'] = 1./(1.-alpha)
        J['BOS', 'power_performance_cost'] = 1./(1.-alpha)
        J['BOS', 'access_road_cost'] = 1./(1.-alpha)
        J['BOS', 'foundation_cost'] = 1./(1.-alpha)
        J['BOS', 'erection_cost'] = 1./(1.-alpha)
        J['BOS', 'electrical_materials_cost'] = 1./(1.-alpha)
        J['BOS', 'electrical_installation_cost'] = 1./(1.-alpha)
        J['BOS', 'insurance_cost'] = 1./(1.-alpha)
        J['BOS', 'markup_cost'] = 1./(1.-alpha)
        # J['BOS', 'alpha_insurance'] =
        # J['BOS', 'alpha_markup'] =
        J['BOS', 'cost'] = -1.
        return J


class BOSgroup(Group):
    """
    Group containing components of BOS
    """
    def __init__(self, nTurbines):

        super(BOSgroup, self).__init__()

        self.add('transportationCost', transportationCost(nTurbines), promotes=['*'])
        self.add('powerPerformanceCost', powerPerformanceCost(nTurbines), promotes=['*'])
        self.add('accessRoadCost', accessRoadCost(nTurbines), promotes=['*'])
        self.add('foundationCost', foundationCost(nTurbines), promotes=['*'])
        self.add('erectionCost', erectionCost(nTurbines), promotes=['*'])
        self.add('electircalMaterialsCost', electircalMaterialsCost(nTurbines), promotes=['*'])
        self.add('electircalInstallationCost', electircalInstallationCost(nTurbines), promotes=['*'])
        self.add('insuranceCost', insuranceCost(nTurbines), promotes=['*'])
        self.add('markupCost', markupCost(nTurbines), promotes=['*'])
        self.add('BOScalc', BOScalc(nTurbines), promotes=['*'])
