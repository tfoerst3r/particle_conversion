import numpy as np
import cantera as ct

class DiffusionCoeff():
    """
    pass
    """
    def __init__(self,mechanism):

        self._gas_constant  = ct.gas_constant/1000 # 8.3144621 J/mol-K
        self._avogadro = ct.avogadro/1000 # (1/mol) / (1/mol)
        self._gas = ct.Solution(mechanism)
        self._oxidizer = 'O2'


    def chapmanAB(self,temp,pressure,gasA,gasB):

        def omega_AB(temp,gasA,gasB):

            # poling 2001, B.1, epsilon/k in K
            coeff = [
                1.06036,
                0.15610,
                0.19300,
                0.47635,
                1.03587,
                1.52996,
                1.76474,
                3.89411
            ]

            # poling 2001, B.1, epsilon/k in K
            epsilondic = {
                'CO':  91.7,
                'CO2': 195.2,
                'H2O': 809.1,
                'O2':  106.7,
                'N2':  71.4,
                'CH4': 148.6,
                'H2':  59.7,
            }

            temp2 = temp / ( np.sqrt( epsilondic[gasA] * epsilondic[gasB] ) )

            value =   coeff[0] / temp2**coeff[1] \
                    + coeff[2] / np.exp(coeff[3] * temp2 )  \
                    + coeff[4] / np.exp(coeff[5] * temp2 )  \
                    + coeff[6] / np.exp(coeff[7] * temp2 )

            return value

        sigmadic = {
            'CO':  3.690e-10,
            'CO2': 3.941e-10,
            'H2O': 2.641e-10,
            'O2':  3.467e-10,
            'N2':  3.798e-10,
            'CH4': 3.758e-10,
            'H2':  2.827e-10
        }

        mw_gasA = self._gas.molecular_weights[self._gas.species_index(gasA)]/1000
        mw_gasB = self._gas.molecular_weights[self._gas.species_index(gasB)]/1000

        sigmaAB = (sigmadic[gasA] + sigmadic[gasB])/2

        #-- bird 2002, (J/mol)**(1/2) * J/K**(3/2)
        constantA = 3/8 * np.sqrt( self._gas_constant**3/(self._avogadro**2 * 2*np.pi ))

        DiffusionCoeff = constantA * \
             temp**(3/2) / (pressure * omega_AB(temp,gasA,gasB) * sigmaAB**2 ) * \
             np.sqrt( 1/mw_gasA + 1/mw_gasB )

        return DiffusionCoeff


    def diffusion_mix(self,comp,temp,pressure):
        """
        @comp: gas composition in mole fractions, first is the diffusing species
        :return:
        """
        elements = []
        value = []
        for i,val in enumerate(comp):
            elements.append(val)
            value.append(comp[val])

        self._oxidizer = elements[0]
        #print('Diffusion species:',self._oxidizer)

        #-- chapman is the diffusion coefficient

        for i,val in enumerate(elements):

            if i == 0:
                dmix = 0
            else:
                diff_coeffAB = self.chapmanAB(temp, pressure, elements[0], elements[i])
                dmix += value[i] / diff_coeffAB

        dmix = (1-value[0])/dmix

        return dmix


#    def diff_limit_const(self,comp,temp0,pressure0,pressure1):
#
#
#        diff0 = self.diffusion_mix(comp,temp0,pressure0)
#
#        M_C  = self._gas.molecular_weights[self._gas.species_index('C')]/1000
#
#        #------------------------#
#        #-- determine nu_C/nu_i
#        #-- it is based on the species, which are considered C + 0.5 O2 -> CO
#        if self._oxidizer == 'O2':
#            print('Used species is:',self._oxidizer)
#            stoichio = 0.5
#        elif self._oxidizer == 'CO2':
#            stoichio = 1
#        elif self._oxidizer == 'H2O':
#            stoichio = rrr.stoichiometry(fuel={'C':1},species=['H2','CO','H2O','N2','SO2'], mw_fuel=M_C)[1]
#        else:
#            raise IOError('Only C+O2 C+CO2 or C+H2O is implemented!')
#
#        for i,val in enumerate(stoichio):
#            if (val != 'fuel' or val != 'C') and stoichio[val] < 0:
#                nu_i = stoichio[val]
#            if (val == 'fuel' or val == 'C') and stoichio[val] < 0:
#                nu_C = stoichio[val]
#        #------------------------#
#
#        print('nu_C/nu_i ratio:',round(nu_C/nu_i,2))
#
#        Sh = 2
#
#        effectivness = 1
#
#        C1 = effectivness * (nu_C/nu_i) * M_C * (Sh * diff0) / (self._gas_constant * temp0**(1.75)) * pressure0/pressure1
#
#        return C1



