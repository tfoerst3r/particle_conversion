# SPDX-FileCopyrightText: 2023 Thomas Förster
#
# SPDX-License-Identifier: MIT

__author__ = 'tf'

import numpy as np
import cantera as ct
from pandas import concat, Series, DataFrame
from scipy.optimize import bisect
from math import pow 

from particle_conversion.models.diffusion import DiffusionCoeff
from particle_conversion.default_parameters import _defaults
from particle_conversion.utils import func_inter, normalize_dict


## CLASS DEFINITION ##
class Particle_class():
    def __init__(self, 
            settings,
            gas_header,
            gas_flags,
            gas_comp,
            others_header,
            others_flags,
            others,
            ptot,
    ):

        # PropDict init
        PropDict = settings['reaction']

        # particle properties
        PropDict['d_p0']          = float(settings['properties']['particle_diameter_init'])
        PropDict['rhop0']         = float(settings['properties']['density_particle_apparent_init'])
        PropDict['rhoa0']         = float(settings['properties']['density_ash_true'])
        self.rho_char_init_true   = \
            PropDict['rho_char_true'] = \
                settings['properties']['density_char_true'] if 'density_char_true' in settings['properties'] else _defaults['particle']['properties']['density_char_true']
        PropDict['sint0']         = 1 # deprecated should be A = s0 * A*; makes it more feasible
        PropDict['mechanism']     = _defaults['numerical']['mechanism']
        
        # proxymate analysis
        PropDict['y_fc0']  = settings['proxymate_analysis']['fixedcarbon_init']
        PropDict['y_ash0'] = settings['proxymate_analysis']['ash_init']
        PropDict['Sdevel'] = settings['conversion_model']['model']

        # model properties
        PropDict['tauf'] = PropDict['tortuosity-factor']


        #.. gas phase properties
        self.gas = ct.Solution(PropDict['mechanism'])
        self.gas.transport_model = 'Multi'

        self._gas_const  = ct.gas_constant/1000   # gas constant 8.3144 J/mol.K
        self.gas_comp    = gas_comp
        self.gas_header  = gas_header
        self.gas_flags   = gas_flags
        self.gas_default = ''
        self.gas_comp_dict = { 'CO':0, 'CO2':0, 'H2O':0, 'O2':0, 'N2':0, 'CH4':0, 'H2':0 }
        self.diffC = DiffusionCoeff(PropDict['mechanism'])   # init class definition
        self.gas_default_diff_partner = { 'CO2':{'CO':0},'O2':{'CO2':0},'H2O':{'CO':0.0,'H2':0.0} }

        'others, Re,Pr,Tg,Tp'
        self.others        = others
        self.others_header = others_header
        self.others_flags  = others_flags

        'pressures in the reactor, Pa'
        self.ptot = ptot # total pressure
        self.pg   = ptot # init value
        self.ps   = ptot # init value

        'particle physical properties'
        self.dp0 = float(PropDict['d_p0'])  # diameter(time=0), m

        'particle thermo properties'
        self.Tp   = 300                     # init temperature, K

        'diffusion constant coeff (only output)'
        self.C1 = {'O2':0.0,'CO2':0.0,'H2O':0.0 }    # init values

        'Surface development model'
        self.Sdevelmodel = PropDict['Sdevel']  # RPM, SPM, SDM

        #--------------------#
        # Proximate Analysis #
        #---------------------------------------------------------------#
        # -- Proxi, kg/kg:
        # -- Output: [ y_wet0, y_vola0, y_Cfix0, y_ash0 ]
        self.proxi = {}

        if 'y_fc0' in PropDict.keys():
            self.proxi['fc'] = float(PropDict['y_fc0'])
        else:
            self.proxi['fc'] = 1.0
        # --
        if 'y_ash0' in PropDict.keys():
            self.proxi['ash'] = float(PropDict['y_ash0'])
        else:
            self.proxi['ash'] = 0.0

        # -- Fixed Carbon Analysis
        self.proxi = normalize_dict(self.proxi)

        self.yc0 = self.proxi['fc']
        self.ya0 = self.proxi['ash']
        #---------------------------------------------------------------#

        # densities, kg/m3
        self.rhop0 = float(PropDict['rhop0'])       # apparent particle density, kg/m3
        self.rhoa0 = float(PropDict['rhoa0'])       # true density ash, kg/m3
        
        # apparent char(daf) density
        self.rhoc0 = (1-self.ya0) / (1/self.rhop0 - self.ya0/self.rhoa0)

        # constants for the apparent density and radius change
        if self.Sdevelmodel == 'RPM':
            self.alpha = 1
        elif self.Sdevelmodel == 'SPM' or self.Sdevelmodel == 'SPMp':
            self.alpha = 0
        elif self.Sdevelmodel == 'SDM' or self.Sdevelmodel == 'SDMp':
            self.alpha = 1
        #else:
        #    self.alpha = float(PropDict['alpha'])

        self.beta  = 1/3 * (1-self.alpha)

        'structural parameter of RPM, -'
        self.etaF_mod = {
            'O2':  1.0,
            'CO2': 1.0,
            'H2O': 1.0
        }

        'particle masses, kg'
        self.mp0 = 1/6*np.pi * pow( self.dp0,3) * self.rhop0
        #-- carbon content
        self.mc0 = self.yc0 * self.mp0
        #-- ash content, constant over time
        self.ma0 = self.ya0 * self.mp0

        'yield'
        self.X = 0.0            # initial value

        'specific intrinsic surface, m2/kg'
        ## will be multiplied with the frequence
        ## factor and will not be used further
        self._sint0 = float(PropDict['sint0'])

        'structural parameter of RPM, -'
        self.psi0 = {
            'O2':  float(PropDict['psi_R1']),
            'CO2': float(PropDict['psi_R2']),
            'H2O': float(PropDict['psi_R3'])
        }
        self.psi1 = 1  # init value

        self.ratep = {'O2':0,'CO2':0,'H2O':0}

        #-------------------#
        # Reaction kinetics #
        #---------------------------------------------------------------#
        #-- [k0] = [A] = kg/(s.Pa^n.m^2)  |  [Ea] = J/kg  |  [n] = -
        #-- BE AWARE THAT sint0 is multiplied here with kin[0] and will not be further used in the code!
        'char kinetic parameters'
        self.kin = [ float(PropDict['A']), float(PropDict['Ea']), float(PropDict['n'])]
        self.kin_R1_available = self.kin_R2_available = self.kin_R3_available = False

        #-- REACTION O2 --#
        if ('A_R1' in PropDict.keys()):
            self.kin_R1 = [ float(PropDict['A_R1'])*self._sint0 , float(PropDict['Ea_R1']), float(PropDict['n_R1'])] # kg/(s.Pa^n.m^2)  |  J/kg  |  Pa
        else:
            raise IOError('no specific reaction coefficients or reaction enthalpy are available for C + 0.5 O2 -> CO ')

        #-- REACTION CO2 --#
        if ('A_R2' in PropDict.keys()):
            self.kin_R2 = [ float(PropDict['A_R2'])*self._sint0 , float(PropDict['Ea_R2']), float(PropDict['n_R2'])] # kg/(s.Pa^n.m^2)  |  J/kg  |  Pa
        else:
            raise IOError('no specific reaction coefficients are available for C + CO2 -> 2 CO ')

        #-- REACTION H2O --#
        if ('A_R3' in PropDict.keys()):
            self.kin_R3 = [ float(PropDict['A_R3'])*self._sint0 , float(PropDict['Ea_R3']), float(PropDict['n_R3'])] # kg/(s.Pa^n.m^2)  |  J/kg  |  Pa
        else:
            raise IOError('no specific reaction coefficients are available for C + H2O -> CO + H2 ')

        #-- Other parameters
        'stoichio metric factor'
        self.stoichio=1     # init value

        'structural char pore parameter = tortuosity*constriction-factor'
        self.tauf = float(PropDict['tauf'])

        #'porosity' #  if assumed as constant value ( transient calculation is possible with self.porosity()
        #self.eps1 = float(PropDict['eps1'])

        'molar weight, carbon'
        self.Mc = 12.011/1000  # kg/mol 0.012011

        self.effF         = {'O2':0,'CO2':0,'H2O':0}   #-- init

        #self._adaption_factor = float(PropDict['adaptionfactor'])

        #-- Init output dataframe --#
        header_time = 'Time'
        header = [
            header_time,    # 01
            'mp',           # 02
            'dmpdt',        # 03
            'X',            # 04
            'dXdt',         # 05
            'dXdt_O2',      # 06
            'dXdt_CO2',     # 07
            'dXdt_H2O',     # 08
            'C1_O2',        # 09
            'C1_CO2',       # 10
            'C1_H2O',       # 11
            'effF_O2',      # 12
            'effF_CO2',     # 13
            'effF_H2O',     # 14
            'Sh',           # 15
            'Bscaling',     # 16
            'porosity',     # 17
            'Tp',           # 18
            'rhop',         # 19
            'dp',           # 20
            'ps',           # 21
        ]

        self.result = DataFrame(columns=header)
        #self.result.set_index(header_time, inplace=True)

        #################
        ## END -- Init ##
        #################


    ###################
    ## Class Methods ##
    ###################
    
    def conv(self,t,y):
        '''
        This function will be solved by the ODE proceedure, conv..conversion
        '''

        X = y  # Carbon conversion, -

        if X >= 1.0:
            self.X = 0.99999999
        else:
            self.X = X

        #--------------------------------------------#
        #-- (A) gas composition at a specific time --#
        for i,val in enumerate(self.gas_header):
                self.gas_comp_dict.update( { val: func_inter(t,self.gas_comp[self.gas_header[val]]) } )

        #-- in case the read data is not normalized!
        self.gas_comp_dict = normalize_dict(self.gas_comp_dict)
        
        #--------------------------------------------#
        #-- (B) check for default conditions --------#
        if self.gas_comp_dict['CO2'] == 1:
            self.gas_default = 'CO2'
        elif self.gas_comp_dict['H2O'] == 1:
            self.gas_default = 'H2O'
        elif self.gas_comp_dict['O2'] == 1:
            self.gas_default = 'O2'


        #--------------------------------------------#
        #-- (C) generate other properties -----------#
        self.Re = func_inter(t,self.others[self.others_header['Re']])  # Reynolds number, -
        self.Tg = func_inter(t,self.others[self.others_header['Tg']])  # Gas temperature, K
        self.Tp = func_inter(t,self.others[self.others_header['Tp']])  # Gas temperature, K
        self.gas.TPX = self.Tg,self.ptot,self.gas_comp_dict
        
        #-- Prandtl number, -
        #[Pr] = [-] = kg/m.s * J/kg.K * s.m.K/J
        #     = [-] = (m2 kg J m K) / (s m3 kg K W) =  m2/s * kg/m3 * J/(kg.K) / W/(m.K)
        #       gas.viscosity, Pa.s (dyn. visco not kin. visco)
        self.Pr = self.gas.viscosity * self.gas.cp/self.gas.thermal_conductivity


        #--------------------------------------------------------------------------#
        #-- (D) reaction rate, 1/s, R ... dXdt = dX1/dt + dX2/dt + dX3/dt ---------#
        self.ratep['O2'],self.ratep['CO2'],self.ratep['H2O'] = self.rate_3Rkt()
        self.dXdt = self.ratep['O2'] + self.ratep['CO2'] + self.ratep['H2O']

        
        #--------------------------------------------#
        #-- (E) output ------------------------------#
        #.. particle mass, kg
        mp          = (1-self.X) * self.mc0 + self.ma0
        #.. conversion particle mass, kg/s
        dmdt        = self.mc0 * self.dXdt
        #.. Sherwood number, Adaption due to Blowing
        Sh,Bscaling = self.Sheerwood(flag=True)

        #.. ps: partial pressure of reactive gas at the surface, defined by Pi_surface()

        data_row = [
            round(float(t),8),                  # time, s
            round(float(mp),13),                # particle mass, kg
            round(float(dmdt),13),              # change in mass, kg/s
            round(float(self.X),6),             # carbon conversion, kg/kg
            round(float(self.dXdt),9),          # total change in carbon conversion, kg/kg 1/s
            round(float(self.ratep['O2']), 7),  # change in carbon conversion for C+O2, kg/kg 1/s 
            round(float(self.ratep['CO2']),7),  # change in carbon conversion for C+CO2, kg/kg 1/s
            round(float(self.ratep['H2O']),7),  # change in carbon conversion for C+H2O, kg/kg 1/s
            round(self.C1['O2'],  17),          # Diffusion Rate Constant for O2
            round(self.C1['CO2'], 17),          # Diffusion Rate Constant for CO2
            round(self.C1['H2O'], 17),          # Diffusion Rate Constant for H2O
            round(self.effF['O2'], 7),          # effectiveness factor for O2
            round(self.effF['CO2'],7),          # effectiveness factor for CO2
            round(self.effF['H2O'],7),          # effectiveness factor for H2O
            round(float(Sh),3),                 # Sheerwood number
            round(float(Bscaling), 3),          # Adaption due to Blowing
            round(float(self.porosity), 5),     # Porosity development, m3/m3
            round(self.Tp,1),                   # particle Temperature, K
            round(float(self.rhop), 2),         # particle density, kg/m3
            round(self.dp, 7),                  # particle diameter, m
            round(self.ps, 1),                  # partial pressure at the surface, only feaseable for one component reactive gas
        ]

        new_data_row = DataFrame([data_row], columns=self.result.columns)
        self.result = concat([self.result, new_data_row], ignore_index=True)

        return self.dXdt

    #============#
    #============#

    def rate_3Rkt(self):
        '''
        todo: check rate definition!
        r = eta * r(Pi) * S/m
        r = eta * k(T) * Pi^n * S/m
        r = eta * k0 * exp(E/R.T) * Pi^n * S/m
        :return:
        '''

        def rate():
             #-- Be aware that sint0 is included in kin[0]
             #-- Definition identical to Kajitani
             #-- dX/dt = eta * r_S * (1-X) * s_r
             #   1/s = 1/s.Pa^n * Pa^n * kg/kg * m2/m2
             #-- is identical to the assumption by DeYoung p.40 dXC/dt -> 1/Sr dmC/dt
             #---------------------------------

             effF = self.effF[self.name_reactive_specie]
             r_S  = self.kinetic_rate_coeff * pow(self.ps, self.kin[2])
             ns_r = self.normalized_Area

             rate = effF * r_S * (1-self.X) * ns_r

             return rate
        
        #=============#
        
        def update_reactive_specie(reactive_specie):
            '''
            Sets global values:
            '''

            #-----------------------------------------------------------------#
            #>> Error handling
            if ((reactive_specie=='O2') or (reactive_specie=='H2O') or (reactive_specie=='CO2')):
                self.name_reactive_specie=reactive_specie
            else:
                raise IOError('WARNING: No update, as no known reactive species is selected ')

            #-----------------------------------------------------------------#
            #>> Set partial pressure of the given species
            self.pg = self.gas_comp_dict[self.name_reactive_specie] * self.ptot

            #-----------------------------------------------------------------#
            #>> update heterogeneous kinetic parameters
            if self.name_reactive_specie == 'O2':

                if self.pg > 0:
                    #>> constants set!
                    self.kin = self.kin_R1 # set Arrhenius expressions A,E,n

                    self.stoichio=0.5
                    self.psi1 = self.psi0[self.name_reactive_specie]

                    #>> dervived values!
                    self.ps = self.Pi_surface()
                    self.effF[self.name_reactive_specie] = self.effectivenessFactor(self.ps)

            elif self.name_reactive_specie == 'CO2':

                if self.pg > 0:
                    #>> constants set!
                    self.kin = self.kin_R2
                    self.stoichio=1.0
                    self.psi1 = self.psi0[self.name_reactive_specie]

                    #>> dervived values!
                    self.ps = self.Pi_surface()
                    self.effF[self.name_reactive_specie] = self.effectivenessFactor(self.ps)


            elif self.name_reactive_specie == 'H2O':

                if self.pg > 0:
                    #>> constants set!
                    self.kin = self.kin_R3
                    self.stoichio=1.0
                    self.psi1 = self.psi0[self.name_reactive_specie]

                    #>> dervived values!
                    self.ps = self.Pi_surface()
                    self.effF[self.name_reactive_specie] = self.effectivenessFactor(self.ps)

            return 0
        
        #=============#

        update_reactive_specie('O2')
        if self.pg > 0:
            R1 = rate()
        elif self.pg == 0:
            R1 = 0
        elif self.pg < 0:
            raise IOError('Negative partial pressure')
        #-----
        update_reactive_specie('CO2')
        if self.pg > 0:
            R2 = rate()
        elif self.pg == 0:
            R2 = 0
        elif self.pg < 0:
            raise IOError('Negative partial pressure')
        #-----
        update_reactive_specie('H2O')
        if self.pg > 0:
            R3 = rate()
        elif self.pg == 0:
            R3 = 0
        elif self.pg < 0:
            raise IOError('Negative partial pressure')
        #-----
        # 1/s * ((J/kg) / (1/s))

        return float(R1),float(R2),float(R3)

    #=============#
    #=============#

    @property
    def kinetic_rate_coeff(self):
        """
        #>> kinetic reaction coefficient
        #>> 1/S * dm/dt , kg/(s m^2 Pa^n)
        #>> 1 /s.Pa^n = 1/s.Pa^n * exp[ (J/mol) / (J/mol.K * K) ]
        """
        kinConst = self.kin[0] * np.exp(-1. * self.kin[1]/(self._gas_const * self.Tp))

        return kinConst

    #=============#
    #=============#

    @property
    def dp(self):
        #>> particle diameter
        #>> Essenhigh correlation
        #>> beta = 1/3 shinking particle , beta = 0 shrinking density
        return self.dp0 * pow( (1.-self.X), self.beta )

    #=============#
    #=============#

    @property
    def rhoc(self):
        # apparent density of char fraction, kg/m^3
        return self.rhoc0 * pow((1.-self.X),self.alpha)

    #=============#
    #=============#

    @property
    def rhop(self):
        # apparent density of particle, particle density
        # rhoc and rhoa0 should be also apperent densities
        # rhop = 1/( (1-self.ya0)/self.rhoc + self.ya0/self.rhoa0 )

        def oneminusXp(X=self.X,wc0=self.yc0):
            return wc0 * ((1-X) + (1/wc0 -1))

        if self.Sdevelmodel == 'SPM' or self.Sdevelmodel == 'SPMp':
            rhop = self.rhop0
        elif self.Sdevelmodel=='SDM' or self.Sdevelmodel=='RPM' or self.Sdevelmodel=='SDMp':
            rhop = self.rhop0 * oneminusXp()
        else:
            KeyError('Model not available!')

        return rhop

    #=============#
    #=============#

    @property
    def rhop_true(self):
        # true density of the particle
        return self.rho_char_init_true

    #=============#
    #=============#

    @property
    def porosity(self):
        """
        Particle porosity.
        The porosity keeps constant with the current implementation.
        """
        return 1 - self.rhop/self.rhop_true

    #=============#
    #=============#

    @property
    def normalized_Area(self):
        # >> s/s0 = (S/m)/(S/m)_0 = m2/m2<<#
        # >> RPM .. Random pore model m2/m2
        # >> SPM .. Shrinking Particle Model m2/m2
        # >> SDM .. Shrinking Density Model, m2/m2
        #------
        def oneminusXp(X=self.X,wc=self.yc0):
            return wc * ((1-X) + (1/wc -1))
        #------

        oneMXp = oneminusXp()

        if self.Sdevelmodel=='RPM':
            if self.X == 1:
                ss0 = 1
            else:
                ss0 = pow(1 - self.psi1 * np.log(1. - self.X), 0.5)

        elif self.Sdevelmodel == 'SPM':
            if self.X == 1:
                ss0 = 1
            else:
                ss0 = pow( (1-self.X) , -1/3)

        elif self.Sdevelmodel == 'SPMp':
            if self.X == 1:
                ss0 = 1
            else:
                ss0 = pow( oneMXp , -1/3)

        elif self.Sdevelmodel == 'SDM':
            if self.X == 1:
                ss0 = 0
            else:
                ss0 = pow( (1-self.X) , -1)

        elif self.Sdevelmodel == 'SDMp':
            if self.X == 1:
                ss0 = 0
            else:
                ss0 = pow( oneMXp , -1)

        else:
            raise IOError('Unkown Methode')

        return ss0

    #=============#
    #=============#

    @property
    def diff_coeff(self):
        # molecular diffusion, multiple estimation methods are given in Prausnitz

        comp = {}
        pressure1 = self.ptot
        Tm = 1/2*self.Tp + 1/2*self.Tg

        #----------------------------------------------
        if self.gas_default:
            comp.update( self.gas_default_diff_partner[self.gas_default] )
            comp.update( {self.gas_default:1} )
        elif self.name_reactive_specie == 'O2':
            comp.update({
                'O2' : self.gas_comp_dict['O2'],
                'CO2': self.gas_comp_dict['CO2'],
                'CH4': self.gas_comp_dict['CH4'],
                'CO' : self.gas_comp_dict['CO'],
                'H2' : self.gas_comp_dict['H2'],
                'H2O': self.gas_comp_dict['CO2'],
                'N2' : self.gas_comp_dict['N2'],
                    })
        elif self.name_reactive_specie == 'CO2':
            comp.update({
                'CO2': self.gas_comp_dict['CO2'],
                'CH4': self.gas_comp_dict['CH4'],
                'CO':  self.gas_comp_dict['CO'],
                'H2':  self.gas_comp_dict['H2'],
                'H2O': self.gas_comp_dict['CO2'],
                'N2':  self.gas_comp_dict['N2'],
                'O2':  self.gas_comp_dict['O2'],
                })
        elif self.name_reactive_specie == 'H2O':
            comp.update({
                'H2O': self.gas_comp_dict['CO2'],
                'CO2': self.gas_comp_dict['CO2'],
                'CH4': self.gas_comp_dict['CH4'],
                'CO':  self.gas_comp_dict['CO'],
                'H2':  self.gas_comp_dict['H2'],
                'N2':  self.gas_comp_dict['N2'],
                'O2':  self.gas_comp_dict['O2'],
            })
        #----------------------------------------------

        diffusionD = self.diffC.diffusion_mix(comp,Tm,pressure1)

        return diffusionD

    #=============#
    #=============#

    def effectivenessFactor(self,ps):
        '''
        Effectivenes factor (Porennutzungsgrad)
        '''

        def D_eff(D=self.diff_coeff, eps1=self.porosity, tauf=self.tauf):
            return eps1 * D * tauf

        def thielemodulus(
            dp       = self.dp,
            n        = self.kin[2],
            stoichio = self.stoichio,
            ns_r     = self.normalized_Area,
            rhoc     = self.rhoc,
            kin      = self.kinetic_rate_coeff,
            pressure_surf = self.ps,
            Ru       = self._gas_const,
            Mc       = self.Mc,
            Tp       = self.Tp,
            Deff     = D_eff()   ):

            '''
            Thiele Modulus for nth Order:
                - Bader 2018, Gärtner 2015 phd, Hla_2006(Eq.(13)), Laurendau_1978(Eq.(97))
            Thiele Modulus for 0th Order:
                - DeYoung or Tremel, simplified version
            '''

            #  molar_rate  = 1/V dn/dt [=] 1/m3 mol/s = mol/kg * kg/m3 * 1/s.Pa^n Pa^n
            molar_rate = stoichio/1 * 1/Mc * rhoc * ns_r * kin * pow( pressure_surf, n )

            if pressure_surf > 0:
                #>> [Th] = m * 1/m = m * (1/m3 mol/s * Pa.m3/K.mol * K 1/Pa s/m2)^0.5:
                Th =  dp/6 * pow ( (n+1)/2 * molar_rate / (pressure_surf/(Ru * Tp) * Deff), 0.5)
            else:
                Th = 2000

            return Th

        #----------------------#
        # effectiveness factor #
        #----------------------#
        Th = thielemodulus(pressure_surf=ps)

        if Th < 0.001:
            effF = 1.0
        elif Th == 2000:
            effF = 0.0
        else:
            effF = 1/Th * ( 1/np.tanh(3*Th) - 1/(3*Th) )

        return effF

    #============#
    #============#

    @property
    def diffusion_rate_coeff(self):
        #>> diffusion constant, kg/s.m2.Pa
        #>>  - Sherwood number = Nusselt number if Lewis = 1
        Tp=self.Tp                  # K
        Tgas=self.gas.T             # K
        dp=self.dp                  # m
        Tm = 1/2 * Tp + 1/2 * Tgas  # K
        DiffC = self.diff_coeff     # m2/s
        # --------------------

        # Sheerwood, -
        Sh = self.Sheerwood(flag=False)[0]

        # kdiff [=] kg/s.m2.Pa = kg/mol * 1/(m * J/mol.K K) * m2/s
        kdiff   = 1/self.stoichio * self.Mc * Sh / (dp * self._gas_const * Tm) * DiffC

        self.C1[self.name_reactive_specie] = kdiff/(pow(Tm,0.75)/dp)

        return kdiff

    #============#
    #============#

    def Sheerwood(self,flag):

        # Definition of the Sheerwood number
        Sh0 = 2.0 + 0.6 * pow(self.Re,0.5) * pow(self.Pr,1/3)

        #----------------------------------
        if flag==True:

            def blowing_number(
                    dXdt = self.dXdt,
                    mc0  = self.mc0,
                    gas  = self.gas,
                    dp   = self.dp,
                    Pr   = self.Pr,
                    Sh0  = 2):
                # [B]org  = [-] = 1/m2 * 1/s s.m2 = 1/kg.m2 * kg/s * kg * s.m/kg * m
                # [B]here = [-] = 1/s * s = 1/s * kg s.m/kg.m
                B = (Pr / Sh0) * dXdt * mc0 / (np.pi * dp * gas.viscosity)
                if (B < 1.e-4):
                    B = 1.e-4
                return B

            B = blowing_number(Sh0 = Sh0)
            scaling = B / (np.exp(B) - 1)
            Sh = Sh0 * scaling
        #----------------------------------

        elif flag==False:
            scaling = 1         #>> needed for the return value
            Sh = Sh0 * scaling

        return Sh,scaling

    #============#
    #============#

    def Pi_surface(self):

        def Func(PS):

            #---------------------------------------------
            #-- Diffusion rate equals Kinetic rate
            k_diff = self.diffusion_rate_coeff
            k_rate = self.kinetic_rate_coeff
            eta    = self.effectivenessFactor(ps=PS)
            #-------------------------------------------
            ptot = self.ptot
            xi_inf = self.pg / ptot  # -- species i mole fraction in the gas phase mol/mol
            xi_s = PS / ptot         # -- species i mole fraction at the surface mol/mol
            rhoc = self.rhoc
            dp   = self.dp
            ns_r = self.normalized_Area
            #-----
            diffRate = k_diff * ptot * np.log((1+xi_inf)/(1+xi_s)) * 6/(rhoc*dp)
            kinRate  = eta * k_rate * pow(PS, self.kin[2]) * ns_r
            #-------------------------------------------
            diff = kinRate - diffRate

            return diff

        #=========#

        if self.pg > 0:
            min_PS=0
            max_PS=self.pg
            ps = bisect(
                Func,
                min_PS,
                max_PS,
                xtol=1e-12, rtol=10e-16, maxiter=100, full_output=False, disp=True)
        else:
            ps = 0

        return ps

##== END of Class Definition ==##
