# SPDX-FileCopyrightText: 2023 Thomas Förster
#
# SPDX-License-Identifier: MIT

__author__ = 'tf'

from particle_conversion.models.solver_zeroDim_reactor import Solver_class

import numpy as np
from particle_conversion.default_parameters import _defaults
from os.path import splitext

class ZEROreactor():
    def __init__(self,settings,output):
                 
        self.settings = settings
        self.output = output
        # Simulation basic numerical boundary conditions
        self.calctime = settings["numerical"]["calctime"] if "calctime" in settings["numerical"].keys() else _defaults["numerical"]["calctime"]
        self.steps    = settings["numerical"]["steps"] if "steps" in settings["numerical"].keys() else _defaults["numerical"]["steps"]
        
        # Simulation basic boundary conditions
        self.default_temp = settings["reactor"]["fallback_temp"] if "fallback_temp" in settings["reactor"].keys() else _defaults["reactor"]["fallback_temp"]
        self.fallback_gas = settings["reactor"]["fallback_gas"] if "fallback_gas" in settings["reactor"].keys() else _defaults["reactor"]["fallback_gas"]

        # Particle transient environment
        self.data_co   = settings["environment"]["co"]  if "co"  in settings["environment"].keys() else _defaults["environment"]["co"]
        self.data_co2  = settings["environment"]["co2"] if "co2" in settings["environment"].keys() else _defaults["environment"]["co2"]
        self.data_h2o  = settings["environment"]["h2o"] if "h2o" in settings["environment"].keys() else _defaults["environment"]["h2o"]
        self.data_o2   = settings["environment"]["o2"]  if "o2"  in settings["environment"].keys() else _defaults["environment"]["o2"]
        self.data_n2   = settings["environment"]["n2"]  if "n2"  in settings["environment"].keys() else _defaults["environment"]["n2"]
        self.data_ch4  = settings["environment"]["ch4"] if "ch4" in settings["environment"].keys() else _defaults["environment"]["ch4"]
        self.data_h2   = settings["environment"]["h2"]  if "h2"  in settings["environment"].keys() else _defaults["environment"]["h2"]
        self.data_Re   = settings["environment"]["Re"]  if "Re"  in settings["environment"].keys() else _defaults["environment"]["Re"]
        self.data_Tg   = settings["environment"]["Tg"]  if "Tg"  in settings["environment"].keys() else _defaults["environment"]["Tg"]
        self.data_Tp   = settings["environment"]["Tp"]  if "Tp"  in settings["environment"].keys() else _defaults["environment"]["Tp"]
        self.data_law  = settings["environment"]["law"] if "law" in settings["environment"].keys() else _defaults["environment"]["law"]


        ################################
        ## (A) Setup the input arrays ##
        ################################

        self._comp_header   = { 'CO':0, 'CO2':1, 'H2O':2, 'O2':3, 'N2':4, 'CH4':5, 'H2':6 }
        self._comp_length   = { 'CO':1, 'CO2':1, 'H2O':1, 'O2':1, 'N2':1, 'CH4':1, 'H2':1 }
        self._comp_flags    = { 'CO':False, 'CO2':False, 'H2O':False, 'O2':False, 'N2':False, 'CH4':False, 'H2':False }
        self._comp_list     = { 'CO':1, 'CO2':1, 'H2O':1, 'O2':1, 'N2':1, 'CH4':1, 'H2':1 }
        self._comp_default  = [[0.0,self.calctime],[0.0,0.0]]
        self.comp_calc      = []

        self._others_header = { 'Re':0, 'Tg':1, 'Tp':2 }
        self._others_length = { 'Re':1, 'Tg':1, 'Tp':1 }
        self._others_flags  = { 'Re':False, 'Tg':False, 'Tp':False, 'Law5':False }
        self._others_list   = { 'Re':1, 'Tg':1, 'Tp':1 }
        self._default_Re    = [[0.0,self.calctime],[0.0,0.0]]
        self._default_Tp    = [[0.0,self.calctime],[self.default_temp,self.default_temp]]
        self._default_Tg    = [[0.0,self.calctime],[self.default_temp,self.default_temp]]
        self._default_law   = [[0.0,self.calctime],[5,5]]
        self.others_calc    = []

        #--------------------------------------------------------------------
        # filter function
        # - data from only law 5 considered
        # - used in self.reading_data()
        if self.data_law:
            self._others_law = self._reading_law(self.data_law)
            self._others_flags['Law5']   = True
        else:
            self._others_law = self._default_law

        #--------------------------------------------------------------------
        # reading data incl. filtering
        if self.data_co:
            self._comp_list['CO']    = self.reading_data(self.data_co)
            self._comp_flags['CO']   = True
            self._comp_length['CO']  = len(self._comp_list['CO'])
        else:
            self._comp_list['CO']    = self._comp_default
        #-----
        if self.data_co2:
            self._comp_list['CO2']   = self.reading_data(self.data_co2)
            self._comp_flags['CO2']  = True
            self._comp_length['CO2'] = len(self._comp_list['CO2'])
        else:
            self._comp_list['CO2']   = self._comp_default
        #-----
        if self.data_h2o:
            self._comp_list['H2O']   = self.reading_data(self.data_h2o)
            self._comp_flags['H2O']  = True
            self._comp_length['H2O'] = len(self._comp_list['H2O'])
        else:
            self._comp_list['H2O']   = self._comp_default
        #-----
        if self.data_o2:
            self._comp_list['O2']    = self.reading_data(self.data_o2)
            self._comp_flags['O2']   = True
            self._comp_length['O2']  = len(self._comp_list['O2'])
        else:
            self._comp_list['O2']    = self._comp_default
        #-----
        if self.data_n2:
            self._comp_list['N2']    = self.reading_data(self.data_n2)
            self._comp_flags['N2']   = True
            self._comp_length['N2']  = len(self._comp_list['N2'])
        else:
            self._comp_list['N2']    = self._comp_default
        #-----
        if self.data_ch4:
            self._comp_list['CH4']   = self.reading_data(self.data_ch4)
            self._comp_flags['CH4']  = True
            self._comp_length['CH4'] = len(self._comp_list['CH4'])
        else:
            self._comp_list['CH4']   = self._comp_default
        #-----
        if self.data_h2:
            self._comp_list['H2']    = self.reading_data(self.data_h2)
            self._comp_flags['H2']   = True
            self._comp_length['H2']  = len(self._comp_list['H2'])
        else:
            self._comp_list['H2']    = self._comp_default
        # -----
        if self.data_Re:
            self._others_list['Re']  = self.reading_data(self.data_Re)
            self._others_flags['Re'] = True
            self._others_length['Re']= len(self._others_list['Re'])
        else:
            self._others_list['Re']  = self._default_Re
        # -----
        if self.data_Tg:
            self._others_list['Tg']  = self.reading_data(self.data_Tg)
            self._others_flags['Tg'] = True
            self._others_length['Tg']= len(self._others_list['Tg'])
        else:
            self._others_list['Tg']  = self._default_Tg
        # -----
        if self.data_Tp:
            self._others_list['Tp']  = self.reading_data(self.data_Tp)
            self._others_flags['Tp'] = True
            self._others_length['Tp']= len(self._others_list['Tp'])
        else:
            self._others_list['Tp']  = self._default_Tp
        # -----


        #--------------------------------------------------------------------
        # Determine the number of cycles
        self._num_length = []
        for i,val in enumerate(self._comp_length):
            self._num_length.append(self._comp_length[val])

        self._num_cycles = np.array(self._num_length).max()


        #--------------------------------------------------------------------
        # When all reactants are zero fallback_gas will be set to 1 mol/mol
        if (self._comp_flags['O2'] == False and self._comp_flags['CO2'] == False and self._comp_flags['H2O'] == False):
            if self.fallback_gas in ['O2','CO2','H2O']:
                self._comp_list[self.fallback_gas] = [[0.0,self.calctime],[1.0,1.0]]
            else:
                raise IOError('No reactive agend are present!')


        #--------------------------------------------------------------------
        # Set the whole dataset
        self._gas_comp   = []
        self._others     = []

        data_tmp = []
        for i,val in enumerate(self._comp_header):
        # co co2 h2o o2 n2 ch4 h2

            if self._comp_flags[val] == False:
                for j in range(self._num_cycles):
                    data_tmp.append(self._comp_list[val])
            elif self._comp_flags[val] == True:
                data_tmp = self._comp_list[val]

            self._gas_comp.append(data_tmp)
            data_tmp = []

        for i, val in enumerate(self._others_header):
            # Re Tg Tp

            if self._others_flags[val] == False:
                for j in range(self._num_cycles):
                    data_tmp.append(self._others_list[val])
            elif self._others_flags[val] == True:
                data_tmp = self._others_list[val]

            self._others.append(data_tmp)
            data_tmp = []

    #=============#
    #== Methods ==#
    #=============#
    def calc(self):

        ################################
        ## (B) Calculating 0D reactor ##
        ################################

        for i in range(self._num_cycles):
            for j,val in enumerate(self._comp_header):
                self.comp_calc.append(self._gas_comp[self._comp_header[val]][i])

            for j,val in enumerate(self._others_header):
                self.others_calc.append(self._others[self._others_header[val]][i])

            #-----------------------------------------------------
            #-- Initialisation of the Solver
            solver_0D_reactor = Solver_class(
                settings                    = self.settings,
                comp_header                 = self._comp_header,
                comp_flags                  = self._comp_flags,
                comp                        = self.comp_calc,
                others_header               = self._others_header,
                others_flags                = self._others_flags,
                others                      = self.others_calc,
            )

            #-----------------------------------------------------
            #-- Solver Call
            #print('=== Round:',i)
            #print('============================================')
            # Start time for calculation, s
            t_start = 0.0
            # Final time, s
            t_end = self.calctime
            # number of time steps
            steps = self.steps
            
            #-- RUN SOLVER (this is the main step of this project)
            solver_0D_reactor.solve_char(t_start,t_end,steps)
            #print('')
            #print('')

            #-- WRITE OUTPUT
            file_output = f'{splitext(self.output)[0]}_run{i:02}.csv'

            solver_0D_reactor.char.result.sort_values('Time',inplace=True)
            solver_0D_reactor.char.result.to_csv(file_output, sep=',')

            self.comp_calc   = [] # reset
            self.others_calc = [] # reset


    def reading_data(self,filename):

        xval = []
        yval = []
        dataset = []
        data = []
        index=0
        time0=0
        time_set=False
        #with open('./pefr/p_co2-source_mod.xy') as f:

        with open(filename) as f:
            for line in f:
                if line == '\n':
                    continue
                elif line == '\t\n':
                    continue
                elif line.startswith('#'):
                    continue
                else:
                    xi,yi=line.split()

                    xi = float(xi)
                    yi = float(yi)

                    if time_set == True and xi < time0:
                        time_set=False

                    # Determines if you either have a dataset or not
                    if self._others_flags['Law5']==True:
                        if func_inter(xi, self._others_law[index]) >= 5:
                            law5_flag = True
                        else:
                            law5_flag = False
                    elif self._others_flags['Law5'] == False:
                        if func_inter(xi, self._others_law) >= 5:
                            law5_flag = True
                        else:
                            law5_flag = False


                    if law5_flag==True:

                        if time_set==False:
                            time0    = xi
                            xi       = 0.0
                            time_set = True

                        if xi == 0.0 and len(xval)>0:
                            dataset = [xval,yval]
                            data.append(dataset)
                            xval = []
                            yval = []
                            index += 1

                        if xi == 0.0:
                            xval.append(float(xi))
                            yval.append(float(yi))
                        else:
                            xval.append(float(xi)-time0)
                            yval.append(float(yi))


        #-- Last Dataset append --#
        dataset = [xval,yval]
        data.append(dataset)

        f.close()

        return data

    #==============#
    #==============#

    def _reading_law(self,filename):

        xval = []
        yval = []
        dataset = []
        data = []

        #with open('./pefr/p_co2-source_mod.xy') as f:
        with open(filename) as f:
            for line in f:
                if line=='\n':
                    continue
                elif line=='\t\n':
                    continue
                elif line.startswith('#'):
                    continue
                else:
                    xi,yi=line.split()

                    xi = float(xi)
                    yi = float(yi)

                    if xi == 0.0 and len(xval)>0:
                        dataset = [xval,yval]
                        data.append(dataset)
                        xval = []
                        yval = []

                    xval.append(float(xi))
                    yval.append(float(yi))

        #-- Last Dataset append --#
        dataset = [xval,yval]
        data.append(dataset)

        f.close()

        return data


def func_inter(timex,data):
    time_array = data[0]
    func_array = data[1]
    funcX = np.interp(timex,time_array,func_array)
    return float(funcX)
