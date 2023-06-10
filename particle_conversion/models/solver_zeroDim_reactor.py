# SPDX-FileCopyrightText: 2023 Thomas Förster
#
# SPDX-License-Identifier: MIT

__author__ = 'tf'

from particle_conversion.models.particle import Particle_class
from matplotlib import pyplot as plt
import numpy as np
from scipy.integrate import ode
import cantera as ct
from particle_conversion.default_parameters import _defaults


class Solver_class():
    def __init__(self,
        settings,
        comp_header,
        comp_flags,
        comp,
        others_header,
        others_flags,
        others,
    ):
        self.comp_header = comp_header
        self.comp_flags = comp_flags
        self.comp = comp
        self.others_header = others_header
        self.others_flags = others_flags
        self.others = others

        #.. numerical settings
        self.X_end = settings["numerical"]["conversion_final"] if "conversion_final" in settings["numerical"].keys() else _defaults["numerical"]["conversion_final"]
        self.X0 = settings["numerical"]["conversion_init"] if "conversion_init" in settings["numerical"].keys() else _defaults["numerical"]["conversion_init"]
        
        #.. reactor boundary conditions
        self.ptot = settings["reactor"]["reactor_pressure"] if "reactor_pressure" in settings["reactor"].keys() else _defaults["reactor"]["reactor_pressure"]

        #---------------------------------------------------------#
        #>> THIS is the INIT for particle conversion EVERYTHING <<#
        #---------------------------------------------------------#
        #-- char particle object initialisation
        self.char = Particle_class(
            settings=settings['particle'],
            gas_header=comp_header,
            gas_flags=comp_flags,
            gas_comp=comp,
            others_header=others_header,
            others_flags=others_flags,
            others=others,
            ptot=self.ptot,
        )
        #-------------------------------------#

        
        # initial value for the solver
        self.y0 = self.X0 
        self.t  = []
        self.number_of_previous_steps = 0

        # -- END Init -- #
        #-------------------------------------#


    def solve_char(self,t_start,t_end,steps):

        # store solver settings globally for later evaluation
        self.t_start = t_start  # t_start = 0.0 # s
        self.t_end = t_end      # t_end = 0.12 # s
        self.steps = steps        # = 100
        dt = (t_end - t_start)/steps
        self.dt = dt

        #======================================
        #== set solver
        # - BDF method suited to stiff systems of ODEs
        # - function which are solved is self.char.conv()
        ode_system = ode(self.char.conv).set_integrator('vode', method='bdf')

        #======================================
        # set initial arrays
        t = np.zeros((steps, 1))
        X = np.zeros((steps, 1))

        #======================================
        #== store initial values of the ODE
        t[0] = t_start
        X[0] = self.y0

        ode_system.set_initial_value(self.y0,t_start)

        #======================================
        #== calculation run
        # - each time step will be evaluated
        k=1
        while ode_system.successful() and k < steps and (X[k-1] < self.X_end):  # << really nice feature!

            test = ode_system.t + dt
            ode_system.integrate(ode_system.t + dt)

            'Store the results of the X and time of the particle'
            X[k] = float(ode_system.y)
            t[k] = ode_system.t
            k += 1
