import numpy as np
import pandas as pd
import re
import sys
import os
import yaml
import ruamel.yaml as rmy

# Rich printing
from rich.console import Console
from rich.table import Table
from rich.markdown import Markdown
from rich import box
from rich.columns import Columns
from rich.panel import Panel
from rich.prompt import Prompt
from rich.prompt import Confirm
from rich.layout import Layout

# Matplotlib related
import matplotlib.pyplot as plt
plt.style.use('science')


class Analysis:
    """Angle-of-Attack Analysis class for Pre- and Postprocessing"""


    def __init__(self, setupFile):
        console = Console()
        # console.print("\n[bold bright_white]Welcome to the Angle-of-Attack Analysis Tool![not bold bright_white]\n\nSelect the mode  [bold yellow]0: Preprocessing    [bold cyan]1: Postprocessing\n")
        self.mode = Prompt.ask("\n[bright_white]Welcome to the Angle-of-Attack Analysis Tool![not bold bright_white]\nSelect the mode  [bold cyan]0: Preprocessing    [bold yellow]1: Postprocessing\n", choices=["0", "1"], default="0")
        with open(setupFile, "r") as f:
            self.setup = yaml.safe_load(f)
    

    def ambientConditions(self):
        """Function to calculate and return ambient condition calculations for free-stream"""
        # 1. Check pressure and temperature. If not given, use standard conditions
        try:
            self.p = float(self.setup["Boundary conditions"]["p"])
        except:
            self.p = 1e5
        try:
            self.T = float(self.setup["Boundary conditions"]["T"])
        except:
            self.T = 298.15
        try:
            self.d = float(self.setup["Boundary conditions"]["d"])
        except:
            self.d = 1.0

        self.rho = 0.029 / 8.314 * self.p / self.T # Density according to material data
        # Calculate viscosity according to Sutherland
        self.mu = 1.716e-5 * (self.T / 273.15)**1.5 * (273.15 + 110.4) / (self.T + 110.4)
        self.nu = self.mu / self.rho
        self.a = np.sqrt(1.4 * 8.314 / 0.029 * self.T)

        self.ambient_dict = {
            "p": float(self.p), "T": float(self.T), "Object width d": float(self.d), "rho": float(self.rho),
            "Re": float(self.Re), "u": float(self.u_inf), "Ma": float(self.Ma), "mu": float(self.mu), "nu":
            float(self.nu), "C_f": float(self.C_f), "tau_w": float(self.tau_w), "u_tau": float(self.u_tau),
            "yplus": float(self.yplus), "y_1": float(self.y_1)
        }


    def objectCalculations(self, obj):
        """Function to calculate and return ambient condition calculations for free-stream"""
        console = Console()
        L = self.setup["Objects"][obj]["L"]

        # TODO: TRANSONIC AND HIGHER MACH NUMBER FLOWS
        def yPlus(self, Re, u_inf):
            C_f = 0.0576 * Re**(-0.2)  # skin friction coefficient
            tau_w = 0.5 * C_f * self.rho * u_inf**2 # wall shear stress
            u_tau = np.sqrt(tau_w / self.rho) # friction velocity
            yplus = 1.0
            y_1 = yplus * self.nu / u_tau # first layer height
            delta = 0.37 * L * Re**(-0.2)
            return y_1, delta

        # 2. Check if the free-stream is well defined
        if self.setup["Boundary conditions"]["Re"]:
            Re = float(self.setup["Boundary conditions"]["Re"])
            u_inf = Re * self.mu / (L * self.rho)
            Ma = u_inf / self.a
            y_1, delta = yPlus(self, Re, u_inf)
        elif self.setup["Boundary conditions"]["u"]:
            u_inf = float(self.setup["Boundary conditions"]["u"])
            Re = u_inf * L * self.rho / self.mu
            Ma = u_inf / self.a
            y_1, delta = yPlus(self, Re, u_inf)
        elif self.setup["Boundary conditions"]["Ma"]:
            Ma = float(self.setup["Boundary conditions"]["Ma"])
            u_inf = Ma * self.a
            Re = u_inf * self.rho * L / self.mu
            y_1, delta = yPlus(self, Re, u_inf)
        else:
            console.print("\n[bold red]Error: [not bold red]No free-stream Reynolds number, velocity or Mach number defined.")
            sys.exit()

        return y_1, Ma, delta


    def exportCalculations(self):
        with open(exportFile, "w") as f:
            # rmy.dump(calculations_dict, f)
            yaml.safe_dump(calculations_dict, f)
        
        explanation_str = """
# C_f:   Skin friction coefficient.......(-)
# Ma:    Mach number.....................(-)
# d:     Object depth....................(m)
# Re:    Reynolds number.................(-)
# T:     Free-stream temperature.........(K)
# mu:    Free-stream dynamic viscosity...(Pa s)
# nu:    Free-stream kinematic viscosity.(m2 s-1)
# p:     Free-stream pressure............(Pa)
# rho:   Free-stream density.............(kg m-3)
# tau_w: Wall shear-stress...............(N m-2)
# u:     Free-stream velocity............(m s-1)
# u_tau: Shear-velocity..................(m s-1)
# y_1:   First layer height for y+.......(m)
# yplus: Dimensionless wall distance.....(-)"""
        # Open a file with access mode 'a'
        with open(exportFile, "a") as f:
            f.write(explanation_str)
        console.print(f"\nPre-calculations written to: \'{exportFile}\'\n")

        # 3. Run Calculations if defined
        try:
            exportFile = self.setup["I/O"]["pre-calculations"]
        except:
            pass
        else:
            exportCalculations(self, exportFile)

    def pre(self, obj):
        """Preprocessing method"""
        pass

        # 4. Create the Run File for the solver
        # def runFile(self):
        #     """Method to create a run """
        #     solver = self.setup["Numerics"]["solver"]

        #     angleList = np.linspace(self.aoaMin, self.aoaMax, self.numberOfIncrements + 1, retstep=True)
        #     delta_alpha = angleList[1]

        #     # Start transcript in order for the console output to be written to a file
        #     fluentString = ';{}\n/file/start-transcript ../{} o\n\n'.format(self.airfoilName, self.transcriptFilename)
        #     fluentString += '/file/replace-mesh {} o\n'.format(self.meshFilename)
        #     fluentString += '/define/bound/set/pres/inlet_farfield () mach n {} ()\n'.format(self.Ma)

        #     # First solution at minimum angle:
        #     fluentString += '''/mesh/rotate {} 0 0 0 0 0 1\n'''.format(-1 * self.aoaMin)  # rotate to initial angle
        #     fluentString += '''/solve/init/hyb\n/solve/it {}\n'''.format(numberOfIterations)  # initial solution
        #     fluentString += '''(display "Angle is: {} deg")\n'''.format(self.aoaMin)
        #     fluentString += '''/report/forces/wall-forces yes 1 0 0 no\n''' # print the drag forces and coefficient
        #     fluentString += '''/report/forces/wall-forces yes 0 1 0 no\n\n''' # print the lift forces and coefficient
        #     # if generateImages:
        #     #     fluentString += '''/display/views/rest view-0\n/display/objects/display logk\n'''
        #     #     fluentString += '''/display/save-picture ../res/img/{}_alpha_{}_logk.png\n'''.format(self.airfoilName, round(self.aoaMin))
        #     #     fluentString += '''/display/views/rest view-0\n/display/objects/display u\n'''
        #     #     fluentString += '''/display/save-picture ../res/img/{}_alpha_{}_u.png\n'''.format(self.airfoilName, round(self.aoaMin))
        #     #     fluentString += '''/display/views/rest view-0\n/display/objects/display p\n'''
        #     #     fluentString += '''/display/save-picture ../res/img/{}_alpha_{}_p.png\n'''.format(self.airfoilName, round(self.aoaMin))

        #     for i in list(range(len(angleList[0]) - 1)):
        #         fluentString += '''/mesh/rotate {} 0 0 0 0 0 1\n'''.format(-1 * delta_alpha)  # rotate by steplength
        #         fluentString += '''/solve/init/hyb o\n/solve/it {}\n'''.format(numberOfIterations)
        #         fluentString += '''(display "Angle is: {} deg")\n'''.format(round(angleList[0][i + 1], 3))
        #         fluentString += '''/report/forces/wall-forces yes 1 0 0 no\n'''
        #         fluentString += '''/report/forces/wall-forces yes 0 1 0 no\n\n'''
        #         # if generateImages:
        #         #     fluentString += '''/display/views/rest view-0\n/display/objects/display logk\n'''
        #         #     fluentString += '''/display/save-picture ../res/img/{}_alpha_{}_logk.png\n'''.format(self.airfoilName, round(angleList[0][i + 1], 3))
        #         #     fluentString += '''/display/views/rest view-0\n/display/objects/display u\n'''
        #         #     fluentString += '''/display/save-picture ../res/img/{}_alpha_{}_u.png\n'''.format(self.airfoilName, round(angleList[0][i + 1], 3))
        #         #     fluentString += '''/display/views/rest view-0\n/display/objects/display p\n'''
        #         #     fluentString += '''/display/save-picture ../res/img/{}_alpha_{}_p.png\n'''.format(self.airfoilName, round(angleList[0][i + 1], 3))

        #     fluentString += '/file/stop-transcript\n\n'
        #     self.fluentString  = fluentString

        #     with open(self.inputFilename, "w") as f:
        #         f.write(fluentString)
        
    
    # For each object, run the pre-processing
    objects = self.setup["Objects"]
    for obj in objects.keys():
        if self.mode == 0:
            Pre(obj)
        elif self.mode == 0:
            Post(obj)



# SCRATCH
        # self, airfoilName, ambientConditions, L, yplus=0.2, aoaMin=-5.0, aoaMax=10.0, numberOfIncrements=10, inputFilename='',
        # transcriptFilename='', resultFilename='', expData=False, plotFilename='', meshFilename='', experimentalDataFilename=''