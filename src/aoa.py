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


    def pre(self, obj):
        """Preprocessing method"""

        console = Console()

        # Check pressure and temperature. If not given, use standard conditions
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
        self.L = self.setup["Objects"][obj]["L"]
        self.a = np.sqrt(1.4 * 8.314 / 0.029 * self.T)

        def yPlus(self):
            self.C_f = 0.0576 * self.Re**(-0.2)  # skin friction coefficient
            self.tau_w = 0.5 * self.C_f * self.rho * self.u_inf**2 # wall shear stress
            self.u_tau = np.sqrt(self.tau_w / self.rho) # friction velocity
            self.yplus = 1.0
            self.y_1 = self.yplus * self.nu / self.u_tau # first layer height

        # First check if the free-stream is well defined
        if self.setup["Boundary conditions"]["Re"]:
            self.Re = float(self.setup["Boundary conditions"]["Re"])
            self.u_inf = self.Re * self.mu / (self.L * self.rho)
            self.Ma = self.u_inf / self.a
            yPlus(self)
        elif self.setup["Boundary conditions"]["u"]:
            self.u_inf = float(self.setup["Boundary conditions"]["u"])
            self.Re = self.u_inf * self.rho * self.rho / self.mu
            self.Ma = self.u_inf / self.a
            yPlus(self)
        elif self.setup["Boundary conditions"]["Ma"]:
            self.Ma = float(self.setup["Boundary conditions"]["Ma"])
            self.u_inf = self.Ma * self.a
            self.Re = self.u_inf * self.rho * self.rho / self.mu
            yPlus(self)
        else:
            console.print("\n[bold red]Error: [not bold red]No free-stream Reynolds number, velocity or Mach number defined.")
            sys.exit()

        def exportCalculations(self, exportFile):
            """Function to export calculations for the mesh and free-stream"""
            calculations_dict = {
                "p": float(self.p), "T": float(self.T), "Object width d": float(self.d), "rho": float(self.rho),
                "Re": float(self.Re), "u": float(self.u_inf), "Ma": float(self.Ma), "mu": float(self.mu), "nu":
                float(self.nu), "C_f": float(self.C_f), "tau_w": float(self.tau_w), "u_tau": float(self.u_tau),
                "yplus": float(self.yplus), "y_1": float(self.y_1)
            }
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
        try:
            exportFile = self.setup["I/O"]["pre-calculations"]
        except:
            pass
        else:
            exportCalculations(self, exportFile)



# SCRATCH
        # self, airfoilName, ambientConditions, L, yplus=0.2, aoaMin=-5.0, aoaMax=10.0, numberOfIncrements=10, inputFilename='',
        # transcriptFilename='', resultFilename='', expData=False, plotFilename='', meshFilename='', experimentalDataFilename=''