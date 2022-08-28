import numpy as np
import pandas as pd
import re
import sys
import os
import yaml

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
            self.u_b = self.Re * self.nu / self.L  # free stream velocity acc. to Re
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

        def exportCalculations(self):
            """Function to export calculations for the mesh and free-stream"""
            try:
                exportFile = self.setup["I/O"]["pre-calculations"]
            except:
                pass
            else:
                calculations_dict = {
                    "p_inf": self.p, "T_inf": self.T, "Object width d": self.d, "rho_inf": self.rho,

                }



# SCRATCH
        # self, airfoilName, ambientConditions, L, yplus=0.2, aoaMin=-5.0, aoaMax=10.0, numberOfIncrements=10, inputFilename='',
        # transcriptFilename='', resultFilename='', expData=False, plotFilename='', meshFilename='', experimentalDataFilename=''