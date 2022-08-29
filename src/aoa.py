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

        self.mode = Prompt.ask("\n[bright_white]Welcome to the Angle-of-Attack Analysis Tool![not bold bright_white]\nSelect the mode  [bold cyan]0: Preprocessing    [bold yellow]1: Postprocessing\n", choices=["0", "1"], default="0")
        with open(setupFile, "r") as f:
            self.setup = yaml.safe_load(f)
        self.objects = self.setup["Objects"]  # Objects
        self.dim = self.setup["Numerics"]["dim"]  # Dimension
        self.solver = self.setup["Numerics"]["solver"].lower()
        self.base_case = self.setup["Numerics"]["base-case"].lower()
        self.run_file = self.setup["Numerics"]["run-file"].lower()
        self.amin = float(self.setup["Parameters"]["amin"])
        self.amax = float(self.setup["Parameters"]["amax"])
        self.inc = int(self.setup["Parameters"]["inc"])
        self.inlet_name = str(self.setup["Boundary conditions"]["inlet"]["name"]).lower()
        self.inlet_type = str(self.setup["Boundary conditions"]["inlet"]["type"]).lower()

        # Preprocessing Mode
        if int(self.mode) == 0:
            self.Pre()
        elif int(self.mode) == 1:
            self.Post()
    

    def ambientConditions(self):
        """Method to calculate and create ambient_dict attribute for free-stream boundary conditions."""
        # 1. Check pressure and temperature. If not given, use standard conditions
        try:
            self.p = float(self.setup["Boundary conditions"]["p"])
        except:
            self.p = 1e5
        try:
            self.T = float(self.setup["Boundary conditions"]["T"])
        except:
            self.T = 298.15
        # TODO FIX depth
        # try:
        #     self.d = float(self.setup["Boundary conditions"]["d"])
        # except:
        #     self.d = 1.0

        self.rho = 0.029 / 8.314 * self.p / self.T # Density according to material data
        # Calculate viscosity according to Sutherland
        self.mu = 1.716e-5 * (self.T / 273.15)**1.5 * (273.15 + 110.4) / (self.T + 110.4)
        self.nu = self.mu / self.rho
        self.a = np.sqrt(1.4 * 8.314 / 0.029 * self.T)

        self.ambientDict = {
            "p": float(self.p), "T": float(self.T), "rho": float(self.rho),
            "mu": float(self.mu), "nu": float(self.nu), "a": float(self.a)
        }


    def objectCalculations(self, obj):
        """Function to calculate and return ambient condition calculations for free-stream"""
        console = Console()
        L = float(self.setup["Objects"][obj]["L"])
        if str(self.dim.lower()) == "2d" or "2":
            d = float(self.setup["Boundary conditions"]["d"])
        elif str(self.dim.lower()) == "3d" or "3":
            d = float(self.setup["Objects"][obj]["d"])
        else:
            d = 1.0


        # TODO: TRANSONIC AND HIGHER MACH NUMBER FLOWS
        def yPlus(self, Re, u_inf):
            C_f = 0.0576 * Re**(-0.2)  # skin friction coefficient
            tau_w = 0.5 * C_f * self.rho * u_inf**2 # wall shear stress
            u_tau = np.sqrt(tau_w / self.rho) # friction velocity
            yplus = 1.0
            y_1 = yplus * self.nu / u_tau # first layer height
            delta = 0.37 * L * Re**(-0.2)
            l = 0.4 * delta
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

        A = d * L

        return y_1, Ma, delta, A, Re, u_inf, L, l


    def exportCalculations(self, exportFile, objDict):
        """Method to export the ambient conditions and object attributes."""

        console = Console()

        # 0. Write header
        with open(exportFile, "w") as f:
            f.write("""
# - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - #
#                      AMBIENT CONDITIONS                    #
# - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - #
""")
        # 1. Write the ambient dict
        with open(exportFile, "a") as f:
            yaml.safe_dump(self.ambientDict, f)
        # 2. Write Obj header
        with open(exportFile, "a") as f:
            f.write("""
# - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - #
#                      OBJECTS TO ANALYZE                    #
# - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - #
""")
        # 3. Write Obj dict
        with open(exportFile, "a") as f:
            yaml.safe_dump(objDict, f)
        
        # 4. Write Explanations
        with open(exportFile, "a") as f:
            f.write("""\n
# - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - #
#                         ABBREVIATIONS                      #
# - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - #
# C_f:   Skin friction coefficient.......(-)
# Ma:    Mach number.....................(-)
# d:     Object depth....................(m)
# A:     Object area.....................(m2)
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
# yplus: Dimensionless wall distance.....(-)""")

        console.print(f"\nPre-calculations written to: \'{exportFile}\'\n")


    def runFileStr(self, objDict):
        """Create the run file based on the objects and solver"""

        console = Console()

        def fluent(self, objDict):
            """Creates a fluent input string"""
            # Initialize all angles as array of angles
            angles = np.linspace(self.amin, self.amax, self.inc + 1, retstep=True)


            # 0. Initialize the input str
            inputStr = f"""\
; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ;
;                        AoA Analysis                       ;
; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ;
; Parameters:
; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ;
; alpha_min: {self.amin}
; alpha_max: {self.amax}
; increments: {self.inc}

; Total no. of simulations: {int(len(objDict) * len(angles))}
; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ;
;                    Beginning Analysis...                  ;
; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ;\n
"""
            # 1. Loop over all objects
            for objName, objAttr in objDict.items():

                # 1.1 Calculate dict of ux and uy components
                angleDict = {}
                for angle in angles:
                    ux = np.cos(np.deg2rad(5)) * objAttr["u_inf"]
                    uy = np.sin(np.deg2rad(5)) * objAttr["u_inf"]
                    angleDict[angle] = {"ux": ux, "uy": uy}

                # 1.2 Print the current object that will be analyzed
                inputStr += f"""
; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ;
;     Analysing object `{objName}`                 ;
; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ;"""

                # 1.3 Read in the mesh of that object
                mesh = self.setup["Objects"][objName]["mesh"]
                inputStr += f"""
; Reading in mesh from {mesh}
/file/replace-mesh {mesh} o
"""

                # 1.4 Loop over angleDict angles
                for angle, components in angleDict.items:
                    # 1.4.1 Set the correct boundary condition according to the calculations for the current angle
                    if self.inlet_type == "velocity-inlet" or "velocity inlet":
                        inputStr += f"""
; Setting 
/define/bound/set/velocity/{self.inlet_name}} () dir-0 no {components["ux"]} no dir-1 {components["uy"]} 
"""
                    # TODO: Farfield
                    # elif self.inlet_type == "pressure-far-field" or "pressure far field":
                    #     pass

            # Start transcript in order for the console output to be written to a file
            # fluentString = ';{}\n/file/start-transcript ../{} o\n\n'.format(self.airfoilName, self.transcriptFilename)
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
        # def openfoam(self, obj):
        # def su2(self, obj):

        if self.solver == "fluent":
            inputStr = fluent(objDict)
        # elif self.solver == "openfoam":
        #     openfoam()
        # elif self.solver == "su2":
        #     su2()
        else:
            console.print("[bold red]Error: [not bold bright_white]Wrong solver specified in \"Numerics -> Solver\".\nAllowed solvers: [fluent, openfoam, su2]")

        return inputStr


    def Pre(self):
        """Preprocessing method"""

        # 1. Determine ambient conditions which return self.ambientDict
        self.ambientConditions()

        # 2. Create objects to analyze, calculate their specific attributes
        #    and generate the runFile String
        objDict = {}
        for obj in self.objects.keys():
            y_1, Ma, delta, A, Re, u_inf, L, l = self.objectCalculations(obj)
            objDict[obj] = {
                "y_1": float(y_1), "delta": float(delta), "A": float(A),
                "Re": float(Re), "Ma": float(Ma), "u_inf": float(u_inf),
                "L": float(L), "l": float(l)
            }

        inputStr = self.runFileStr(objDict)
        
        # 3. Export calculations if specified
        try:
            exportFile = self.setup["I/O"]["pre-calculations"]
        except:
            pass
        else:
            if exportFile:
                self.exportCalculations(exportFile, objDict)

        # 4. Export input str
        with open(self.run_file, "w") as f:
            f.write(inputStr)
        
    


# SCRATCH
        # self, airfoilName, ambientConditions, L, yplus=0.2, aoaMin=-5.0, aoaMax=10.0, numberOfIncrements=10, inputFilename='',
        # transcriptFilename='', resultFilename='', expData=False, plotFilename='', meshFilename='', experimentalDataFilename=''