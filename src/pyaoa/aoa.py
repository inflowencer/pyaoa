import numpy as np
import pandas as pd
import re
import sys
import os
import yaml
import subprocess
import glob

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
        # Read in yaml setup
        with open(setupFile, "r") as f:
            self.setup = yaml.safe_load(f)
        # Folder related stuff
        self.pre_folder = str(self.setup["I/O"]["pre-folder"])
        self.post_folder = str(self.setup["I/O"]["post-folder"])
        self.run_folder = str(self.setup["I/O"]["run-folder"])
        self.working_dir = self.setup["I/O"]["working-dir"]
        os.chdir(self.working_dir)

        def folderSlash(folder):
            """Appends '/' if folder isn't specified like that"""
            if not folder[-1] == "/":
                folder += "/"
                return folder
            else:
                return folder

        self.pre_folder = folderSlash(self.pre_folder)
        self.post_folder = folderSlash(self.post_folder)
        self.run_folder = folderSlash(self.run_folder)
        self.working_dir = folderSlash(self.working_dir)
        self.operating_system = str(self.setup["I/O"]["OS"]).lower()
        
        # Numerics related
        self.dim = self.setup["Numerics"]["dim"]  # Dimension
        self.solver = self.setup["Numerics"]["solver"].lower()
        self.n_iter = int(str(self.setup["Numerics"]["iter"]))
        self.np = int(str(self.setup["Numerics"]["np"]))
        # I/O related
        self.base_case = self.run_folder + self.setup["I/O"]["base-case"]
        self.run_file = self.working_dir + self.run_folder + self.setup["I/O"]["run-file"]
        self.run_script = self.setup["I/O"]["run-script"]
        self.export_csv = self.setup["I/O"]["export-csv"]
        # Objects and parameter related
        self.objects = self.setup["Objects"]  # Objects
        self.amin = float(self.setup["Parameters"]["amin"])
        self.amax = float(self.setup["Parameters"]["amax"])
        self.inc = int(self.setup["Parameters"]["inc"])
        self.avg = int(self.setup["Parameters"]["avg"])
        # Boundary condition related
        self.inlet_name = str(self.setup["Boundary conditions"]["inlet"]["name"]).lower()
        self.inlet_type = str(self.setup["Boundary conditions"]["inlet"]["type"]).lower()
        try:
            self.I = float(self.setup["Boundary conditions"]["I"])
        except:
            self.I = 0.3

        # Execute the pre- or postprocessing mode
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
            return y_1, delta, l

        # 2. Check if the free-stream is well defined
        if self.setup["Boundary conditions"]["Re"]:
            Re = float(self.setup["Boundary conditions"]["Re"])
            u_inf = Re * self.mu / (L * self.rho)
            Ma = u_inf / self.a
            y_1, delta, l = yPlus(self, Re, u_inf)
        elif self.setup["Boundary conditions"]["u"]:
            u_inf = float(self.setup["Boundary conditions"]["u"])
            Re = u_inf * L * self.rho / self.mu
            Ma = u_inf / self.a
            y_1, delta, l = yPlus(self, Re, u_inf)
        elif self.setup["Boundary conditions"]["Ma"]:
            Ma = float(self.setup["Boundary conditions"]["Ma"])
            u_inf = Ma * self.a
            Re = u_inf * self.rho * L / self.mu
            y_1, delta, l = yPlus(self, Re, u_inf)
        else:
            console.print("\n[bold red]Error: [not bold red]No free-stream Reynolds number, velocity or Mach number defined.")
            sys.exit()

        A = d * L

        return y_1, round(Ma, 3), round(delta, 5), A, round(Re, 1), round(u_inf, 4), L, round(l, 5)


    def ObjAttr(self):
        objDict = {}
        for obj in self.objects.keys():
            y_1, Ma, delta, A, Re, u_inf, L, l = self.objectCalculations(obj)
            objDict[obj] = {
                "y_1": float(y_1), "delta": float(delta), "A": float(A),
                "Re": float(Re), "Ma": float(Ma), "u_inf": float(u_inf),
                "L": float(L), "l": float(l)
            }
        return objDict


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
# l:     Turbulent length scale d99*0.4..(m)
# L:     Object chord length.............(m)
# delta: Boundary layer thickness d99....(m)
# mu:    Free-stream dynamic viscosity...(Pa s)
# nu:    Free-stream kinematic viscosity.(m2 s-1)
# p:     Free-stream pressure............(Pa)
# rho:   Free-stream density.............(kg m-3)
# tau_w: Wall shear-stress...............(N m-2)
# u:     Free-stream velocity............(m s-1)
# u_tau: Shear-velocity..................(m s-1)
# y_1:   First layer height for y+.......(m)
# yplus: Dimensionless wall distance.....(-)""")

        console.print(f"\nPre-calculations written to: \'{exportFile}\'")


    def runFileStr(self, objDict):
        """Create the run file based on the objects and solver"""

        console = Console()

        def fluent(self, objDict):
            """Creates a fluent input string"""
            # Initialize all angles as array of angles
            angles = np.linspace(self.amin, self.amax, self.inc + 1, retstep=False)


            # 0. Initialize the input str and start the transcript
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
; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ;
; Reading in Base Case
/file/read-case {self.base_case}
""" # sync-chdir ../
            # 1. Loop over all objects
            for objName, objAttr in objDict.items():

                # 1.1 Make sure folder exists and set the correct file-name for the output of lift and drag
                os.system(f"mkdir -p {self.post_folder + objName}")

                # 1.2 Calculate dict of ux and uy components
                angleDict = {}
                for angle in angles:
                    angle = round(float(angle), 4)
                    # x-component
                    ux = round(np.cos(np.deg2rad(angle)) * objAttr["u_inf"], 4)
                    Fx = round(np.cos(np.deg2rad(angle)), 4)  # Force components
                    # y-component
                    uy = round(np.sin(np.deg2rad(angle)) * objAttr["u_inf"], 4)
                    Fy = round(np.sin(np.deg2rad(angle)), 4)  # Force components
                    angleDict[str(angle)] = {"ux": ux, "Fx": Fx, "uy": uy, "Fy": Fy}

                # 1.2 Print the current object that will be analyzed
                inputStr += f"""\n
; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ;
;              Analysing object `{objName}`
; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ;
; Opening transcript...
/file/start-transcript run/{objName}.out"""

                # 1.3 Read in the mesh of that object
                mesh = self.setup["Objects"][objName]["mesh"]
                inputStr += f"""
; Reading in mesh from {mesh}
/file/replace-mesh {mesh} o ()
"""

                # 1.4 Loop over angleDict angles
                for angle, components in angleDict.items():
                    # 1.4.1 Set the correct boundary condition according to the calculations for the current angle
                    if self.inlet_type == "velocity-inlet" or "velocity inlet":
# ; Setting the velocity components for angle {angle} deg to [ux: {components["ux"]}, uy: {components["uy"]}]
                        inputStr += f"""
; Angle: {angle} deg
(display "Angle is: {angle} deg")
/define/bound/set/velocity/{self.inlet_name} () dir-0 no {components["ux"]} dir-1 no {components["uy"]} ke-spec no yes turb-len {objAttr["l"]} turb-int {self.I} ()
"""
                    # TODO: Farfield
                    # elif self.inlet_type == "pressure-far-field" or "pressure far field":

# ; Setting the output folder name for angle {angle} deg to "{self.post_folder}{objName}"
# ; Setting the correct velocity vector components
                    inputStr += f"""\
/solve/report-files/edit drag file {self.post_folder}{objName}/alpha_{angle}_drag.out ()
/solve/report-definitions/edit/drag force-vector {components["Fx"]} {components["Fy"]} ()
/solve/report-files/edit lift file {self.post_folder}{objName}/alpha_{angle}_lift.out ()
/solve/report-definitions/edit/lift force-vector {components["Fy"]} {components["Fx"]} ()
"""
                    # 1.5 Initialize the solution
                    inputStr += f"""\
/solve/init/hyb o ()\n/solve/it {self.n_iter}\n"""

                # 1.6 Stop transcript for this object
                inputStr += f"""
; Stopping transcript...
/file/stop-transcript\n"""
            return inputStr

        # def openfoam(self, obj):
        # def su2(self, obj):

        if self.solver == "fluent":
            inputStr = fluent(self, objDict)
        # elif self.solver == "openfoam":
        #     openfoam()
        # elif self.solver == "su2":
        #     su2()
        else:
            console.print("[bold red]Error: [not bold bright_white]Wrong solver specified in \"Numerics -> Solver\".\nAllowed solvers: [fluent, openfoam, su2]")

        return inputStr




    def Pre(self):
        """Preprocessing method"""

        console = Console()

        # 1. Determine ambient conditions which return self.ambientDict
        self.ambientConditions()

        # 2. Create objects to analyze, calculate their specific attributes
        #    and generate the runFile String

        objDict = self.ObjAttr()

        inputStr = self.runFileStr(objDict)
        
        # 3. Export calculations if specified
        try:
            exportFolder = self.pre_folder
        except:
            exportFile = f"pre/calculations.yaml"
            self.exportCalculations(exportFile, objDict)
        else:
            if exportFolder:
                exportFile = exportFolder + "calculations.yaml"
                self.exportCalculations(exportFile, objDict)

        # 4. Export input str
        with open(self.run_file, "w") as f:
            f.write(inputStr)
            console.print(f"Input file written to:       \'{self.run_file}\'")

        # 5. Input File
        if self.run_script:
            self.createRunScript()

        
    def clearFluentFiles(self):
        """Clears all Fluent files that have been written"""
        # get a recursive list of file paths that matches pattern including sub directories
        trn_files = glob.glob(f"{self.working_dir}/**/*.trn", recursive=True)
        log_files = glob.glob(f"{self.working_dir}/**/*.log", recursive=True)
        hist_files = glob.glob(f"{self.working_dir}/**/*.hist", recursive=True)

        fileList = trn_files + log_files + hist_files
        # Iterate over the list of filepaths & remove each file.
        for filePath in fileList:
            try:
                os.remove(filePath)
            except:
                pass
        

    def clearResults(self):
        previousResults = False
        for obj in self.objects.keys():
            # First check if the path already exists and if yes check if its empty
            objResFolder = f"{self.post_folder}{obj}"
            if os.path.exists(objResFolder):
                contents = os.listdir(objResFolder)
                if contents:
                    previousResults = True
                    break
        # If results have been found, clear them 
        if contents:
            self.clearRes = Prompt.ask("[bright_white]Previous results have been found, would you like to clear them?", choices=["y", "n"], default="y")
            if self.clearRes:
                for obj in self.objects.keys():
                    try:
                        os.rmdir(f"post/{obj}")
                    except:
                        pass
                    try:
                        os.remove(f"run/{obj}.hist")
                    except:
                        pass
        

    def createRunScript(self):
        """Creates a Run script to automatically run the simulation depending on the platform"""
        console = Console()
        if self.operating_system == "windows":
            result = subprocess.run(['wslpath', '-w', self.working_dir], stdout=subprocess.PIPE)
            win_path = result.stdout[:-1]
            win_path = win_path.decode("UTF-8")
            self.script_path = f"{self.run_folder}run_script.bat"
            runScript = f"cd \"{win_path}\"\n"
            jouFile = self.run_folder + self.setup["I/O"]["run-file"]
            if str(self.dim.lower()) == "2d" or "2":
                runScript += f"fluent 2ddp -g -t{self.np} < {jouFile}"
            else:
                runScript += f"fluent 3ddp -g -t{self.np} < {jouFile}"
            with open(self.script_path, "w") as f:
                f.write(runScript)
            console.print(f"Run script written to:       \'{self.script_path}\'")


    def Plot(self):
        """Method to plot depending on the setup.yaml"""
        console = Console()
        # 1. Try to read the objects which shall be plotted otherwise use the objects to analyze
        try:
            self.plotObjList = list(self.setup["Plot"].keys())
        except:
            self.plotObjList = list(self.setup["Objects"].keys())
        else:
            if len(self.plotObjList) == 0:
                self.plotObjList = list(self.setup["Objects"].keys())

        # 2. Try to determine the type of plot otherwise set default
        try:
            self.plot_type = self.setup["Plot"]["type"]
        except:
            self.plot_type = "both"

        # 3. Determine the layout
        try:
            self.plot_layout = self.setup["Plot"]["layout"]
        except:
            self.plot_layout = "side-to-side"

        # 4. Determine the export file type
        try:
            self.plot_ftype = self.setup["Plot"]["export"]
        except:
            self.plot_ftype = "png"

        # 5. Functions for plotting
        def alpha_cl(self):
            """Plotting AOA-Cl"""


        # 8. Initialize the plot based on the settings
        fig, axs = plt.subplots(1, 2, figsize=(12, 3), sharey=True)
        # if self.plot_type.lower() == "both":
        #     if self.plot_layout.lower() == "side-to-side":
        #     else:


        # 10. Loop over all objects
        for obj in self.plotObjList:
            # 10.1 Determine if there's reference data available
            try:
                ref_data = self.setup["Plot"][obj]["ref-data"]
            except:
                ref_data = False
            else:
                df_ref = pd.read_csv(ref_data, sep=",")
            objResFolder = f"{self.post_folder}{obj}"
            if not objResFolder[-1] == "/":
                objResFolder += "/"
            if not os.path.exists(objResFolder):
                console.print(f"[bold red]Error:[not bold bright_white] Result folder for {obj} not found...")
                sys.exit()

            # Read in CSV
            df = pd.read_csv(f"{objResFolder}{obj}.csv", sep=",")

            # Cl/AOA polar Layout
            axs[0].axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
            axs[0].axvline(x=0, color='gray', linestyle='-', linewidth=0.5)
            axs[0].set_xlabel(r'Angle of attack $\alpha$')
            axs[0].set_ylabel(r'Lift coeffient $c_l$')
            axs[0].grid(True, which='both', axis='both', linewidth=0.1, color='grey')
            axs[0].plot(df['alpha'], df['cl'], label=f"{obj}")
            if ref_data:
                axs[0].scatter(df_ref['alpha'], df_ref['cl'], color='red', marker='+', label=f"{obj} ref.")

            # Cl/Cd (Lilienthal) Layout
            axs[1].axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
            axs[1].set_xlabel(r'Drag coefficient $c_d$')
            axs[1].grid(True, which='both', axis='both', linewidth=0.1, color='grey')

        leg = plt.legend()
    

#     # leg.get_lines()[0].set_linewidth(0.2)
#     for airfoil in airfoils:
#         # If not data has been read in already
#         # if airfoil.df_sim.empty():
#         #     airfoil.df_sim = pd.read_csv(f'airfoil-data/{airfoil.airfoilName}_sim.csv', sep=',')
#         if airfoil.expData:
#             airfoil.df_exp = pd.read_csv("airfoil-data/" + airfoil.airfoilName + "_exp.csv", sep=',')
#         if plotType == 'polar-and-lilienthal':
#             # Polar
#             if airfoil.expData:
#             # Lilienthal
#             axs[1].plot(airfoil.df_sim['cd'], airfoil.df_sim['cl'], label='{}'.format(airfoil.airfoilName)) # Simulation data
#             if airfoil.expData:
#                 axs[1].scatter(airfoil.df_exp['cd'], airfoil.df_exp['cl'], color='red', marker='+', label=f'{airfoil.airfoilName} Experiment')  # Experimental data
#     # leg = axs[1].legend(facecolor='white', frameon=True, framealpha=1)
#     leg = plt.legend(bbox_to_anchor=(-1, -0.3), loc="lower left", ncol=3, prop={'size': 8})  # bbox_transform=fig.transFigure
#     fig.savefig(plotFilename, dpi=300)
        # plt.plot(clean_df["cd"], clean_df["cl"])
        # plt.grid()
        # plt.savefig("test.png")
    
            

    def Post(self):
        """Post-processing method"""
        console = Console()
        objResDict = {}
        # First prepare the raw output data depending on the solver
        if self.solver == "fluent":
            for obj in self.objects:
                objResFolder = f"{self.post_folder}{obj}"
                if not objResFolder[-1] == "/":
                    objResFolder += "/"
                if not os.path.exists(objResFolder):
                    console.print(f"[bold red]Error:[not bold bright_white] Result folder for {obj} not found...")
                    sys.exit()
                files = os.listdir(objResFolder)
                if not files:
                    console.print(f"[bold red]Error:[not bold bright_white] Result files for {obj} not found...")
                    sys.exit()

                def fluentToDataFrame(self, raw_file, quantity):
                    """Converts raw fluent .out files to a dataframe"""

                    df = pd.read_csv(f"{objResFolder}{raw_file}", sep=" ", skiprows=2)

                    if len(df.columns) == 2:
                        df.columns = ["iter", f"{quantity}"]
                    elif len(df.columns) > 2:
                        df.columns = ["iter", f"{quantity}", f"{quantity}_inst"]

                    avg = np.mean(df[f"{quantity}"].iloc[-self.avg:].values)

                    return df, avg

                # Initialize postDict for all objects
                # Each object in the objResDict() has a key for each AOA and a corresponding df for lift and drag
                AOA_dict = {}

                # 1. Get all AOAs
                AOA_list = []
                for f in files:
                    num_list = re.findall('-?\d+\.?\d*', f)
                    aoa = num_list[0]
                    if not float(aoa) in AOA_list:
                        AOA_list.append(float(aoa))
                AOA_list = sorted(AOA_list)

                # 2. Loop over each angle and append to AOA_dict of this object
                for angle in AOA_list:
                    for f in files:
                        if str(angle) in f:
                            if not angle in AOA_dict.keys():
                                AOA_dict[angle] = {"drag": {"df": 0, "avg": 0}, "lift": {"df": 0, "avg": 0}}
                            if "drag" in f.lower():
                                df_drag, avg_drag = fluentToDataFrame(self, f, "drag")
                                AOA_dict[angle]["drag"] = {"df": df_drag, "avg": avg_drag}
                            elif "lift" in f.lower():
                                df_lift, avg_lift = fluentToDataFrame(self, f, "lift")
                                AOA_dict[angle]["lift"] = {"df": df_lift, "avg": avg_lift}


                # 3. Initialize a clean df which we can later export as a csv
                clean_drag = []
                clean_lift = []
                for angle in AOA_list:
                    clean_drag.append(AOA_dict[angle]["drag"]["avg"])
                    clean_lift.append(AOA_dict[angle]["lift"]["avg"])
                clean_df = pd.DataFrame({"alpha": AOA_list, "drag_force": clean_drag, "lift_force": clean_lift})

                # 4. Calculate cd and cl based on the objects' properties
                self.ambientConditions()
                objDict = self.ObjAttr()
                clean_df["cd"] = clean_df["drag_force"] / (0.5 * self.ambientDict["rho"] * objDict[obj]["u_inf"]**2 * objDict[obj]["A"])
                clean_df["cl"] = clean_df["lift_force"] / (0.5 * self.ambientDict["rho"] * objDict[obj]["u_inf"]**2 * objDict[obj]["A"])

                if self.export_csv == True:
                    clean_df.to_csv(f"{objResFolder}{obj}.csv", sep=',', header=True, index=False)
                # console.print(clean_df)
                # 4. Append DF objResDict
                # objResDict[obj] = clean_df
                # console.print


# SCRATCH
        # self, airfoilName, ambientConditions, L, yplus=0.2, aoaMin=-5.0, aoaMax=10.0, numberOfIncrements=10, inputFilename='',
        # transcriptFilename='', resultFilename='', expData=False, plotFilename='', meshFilename='', experimentalDataFilename=''

# else:
#     quantity = Prompt.ask("[bright_white]Couldn't determine which quantity should be analyzed.\n[green]Please enter the quantity name")
#     df = fluentToDataFrame(self, f, quantity)
#     files.remove(_)

# def multiAirfoilPlot(airfoils, plotType='polar-and-lilienthal', plotLimits={}, plotFilename='plots/multiAirfoilComparison.pdf'):
#     '''Generate plots from airfoil data'''
#     '''Generates and exports LaTeX figures from the airfoil data.

#     Parameters
#     ----------
#     plotType : \'polar-and-lilienthal\', \'polar\', \'lilienthal\'
#     exportPlot : True
#     expDataAvailable : True
#     plotLimits : dict(xmin_aoa=float, xmax_aoa=float, xmin_cd=float, xmax_cd=float, ymin_cl=float, ymax_cl=float)
#     '''
#     # Initialize plot
#     fig, axs = plt.subplots(1, 2, figsize=(12, 3), sharey=True)
