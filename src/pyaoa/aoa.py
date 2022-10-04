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
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
plt.rcParams.update({'figure.max_open_warning': 0})

# import importlib.resources # For importing matplotlib settings

# with importlib.resources.open_text("pyaoa", "matplotlib.yml") as f:
#     data = yaml.safe_load(f)  


class Analysis:
    """Angle-of-Attack Analysis class for Pre- and Postprocessing"""


    def __init__(self, setupFile):

        console = Console()

        self.mode = Prompt.ask("""\n[bright_white]Welcome to the Angle-of-Attack Analysis Tool![not bold bright_white]
Select the mode  [bold cyan]0: Preprocess   [bold yellow]1: Postprocess   [bold green]2: Plot\n""", choices=["0", "1", "2"], default="0")
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
        elif int(self.mode) == 2:
            self.Plot()

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
/file/start-transcript run/{objName}.hist o ()"""

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


    def clearResults(self):
        previousResults = False
        for obj in self.objects.keys():
            # First check if the path already exists and if yes check if its empty
            objResFolder = f"{self.post_folder}{obj}"
            if os.path.exists(objResFolder):
                contents = os.listdir(objResFolder)
                if len(contents) > 0:
                    self.clearRes = Prompt.ask("[bright_white]Previous results have been found, would you like to clear them?", choices=["y", "n"], default="n")
                    if self.clearRes == "y":
                        for item in contents:
                            try:
                                check_aoa = re.findall('-?\d+\.?\d*', item)[0]
                            except:
                                check_aoa = False
                            if "alpha" in item and check_aoa and not ".csv" in item and not ".pdf" in item:
                                os.remove(os.path.join(objResFolder, item))
                        try:
                            os.remove(f"{self.run_folder}{obj}.out")
                        except:
                            pass

        
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
                    if ".out" in f:
                        num_list = re.findall('-?\d+\.?\d*', f)
                        aoa = num_list[0]
                        if not float(aoa) in AOA_list:
                            AOA_list.append(float(aoa))
                AOA_list = sorted(AOA_list)

                # 2. Loop over each angle and append to AOA_dict of this object
                for angle in AOA_list:
                    for f in files:
                        try:
                            check_aoa = re.findall('-?\d+\.?\d*', f)[0]
                        except:
                            check_aoa = 1e10
                        if str(angle) in f and str(angle) == str(check_aoa) and "alpha" in f:
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


    def Plot(self):
        """Method to plot depending on the setup.yaml"""
        console = Console()

        if self.mode == 0:
            self.plot = Prompt.ask("[yellow]Current mode was set as `pre`.[bright_white]Would you still like to plot?", choices=["y", "n"], default="y")
            if self.clearRes == "y":
                pass
            else:
                sys.exit()

        # 1. Try to read the objects which shall be plotted otherwise use the objects to analyze
        self.plotObjList = list()
        for obj in self.setup["Objects"].keys():
            try:
                enable_plot = self.setup["Objects"][obj]["plot"]["enable"]
            except:
                pass
            else:
                if enable_plot:
                    self.plotObjList.append(obj)
        
        # 2. Read the general Plot settings
        try:
            self.plot_file_type = self.setup["Plot"]["export"]
        except:
            self.plot_file_type = "pdf"
        try:
            self.plot_grid = self.setup["Plot"]["grid"]["enable"]
        except:
            self.plot_grid = False
        try:
            self.legend = self.setup["Plot"]["legend"]
        except:
            self.legend = False
        

        # - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - #
        #                          D R A G                           #
        # - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - #
        # 3. Read the drag plot settings
        try:
            self.plot_drag_x_label = self.setup["Plot"]["drag"]["x-label"]
        except:
            self.plot_drag_x_label = r"Angle of attack $\alpha$"
        try:
            self.plot_drag_y_label = self.setup["Plot"]["drag"]["y-label"]
        except:
            self.plot_drag_y_label = r"Drag coefficient $c_d$"
        try:
            self.plot_drag_xlims = self.setup["Plot"]["drag"]["x-lims"]
        except:
            self.plot_drag_xlims = False
        try:
            self.plot_drag_ylims = self.setup["Plot"]["drag"]["y-lims"]
        except:
            self.plot_drag_ylims = False
        try:
            self.plot_drag_x_tick_major = self.setup["Plot"]["drag"]["x-ticks"]["major"]
        except:
            self.plot_drag_x_tick_major = False
        try:
            self.plot_drag_x_tick_minor = self.setup["Plot"]["drag"]["x-ticks"]["minor"]
        except:
            self.plot_drag_x_tick_minor = False
        try:
            self.plot_drag_y_tick_major = self.setup["Plot"]["drag"]["y-ticks"]["major"]
        except:
            self.plot_drag_y_tick_major = False
        try:
            self.plot_drag_y_tick_minor = self.setup["Plot"]["drag"]["y-ticks"]["minor"]
        except:
            self.plot_drag_y_tick_minor = False
        # 3.1 Initialize Drag plot for each obj and a global one
        # Obj
        def dragPlotLayout(self, mode):
            fig_cd, axs_cd = plt.subplots(1, 1, figsize=(3, 3))
            axs_cd.axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
            axs_cd.axvline(x=0, color='gray', linestyle='-', linewidth=0.5)
            axs_cd.set_xlabel(self.plot_drag_x_label)
            axs_cd.set_ylabel(self.plot_drag_y_label)
            # Global
            fig_cd_global, axs_cd_global = plt.subplots(1, 1, figsize=(3, 3))
            axs_cd_global.axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
            axs_cd_global.axvline(x=0, color='gray', linestyle='-', linewidth=0.5)
            axs_cd_global.set_xlabel(self.plot_drag_x_label)
            axs_cd_global.set_ylabel(self.plot_drag_y_label)
            if self.plot_grid:
                axs_cd.grid(True, which='major', axis='both', linewidth=0.1, color='grey')
                axs_cd_global.grid(True, which='major', axis='both', linewidth=0.1, color='grey')
            if self.plot_drag_xlims:
                axs_cd.set_xlim(self.plot_drag_xlims[0], self.plot_drag_xlims[1])
                axs_cd_global.set_xlim(self.plot_drag_xlims[0], self.plot_drag_xlims[1])
            if self.plot_drag_ylims:
                axs_cd.set_ylim(self.plot_drag_ylims[0], self.plot_drag_ylims[1])
                axs_cd_global.set_ylim(self.plot_drag_ylims[0], self.plot_drag_ylims[1])
            if self.plot_drag_x_tick_major:
                axs_cd.xaxis.set_major_locator(MultipleLocator(self.plot_drag_x_tick_major))
                axs_cd_global.xaxis.set_major_locator(MultipleLocator(self.plot_drag_x_tick_major))
                # ax.xaxis.set_major_formatter('{x:.0f}')
            if self.plot_drag_x_tick_minor:
                axs_cd.xaxis.set_minor_locator(MultipleLocator(self.plot_drag_x_tick_minor))
                axs_cd_global.xaxis.set_minor_locator(MultipleLocator(self.plot_drag_x_tick_minor))
            if self.plot_drag_y_tick_major:
                axs_cd.yaxis.set_major_locator(MultipleLocator(self.plot_drag_y_tick_major))
                axs_cd_global.yaxis.set_major_locator(MultipleLocator(self.plot_drag_y_tick_major))
                # ax.xaxis.set_major_formatter('{x:.0f}')
            if self.plot_drag_y_tick_minor:
                axs_cd.yaxis.set_minor_locator(MultipleLocator(self.plot_drag_y_tick_minor))
                axs_cd_global.yaxis.set_minor_locator(MultipleLocator(self.plot_drag_y_tick_minor))

            if mode == "local":
                return fig_cd, axs_cd
            elif mode == "global":
                return fig_cd_global, axs_cd_global

        # - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - #
        #                          L I F T                           #
        # - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - #
        # 3. Read the lift plot settings
        try:
            self.plot_lift_x_label = self.setup["Plot"]["lift"]["x-label"]
        except:
            self.plot_lift_x_label = r"Angle of attack $\alpha$"
        try:
            self.plot_lift_y_label = self.setup["Plot"]["lift"]["y-label"]
        except:
            self.plot_lift_y_label = r"lift coefficient $c_d$"
        try:
            self.plot_lift_xlims = self.setup["Plot"]["lift"]["x-lims"]
        except:
            self.plot_lift_xlims = False
        try:
            self.plot_lift_ylims = self.setup["Plot"]["lift"]["y-lims"]
        except:
            self.plot_lift_ylims = False
        try:
            self.plot_lift_x_tick_major = self.setup["Plot"]["lift"]["x-ticks"]["major"]
        except:
            self.plot_lift_x_tick_major = False
        try:
            self.plot_lift_x_tick_minor = self.setup["Plot"]["lift"]["x-ticks"]["minor"]
        except:
            self.plot_lift_x_tick_minor = False
        try:
            self.plot_lift_y_tick_major = self.setup["Plot"]["lift"]["y-ticks"]["major"]
        except:
            self.plot_lift_y_tick_major = False
        try:
            self.plot_lift_y_tick_minor = self.setup["Plot"]["lift"]["y-ticks"]["minor"]
        except:
            self.plot_lift_y_tick_minor = False
        # 3.1 Initialize lift plot for each obj and a global one
        def liftPlotLayout(mode):
            # Obj
            fig_cl, axs_cl = plt.subplots(1, 1, figsize=(3, 3))
            axs_cl.axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
            axs_cl.axvline(x=0, color='gray', linestyle='-', linewidth=0.5)
            axs_cl.set_xlabel(self.plot_lift_x_label)
            axs_cl.set_ylabel(self.plot_lift_y_label)
            # Global
            fig_cl_global, axs_cl_global = plt.subplots(1, 1, figsize=(3, 3))
            axs_cl_global.axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
            axs_cl_global.axvline(x=0, color='gray', linestyle='-', linewidth=0.5)
            axs_cl_global.set_xlabel(self.plot_lift_x_label)
            axs_cl_global.set_ylabel(self.plot_lift_y_label)
            if self.plot_grid:
                axs_cl.grid(True, which='major', axis='both', linewidth=0.1, color='grey')
                axs_cl_global.grid(True, which='major', axis='both', linewidth=0.1, color='grey')
            if self.plot_lift_xlims:
                axs_cl.set_xlim(self.plot_lift_xlims[0], self.plot_lift_xlims[1])
                axs_cl_global.set_xlim(self.plot_lift_xlims[0], self.plot_lift_xlims[1])
            if self.plot_lift_ylims:
                axs_cl.set_ylim(self.plot_lift_ylims[0], self.plot_lift_ylims[1])
                axs_cl_global.set_ylim(self.plot_lift_ylims[0], self.plot_lift_ylims[1])
            if self.plot_lift_x_tick_major:
                axs_cl.xaxis.set_major_locator(MultipleLocator(self.plot_lift_x_tick_major))
                axs_cl_global.xaxis.set_major_locator(MultipleLocator(self.plot_lift_x_tick_major))
                # ax.xaxis.set_major_formatter('{x:.0f}')
            if self.plot_lift_x_tick_minor:
                axs_cl.xaxis.set_minor_locator(MultipleLocator(self.plot_lift_x_tick_minor))
                axs_cl_global.xaxis.set_minor_locator(MultipleLocator(self.plot_lift_x_tick_minor))
            if self.plot_lift_y_tick_major:
                axs_cl.yaxis.set_major_locator(MultipleLocator(self.plot_lift_y_tick_major))
                axs_cl_global.yaxis.set_major_locator(MultipleLocator(self.plot_lift_y_tick_major))
                # ax.xaxis.set_major_formatter('{x:.0f}')
            if self.plot_lift_y_tick_minor:
                axs_cl.yaxis.set_minor_locator(MultipleLocator(self.plot_lift_y_tick_minor))
                axs_cl_global.yaxis.set_minor_locator(MultipleLocator(self.plot_lift_y_tick_minor))

            if mode == "local":
                return fig_cl, axs_cl
            elif mode == "global":
                return fig_cl_global, axs_cl_global

        # - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - #
        #                    L I L I E N T H A L                     #
        # - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - #
        def LLPlotLayout(mode):
            # Obj
            fig_ll, axs_ll = plt.subplots(1, 1, figsize=(3, 3))
            axs_ll.axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
            axs_ll.axvline(x=0, color='gray', linestyle='-', linewidth=0.5)
            axs_ll.set_xlabel(self.plot_drag_y_label)
            axs_ll.set_ylabel(self.plot_lift_y_label)
            # Global
            fig_ll_global, axs_ll_global = plt.subplots(1, 1, figsize=(3, 3))
            axs_ll_global.axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
            axs_ll_global.axvline(x=0, color='gray', linestyle='-', linewidth=0.5)
            axs_ll_global.set_xlabel(self.plot_drag_y_label)
            axs_ll_global.set_ylabel(self.plot_lift_y_label)
            if self.plot_grid:
                axs_ll.grid(True, which='major', axis='both', linewidth=0.1, color='grey')
                axs_ll_global.grid(True, which='major', axis='both', linewidth=0.1, color='grey')
            if self.plot_lift_xlims:
                axs_ll.set_xlim(self.plot_drag_ylims[0], self.plot_drag_ylims[1])
                axs_ll_global.set_xlim(self.plot_drag_ylims[0], self.plot_drag_ylims[1])
            if self.plot_lift_ylims:
                axs_ll.set_ylim(self.plot_lift_ylims[0], self.plot_lift_ylims[1])
                axs_ll_global.set_ylim(self.plot_lift_ylims[0], self.plot_lift_ylims[1])
            if self.plot_lift_x_tick_major:
                axs_ll.xaxis.set_major_locator(MultipleLocator(self.plot_drag_y_tick_major))
                axs_ll_global.xaxis.set_major_locator(MultipleLocator(self.plot_drag_y_tick_major))
                # ax.xaxis.set_major_formatter('{x:.0f}')
            if self.plot_lift_x_tick_minor:
                axs_ll.xaxis.set_minor_locator(MultipleLocator(self.plot_drag_y_tick_minor))
                axs_ll_global.xaxis.set_minor_locator(MultipleLocator(self.plot_drag_y_tick_minor))
            if self.plot_lift_y_tick_major:
                axs_ll.yaxis.set_major_locator(MultipleLocator(self.plot_lift_y_tick_major))
                axs_ll_global.yaxis.set_major_locator(MultipleLocator(self.plot_lift_y_tick_major))
                # ax.xaxis.set_major_formatter('{x:.0f}')
            if self.plot_lift_y_tick_minor:
                axs_ll.yaxis.set_minor_locator(MultipleLocator(self.plot_lift_y_tick_minor))
                axs_ll_global.yaxis.set_minor_locator(MultipleLocator(self.plot_lift_y_tick_minor))

            if mode == "local":
                return fig_ll, axs_ll
            elif mode == "global":
                return fig_ll_global, axs_ll_global

        # - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - #
        #          L I F T   -   S I D E - b y - S I D E              #
        # - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - #
        def sbsPlotLayout(mode):
            # Obj
            fig_sbs, axs_sbs = plt.subplots(1, 2, figsize=(7, 3))
            fig_sbs.subplots_adjust(wspace=0.05)
            # Left plot (lift)
            axs_sbs[0].axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
            axs_sbs[0].axvline(x=0, color='gray', linestyle='-', linewidth=0.5)
            axs_sbs[0].set_xlabel(self.plot_lift_x_label)
            axs_sbs[0].set_ylabel(self.plot_lift_y_label)
            # Right plot (Lilienthal)
            axs_sbs[1].axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
            axs_sbs[1].axvline(x=0, color='gray', linestyle='-', linewidth=0.5)
            axs_sbs[1].set_xlabel(self.plot_drag_y_label)
            axs_sbs[1].yaxis.set_ticklabels([])

            # Global
            # Left plot (lift)
            fig_sbs_global, axs_sbs_global = plt.subplots(1, 2, figsize=(7, 3), sharey=True)
            fig_sbs_global.subplots_adjust(wspace=0.05)
            axs_sbs_global[0].axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
            axs_sbs_global[0].axvline(x=0, color='gray', linestyle='-', linewidth=0.5)
            axs_sbs_global[0].set_xlabel(self.plot_lift_x_label)
            axs_sbs_global[0].set_ylabel(self.plot_lift_y_label)
            # Lilienthal
            axs_sbs_global[1].axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
            axs_sbs_global[1].axvline(x=0, color='gray', linestyle='-', linewidth=0.5)
            axs_sbs_global[1].set_xlabel(self.plot_drag_y_label)
            axs_sbs_global[1].yaxis.set_ticklabels([])
            if self.plot_grid:
                axs_sbs[0].grid(True, which='major', axis='both', linewidth=0.1, color='grey')
                axs_sbs[1].grid(True, which='major', axis='both', linewidth=0.1, color='grey')
                axs_sbs_global[0].grid(True, which='major', axis='both', linewidth=0.1, color='grey')
                axs_sbs_global[1].grid(True, which='major', axis='both', linewidth=0.1, color='grey')
            if self.plot_lift_xlims:
                axs_sbs[0].set_xlim(self.plot_lift_xlims[0], self.plot_lift_xlims[1])
                axs_sbs[1].set_xlim(self.plot_drag_ylims[0], self.plot_drag_ylims[1])
                axs_sbs_global[0].set_xlim(self.plot_lift_xlims[0], self.plot_lift_xlims[1])
                axs_sbs_global[1].set_xlim(self.plot_drag_ylims[0], self.plot_drag_ylims[1])
            if self.plot_lift_ylims:
                axs_sbs[0].set_ylim(self.plot_lift_ylims[0], self.plot_lift_ylims[1])
                axs_sbs_global[0].set_ylim(self.plot_lift_ylims[0], self.plot_lift_ylims[1])
            if self.plot_lift_x_tick_major:
                axs_sbs[0].xaxis.set_major_locator(MultipleLocator(self.plot_lift_x_tick_major))
                axs_sbs[1].xaxis.set_major_locator(MultipleLocator(self.plot_drag_y_tick_major))
                axs_sbs_global[0].xaxis.set_major_locator(MultipleLocator(self.plot_lift_x_tick_major))
                axs_sbs_global[1].xaxis.set_major_locator(MultipleLocator(self.plot_drag_y_tick_major))
                # ax.xaxis.set_major_formatter('{x:.0f}')
            if self.plot_lift_x_tick_minor:
                axs_sbs[0].xaxis.set_minor_locator(MultipleLocator(self.plot_lift_x_tick_minor))
                axs_sbs[1].xaxis.set_minor_locator(MultipleLocator(self.plot_drag_y_tick_minor))
                axs_sbs_global[0].xaxis.set_minor_locator(MultipleLocator(self.plot_lift_x_tick_minor))
                axs_sbs_global[1].xaxis.set_minor_locator(MultipleLocator(self.plot_drag_y_tick_minor))
            if self.plot_lift_y_tick_major:
                axs_sbs[0].yaxis.set_major_locator(MultipleLocator(self.plot_lift_y_tick_major))
                axs_sbs[1].yaxis.set_major_locator(MultipleLocator(self.plot_drag_y_tick_major))
                axs_sbs_global[0].yaxis.set_major_locator(MultipleLocator(self.plot_lift_y_tick_major))
                axs_sbs_global[1].yaxis.set_major_locator(MultipleLocator(self.plot_drag_y_tick_major))
            if self.plot_lift_y_tick_minor:
                axs_sbs[0].yaxis.set_minor_locator(MultipleLocator(self.plot_lift_y_tick_minor))
                axs_sbs[1].yaxis.set_minor_locator(MultipleLocator(self.plot_drag_y_tick_minor))
                axs_sbs_global[0].yaxis.set_minor_locator(MultipleLocator(self.plot_lift_y_tick_minor))
                axs_sbs_global[1].yaxis.set_minor_locator(MultipleLocator(self.plot_drag_y_tick_minor))

            if mode == "local":
                return fig_sbs, axs_sbs
            elif mode == "global":
                return fig_sbs_global, axs_sbs_global

        # Initialize global figures
        fig_cd_global, axs_cd_global = liftPlotLayout("global")
        fig_cl_global, axs_cl_global = liftPlotLayout("global")
        fig_ll_global, axs_ll_global = LLPlotLayout("global")
        fig_sbs_global, axs_sbs_global = sbsPlotLayout("global")

        # Initialize DataFrame for 3D plotting
        df = pd.read_csv(f"{self.post_folder}{obj}/{obj}.csv", sep=",", index_col="alpha")
        df_cl = df.drop(columns=["cd", "cl", "drag_force", "lift_force"])

        # 4. Loop over all Objects and generate drag, lift, LL and side-by-side plot
        counter = 0
        for obj in self.plotObjList:
            tmp_df = pd.read_csv(f"{self.post_folder}{obj}/{obj}.csv", sep=",", index_col="alpha")
            df_cl[counter] = tmp_df["cl"]
            print(f"object: {obj}", tmp_df)
            counter += 1
            try:
                ref_data = self.setup["Objects"][obj]["ref-data"]
            except:
                ref_data = False
            else:
                if ref_data:
                    df_ref = pd.read_csv(ref_data, sep=",")
            try:
                legend_name = self.setup["Objects"][obj]["plot"]["legend-name"]
            except:
                legend_name = obj
            try:
                line_width = self.setup["Objects"][obj]["plot"]["line-width"]
            except:
                line_width = None

            objResFolder = f"{self.post_folder}{obj}"
            if not objResFolder[-1] == "/":
                objResFolder += "/"
            if not os.path.exists(objResFolder):
                console.print(f"[bold red]Error:[not bold bright_white] Result folder for {obj} not found...")
                sys.exit()

            # Read in CSV
            df = pd.read_csv(f"{objResFolder}{obj}.csv", sep=",")

            # Initiliaze Obj figures
            fig_cd, axs_cd = liftPlotLayout("local")
            fig_cl, axs_cl = liftPlotLayout("local")
            fig_ll, axs_ll = LLPlotLayout("local")
            fig_sbs, axs_sbs = sbsPlotLayout("local")
            # Local plot
            if self.legend:
                axs_cd.plot(df['alpha'], df['cd'], linewidth=line_width, label=f"{legend_name}")
                axs_cl.plot(df['alpha'], df['cl'], linewidth=line_width, label=f"{legend_name}")
                axs_ll.plot(df['cd'], df['cl'], linewidth=line_width, label=f"{legend_name}")
                axs_sbs[0].plot(df['alpha'], df['cl'], linewidth=line_width)
                axs_sbs[1].plot(df['cd'], df['cl'], linewidth=line_width, label=f"{legend_name}")
            else:
                axs_cd.plot(df['alpha'], df['cd'], linewidth=line_width)
                axs_cl.plot(df['alpha'], df['cl'], linewidth=line_width)
                axs_ll.plot(df['cd'], df['cl'], linewidth=line_width)
                axs_sbs[0].plot(df['alpha'], df['cl'], linewidth=line_width)
                axs_sbs[1].plot(df['cd'], df['cl'], linewidth=line_width)
            # Global plot
            axs_cd_global.plot(df['alpha'], df['cd'], linewidth=line_width, label=f"{legend_name}")
            axs_cl_global.plot(df['alpha'], df['cl'], linewidth=line_width, label=f"{legend_name}")
            axs_ll_global.plot(df['cd'], df['cl'], linewidth=line_width, label=f"{legend_name}")
            axs_sbs_global[0].plot(df['alpha'], df['cl'], linewidth=line_width)
            axs_sbs_global[1].plot(df['cd'], df['cl'], linewidth=line_width, label=f"{legend_name}")
            if ref_data:
                if self.legend:
                    axs_cd.scatter(df_ref['alpha'], df_ref['cd'], color='red', marker='+', label=f"{legend_name} Ref.")
                    axs_cl.scatter(df_ref['alpha'], df_ref['cl'], color='red', marker='+', label=f"{legend_name} Ref.")
                    axs_ll.scatter(df_ref['cd'], df_ref['cl'], color='red', marker='+', label=f"{legend_name} Ref.")
                    axs_sbs[0].scatter(df_ref['alpha'], df_ref['cl'], label=f"{obj} ref")
                    axs_sbs[1].scatter(df_ref['cd'], df_ref['cl'], label=f"{obj} ref")
                else:
                    axs_cd.scatter(df_ref['alpha'], df_ref['cd'], color='red', marker='+')
                    axs_cl.scatter(df_ref['alpha'], df_ref['cl'], color='red', marker='+')
                    axs_ll.scatter(df_ref['cd'], df_ref['cl'], color='red', marker='+')
                    axs_sbs[0].scatter(df_ref['alpha'], df_ref['cl'])
                    axs_sbs[1].scatter(df_ref['cd'], df_ref['cl'])
                # Global plot with legend
                axs_cd_global.scatter(df_ref['alpha'], df_ref['cd'], marker='+', label=f"{legend_name} Ref.")
                axs_cl_global.scatter(df_ref['alpha'], df_ref['cl'], marker='+', label=f"{legend_name} Ref.")
                axs_ll_global.scatter(df_ref['cd'], df_ref['cl'], marker='+', label=f"{legend_name} Ref.")
                axs_sbs_global[0].scatter(df_ref['alpha'], df_ref['cl'], marker='+', label=f"{legend_name} Ref.")
                axs_sbs_global[1].scatter(df_ref['cd'], df_ref['cl'], marker='+')

            if self.legend:
                axs_cd.legend()
                axs_cl.legend()
                axs_ll.legend()
                axs_sbs.legend(bbox_to_anchor=(-1, -0.3), loc="lower left", ncol=4, prop={'size': 6})  # bbox_transform=fig.transFigure
            fig_cd.savefig(f"{self.post_folder}{obj}/{obj}_cd.pdf")
            fig_cl.savefig(f"{self.post_folder}{obj}/{obj}_cl.pdf")
            fig_ll.savefig(f"{self.post_folder}{obj}/{obj}_LL.pdf")
            fig_sbs.savefig(f"{self.post_folder}{obj}/{obj}_SBS.pdf")
            # Clear figures for next Obj
            plt.close(fig_cd)
            plt.close(fig_cl)
            plt.close(fig_ll)
            plt.close(fig_sbs)

        # Enable legend for the global plots
        # Drag
        cd_leg = axs_cd_global.legend(ncol=2, bbox_to_anchor=(1.05, 1.00), prop={'size': 6}, frameon=True, fancybox=False, shadow=True) # loc="lower left", 
        cd_leg.get_frame().set_edgecolor('black')
        cd_leg.get_frame().set_linewidth(0.2)
        # Lift
        cl_leg = axs_cl_global.legend(ncol=2, bbox_to_anchor=(1.05, 1.00), prop={'size': 6}, frameon=True, fancybox=False, shadow=True) # loc="lower left", 
        cl_leg.get_frame().set_edgecolor('black')
        cl_leg.get_frame().set_linewidth(0.2)
        # Lilienthal
        ll_leg = axs_ll_global.legend(ncol=2, bbox_to_anchor=(1.05, 1.00), prop={'size': 6}, frameon=True, fancybox=False, shadow=True) # loc="lower left", 
        ll_leg.get_frame().set_edgecolor('black')
        ll_leg.get_frame().set_linewidth(0.2)
        # Side-by-Side
        sbs_leg = axs_sbs_global[1].legend(ncol=8, bbox_to_anchor=(0.95, 1.2), prop={'size': 6}, frameon=True, fancybox=False, shadow=True) # loc="lower left", 
        sbs_leg.get_frame().set_edgecolor('black')
        sbs_leg.get_frame().set_linewidth(0.2)

        # Save global plots
        fig_cd_global.savefig(f"{self.post_folder}cd.pdf")
        fig_cl_global.savefig(f"{self.post_folder}cl.pdf")
        fig_ll_global.savefig(f"{self.post_folder}LL.pdf")
        fig_sbs_global.savefig(f"{self.post_folder}SBS.pdf")

        # Export dataframe
        df_cl.to_csv(f"{self.post_folder}3D_cl.csv", index=True, sep=",")