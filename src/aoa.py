#!/usr/bin/env python3

import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt
import os

plt.style.use('science')


class aoaAnalysis:
    '''Class for generating airfoil objects to carry out angle-of-attack analysis

    Parameters
    ----------
    airfoilName : string
    aoaMin : float
    aoaMax : float
    numberOfIncrements : int
    '''

    def __init__(
        self, airfoilName, ambientConditions, L, yplus=0.2, aoaMin=-5.0, aoaMax=10.0, numberOfIncrements=10, inputFilename='',
        transcriptFilename='', resultFilename='', expData=False, plotFilename='', meshFilename='', experimentalDataFilename=''
        ) -> None:
        

        self.airfoilName = str(airfoilName)
        self.aoaMin = float(aoaMin)
        self.aoaMax = float(aoaMax)
        self.numberOfIncrements = int(numberOfIncrements)
        self.expData = expData
        self.Re = ambientConditions['Re']
        self.p = ambientConditions['p']
        self.T = ambientConditions['T']
        self.L = L
        self.C_f = 0.0576 * self.Re**(-0.2)  # skin friction coefficient
        self.rho = 0.029 / 8.314 * self.p / self.T                       # density according to material data
        self.mu = 1.9e-5
        self.nu = self.mu / self.rho
        self.u_b = self.Re * self.nu / self.L                 # free stream velocity acc. to Re
        self.tau_w = 0.5 * self.C_f * self.rho * self.u_b**2  # wall shear stress
        self.u_tau = np.sqrt(self.tau_w / self.rho)           # friction velocity
        self.yplus = yplus
        self.y_1 = self.yplus * self.nu / self.u_tau          # first layer height
        self.a = np.sqrt(1.4 * 8.314 / 0.029 * self.T)
        self.Ma = self.u_b / self.a
        # # print(self.nu)
        # print(self.Ma)
        # print(self.u_b)

        # Inputfile for solver
        if not inputFilename:
            self.inputFilename = "input/" + self.airfoilName + ".jou"
        else:
            self.inputFilename = inputFilename

        # Transciptfile where console output should be written to
        if not transcriptFilename:
            self.transcriptFilename = "transcripts/" + self.airfoilName + ".out"
        else:
            self.transcriptFilename = transcriptFilename

        # Resultfile should be same transcriptfile
        if not resultFilename:
            self.resultFilename = "transcripts/" + self.airfoilName + ".out"
        else:
            self.resultFilename = resultFilename
        
        # Experimental data
        if not experimentalDataFilename:
            self.experimentalDataFilename = "airfoil-data/" + self.airfoilName + "_exp.csv"
        else:
            self.experimentalDataFilename = experimentalDataFilename
        
        # Plot file path
        if not plotFilename:
            self.plotFilename = "plots/" + self.airfoilName + ".pdf"
        else:
            self.plotFilename = plotFilename

        # Mesh file path
        if not meshFilename:
            # self.meshFilename = "mesh/" + self.airfoilName + ".msh"
            self.meshFilename = "mesh/" + self.airfoilName + ".cas"
        else:
            self.meshFilename = meshFilename

        self.df_sim = None
        self.df_exp = None
        self.fluentString = None


    def generateInputFile(self, numberOfIterations=500, generateImages=True):
        '''Function to export the input string for your CFD solver'''

        angleList = np.linspace(self.aoaMin, self.aoaMax, self.numberOfIncrements + 1, retstep=True)
        delta_alpha = angleList[1]

        # Start transcript in order for the console output to be written to a file
        fluentString = ';{}\n/file/start-transcript ../{} o\n\n'.format(self.airfoilName, self.transcriptFilename)
        fluentString += '/file/replace-mesh {} o\n'.format(self.meshFilename)
        fluentString += '/define/bound/set/pres/inlet_farfield () mach n {} ()\n'.format(self.Ma)

        # First solution at minimum angle:
        fluentString += '''/mesh/rotate {} 0 0 0 0 0 1\n'''.format(-1 * self.aoaMin)  # rotate to initial angle
        fluentString += '''/solve/init/hyb\n/solve/it {}\n'''.format(numberOfIterations)  # initial solution
        fluentString += '''(display "Angle is: {} deg")\n'''.format(self.aoaMin)
        fluentString += '''/report/forces/wall-forces yes 1 0 0 no\n''' # print the drag forces and coefficient
        fluentString += '''/report/forces/wall-forces yes 0 1 0 no\n\n''' # print the lift forces and coefficient
        # if generateImages:
        #     fluentString += '''/display/views/rest view-0\n/display/objects/display logk\n'''
        #     fluentString += '''/display/save-picture ../res/img/{}_alpha_{}_logk.png\n'''.format(self.airfoilName, round(self.aoaMin))
        #     fluentString += '''/display/views/rest view-0\n/display/objects/display u\n'''
        #     fluentString += '''/display/save-picture ../res/img/{}_alpha_{}_u.png\n'''.format(self.airfoilName, round(self.aoaMin))
        #     fluentString += '''/display/views/rest view-0\n/display/objects/display p\n'''
        #     fluentString += '''/display/save-picture ../res/img/{}_alpha_{}_p.png\n'''.format(self.airfoilName, round(self.aoaMin))

        for i in list(range(len(angleList[0]) - 1)):
            fluentString += '''/mesh/rotate {} 0 0 0 0 0 1\n'''.format(-1 * delta_alpha)  # rotate by steplength
            fluentString += '''/solve/init/hyb o\n/solve/it {}\n'''.format(numberOfIterations)
            fluentString += '''(display "Angle is: {} deg")\n'''.format(round(angleList[0][i + 1], 3))
            fluentString += '''/report/forces/wall-forces yes 1 0 0 no\n'''
            fluentString += '''/report/forces/wall-forces yes 0 1 0 no\n\n'''
            # if generateImages:
            #     fluentString += '''/display/views/rest view-0\n/display/objects/display logk\n'''
            #     fluentString += '''/display/save-picture ../res/img/{}_alpha_{}_logk.png\n'''.format(self.airfoilName, round(angleList[0][i + 1], 3))
            #     fluentString += '''/display/views/rest view-0\n/display/objects/display u\n'''
            #     fluentString += '''/display/save-picture ../res/img/{}_alpha_{}_u.png\n'''.format(self.airfoilName, round(angleList[0][i + 1], 3))
            #     fluentString += '''/display/views/rest view-0\n/display/objects/display p\n'''
            #     fluentString += '''/display/save-picture ../res/img/{}_alpha_{}_p.png\n'''.format(self.airfoilName, round(angleList[0][i + 1], 3))

        fluentString += '/file/stop-transcript\n\n'
        self.fluentString  = fluentString

        with open(self.inputFilename, "w") as f:
            f.write(fluentString)


    def generateAerodynamicData(self, aoa=True, cd=True, cl=True, exportData=True, exportFilename=''):
        '''Converts the information from the output file from your CFD solver into tabular data (aoa, cd and cl)'''

        # Create airfoil-data directory if it doesn't exist already
        try:
            os.makedirs('airfoil-data')
        except:
            pass

        # Create export file name if 
        if exportData:
            if not exportFilename:
                exportFileName = f'airfoil-data/{self.airfoilName}_sim.csv'

        outputFileFluent = self.resultFilename

        with open(outputFileFluent, "r") as f:
            lines = f.readlines()

        angleList = []
        cdList = []
        clList = []

        for i, line in enumerate(lines):
            if line.startswith("Angle is"):
                # angle = [int(d) for d in re.findall(r'-?\d+', line)][0]
                angle = [float(d) for d in re.findall(r"[-+]?(?:\d*\.\d+|\d+)", line)][0]
                angleList.append(angle)

            if line.startswith("Forces - Direction Vector"):
                # Find drag
                if line[-7] == '1':
                    coeffs = re.findall(r"[-+]?(?:\d*\.\d+|\d+)", lines[i+3])
                    cdList.append(float(coeffs[-1]))
                # Find lift
                elif line[-5] == '1':
                    coeffs = re.findall(r"[-+]?(?:\d*\.\d+|\d+)", lines[i+3])
                    clList.append(float(coeffs[-1]))

        df_sim = pd.DataFrame({'alpha': angleList, 'cd': cdList, 'cl': clList})

        # Export data to CSV if true
        if exportData:
            df_sim.to_csv(exportFileName, sep=',', index=False)

        self.df_sim = df_sim


    def generatePlot(self, plotType='polar-and-lilienthal', exportPlot=True, expDataAvailable=True, plotLimits={}):
        '''Generates and exports LaTeX figures from the airfoil data.

        Parameters
        ----------
        plotType : \'polar-and-lilienthal\', \'polar\', \'lilienthal\'
        exportPlot : True
        expDataAvailable : True
        plotLimits : dict(xmin_aoa=float, xmax_aoa=float, xmin_cd=float, xmax_cd=float, ymin_cl=float, ymax_cl=float)
        '''

        # If not data has been read in already
        if not self.df_sim:
            self.df_sim = pd.read_csv(f'airfoil-data/{self.airfoilName}_sim.csv', sep=',')

        if expDataAvailable:
            self.df_exp = pd.read_csv(self.experimentalDataFilename, sep=',')

        # Check if plotlimits are defined and if yes which ones
        if plotLimits:
            activeLimits = dict.fromkeys( ['xmin_aoa', 'xmax_aoa', 'xmin_cd', 'xmax_cd', 'ymin_cl', 'ymax_cl'] )
            for key, value in activeLimits.items():
                if key in plotLimits:
                    activeLimits[key] = plotLimits[key]


        if plotType == 'polar-and-lilienthal':

            # Initialize 
            fig, axs = plt.subplots(1, 2, figsize=(6, 3), sharey=True)

            # aoa-Lift polar
            axs[0].plot(self.df_sim['aoa'], self.df_sim['aoa'], color='black',)# label='Simulation')
            axs[0].scatter(self.df_exp['alpha'], self.df_exp['cl'], color='red', marker='+',)# label='Experiment')
            axs[0].axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
            axs[0].axvline(x=0, color='gray', linestyle='-', linewidth=0.5)
            axs[0].set_xlabel(r'Angle of attack $\alpha$')
            axs[0].grid(True, which='both', axis='both', linewidth=0.1, color='grey')
            # Set aoa xlimits if defined
            if activeLimits['xmin_aoa'] and activeLimits['xmax_aoa']:
                axs[0].set_xlim(activeLimits['xmin_aoa'], activeLimits['xmax_aoa'])
            # Set cd ylimits if defined
            if activeLimits['ymin_cl'] and activeLimits['ymax_cl']:
                axs[0].set_ylim(activeLimits['ymin_cl'], activeLimits['ymax_cl'])
            axs[0].set_ylabel(r'Lift coefficient $c_l$')

            # cd-cl (Lilienthal)
            axs[1].plot(self.df_sim['aoa'], self.df_sim['cl'], color='black', label='Simulation')
            axs[1].scatter(self.df_sim['cd'], self.df_sim['cl'], color='red', marker='+', label='Experiment')
            axs[1].axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
            axs[1].set_xlabel(r'Drag coefficient $c_d$')
            axs[1].grid(True, which='both', axis='both', linewidth=0.1, color='grey')
            # Set cd xlimits if defined
            if activeLimits['xmin_cd'] and activeLimits['xmax_cd']:
                axs[1].set_xlim(activeLimits['xmin_cd'], activeLimits['xmax_cd'])
            leg = axs[1].legend(facecolor='white', frameon=True, framealpha=1)
            leg.get_lines()[0].set_linewidth(0.2)

            fig.savefig(self.plotFilename, dpi=300)

            # Clear the plot
            # plt.clf()
            # plt.cla()
            # plt.close()

def multiAirfoilAnalysis(airfoils, baseCaseFilename='base', inputFilename='', numberOfIterations=500, generateImages=True):
        '''Generate an input file .jou and export it.

        Parameters
        ----------
        airfoils : list() of airfoil objects generated by the aoaAnalysis() class
        baseCaseFilename : str() path to the base case \'.cas.h5\'
        inputFilename : str() path to the file where the .jou input file shall be written to
        '''

        fluentString = f''';; /file/read-case {baseCaseFilename}\n'''

        for airfoil in airfoils:
            airfoil.generateInputFile(numberOfIterations, generateImages)
            if not airfoil.fluentString:
                with open("input/" + airfoil.airfoilName + ".jou", 'r') as f:
                    inputString = f.read()
                    fluentString += inputString
            else:
                fluentString += airfoil.fluentString

        if not inputFilename:
            with open("input/multiAirfoilAnalysis.jou", "w") as f:
                f.write(fluentString)
        else:
            with open(inputFilename, "w") as f:
                f.write(fluentString)


def multiAirfoilData(airfoils):
    for airfoil in airfoils:
        airfoil.generateAerodynamicData()



def multiAirfoilPlot(airfoils, plotType='polar-and-lilienthal', plotLimits={}, plotFilename='plots/multiAirfoilComparison.pdf'):
    '''Generate plots from airfoil data'''
    '''Generates and exports LaTeX figures from the airfoil data.

    Parameters
    ----------
    plotType : \'polar-and-lilienthal\', \'polar\', \'lilienthal\'
    exportPlot : True
    expDataAvailable : True
    plotLimits : dict(xmin_aoa=float, xmax_aoa=float, xmin_cd=float, xmax_cd=float, ymin_cl=float, ymax_cl=float)
    '''
    # Initialize plot
    fig, axs = plt.subplots(1, 2, figsize=(12, 3), sharey=True)

    # Check if plotlimits are defined and if yes which ones
    if plotLimits:
        activeLimits = dict.fromkeys( ['xmin_aoa', 'xmax_aoa', 'xmin_cd', 'xmax_cd', 'ymin_cl', 'ymax_cl'] )
        for key, value in activeLimits.items():
            if key in plotLimits:
                activeLimits[key] = plotLimits[key]
        # Set aoa xlimits if defined
        if activeLimits['xmin_aoa'] and activeLimits['xmax_aoa']:
            axs[0].set_xlim(activeLimits['xmin_aoa'], activeLimits['xmax_aoa'])
        # Set cd ylimits if defined
        if activeLimits['ymin_cl'] and activeLimits['ymax_cl']:
            axs[0].set_ylim(activeLimits['ymin_cl'], activeLimits['ymax_cl'])
        axs[0].set_ylabel(r'Lift coefficient $c_l$')
        # Set cd xlimits if defined
        if activeLimits['xmin_cd'] and activeLimits['xmax_cd']:
            axs[1].set_xlim(activeLimits['xmin_cd'], activeLimits['xmax_cd'])

    # aoa-Lift polar
    axs[0].axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
    axs[0].axvline(x=0, color='gray', linestyle='-', linewidth=0.5)
    axs[0].set_xlabel(r'Angle of attack $\alpha$')
    axs[0].set_ylabel(r'Lift coeffient $c_l$')
    axs[0].grid(True, which='both', axis='both', linewidth=0.1, color='grey')

    # cd-cl (Lilienthal)
    axs[1].axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
    axs[1].set_xlabel(r'Drag coefficient $c_d$')
    axs[1].grid(True, which='both', axis='both', linewidth=0.1, color='grey')

    # leg.get_lines()[0].set_linewidth(0.2)

    for airfoil in airfoils:
        # If not data has been read in already
        # if airfoil.df_sim.empty():
        #     airfoil.df_sim = pd.read_csv(f'airfoil-data/{airfoil.airfoilName}_sim.csv', sep=',')

        if airfoil.expData:
            airfoil.df_exp = pd.read_csv("airfoil-data/" + airfoil.airfoilName + "_exp.csv", sep=',')

        if plotType == 'polar-and-lilienthal':
            # Polar
            axs[0].plot(airfoil.df_sim['alpha'], airfoil.df_sim['cl'], label='{}'.format(airfoil.airfoilName)) # Simulation data
            if airfoil.expData:
                axs[0].scatter(airfoil.df_exp['alpha'], airfoil.df_exp['cl'], color='red', marker='+', label=f'{airfoil.airfoilName} Experiment')   # Experimental data

            # Lilienthal
            axs[1].plot(airfoil.df_sim['cd'], airfoil.df_sim['cl'], label='{}'.format(airfoil.airfoilName)) # Simulation data
            if airfoil.expData:
                axs[1].scatter(airfoil.df_exp['cd'], airfoil.df_exp['cl'], color='red', marker='+', label=f'{airfoil.airfoilName} Experiment')  # Experimental data


    # leg = axs[1].legend(facecolor='white', frameon=True, framealpha=1)
    leg = plt.legend(bbox_to_anchor=(-1, -0.3), loc="lower left", ncol=3, prop={'size': 8})  # bbox_transform=fig.transFigure
    fig.savefig(plotFilename, dpi=300)