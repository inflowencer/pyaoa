# Helper to import non packaged Python source file
import sys
sys.path.insert(1, '/home/joe/airfoil_lab/project_db/projects/software/iikj/src')

import aoa

ambientConditions = {
    'Re' : 8e5,
    'p' : 9.8e4,
    'T' : 300.
}

aoaMin = -7.0
aoaMax = 17.0

# basis = aoaAnalyze.aoaAnalysis('TEG2618_b', ambientConditions, 0.35, aoaMin=aoaMin, aoaMax=aoaMax, numberOfIncrements=30, expData=True)
var2 = aoa.aoaAnalysis('TEG2618_var2', ambientConditions, 0.355, aoaMin=aoaMin, aoaMax=aoaMax, numberOfIncrements=30)
# var3 = aoaAnalyze.aoaAnalysis('TEG2618_var3', ambientConditions, 0.355, aoaMin=aoaMin, aoaMax=aoaMax, numberOfIncrements=30)

airfoils = [var2]
iter = 1000

# aoaAnalyze.multiAirfoilAnalysis(airfoils, numberOfIterations=iter, generateImages=False)
# aoaAnalyze.multiAirfoilData(airfoils)

# Plot
# aoaAnalyze.multiAirfoilPlot(airfoils)

