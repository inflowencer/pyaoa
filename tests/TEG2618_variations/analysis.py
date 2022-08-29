# Helper to import non packaged Python source file
import sys
from pathlib import Path
home = str(Path.home())
sys.path.insert(1, f"{home}/airfoil_lab/project_db/projects/software/iikj/src")

import aoa

simulation = aoa.Analysis("setup.yaml")

# aoaAnalyze.multiAirfoilAnalysis(airfoils, numberOfIterations=iter, generateImages=False)
# aoaAnalyze.multiAirfoilData(airfoils)

# Plot
# aoaAnalyze.multiAirfoilPlot(airfoils)

