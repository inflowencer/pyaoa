# Helper to import non packaged Python source file
import sys
from pathlib import Path
home = str(Path.home())
sys.path.insert(1, f"{home}/airfoil_lab/project_db/projects/software/iikj/src")

import pyaoa.aoa as aoa

simulation = aoa.Analysis("setup.yml")