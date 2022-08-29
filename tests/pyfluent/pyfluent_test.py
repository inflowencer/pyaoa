import ansys.fluent.core as pyfluent
session = pyfluent.launch_fluent(mode="solver")
session.check_health()

session.solver.tui.file.read_case()