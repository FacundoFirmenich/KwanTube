import qutip
from qutip.solver.heom import HEOMSolver, DrudeLorentzPadeBath
import numpy as np

N = 2
H = qutip.qeye(N)
Q = qutip.sigmaz()
bath = DrudeLorentzPadeBath(Q, lam=0.1, gamma=0.1, T=300, Nk=1)
solver = HEOMSolver(H, [bath], max_depth=1)

print("--- HEOMSolver attributes ---")
for attr in sorted(dir(solver)):
    if not attr.startswith("_"):
        print(attr)

rho0 = qutip.basis(N, 0) * qutip.basis(N, 0).dag()
try:
    state = solver.create_state(rho0)
    print("\nSUCCESS: solver.create_state exists.")
except AttributeError:
    print("\nFAILED: solver.create_state does not exist.")
