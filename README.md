# ci-GSPT-reductions-of-chemical-reaction-networks
Code for the numerical simulations in the paper titled 'Coordinate-independent model reductions of chemical reaction networks based on geometric singular perturbation theory' by T.E.F. Lapuz and M. Wechselberger (2026)

In particular,
  - MM_solutions.m reproduces Figure 2b in the paper. It accesses MM_ODE.m and MM_ODE_GSPT.m in order to numerically solve the relevant ODEs.
- KF_solutions.m reproduces Figure 3b in the paper. It accesses KF_ODE.m and KF_ODE_GSPT.m in order to numerically solve the relevant ODEs. There is an instruction in the code to reproduce Figure 3a.
