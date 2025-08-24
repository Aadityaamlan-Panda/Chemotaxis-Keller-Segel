# 🌊 Chemotaxis–Keller–Segel MATLAB Simulation Suite

A MATLAB toolkit for simulating bacterial chemotaxis and reaction–diffusion systems using the Keller–Segel framework and related PDE/ODE/BVP formulations. Explore 38+ scripts covering explicit/implicit finite differences, Method of Lines (MoL), similarity transforms, PDEPE-based solvers, ODE solvers (Euler, RK2, RK4), and nonlinear least-squares formulations.

- 🧫 1D chemotaxis in a pore with ammonia gradients  
- 🔗 Coupled bacteria–chemoattractant dynamics  
- 🧮 Boundary-value formulations for density, fuel, and sensing variables  
- 🧪 Grid independence studies  
- ✅ Analytical checks and verification utilities

***

## ✨ Features

- Chemotaxis PDEs with drift from concentration gradients  
- Multiple discretization schemes: central, upwind-like, higher-order  
- Time integration via explicit schemes and MoL with ODE integrators  
- Eigenfunction-based closed-form building blocks  
- BVP solutions using bvp4c and finite-difference linear systems  
- Rich visualization: space–time surfaces and profile plots

***

## 📦 Repository Contents

Each item below corresponds to a section in the combined file (look for headers like “File: … end …”):

- PDEPE-based chemotaxis and diffusion  
  - Bacteria_ammonia_pore_diffusion.m  
  - Chemotaxis_Ecoli.m  
  - Coupled_chemotaxis.m  
  - fuel.m + fuel_CD.m, fuel_CD_ic.m, fuel_CD_bc.m

- Explicit/implicit finite-difference solvers  
  - explicit_FD_resurrected.m  
  - explicit_FD_revised_model.m  
  - implicit_FD.m  
  - implicit_FD_revised_model.m  
  - implicit_FD_similarity_transform.m  
  - modified_FD_Bacteria_ammonia_pore_diffusion.m  
  - FD_Bacteria_ammonia_pore_diffusion.m

- Method of Lines and ODE-based updates  
  - bact_method_of_liiines.m  
  - grid_variate_MoL.m  
  - grid_variate.m  
  - bact_odefunc.m

- Grid-independence and verification  
  - grid_dep_bacteris_conc.m  
  - grid_dep_check_bact_method_of_lines.m  
  - c_values.m  
  - myfunc.m (diffusion residual)

- Analytical/alternative formulations  
  - similarity_transform_Bacteria_ammonia_pore_diffusion.m  
  - odefcn_BACD.m (similarity ODE system)  
  - bvp_fuel.m (BVP via bvp4c)  
  - matrix_fuel.m (FD linear system)  
  - fuel_euler.m, fuel_RK2.m, fuel_RK4.m (time integrators)  
  - root_BACD.m, root3_BACD.m, root_sim_BACD.m (nonlinear residuals)  
  - ode_fuel.m (symbolic setup)  
  - test.m (quick analytic comparison)

A full “Table of Contents” is embedded at the top of the combined file.

***

## 🚀 Getting Started

### Requirements
- MATLAB R2018a or newer (recommended)  
- Optimization Toolbox (for lsqnonlin/fsolve used by implicit/residual solvers)  
- Symbolic Math Toolbox (optional, for ode_fuel.m)  
- No external dependencies beyond standard toolboxes

### Quick Start

Option A — Use the combined file as an index:
- Open combined_matlab_code.m (or the generated combined file).
- Copy any section (between “File:” and “End of …”) into its own .m file using the same name.
- Ensure helper functions are in the same folder or added to the MATLAB path.

Option B — Run individual scripts:
- Put all .m files in one folder.
- In MATLAB: set that folder as the working directory.
- Run primary scripts:
  - Bacteria_ammonia_pore_diffusion.m
  - Chemotaxis_Ecoli.m
  - explicit_FD_revised_model.m
  - implicit_FD_revised_model.m
  - bact_method_of_liiines.m
  - fuel.m

Most scripts auto-generate figures (surf, plot).

***

## 🧭 Typical Workflows

- PDEPE chemotaxis in a pore  
  - Run: Bacteria_ammonia_pore_diffusion.m  
  - Needs: coupledpde_BACD.m, coupledic_BACD.m, coupledbc_BACD.m

- Explicit FD with gradient-induced velocity  
  - Run: explicit_FD_revised_model.m  
  - Tune: d_ratio, v (modes), N/M, points-of-interest

- Method of Lines (chemotactic sensitivity law)  
  - Run: bact_method_of_liiines.m or grid_variate_MoL.m  
  - Uses: bact_odefunc.m internally (semi-discrete advection–diffusion)

- Grid independence  
  - Time refinement, space coarsening: grid_dep_bacteris_conc.m  
  - Space refinement for MoL: grid_dep_check_bact_method_of_lines.m

- Boundary-value “fuel” system  
  - PDEPE version: fuel.m (+ fuel_CD.m, fuel_CD_ic.m, fuel_CD_bc.m)  
  - Alternatives: bvp_fuel.m (bvp4c), matrix_fuel.m (FD linear system)  
  - Time integrators: fuel_euler.m, fuel_RK2.m, fuel_RK4.m

- Similarity transform check  
  - Run: similarity_transform_Bacteria_ammonia_pore_diffusion.m  
  - Uses: odefcn_BACD.m; compare qualitative profiles

***

## ⚙️ Parameters and Tuning

- Domain & grids: L, T, N (space), M (time)  
- Transport ratios: Db, Dc, d_ratio = Db/Dc  
- Chemotactic constants: X0, k1, k2, Ds, Cch  
- Series accuracy: v (eigenmodes)  
- Stability: for explicit schemes, ensure Db*dt/dx^2 <= 0.5

Tip: Most scripts mark “DEFINE CONSTANTS” / “ENTER …” blocks for quick edits.

***

## 📝 Notes and Known Issues

- Keep helper functions in the same folder or addpath the directory.  
- test.m typo: change plot(z,f _anly) → plot(z,f_anly).  
- edited_Bacteria_ammonia_pore_diffusion.m: replace xdom/tdom with x/t when plotting.  
- Optimization-based solvers (implicit_FD, implicit_FD_revised_model) require Optimization Toolbox.  
- ode_fuel.m (symbolic) is optional and may be slow without assumptions.

***

## 🤝 Contributing

- Open issues for bugs, numerical instability, or enhancements  
- PRs welcome with:
  - Clear description of changes  
  - Comments for scheme choices and assumptions  
  - Before/after plots where helpful

***

## 📄 License

MIT License.
***

## 🏷️ Credits

- Author = Aaditya Amlan Panda,
- Mentor  = Prof. Akash Choudhary, Department of Chemical Engineering, IIT Kanpur
- Peer Support  = Punam Singh, Kushagra Tiwari
- Year   = 2024


***

## 🙏 Acknowledgements

Developed in the context of project on chemotaxis, Keller–Segel modeling, and numerical PDEs/ODEs. Thanks to mentors and peers for insights on stability, discretization, and verification.

***


