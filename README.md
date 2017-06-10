# barotropic_quasigeostrophic_flow_solver
This collection of Julia scripts solves for nonlinear barotropic quasi-geostrophic ocean circulation in a basin. 
There are options to solve for a single gyre, symmetric double gyres, or asymmetric double gyres. The solver is 
pseudo-spectral (advection is de-aliased by the 2/3 rule), includes explicit 4th-order Runge-Kutta time integration, and 
the boundary conditions at the basin walls are set by the spectral bases (discrete sine transforms a la FFTW) are free 
slip (no stress in tangential velocity, impermiable).

To run the code, make sure you have the directories /psi and /q if you have plot_flag = 1 (output plots) and the directory 
/output if you have write_flag = 1 (output .h5 files). The execute command is "$ julia qg1l_basin.jl" and the code can be 
used as a black box by merely editing the physical and grid parameters in params.jl
