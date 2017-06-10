# barotropic_quasigeostrophic_flow_solver
This collection of Julia scripts solves for nonlinear barotropic quasi-geostrophic ocean circulation in a basin. 
There are options to solve for a single gyre, symmetric double gyres, or asymmetric double gyres. The solver is 
pseudo-spectral (advection is de-aliased by the 2/3 rule), includes explicit 4th-order Runge-Kutta time integration, 
and the boundary conditions at the basin walls are free slip (no stress in tangential velocity, impermiable).
