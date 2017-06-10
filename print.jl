# Initial set up print out to terminal window, qg1l params
# Bryan Kaiser
# 6/2/17

# v.6

# ==========================================================================
# simulation description (printout of important variables)

fileIO = open("./simulation_parameters.txt","w")

Ui = taux/(rho*H*beta*Lx) # m/s, interior velocity from wind-stress - beta term balance
nu = 1.5e-6 # m^2/s, kinematic viscosity of water at 5 degrees C
Roi = Ui/(f*Lx); # interior Rossby number (hypothesis)
dI = sqrt(Ui/beta); # m, inertial thickness
dI0 = dI/Lx
dM = (kappa0/beta)^(1.0/3.0) # m, horizontal friction thickness
dM0 = dM/Lx
dS = r/(2.0*beta); # m
dS0 = dS/Lx;
dBL = max(dI,dM,dS); # km, the predicted boundary layer thickness for the flow
dBL0 = dBL/Lx; # dimensionless, predicted boundary layer thickness for the flow
UBL = Ui*(Lx-dBL)/dBL; # m/s, predicted boundary layer velocity
dx0 = dx/Lx; # dimensionless grid spacing
dy0 = dy/Ly; # dimensionless grid spacing
CoX = UBL*(dt/dx); # zonal Courant number estimate
CoY = UBL*(dt/dy); # meridonal Courant number estimate
DX = kappa0*dt/(dx^2.0); # zonal diffusion stability constant
DY = kappa0*dt/(dy^2.0); # meridonal diffusion stability constant

Re_T = ((dI/dM)^3.0)*sqrt(beta/Ui)*Lx; # the qg Peclet number = UL/kappa
dIdM3 = (dI/dM)^3.0; # Fox-Kemper (2005) definition: "boundary layer Reynolds number"
Re = Ui*Lx/nu; # true Reynolds number assiming wind / beta balance
ReBL = UBL*dBL/nu; # true BL Reynolds number 
dxdt = dx/(dt_days*3600.0*24.0);
dydt = dy/(dt_days*3600.0*24.0);
dEb = r*H/f; # m, benthic BL thickness

write(fileIO,"\n********** dimensional parameters **********\n 
$(round(nm*dt_days,3)) days = time at which averaging begins
$(round(ns*dt_days,3)) days = time at which statistics are collected
dt = $dt_days days, time step
kappa = $kappa0 m^2/s, eddy vorticity diffusivity (interior)
tau = $taux Pa, zonal wind stress magnitude
U = $Ui m/s, interior velocity
dI = $(dI/1000.0) km, inertial boundary layer thickness
dM = $(dM/1000.0) km, horizontal friction boundary layer thickness
dS = $(dS/1000.0) km, bottom friction boundary layer thickness
dBL = $(dBL/1000.0) km, apparent boundary layer thickness
dx = $(dx/1000.0) km, zonal grid spacing
dy = $(dy/1000.0) km, meridonal grid spacing
Lx = $(Lx/1000.0) km, zonal grid spacing
Ly = $(Ly/1000.0) km, meridonal grid spacing
H = $(H/1000.0) km, layer depth
dE = $dEb m, benthic Ekman boundary layer thickness
UB = $UBL m/s, boundary layer velocity
u,grid = $dxdt m/s, fastest zonal resolved velocity
v,grid = $dydt m/s, fastest meridonal resolved velocity
\n********** non-dimensional parameters **********\n
Ro = $(round(Roi,7)) interior Rossby number
Re = $Re true interior Reynolds number (Ui*L/nu)
Re_T = $Re_T the qgpv Peclet number, Ui*Lx/kappa=(dI/dM)^3*sqrt(beta/Ui)*Lx
(dI/dM)^3 = $dIdM3 Fox-Kemper (2005) boundary layer Reynolds number definition
dI/Lx = $dI0 dimensionless inertial thickness
dM/Lx = $dM0 dimensionless horizontal friction boundary layer thickness
dS/Lx = $dS0 dimensionless bottom friction boundary layer thickness
dBL/Lx = $dBL0 dimensionless apparent boundary layer thickness
dx/Lx = $dx0 dimensionless zonal grid spacing
dy/Ly = $dy0 dimensionless meridonal grid spacing
Co,x = $CoX estimated zonal Courant number (CFL condition: must be <1)
Co,y = $CoY estimated meridonal Courant number (CFL condition: must be <1)
D,x = $DX estimated zonal diffusion stability constant (must be <0.5)
D,y = $DY estimated meridonal diffusion stability constant (must be <0.5)")

close(fileIO)




