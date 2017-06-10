# Initialization of fields and parameters for simulation
# Bryan Kaiser
# 6/2/17

# v.6

# =============================================================================
# computational grid

if gyre_config == 1 || gyre_config == 2 # single or double gyre

  Lxcenter = 0.0; Lycenter = 0.0; # the center of the grid
  dx = Lx/Float64(Nx); dy = Ly/Float64(Ny); # m, uniform grid spacing
  x = collect(0.5*dx:dx:dx*Nx)-(Lx/2.0-Lxcenter); # m
  y = collect(0.5*dy:dy:dy*Ny)-(Ly/2.0-Lycenter); # m
  const X,Y = meshgrid( x , y ) # m

  # sine transform wavenumbers for derivatives/inversion
  ks = collect(1:Nx).*(pi/Lx); ls = collect(1:Ny).*(pi/Ly) # rad/m, for DST-II
  const Ks,Ls = meshgrid(ks,ls) 
  const Ks_2 = Ks.^2.0 # rad^2/m^2
  const Ls_2 = Ls.^2.0 # rad^2/m^2
  Kmags = (Ks.^2.0 + Ls.^2.0).^(1.0/2.0); # gridded wavenumber magnitudes

  # cosine transform wavenumbers for derivatives/inversion
  kc = collect(0:Nx-1).*(pi/Lx) # rad/m, for DCT-II
  lc = collect(0:Ny-1).*(pi/Ly) # rad/m, for DCT-II
  const Kc,Lc = meshgrid( kc , lc )

else # gyre_config == 3, two-gyre

  A = ((11.0+sqrt(21.0))/10.0) # ~1.56/2.0 = Ly2/Ly
  Ly2 = A*Ly/2.0;
  Lxcenter = 0.0; Lycenter = 0.0; # the center of the grid
  dx = Lx/Float64(Nx); dy = Ly2/Float64(Ny); # m, uniform grid spacing
  x = collect(0.5*dx:dx:dx*Nx)-(Lx/2.0-Lxcenter); # m
  y = collect(0.5*dy:dy:dy*Ny)-(Ly/2.0-Lycenter); # m
  const X,Y = meshgrid( x , y ) # m

  # sine transform wavenumbers for derivatives/inversion
  ks = collect(1:Nx).*(pi/Lx); ls = collect(1:Ny).*(pi/Ly) # rad/m, for DST-II
  const Ks,Ls = meshgrid(ks,ls) 
  const Ks_2 = Ks.^2.0 # rad^2/m^2
  const Ls_2 = Ls.^2.0 # rad^2/m^2
  Kmags = (Ks.^2.0 + Ls.^2.0).^(1.0/2.0); # gridded wavenumber magnitudes

  # cosine transform wavenumbers for derivatives/inversion
  kc = collect(0:Nx-1).*(pi/Lx) # rad/m, for DCT-II
  lc = collect(0:Ny-1).*(pi/Ly) # rad/m, for DCT-II
  const Kc,Lc = meshgrid( kc , lc )

end

# mask for 2/3 rule padding for de-aliasing a quadratic signal via fft
const mask = makemask( Ks , Ls , Kmags )


# =============================================================================
# initialization of physical constants

# wind stress curl
rho = 1000.0; # kg/m^2
taux_rho = taux/rho; # m^2/s^2
curltau = taux/(rho*H)*(2.0*pi/Ly) # 1/s^2 
if gyre_config == 1; # single gyre
  wind_profile = sin((Y-0.5*Ly).*(pi/Ly)) 
  const plotsize = [9.0,7.0]; # size of output plots 
elseif gyre_config == 2; # double gyre
  wind_profile = sin(Y.*2.0*pi/Ly) # 0<=Y<=Ly 
  const plotsize = [7.0,10.0]; # size of output plots
elseif gyre_config == 3; # two-gyre
  Yp = (Y+Ly/2.0).*(2.0/Ly) 
  wind_profile = -(1.0+Yp.*(4.0/A)).*(Yp.^4.0*(5.0/A^4.0)-Yp.^3.0*(16.0/A^3.0)+Yp.^2.0*(16.0/A^2.0)-Yp.*(5.0/A))./(10.0/(3.0*A^5.0)-59.0/(5.0*A^4.0)+12.0/A^3.0-4.0/(3.0*A^2.0)-5.0/(2.0*A)).*2.0/pi;
  const plotsize = [7.0,9.0]; # size of output plots
end
const wind = wind_profile.*curltau # 1/s^2

# deformation radius
if LR == 0.0; 
  const LR_m2 = 0.0; #  (a completely barotropic gyre).
else 
  const LR_m2 = LR^(-2.0); # m^-2 
end

const scalars = [ LR_m2 , beta , r ]; 

if kappa_profile == 1 # constant kappa
  const kappa = kappa0
elseif kappa_profile == 2 # boundary enhanced
  Ui = taux/(rho*H*beta*Lx) # m/s, interior velocity from wind-stress - beta term balance
  dI = sqrt(Ui/beta); # m, inertial thickness
  dI0 = dI/Lx
  dI03 = (dI/Lx)^3.0
  dM = (kappa0/beta)^(1.0/3.0) # m, horizontal friction thickness
  dM0 = dM/Lx
  dIdM3 = (dI/dM)^3.0; # Fox-Kemper (2005) definition: "boundary layer Reynolds number"
  d_k = dI0/sqrt(dIdM3) # d_kappa = dI/Lx/sqrt((dI/dM)^3)
  Xd = (X+Lx/2.0)./Lx;
  #print("$(minimum(Xd)) $(maximum(Xd))\n")
  #const kappa = ( dI03/dIdM3+( exp(-Xd./d_d)+exp(-(1.0-Xd)./d_d) ).*(dI03/Re_b-dI03/dIdM3) ).*(beta*Lx^3.0)
  const kappa = ( kappa0+( exp(-Xd./d_k)+exp(-(1.0-Xd)./d_k) ).*(kappa0*10.0-kappa0) )
end

#=
fig1 = figure(figsize=(plotsize)) # works!
CP1 = contourf(X./1000.0,Y./1000.0,kappa,Nres,cmap="gray") 
xlabel("x (km)"); ylabel("y (km)"); title("kappa")
colorbar(CP1); savefig("./kappa.png",format="png"); 
=#	


# =============================================================================
# initialization of qgpv, mean qgpv, and statistics fields

if n0 == 0 # quiescent initial fields
  qn = zeros( Ny , Nx ); # 1/s, qgpv initial field
  time_elapsed_start = 0.0;
  mean_convergence = 0.0;
  if nm < n0 # initialization for output statistics 
    qm = zeros( Ny , Nx ); # 1/s, mean qgpv initial field
    #q0 = zeros( Ny , Nx );
    um = zeros( Ny , Nx ); 
    vm = zeros( Ny , Nx ); 
    psim = zeros( Ny , Nx ); 
    zetam = zeros( Ny , Nx ); 
    uzm = zeros( Ny , Nx ); 
    vzm = zeros( Ny , Nx ); 
    #Rqq = zeros( Ny , Nx ); 
  end
else # input from .h5 file	
  if n0 <10 # filename
    readname = "$(output_path)/output/qg1l_00000$(convert(Int64,n0)).h5"; elseif n0 <100 
    readname = "$(output_path)/output/qg1l_0000$(convert(Int64,n0)).h5"; elseif n0 <1000
    readname = "$(output_path)/output/qg1l_000$(convert(Int64,n0)).h5"; elseif n0 <10000
    readname = "$(output_path)/output/qg1l_00$(convert(Int64,n0)).h5"; elseif n0 <100000
    readname = "$(output_path)/output/qg1l_0$(convert(Int64,n0)).h5"; elseif n0 <1000000
    readname = "$(output_path)/output/qg1l_$(convert(Int64,n0)).h5"
  end
  qn = h5read( readname , "q" ); # 1/s, imported qgpv initial field
  time_elapsed_start = h5read( readname , "time_elapsed" ); # days
  if n0 < nm
    qm = zeros( Ny , Nx ); # 1/s, mean qgpv initial field
    um = zeros( Ny , Nx ); 
    vm = zeros( Ny , Nx ); 
    psim = zeros( Ny , Nx ); 
    zetam = zeros( Ny , Nx ); 
    uzm = zeros( Ny , Nx ); 
    vzm = zeros( Ny , Nx ); 
  elseif n0 >= ns # import mean and statistics
    qm = h5read( readname , "qm" ); # 1/s, mean qgpv initial field
    um = h5read( readname , "um" );
    vm = h5read( readname , "vm" ); 
    psim = h5read( readname , "psim" ); 
    zetam = h5read( readname , "zetam" );
    uzm = h5read( readname , "uzm" ); 
    vzm = h5read( readname , "vzm" ); 
    mean_convergence = h5read( readname , "mean_convergence" );
  else # n0 >= nm, n0 < ns # import mean, but the other statistics are zero
    qm = h5read( readname , "qm" ); # 1/s, mean qgpv initial field
    um = zeros( Ny , Nx ); 
    vm = zeros( Ny , Nx ); 
    psim = zeros( Ny , Nx ); 
    zetam = zeros( Ny , Nx ); 
    uzm = zeros( Ny , Nx ); 
    vzm = zeros( Ny , Nx ); 
    mean_convergence = h5read( readname , "mean_convergence" );
  end 
  if nudge == 1 # Gaussian perturbations to kick off nonlinearity (add perturbation at 12 days)
  sigma = Lx/30.0;
  qperturb = exp(-((X-Lxcenter+10E5).^2.0/(1.0*sigma)^2.0+(Y+Lycenter+2E5).^2.0./(2.0*sigma^2.0))).*randn(size(qn)).*(maximum(qn)*1.25)-exp(-((X-Lxcenter+10.5E5).^2.0/(1.5*sigma)^2.0+(Y+Lycenter-1.75E5).^2.0./(2.0*sigma^2.0))).*randn(size(qn)).*(maximum(qn)*2.0)+exp(-((X-Lxcenter+10E5).^2.0/(2.0*sigma)^2.0+(Y-Lycenter+22E5).^2.0./(2.0*sigma^2.0))).*randn(size(qn)).*(maximum(qn)*2.0)-exp(-((X-Lxcenter+12E5).^2.0/(2.0*sigma)^2.0+(Y-Lycenter-22E5).^2.0./(2.0*sigma^2.0))).*randn(size(qn)).*(maximum(qn)*2.0);
  qn = qn + qperturb;
  end
  println( "\ninput .h5 file: "readname ); 
end # initialization set up


# =============================================================================
# time step initializations

# counters
nplot = 0; # counter for output plots
nwrite = 0; # counter for output files
n10step = 0; # counter for time steps
avg_step_time = 0.0; # average time per time step initialization

# time step intervals (4th-order Runge-Kutta) 
const dt = dt_days*3600.0*24.0; # s
const dt_2 = dt/2.0; # s
const dt_6 = dt/6.0; # s


# =============================================================================
# output plot dimensions 

# field initialization for time advancement
k = zeros( Ny , Nx ); 
Qn = zeros( Ny , Nx ); 
psin = zeros( Ny , Nx ); 
PSIn = zeros( Ny , Nx ); 


# =============================================================================
# FFTW plan initialization 

# dst plan initialization
qs = copy( qn ); 
FFTW.set_num_threads( Np );
Pcos = plan_r2r( qs , FFTW.REDFT10 , 1:2 , flags=FFTW.MEASURE ); # 2D cosine transform
Psin = plan_r2r( qs , FFTW.RODFT10 , 1:2 , flags=FFTW.MEASURE ); # 2D sine transform
Psinx = plan_r2r( qs , FFTW.RODFT10 , 2 , flags=FFTW.MEASURE ); # 1D sine transform in x
Psiny = plan_r2r( qs , FFTW.RODFT10 , 1 , flags=FFTW.MEASURE ); # 1D sine transform in y 
Psininv = plan_r2r( qs , FFTW.RODFT01 , 1:2 , flags=FFTW.MEASURE ); # 2D inverse sine transform
Pcosinv = plan_r2r( qs , FFTW.REDFT01 , 1:2 , flags=FFTW.MEASURE ); # 2D inverse cosine transform
Pcosxinv = plan_r2r( qs , FFTW.REDFT01 , 2 , flags=FFTW.MEASURE ); # 1D inverse cosine transform in x
Pcosyinv = plan_r2r( qs , FFTW.REDFT01 , 1 , flags=FFTW.MEASURE ); # 1D inverse cosine transform in y 

# is this part necessary:
Ux = zeros( size( qs ) ); Uy = zeros( size( qs ) ); # example of shift for derivative
Ux[:,2:Nx] = (Psinx*qs)[:,1:Nx-1]; Ux = Ux./(2.0*Float64(Nx)); # example of shift for derivative
Uy[2:Ny,:] = (Psiny*qs)[1:Ny-1,:]; Uy = Uy./(2.0*Float64(Ny)); # example of shift for derivative 
Pdycosinv = plan_r2r(Uy.*Lc,FFTW.REDFT01,1,flags=FFTW.MEASURE); # 1D inverse cosine transform in x
Pdxcosinv = plan_r2r(Ux.*Kc,FFTW.REDFT01,2,flags=FFTW.MEASURE); # 1D inverse cosine transform in y
Pdx2dy2sininv = plan_r2r( -(Psin*qs).*( Ks.^2 + Ls.^2 ) , FFTW.RODFT01,1:2 , flags=FFTW.MEASURE ); # 2D inverse sine transform
