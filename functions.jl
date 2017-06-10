# Functions for qg1l_basin.jl, including plan_dst and de-aliasing
# Bryan Kaiser
# 6/2/17

# v.6

# =============================================================================
# function declaration 

function meshgrid{T}( vx::AbstractVector{T} , vy::AbstractVector{T} )
    m, n = length(vy), length(vx)
    vx = reshape(vx, 1, n); vy = reshape(vy, m, 1);
  return (repmat(vx, m, 1), repmat(vy, 1, n))
end

function bandavg( U::Array{Float64,1} , nb::Int64 )
# nb = number of bands to average over (window size)
# U = spectral variable to be band-averaged, low to high frequency.
    N = length(U);
    Nb = Int64(floor(N/nb)); Ub = zeros(Nb);
    r = Int64(N-Nb*nb); # remainder 
    for j = 1:Nb
      Ub[j] = sum(U[(r+1+(j-1)*nb):(r+j*nb)])/Float64(nb);
    end
  return Ub
end

function makemask( Ks::Array{Float64,2} , Ls::Array{Float64,2} , Kmags::Array{Float64,2} )
# mask for 2/3 rule padding for de-aliasing a quadratic signal via fft
  mask0 = ones( size( Kmags ) );
  for j = 1:length( Kmags )
    if abs( Ks[j] ) >= 2.0/3.0*(Float64(Nx)*pi/Lx);
      mask0[j] = 0.0;
    end
    if abs( Ls[j] ) >= 2.0/3.0*(Float64(Ny)*pi/Ly);
      mask0[j] = 0.0;
    end
  end
  return mask0
end

function power_spectrum_2D(S::Array{Float64,2},Kmag::Array{Float64,2})
# takes 2D spectra and generates a 1D power spectra for plotting. 
# convert the gridded wavenumber magnitudes to a vector, remove 
# repeated values, and sort:
    Kmag_vec_union = sort(union(vec(Kmag))); 
    S_vec = zeros(size(Kmag_vec_union)); # power spectrum
    for j = 1:length(Kmag_vec_union) # each wavenumber magnitude
      s = 0.0; count = 0.0;
      for n = 1:(size(Kmag,1)*size(Kmag,2)) # loop over Kmag, S
        if Kmag[n] == Kmag_vec_union[j];
          s = s + S[n];
          count = count+1.0;
        end
      end 
      S_vec[j] = s/count; # averaged magnitude
    end
  return S_vec, Kmag_vec_union;
end


function dealias( U::Array{Float64,2} , V::Array{Float64,2} , mask::Array{Float64,2} )
# 2/3 rule padding for de-aliasing a quadratic signal via dst
  return r2r(U.*mask,FFTW.RODFT01,1:2).*r2r(V.*mask,FFTW.RODFT01,1:2);
end

function inversion( Q::Array{Float64,2} ) 
    PSI = -Q.*((Ks_2+Ls_2+LR_m2).^(-1));
    psi = Pdx2dy2sininv*PSI;
  return psi,PSI
end

function advection( psi::Array{Float64,2} , Q::Array{Float64,2} )
    # relative and planetary advection of vorticity
    PSIx = zeros(size(Q)); PSIy = zeros(size(Q)); 
    U = zeros(size(Q)); V = zeros(size(Q)); 
    UQ = zeros(size(Q)); VQ = zeros(size(Q));
    # spectral velocity  
    PSIx[:,1] = zeros(Ny); PSIx[:,2:Nx] = (Psinx*psi)[:,1:Nx-1]; PSIx = PSIx./(2.0*Float64(Nx));
    PSIy[1,:] = zeros(Nx); PSIy[2:Ny,:] = (Psiny*psi)[1:Ny-1,:]; PSIy = PSIy./(2.0*Float64(Ny));  
    U[:,:] = -(Psin*(Pdycosinv*(PSIy.*Lc)))./(4.0*Float64(Nx)*Float64(Ny)); # dst spectral u
    V[:,:] = (Psin*(Pdxcosinv*(PSIx.*Kc)))./(4.0*Float64(Nx)*Float64(Ny)); # dst spectral v
    # spectral (de-aliased) advection:
    UQ[:,1] = zeros(Ny); UQ[:,2:Nx] = (Psinx*dealias(U,Q,mask))[:,1:Nx-1]; UQ = UQ./(2.0*Float64(Nx));
    VQ[1,:] = zeros(Nx); VQ[2:Ny,:] = (Psiny*dealias(V,Q,mask))[1:Ny-1,:]; VQ = VQ./(2.0*Float64(Ny)); 
    # flux conserving advection = Jacobian(psi,q), including planetary advection
    #adv = -(Pdxcosinv*(UQ.*Kc)+Pdycosinv*(VQ.*Lc)+(Pdxcosinv*(PSIx.*Kc)).*beta);
    #adv = -Pdxcosinv*(UQ.*Kc+PSIx.*Kc.*beta)-Pdycosinv*(VQ.*Lc); 
  return -(Pdxcosinv*(UQ.*Kc+PSIx.*Kc.*beta)+Pdycosinv*(VQ.*Lc)) # sign for RHS of equation
end

function diffusion( Q::Array{Float64,2} )
  return (Pdx2dy2sininv*(-Q.*(Ks_2+Ls_2))).*kappa  # sign for RHS of equation
end

function bottom_friction( PSI::Array{Float64,2} ) 
	return -(Pdx2dy2sininv*(-PSI.*(Ks_2+Ls_2))).*r # sign for RHS of equation
end

function RK4( q::Array{Float64,2} , Q::Array{Float64,2} , psi::Array{Float64,2} , PSI::Array{Float64,2} )
    Q[:,:] = (Psin*q)./(4.0*Float64(Nx)*Float64(Ny));
    psi[:,:],PSI[:,:] = inversion( Q ) # inversion, scalars = [LR_m2,beta,kappa,r,phi];   
  return advection( psi , Q ) + diffusion( Q ) + bottom_friction( PSI ) + wind; # Runge-Kutta coefficient
end

function outplot( q::Array{Float64,2} , n::Int64 , time_elapsed::Float64 )
  # output plots
  # psi
  psip = Pdx2dy2sininv*(-((Psin*q)./(4.0*Float64(Nx)*Float64(Ny))).*((Ks_2+Ls_2+LR_m2).^(-1.0)));
  fig1 = figure(figsize=(plotsize)) 
  CP1 = contourf(X./1000.0,Y./1000.0,psip,Nres,cmap="seismic") 
  xlabel("x (km)"); ylabel("y (km)"); title("psi"); colorbar(CP1)  
  if n <10 
    plotname = "$(output_path)/psi/00000$(convert(Int64,n)).png"; elseif n <100
    plotname = "$(output_path)/psi/0000$(convert(Int64,n)).png"; elseif n <1000
    plotname = "$(output_path)/psi/000$(convert(Int64,n)).png"; elseif n <10000
    plotname = "$(output_path)/psi/00$(convert(Int64,n)).png"; elseif n <100000
    plotname = "$(output_path)/psi/0$(convert(Int64,n)).png"; elseif n <1000000
    plotname = "$(output_path)/psi/$(convert(Int64,n)).png"
  end
  savefig(plotname,format="png"); close(fig1); 
  println("output plot "plotname);
  # q
  fig2 = figure(figsize=(plotsize)) # q 
  CP2 = contourf(X./1000.0,Y./1000.0,q,Nres,cmap="inferno") 
  xlabel("x (km)"); ylabel("y (km)"); title("q(t), t=$(round(time_elapsed,1)) days"); colorbar(CP2) 
  if n <10 
    plotname = "$(output_path)/q/00000$(convert(Int64,n)).png"; elseif n <100 
    plotname = "$(output_path)/q/0000$(convert(Int64,n)).png"; elseif n <1000
    plotname = "$(output_path)/q/000$(convert(Int64,n)).png"; elseif n <10000
    plotname = "$(output_path)/q/00$(convert(Int64,n)).png"; elseif n <100000
    plotname = "$(output_path)/q/0$(convert(Int64,n)).png"; elseif n <1000000
    plotname = "$(output_path)/q/$(convert(Int64,n)).png"
  end
  savefig(plotname,format="png"); close(fig2);
  println("output plot "plotname);	
end 
 
function outwrite( q::Array{Float64,2} , n::Int64 , time_elapsed::Float64 )
  # output hdf5 files
  if n <10 # filename
    writename = "$(output_path)/output/qg1l_00000$(convert(Int64,n)).h5"; elseif n <100 
    writename = "$(output_path)/output/qg1l_0000$(convert(Int64,n)).h5"; elseif n <1000
    writename = "$(output_path)/output/qg1l_000$(convert(Int64,n)).h5"; elseif n <10000
    writename = "$(output_path)/output/qg1l_00$(convert(Int64,n)).h5"; elseif n <100000
    writename = "$(output_path)/output/qg1l_0$(convert(Int64,n)).h5"; elseif n <1000000
    writename = "$(output_path)/output/qg1l_$(convert(Int64,n)).h5"
  end
  h5open(writename, "w") do file # alternatively, say "@write file A"
    write( file, "q" , qn );   
    write( file, "Nx" , Nx ); 
    write( file , "Ny" , Ny );
    write( file, "Lx" , Lx ); 
    write( file ,"Ly" , Ly );  
    write( file, "n" , n ); 
    write( file , "dt" , dt );
    write( file, "time_elapsed" , time_elapsed );
    write( file, "curltau" , curltau );
    write( file, "gyre_config" , gyre_config );
    write( file, "scalars" , scalars ); # = [ LR_m2 , beta , r ]; 
    write( file, "kappa_profile" , kappa_profile );
    write( file, "kappa" , kappa );
    write( file, "avg_step_time" , avg_step_time );
  end # h5open
  println("\noutput .h5 file "writename); 

  # add: scalars = [ LR_m2 , beta , kappa , r ];
end

function outwrite_stat( q::Array{Float64,2} , qm::Array{Float64,2} , n::Int64 , time_elapsed::Float64 , um::Array{Float64,2} , vm::Array{Float64,2} , psim::Array{Float64,2} , zetam::Array{Float64,2} , uzm::Array{Float64,2} , vzm::Array{Float64,2} )
  # output hdf5 files
  if n <10 # filename
    writename = "$(output_path)/output/qg1l_00000$(convert(Int64,n)).h5"; elseif n <100 
    writename = "$(output_path)/output/qg1l_0000$(convert(Int64,n)).h5"; elseif n <1000
    writename = "$(output_path)/output/qg1l_000$(convert(Int64,n)).h5"; elseif n <10000
    writename = "$(output_path)/output/qg1l_00$(convert(Int64,n)).h5"; elseif n <100000
    writename = "$(output_path)/output/qg1l_0$(convert(Int64,n)).h5"; elseif n <1000000
    writename = "$(output_path)/output/qg1l_$(convert(Int64,n)).h5"
  end
  h5open(writename, "w") do file # alternatively, say "@write file A"
    write( file, "q" , qn ); 
    write( file , "qm" , qm );  
    write( file, "Nx" , Nx ); 
    write( file , "Ny" , Ny );
    write( file, "Lx" , Lx ); 
    write( file ,"Ly" , Ly );  
    write( file, "n" , n ); 
    write( file , "dt" , dt );
    write( file, "time_elapsed" , time_elapsed );
    write( file, "um" , um ); 
    write( file , "vm" , vm );  
    write( file, "psim" , psim ); 
    write( file , "zetam" , zetam );
    write( file, "uzm" , uzm ); 
    write( file ,"vzm" , vzm ); 
    write( file, "curltau" , curltau );
    write( file, "gyre_config" , gyre_config );
    write( file, "scalars" , scalars ); # = [ LR_m2 , beta , r ]; 
    write( file, "kappa_profile" , kappa_profile );
    write( file, "kappa" , kappa );
    write( file, "mean_convergence" , mean_convergence );
    write( file, "avg_step_time" , avg_step_time );
  end # h5open
  println("\noutput .h5 file "writename); 
end

function outwrite_mean( q::Array{Float64,2} , qm::Array{Float64,2} , n::Int64 , time_elapsed::Float64 )
  # output hdf5 files
  if n <10 # filename
    writename = "$(output_path)/output/qg1l_00000$(convert(Int64,n)).h5"; elseif n <100 
    writename = "$(output_path)/output/qg1l_0000$(convert(Int64,n)).h5"; elseif n <1000
    writename = "$(output_path)/output/qg1l_000$(convert(Int64,n)).h5"; elseif n <10000
    writename = "$(output_path)/output/qg1l_00$(convert(Int64,n)).h5"; elseif n <100000
    writename = "$(output_path)/output/qg1l_0$(convert(Int64,n)).h5"; elseif n <1000000
    writename = "$(output_path)/output/qg1l_$(convert(Int64,n)).h5"
  end
  h5open(writename, "w") do file # alternatively, say "@write file A"
    write( file, "q" , qn ); 
    write( file , "qm" , qm );  
    write( file, "Nx" , Nx ); 
    write( file , "Ny" , Ny );
    write( file, "Lx" , Lx ); 
    write( file ,"Ly" , Ly );  
    write( file, "n" , n ); 
    write( file , "dt" , dt );
    write( file, "time_elapsed" , time_elapsed ); 
    write( file, "curltau" , curltau );
    write( file, "gyre_config" , gyre_config );
    write( file, "scalars" , scalars ); # = [ LR_m2 , beta , r ]; 
    write( file, "kappa_profile" , kappa_profile );
    write( file, "kappa" , kappa );
    write( file, "mean_convergence" , mean_convergence );
    write( file, "avg_step_time" , avg_step_time );
  end # h5open
  println("\noutput .h5 file "writename); 
end

function CFL( psi::Array{Float64,2} )
#psin = Pdx2dy2sininv*(-((Psin*qn)./(4.0*Float64(Nx)*Float64(Ny))).*((Ks_2+Ls_2+LR_m2).^(-1.0)));
# computes the velocities and checks the CFL criterion
  PSIx = zeros( Ny , Nx ); #PSIy[1,:] = zeros(Nx);
  PSIy = zeros( Ny , Nx ); #PSIx[:,1] = zeros(Ny); 
  PSIx[:,2:Nx] = (Psinx*psi)[:,1:Nx-1]; PSIx = PSIx./(2.0*Float64(Nx));
  PSIy[2:Ny,:] = (Psiny*psi)[1:Ny-1,:]; PSIy = PSIy./(2.0*Float64(Ny));  
  u = r2r(-(Psin*(Pdycosinv*(PSIy.*Lc)))./(4.0*Float64(Nx)*Float64(Ny)),FFTW.RODFT01,1:2);
  v = r2r((Psin*(Pdxcosinv*(PSIx.*Kc)))./(4.0*Float64(Nx)*Float64(Ny)),FFTW.RODFT01,1:2);
  CFL = maximum( [(abs(u).*(dt*Float64(Nx)/Lx)) ; (abs(v).*(dt*Float64(Ny)/Ly)) ] );
  if CFL >= 0.8
    println("WARNING: Courant number approaching unity, Co = $CFL");
  end
  return u,v
end

function statcalc( psi::Array{Float64,2} , um::Array{Float64,2} , vm::Array{Float64,2} , psim::Array{Float64,2} , zetam::Array{Float64,2} , uzm::Array{Float64,2} , vzm::Array{Float64,2} , n::Int64 )
# must turn CFL on to run statcalc()
# To compute the full vorticity budget, including the mean flux, eddy flux, beta flux, 
# lateral friction flux, and bottom friction flux, one must compute: 
# <u>, <v>, <psi>, <zeta>=nabla^2<psi>, <(u-<u>)*(zeta-<zeta>)>, <(v-<v>)*(zeta-<zeta>)>. 
  u,v = CFL( psi ); # zonal and meridonal velocities
  um = u./Float64(n-ns+1) + Float64(n-ns)/Float64(n-ns+1).*um; # mean zonal velocity
  vm = v./Float64(n-ns+1) + Float64(n-ns)/Float64(n-ns+1).*vm; # mean meridonal velocity
  psim = psi./Float64(n-ns+1) + Float64(n-ns)/Float64(n-ns+1).*psim; # mean meridonal velocity
  zeta = (Pdx2dy2sininv*(-(Psin*psi)./(4.0*Float64(Nx)*Float64(Ny)).*(Ks_2+Ls_2))); # relative vorticity
  zetam = zeta./Float64(n-ns+1) + Float64(n-ns)/Float64(n-ns+1).*zetam; # mean relative vorticity
  #uz = (u-um).*(zeta-zetam); # aliased
  UU = (Psin*(u-um))./(4.0*Float64(Nx)*Float64(Ny));
  ZZ = (Psin*(zeta-zetam))./(4.0*Float64(Nx)*Float64(Ny));
  uz = dealias(UU,ZZ,mask);
  uzm = uz./Float64(n-ns+1) + Float64(n-ns)/Float64(n-ns+1).*uzm; # zonal eddy vorticity flux
  #vz = (v-vm).*(zeta-zetam); # aliased
  VV = (Psin*(v-vm))./(4.0*Float64(Nx)*Float64(Ny));
  vz = dealias(VV,ZZ,mask);
  vzm = vz./Float64(n-ns+1) + Float64(n-ns)/Float64(n-ns+1).*vzm; # zonal eddy vorticity flux
  return um, vm, psim, zetam, uzm, vzm
end




