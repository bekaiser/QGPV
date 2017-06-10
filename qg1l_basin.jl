# A quasi-gestrophic one layer ocean flow solver for flow in a closed basin.
# Bryan Kaiser
# 6/2/17

# v.6

# <u'>=0 if 1/T int_0^T b_u'u'(tau) dtau = 0 as T -> infinty Slutsky (1938).
# this is the "law of large numbers" or "ergodic theorem"
# need the mean to be the true mean...maybe import it from an initial run.

# add choice of colorbar bounds

# add nudge at n=1200
# add perturbation at 12 days 

# this package is needed for computation of derivatives:
using Base.FFTW 

# this package is only required for writing .h5 output files:
using HDF5 

# import simulation parameters & functions:
include("./params.jl")
include("./functions.jl")

# initialize the grid, variables, parameters, and FFTW plan: 
include("./initialization.jl")

# display parameters in terminal window:
if print_flag == 1
  include("./print.jl")
end

# the following packages are necessary for result plots:
if plot_flag == 1
  using PyPlot 
  using PyCall 
end

# main loop
for n = n0:Nt+n0
  tic()

  # terminal output:
  if ( ( n - n0 ) / 10 ) == n10step 
    time_elapsed = ( n - n0 )*dt_days + time_elapsed_start;
    println("\ntime: $( time_elapsed ) days") 
    println("$( round( (avg_step_time) , 3 ) ) s per time step")
    if n >= (nm + 1)
      println("max [ <q>_{n+1} - <q>_{n} ] = $mean_convergence")
    end
    n10step = n10step+1;
  end 

  # output plots to .png for a movie:
  if plot_flag == 1 
    if ( n - n0 ) >= nplot_start # start output
      if ( ( n - nplot_start - n0 ) / nplot_interval ) == nplot # plot interval
        outplot( qn , Int64(n), time_elapsed );
        nplot = nplot + 1;
      end 
    end 
  end 

  # output .h5 files:
  if write_flag == 1 
    if ( n - n0 ) >= nwrite_start # start output
      if ( ( n - nwrite_start - n0 ) / nwrite_interval ) == nwrite # output interval
        if n >= ns # begin computing statistics
          outwrite_stat( qn , qm , Int64(n) , time_elapsed , um , vm , psim , zetam , uzm , vzm )
        elseif n >= nm # begin averaging  
          outwrite_mean( qn , qm , Int64(n) , time_elapsed )
        else # just output the field
          outwrite( qn , Int64(n) , time_elapsed )
        end 
        nwrite = nwrite + 1; 
      end 
    end
  end 

  # time advancement:
  q1 = qn; k1 = RK4( q1 , Qn, psin, PSIn ) # first RK4 coefficient
  q2 = qn + k1.*dt_2; k2 = RK4( q2 , Qn, psin, PSIn ) # second RK4 coefficient
  q3 = qn + k2.*dt_2; k3 = RK4( q3 , Qn, psin, PSIn ) # third RK4 coefficient
  q4 = qn + k3.*dt; k4 = RK4( q4 , Qn, psin, PSIn ) # fourth RK4 coefficient
  qn = qn + ( k1 + k2.*2.0 + k3.*2.0 + k4 ).*dt_6; # 1/s

  # wall time per time step:
  avg_step_time = toq()/Float64( n - n0 + 1 ) + Float64( n - n0 )/Float64( n - n0 + 1 )*avg_step_time; 

  # compute statistics:
  if n >= nm 
    qm0 = copy(qm);
    qm = qn/Float64( n - nm + 1 ) + Float64( n - nm )/Float64( n - nm + 1 )*qm; # mean qgpv
    mean_convergence = maximum(abs(qm-qm0));
    if n >= ns
      psi,PSI = inversion( (Psin*qn)./(4.0*Float64(Nx)*Float64(Ny)) )
      um, vm, psim, zetam, uzm, vzm = statcalc( psi , um , vm , psim , zetam , uzm , vzm , Int64(n) );
    end
  end 

  # check for NaNs:
  if sum(isnan(qn)) >= 1 
    println("Simulation stopped: NaNs detected!"); break
  end  

end # time advancement loop




