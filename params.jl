# Initial conditions for a quasi-gestrophic one layer ocean flow solver.
# Bryan Kaiser
# 6/4/17

# v.6

# ============================================================================
# physical parameters

# beta plane parameters
Lx = 3e6 # m, longitudinal domain size
Ly = 2.0*Lx # m, latitudinal domain size
const Nx = 2^11 # longitudinal grid resolution (must be an even integer)
const Ny = 2^12 # latitudinal grid resolution (must be an even integer)
const f = 1e-4 # 1/s
const beta = 2e-11 # 1/(m*s), beta (gradient of planetary vorticity)
H = 1e3 # m, layer depth 
LR = 0.0 # m, Rossby radius of deformation

# forcing parameters
taux = 4.5 # N/m^2 or Pa, longitudinal wind stress (0.1 as reported by 1968-1996 NCEP 
# Reanalysis, Kalnay et al. 1996)
gyre_config = 2 # 1 = single gyre longitudinal wind stress, 2 = double gyre 
# 3 = two gyre

# dissipative parameters
const r = 0.0; # 1/s, Rayleigh bottom friction (typically ~5E-7 1/s for the ocean)
const kappa0 = 920.0 # m/s^2, turbulent viscosity 
kappa_profile = 1 # 1 is constant, 2 is enhanced boundary
# d_kappa = dI/Lx/sqrt((dI/dM)^3) only matters for enhanced boundary
# Estimated to be between 10^3 and 10^4 m^2/s, with large variability
# Zhurbos & Oh (2004), from L.Talley "Descriptive Oceanography" book p.193. 

# time steppin'
dt_days = round(1.0/800.0,6)
const Nt = 1e9 # number of time steps beyond n0
const n0 = 6400 # initial time step (NOTE: if this is not zero the corresponding 
# .h5 input file to the time step specified must be in the /output folder)

nudge = 0

# =============================================================================
# computation parameters

# set the number of threads equal to the number of CPU cores on your machine:
Np = 4

# set the output path for plots and h5 files:
output_path = "/home/bryan/Documents/data/gyre/gyre_sim5_Re5_double_hires"
#output_path = "/nobackup1/bkaiser/gyre/gyre_sim5_Re5_double_hires"

# print out document of chosen parameters
print_flag = 0

# compute the running mean q, psi
const nm = 80000 # once past this time step, begin taking the mean / testing convergence 
const ns = 320000 # once past this time step, begin computing statistics

# output plots
const plot_flag = 0 # enter 1 for .png output plots or 0 for no output plots
const nplot_interval = 10 # number of time steps between output plots
const nplot_start = 10 # number of time steps AFTER n0 before starting output plots
const Nres = 200 # number of contour levels for output plots (resolution)

# output .h5 files
const write_flag = 1 # enter 1 for .h5 file output or anything else for no output files
const nwrite_interval = 800 # number of time steps between output files
const nwrite_start = 800 # number of time steps AFTER n0 before starting output plots

