#!/bin/bash

# Constants

# M_PI=3.14159265358979323846     # Pi

# Parameters have to be written in command line in the order specified in the "varlist" of the .cpp file

# System
N=32                   # number of lattice points in each direction
DIM=3                  # dimension of lattice (arbitrary)
Nx=1                   # number of lattice points in x-direction
Ny=1                   # number of lattice points in y-direction
Nz=1                    # number of lattice points in z-direction (use by default Nx,Ny,Nz unless all set to 1)
ls=1.0            # Lattice spacing

# GPE parameters
m=0.5                # parameter of the sigma_x term in Hamiltonian
g=1            # coupling constant of the phi^4 term
gn=0                    # Coupling constant of the phi^n term
glong=0                 # Prefactor of long-range interactions
expphi=6                # Power n of phi^n term in action
explong=6               # Power of power-law decaying interactions
h=0.01                 # Radial part of external field, phase arbitrary.
hDP=0.001				# Variance of external perturbation for density-phase commutators

# Initial Conditions
IC=1                    # Type of initial condition. 1: box, 2: homog Density
amp=50                # Amplitude of box initial conditions
Qmax=1                # Maximal momentum of the box initial conditions
Qmin=-0.1               # Minimum momentum for box IC

# Dynamics
tmax=100               # differential equation will be solved until this maximal time    # $( echo "20.0/10000.0" | bc )
dt=0.1                  # constant step size for numerical solver of differential equation
iter=10                  # Average over "iter" number of runs within program
tw=-50                   # "Waiting time": after time tw system is doubled and one part of it evolves with external field h. NEGATIVE -> no doubling
twF=-200                  # At time tw save state of the system and measure F(t,twF) for t>twF
twDPrho=50				# At time twDPrho system is tripled, and the two copies are perturbed with hDP to compute rho_DP(t,twDPF) for t>twDPF
twDPF=-50				# At time twDPF save phase and density of system and measure F(t,twDPF) for t>twDPF

# Sonstiges
bsize=1024              # size of buffer to store name of output file
tnum=2                  # Number of threads to be used in parallel
seed=0                 # seed for random number generator. If negative seed -> extract from current time.

# Output
nfile=8                 # This number will be appended to the name of the output files
Nbins=250               # Number of bins for the momenta

spstep=-10              # print spectrum after each spstep time interval. NEGATIVE -> no output
bogstep=-100             # print Bogoliubov spectrum after each bogstep time interval. NEGATIVE -> no output. Only for phi^4
opstep=-2               # print op file - containing magnetization, total N, etc. - after each opstep time interval. NEGATIVE -> no output
npstep=-10               # print np file - containing phi(x)^2, phi(x)^4, etc. - after each npstep time interval. NEGATIVE -> no output
phi0step=-10             # print phi(t,p=0) after each phi0step time interval. NEGATIVE -> no output
Fstep=-0.2                # output F(t,twF) every Fstep
CRstep=-0.2                # Autocorrelation and response functions output. NEGATIVE -> no output

DPspstep=-1				# Output density-phase spectrum in intervals of DPspstep
DPFstep=-0.02				# Output density-phase anticommutator F in intervals of DPFstep
DPrhostep=0.02				# Output density-phase commutator rho in intervals of DPrhostep

# Settings linear response
CRmodeResp=2            # Choose which mode to perturb and output for linear response. 0: p=0 mode, 1: integrated autoresponse, 2: full p-dep, 3: full p-dep with uniform field
cphysMom=0               # Output list of physical momenta. 0: No, 1: Yes
Pcutoff=4              # Cutoff in momentum space for perturbation field

# Folders
folder=data             # Folder to save the data

# Note: careful with bash arithmetics. Bash does not perform double divisions automatically. For this use "bc" or use a C++ expression parser.


./bosons_F $N $DIM $Nx $Ny $Nz $ls $m $g $gn $glong $expphi $explong $h $hDP $IC $amp $Qmax $Qmin $tmax $dt $iter $tw $twF $twDPrho $twDPF $bsize $tnum $seed $nfile $Nbins $spstep $bogstep $opstep $npstep $phi0step $Fstep $CRstep $DPspstep $DPFstep $DPrhostep $CRmodeResp $cphysMom $Pcutoff $folder








