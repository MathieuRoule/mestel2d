#!/bin/bash

##############################
# Dump all files in the same
# directory
##############################

# Move to the directory where the script is located
script_dir="$(dirname "${BASH_SOURCE[0]}")"
cd ${script_dir}

# Initial conditions, simulation outputs
outpath="./data/"

# Plots
outplots="${outpath}plots/"
mkdir -p ${outplots}

##############################
# Generate initial conditions
#############################

icfile="${outpath}ics.h5"

# Activate your favorite python environment if you want to
#conda activate myenv

# Install required packages
pip -q install -r ics_generation/requirements.txt

# Run IC generator

python ics_generation/genics.py -N=100000 -q=6 -seed=1 ${icfile}

##############################
# Simulation
##############################

# Parameters
dt=0.01		# Time step
nsteps=2000 # Total nb of time steps
xmax=20.0	# Box size (2*xmax)^2
nx=128		# Cells size (2*xmax/nx)^2
kernel="kuzmin" # Gravity softening kernel(plummer/kuzmin)
eps=0.16	# Softening length
nring=2048	# Radial grid (if -1, default to 20*nx)
nphi=1024	# Azimuthal grid
fractive=1.0 # Active fraction
# Output frequencies (number of step between 2 outputs)
# Initial state is allways dumped (min is 1 and don't forget to count it)
fulldumpevery=100 # Full positions/velocities
partdumpevery=10 # Summary statistics

# Coefficient computation parameters (summary stat)
basis="gate" # Radial function shape (gate/lognorm)
rb=3.0
sigmar=2.0
lmax=2

# Others
verbose=1

# Files and directories
outfile="${outpath}out.h5"
rm -f ${outfile} # Preemptive deletion

# Run (m=2 only -- m2 option)
../mestel2d -xmax=${xmax} -nx=${nx} \
	        -dt=${dt} -partdumpevery=${partdumpevery} -fulldumpevery=${fulldumpevery} \
			-nstep=${nsteps} \
	        -nring=${nring} -nphi=${nphi} -kernel=${kernel} -eps=${eps} \
	        -fractive=${fractive} -m2 \
            -basis=${basis} -rb=${rb} -sigmar=${sigmar} -lmax=${lmax} \
            -verbose=${verbose} ${icfile} ${outfile}


##############################
# Post-processing
##############################

# # Install required packages
pip -q install -r postprocess/requirements.txt

# Plot density
vmax=3 # Maximum density value for the plot
python postprocess/plot_density.py -vmax=${vmax} ${outfile} ${outplots}

# If ffmpeg is installed, make a movie of the frames
if command -v ffmpeg &> /dev/null; then
	# Store location/name
	movie="${outpath}Zang4_instability.mp4"
    # Create movie
    ffmpeg -loglevel quiet -r 10 -i ${outplots}frame_%03d.png -q:v 0 -pix_fmt yuv420p ${movie}

	# Play movie
	if command -v xdg-open &> /dev/null; then
		xdg-open ${movie}
	elif command -v open &> /dev/null; then
		open ${movie}
	else
		echo "Movie created at ${movie}"
	fi
else
    echo "ffmpeg is not installed, install it if you want to create a movie." 
	echo "You can still find the frames in /data/plots/."
fi