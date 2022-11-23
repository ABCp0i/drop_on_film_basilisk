#!/bin/sh
#SBATCH -J DropTest2_t31        # Job name
#SBATCH -o drop.o%j        # Name of stdout output file
#SBATCH -e drop.e%j        # Name of stderr error file
##SBATCH -p skx-dev      # Queue (partition) name
#SBATCH -p skx-normal      # Queue (partition) name
#SBATCH -N 1              # Total # of nodes 
#SBATCH -n 1             # Total # of mpi tasks n=N*24
#SBATCH -t 08:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=jiny9@rpi.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A Evap-Drop-Breakup   # Allocation name (req'd if you have more than 1)

# Other commands must follow all #SBATCH directives...

module list
pwd
date

max_level=11
L=8

t_out=0.002
t_end=4.990
h=4.

rhog=0.001227
mug=4.601463e-5
mul=.002487
u0=1.826726
e=.013453

femax=0.001
uemax=0.001

spherewidth=.912871
sphereheight=1.2

export OMP_NUM_THREADS=48

#Please change the path to the executable file!
$HOME/local_compile_ellipse/impact $max_level $L $u0 $t_out $t_end $e $h $rhog $mul $mug $femax $uemax $spherewidth $sphereheight  2>log 1>out

echo
echo "Job finished at `date`"
