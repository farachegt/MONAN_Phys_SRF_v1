#!/bin/bash
######################################################
##PARTITION AVAIL  TIMELIMIT  NODES  STATE NODELIST
##batch        up    4:00:00     33   idle n[01-33]
##proc         up    4:00:00      3   idle proc[01-03]
##OPER         up    3:00:00     15   idle n[01-15]
##PESQ1        up   12:00:00      8   idle n[16-23]
##PESQ2        up   12:00:00      8   idle n[24-31]
##PESQ3        up    2:00:00      2   idle n[32-33]
##DEBUG        up      15:00      1   idle proc01
#####################################################

cat << Eof > submit.bash
#!/bin/bash
#SBATCH --job-name=MONANPDv
#SBATCH --nodes=${1}
#SBATCH --partition=${2}
#SBATCH --tasks-per-node=${3}
#SBATCH --ntasks=${4}
##SBATCH --time=02:00:00
#SBATCH --time=04:00:00
##SBATCH --mem=64000M
#SBATCH --output=/home/guilherme.farache/logs/%j.log
####SBATCH --exclusive

cd \$SLURM_SUBMIT_DIR
echo \$SLURM_SUBMIT_DIR

module purge
module load ohpc
module unload openmpi4
module load phdf5
module load netcdf
module load netcdf-fortran
module load mpich-4.0.2-gcc-9.4.0-gpof2pv
module load hwloc
module list

export OMP_NUM_THREADS=1
export OMPI_MCA_btl_openib_allow_ib=1
export OMPI_MCA_btl_openib_if_include="mlx5_0:1"
export PMIX_MCA_gds=hash

export NETCDF=/mnt/beegfs/monan/libs/netcdf
export PNETCDF=/mnt/beegfs/monan/libs/PnetCDF

MPI_PARAMS="-iface ib0 -bind-to core -map-by core"
export MKL_NUM_THREADS=1
export I_MPI_DEBUG=5
export MKL_DEBUG_CPU_TYPE=5
export I_MPI_ADJUST_BCAST=12 ## NUMA aware SHM-Based (AVX512)

export NETCDF=/mnt/beegfs/monan/libs/netcdf
export PNETCDF=/mnt/beegfs/monan/libs/PnetCDF
# PIO is not necessary for version 8.* If PIO is empty, MPAS Will use SMIOL
export PIO=
module list
echo "============================="

ulimit -s unlimited
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export I_MPI_DEBUG=5
export MKL_DEBUG_CPU_TYPE=5
#export I_MPI_ADJUST_BCAST=12 ## NUMA aware SHM-Based (AVX512)
export I_MPI_FABRICS=shm:ofi

mpirun -genvall ${5} &> ${5}.log 
#mpirun -env MKL_DEBUG_CPU_TYPE=5 -env UCX_NET_DEVICES=mlx5_0:1 -genvall ${5} &> ${5}.log

Eof


