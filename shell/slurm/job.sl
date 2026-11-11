#!/bin/bash

#####################################################################
##  Slurm job for a UKRmol-Scripts run.                            ##
##  The intention for it to be submitted by the submission script. ##
#####################################################################

#SBATCH --output sbatch_out/%x.o%j
#SBATCH --error sbatch_out/%x.e%j
#SBATCH --partition tcourt
#SBATCH --nodes 1
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user EMAIL@DOMAIN

die() { echo "ERROR: $*" >&2 ; exit 1 ; }

for var in UKR_GEOMS UKR_SYM UKR_RANKS UKR_THREADS UKR_MODEL UKR_CONFIG UKR_GEOM UKR_DIRS UKR_MAIN ; do
  [[ -n "${!var}" ]] || die  "$var is not set"
done

for var in UKR_MODEL UKR_CONFIG UKR_GEOM UKR_DIRS UKR_MAIN ; do
 [[ -e "${!var}" ]] || die "$var: ${!var}: file not found"
done

name="$SLURM_JOB_NAME"

# -- ensure allocation match
nproc=$(( UKR_GEOMS * ( UKR_SYM> 0 ?UKR_SYM : 1 ) ))
(( SLURM_NTASKS == nproc * UKR_RANKS && SLURM_CPUS_PER_TASK == UKR_THREADS )) ||
  die "allocation (ntasks=$SLURM_NTASKS, cpus-per-task=$SLURM_CPUS_PER_TASK) mismatch" \
      "geoms=$UKR_GEOMS sym=$UKR_SYM ranks=$UKR_RANKS threads=$UKR_THREADS"

# -- threads for non-MPI programs
export OMP_NUM_THREADS=$(( SLURM_NTASKS * SLURM_CPUS_PER_TASK / nproc ))
export OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS" MKL_NUM_THREADS="$OMP_NUM_THREADS"

#######################################################################################
# -- the details here may need to be modified on a different system
runner=srun
plugin=cray_shasta # <-- don't let molpro boss us around
#######################################################################################

if ((  UKR_RANKS > 1 || UKR_SYM == 0 )) ; then
  luncher="${runner} --exact --mpi=$plugin"
  luncher+=" --export=ALL,OMP_NUM_THREADS=${UKR_THREADS},OPENBLAS_NUM_THREADS=${UKR_THREADS},MKL_NUM_THREADS=${UKR_THREADS}"
  luncher+=" -n ${UKR_RANKS} -c ${UKR_THREADS}"
else
  luncher="" # forlegacy scatci and serial rsolve
fi

#######################################################################################
# -- also critically system-dependent are the module loads/unloads, etc
#######################################################################################
# CRIANN:
source /opt/cray/pe/cpe/24.07/restore_lmod_system_defaults.sh
module load PrgEnv-gnu cray-mpich
module unload cray-libsci
source /soft/RH9/2025032/lomc/ukrmolenv
module load atomic_simu/lomc-molpro/2025.4
module list
# -- greedy molpro wants to start itself on all tasks >:c
#    Make all `srun`s take `exact`ly the cores it asks for
export MOLPRO_OPTIONS="-n ${UKR_RANKS}"
export SLURM_EXACT=1

export TMPDIR="$LOCAL_WORK_DIR"/tmp
mkdir -p "$TMPDIR"
cd "$LOCAL_WORK_DIR" || die "couldn't change directory to $LOCAL_WORK_DIR"

# -- copy perl files to cwd
cp "$UKR_MAIN" "$UKR_DIRS" "$UKR_MODEL" "$UKR_GEOM" . || die "failed to copy the perl files"
main_file="${UKR_MAIN##*/}"
dirs_file="${UKR_DIRS##*/}"
model_file="${UKR_MODEL##*/}"
geom_file="${UKR_GEOM##*/}"

# -- fill config template
config_file="config.${name}.pl"
sed -e "s|<<SUFFIX>>|${name}|g" \
    -e "s|<<MPI_LUNCHER>>|${luncher}|g" \
    -e "s|<<PARALLEL_GEOM>>|${UKR_GEOMS}|g" \
    -e "s|<<PARALLEL_SYMM>>|${UKR_SYM}|g" \
    "$UKR_CONFIG" > "$config_file"
grep -n '<<[A-Z_]*>>' "$config_file" && die "detected unfilled template items in $config_file"

echo "=== $name ==="
echo "    $UKR_GEOMS geometries in parallel, parallel_symm=$UKR_SYM, $UKR_RANKS ranks x $UKR_THREADS threads each"
echo "    OMP_NUM_THEADS for non-MPI programs: $OMP_NUM_THREADS"
echo "    mpi launcher: ${luncher:-<none; legacy scatci/rsolve>}"
grep -nE "^\s*'(suffix|scattering|mpi_scatci|mpi_rsolve|parallel_geom|parallel_symm)',"  "$config_file"

start="$(date +%s)"
perl "$main_file" "$dirs_file" "$config_file" "$model_file" "$geom_file"
exitstatus="$?"
end="$(date +%s)"
echo "${name} wall time: $(( end - start )) s, exit status $exitstatus"

exit "$exitstatus"
