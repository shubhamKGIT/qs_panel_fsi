#!/bin/bash
#SBATCH --account=rdeshmukh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=4G
#SBATCH --job-name=qs_sweep
#SBATCH --output=logs/qs_sweep_%j.out
#SBATCH --error=logs/qs_sweep_%j.err
# =====================================================================
#  One slice of a qs_framework sweep.
#
#  Submitted by submit_sweep.sh as N independent sbatch calls (NOT a job
#  array). That is deliberate: the array route lost whole blocks of jobs on
#  this cluster, and N plain submissions with the slice number exported has
#  been verified to dispatch every task on the first try.
#
#  Required in the environment (submit_sweep.sh exports them):
#      QS_ROOT   framework root directory
#      QS_CASE   case id
#      QS_STUDY  study id
#      SLICE     this task's slice number, 1..NSLICE
#      NSLICE    total number of slices
# =====================================================================

module load matlab

cd "${QS_ROOT:?QS_ROOT not set}"
mkdir -p logs

echo "host   : $(hostname)"
echo "case   : ${QS_CASE}"
echo "study  : ${QS_STUDY}"
echo "slice  : ${SLICE} of ${NSLICE}"

# -singleCompThread: one core per task, so tasks sharing a node do not
# oversubscribe it.
matlab -nodisplay -singleCompThread -r \
  "run('${QS_ROOT}/qs_startup.m'); qs_run_sweep('${QS_CASE}','${QS_STUDY}',${SLICE},${NSLICE}); exit;"
