#!/bin/bash
#SBATCH --account=rdeshmukh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --job-name=qs_long
#SBATCH --output=logs/qs_long_%j.out
#SBATCH --error=logs/qs_long_%j.err
# =====================================================================
#  A long single-point run (qs_run_long). One job, not a fan-out.
#
#  Required in the environment:
#      QS_ROOT   framework root directory
#      QS_CASE   case id
#      QS_STUDY  study id  (must have "mode": "long")
#  Optional:
#      QS_POINT  index of a single point in the study's list; omit for all
#
#  Submit with:
#      export QS_ROOT=$PWD QS_CASE=<case> QS_STUDY=<study>
#      sbatch --exclude=ec77 --export=ALL,QS_ROOT,QS_CASE,QS_STUDY slurm/qs_long.slurm
#
#  Memory and time are larger than the sweep's because a long run keeps the
#  full modal history in memory and may stream instantaneous fields to disk.
#  Fields are written one file per simulated second, so a job that dies part
#  way through still leaves every completed second behind.
# =====================================================================

module load matlab

cd "${QS_ROOT:?QS_ROOT not set}"
mkdir -p logs

echo "host  : $(hostname)"
echo "case  : ${QS_CASE}"
echo "study : ${QS_STUDY}"
echo "point : ${QS_POINT:-all}"

if [ -n "${QS_POINT}" ]; then
  ARG=",${QS_POINT}"
else
  ARG=""
fi

matlab -nodisplay -singleCompThread -r \
  "run('${QS_ROOT}/qs_startup.m'); qs_run_long('${QS_CASE}','${QS_STUDY}'${ARG}); exit;"
