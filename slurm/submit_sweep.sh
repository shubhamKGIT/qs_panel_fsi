#!/bin/bash
# =====================================================================
#  Submit a qs_framework sweep as N independent SLURM jobs.
#
#    bash slurm/submit_sweep.sh <case_id> <study_id> [nslice]
#
#  e.g.
#    bash slurm/submit_sweep.sh heated_c1_NoSBLI_Periodic sweep_coarse_v1 48
#
#  What it does, in order:
#    1. initialises the study (once, on the login node) so every slice reads
#       the same point list instead of racing to create one
#    2. submits N plain sbatch jobs with SLICE/NSLICE exported
#    3. prints the job ids and the sacct command to check where they landed
#
#  --exclude=ec77 is passed on the command line, not baked into the .slurm
#  file: node ec77 kills every job it receives at launch (1 s, "CANCELLED by
#  0"), and hardcoding node names into files makes them go stale. Drop the
#  flag once an admin has fixed the node.
# =====================================================================
set -euo pipefail

CASE="${1:?usage: submit_sweep.sh <case_id> <study_id> [nslice]}"
STUDY="${2:?usage: submit_sweep.sh <case_id> <study_id> [nslice]}"
NSLICE="${3:-48}"

QS_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXCLUDE="${QS_EXCLUDE:---exclude=ec77}"

cd "$QS_ROOT"
mkdir -p logs

echo "=== initialising study (login node) ==="
module load matlab
matlab -nodisplay -nosplash -r \
  "run('${QS_ROOT}/qs_startup.m'); qs_init_study('${CASE}','${STUDY}'); exit;" \
  | tee "logs/init_${CASE}_${STUDY}.log"

echo
echo "=== submitting ${NSLICE} slices ==="
IDS=()
for k in $(seq 1 "$NSLICE"); do
  JID=$(sbatch --parsable $EXCLUDE \
        --job-name="qs_${STUDY}_${k}" \
        --export=ALL,QS_ROOT="$QS_ROOT",QS_CASE="$CASE",QS_STUDY="$STUDY",SLICE="$k",NSLICE="$NSLICE" \
        "$QS_ROOT/slurm/qs_sweep.slurm")
  IDS+=("$JID")
  echo "  slice ${k}/${NSLICE} -> job ${JID}"
done

JOBLIST=$(IFS=,; echo "${IDS[*]}")
echo
echo "submitted: $JOBLIST"
echo
echo "ALWAYS check where they landed before trusting the run:"
echo "  sacct -j $JOBLIST -X --format=JobID%12,State%20,Elapsed,NodeList%20"
echo "  sinfo -R    # drained/down nodes and why"
echo
echo "A job that shows Elapsed 00:00:01 and 'CANCELLED by 0' never started."
echo "Re-submit just those slices; finished points are skipped on resume."
echo
echo "When they are done:"
echo "  matlab -nodisplay -r \"run('$QS_ROOT/qs_startup.m'); qs_merge('$CASE','$STUDY'); qs_plot_map('$CASE','$STUDY'); exit\""
