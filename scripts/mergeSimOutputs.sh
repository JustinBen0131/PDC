#!/usr/bin/env bash
###############################################################################
# mergeSimOutputs.sh
#
# Two-step merging workflow for many “PositionDep_sim_chunk*.root” files.
#
#   1)  ./mergeSimOutputs.sh condor [test|firstHalf|asManyAsCan]
#       ▸ Split the chunk list into groups of N files and launch Condor jobs
#         that run ‘hadd’ in parallel.
#       ▸ “test”         – submit **one** job (first 100 files only).
#       ▸ “firstHalf”    – submit **half** the jobs (first 30 000 files ⇒
#                          300 partial outputs rather than 600).
#       ▸ “asManyAsCan”  – split the list into chunks of 1 000 files and submit
#                          every possible job; the final sub-list may contain
#                          <1 000 files.
#       ▸ (no sub-mode)  – submit all jobs with the default 100-file grouping.
#       ▸ Before any submission, previously produced partial files
#         (chunkMerge_*.root) **and** the final PositionDep_sim_ALL.root
#         in OUTPUT_DIR are deleted so the run starts 100 % fresh.
#
#   2)  ./mergeSimOutputs.sh addChunks [condor]
#       ▸ Locally (default) or **via one Condor job** (“condor” sub-mode)
#         merges *whatever* partial outputs exist in OUTPUT_DIR into
#         the final PositionDep_sim_ALL.root.
#
# USAGE:
#       ./mergeSimOutputs.sh condor
#       ./mergeSimOutputs.sh condor test
#       ./mergeSimOutputs.sh condor firstHalf
#       ./mergeSimOutputs.sh condor asManyAsCan
#       ./mergeSimOutputs.sh addChunks
#       ./mergeSimOutputs.sh addChunks condor
###############################################################################
# Uncomment for script-level debugging
# set -x
###############################################################################
# 0) Configuration
###############################################################################
SIM_CHUNK_DIR="/sphenix/tg/tg01/bulk/jbennett/PDC/SimOut"
OUTPUT_DIR="/sphenix/u/patsfan753/scratch/PDCrun24pp/output/simOutput"
MERGED_FILE="PositionDep_sim_ALL.root"
LISTFILE="sim_chunks_fulllist.txt"
PARTIAL_PREFIX="chunkMerge"
FILES_PER_GROUP=100          # default; may be overwritten for asManyAsCan
CONDOR_OUTDIR="/sphenix/u/patsfan753/scratch/PDCrun24pp/stdout"
CONDOR_ERRDIR="/sphenix/u/patsfan753/scratch/PDCrun24pp/error"
CONDOR_LOGDIR="/sphenix/u/patsfan753/scratch/PDCrun24pp/log"
CONDOR_EXEC="hadd_condor.sh"
TMP_SUBDIR="/sphenix/u/patsfan753/scratch/PDCrun24pp/tmp_sublists_for_merging"

usage() {
  echo "Usage:"
  echo "  $0 condor [test|firstHalf|asManyAsCan]"
  echo "  $0 addChunks [condor]"
  exit 1
}

[[ $# -eq 0 || $# -gt 2 ]] && usage
MODE="$1"; SUBMODE="${2:-}"
[[ "$MODE" != "condor" && "$MODE" != "addChunks" ]] && usage

# asManyAsCan ⇒ 1 000-file groups
[[ "$MODE" == "condor" && "$SUBMODE" == "asManyAsCan" ]] && FILES_PER_GROUP=1000

echo "============================================================================"
echo "[mergeSimOutputs.sh] START  $(date)"
echo "[INFO] Mode = $MODE $SUBMODE  |  FILES_PER_GROUP=$FILES_PER_GROUP"
echo "============================================================================"
mkdir -p "$OUTPUT_DIR" || { echo "[ERROR] Cannot create $OUTPUT_DIR"; exit 1; }

###############################################################################
# 1) CONDOR   -----------------------------------------------------------------
###############################################################################
if [[ "$MODE" == "condor" ]]; then
  testMode=false; firstHalf=false; asMany=false
  case "$SUBMODE" in
    test)          testMode=true  ;;
    firstHalf)     firstHalf=true ;;
    asManyAsCan)   asMany=true    ;;
    "")            ;;
    *) usage ;;
  esac

  # clean old outputs
  echo "[CLEAN] Removing old chunkMerge_*.root and $MERGED_FILE from $OUTPUT_DIR"
  find "$OUTPUT_DIR" -maxdepth 1 -type f \( -name "${PARTIAL_PREFIX}_*.root" \
                                         -o -name "$MERGED_FILE" \) -delete

  # build full list of chunks
  find "$SIM_CHUNK_DIR" -maxdepth 1 -type f -name "PositionDep_sim_chunk*.root" \
       | sort > "$LISTFILE"
  [[ ! -s "$LISTFILE" ]] && { echo "[ERROR] No chunk files found"; exit 1; }

  # split into sub-lists
  rm -rf "$TMP_SUBDIR"; mkdir -p "$TMP_SUBDIR"
  if $testMode; then
    head -n "$FILES_PER_GROUP" "$LISTFILE" > "$TMP_SUBDIR/sublist_00"
    totalSublists=1
  else
    split -d -l "$FILES_PER_GROUP" "$LISTFILE" "$TMP_SUBDIR/sublist_"
    totalSublists=$(ls "$TMP_SUBDIR"/sublist_* | wc -l)
  fi
  queueN=$totalSublists
  $firstHalf && queueN=$(( totalSublists / 2 ))

  # ---------------------------------------------------------------------------
  # helper script **with nounset off while sourcing env**
  # ---------------------------------------------------------------------------
  cat > "$CONDOR_EXEC" <<'EOS'
#!/usr/bin/env bash
set -eo pipefail
set +u                      # ------------ disable nounset ------------
export USER="$(id -un)"
export LOGNAME="$USER"
export HOME="/sphenix/u/$USER"
MYINSTALL="/sphenix/user/$USER/install"
source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh "$MYINSTALL"
set -u                      # ------------ re-enable nounset -----------
[[ $# -eq 2 ]] || { echo "[ERROR] Usage: $0 <list> <out>"; exit 1; }
LIST="$1"; OUT="$2"
echo "[hadd_condor] merging $(wc -l < "$LIST") files -> $OUT"
hadd -v 3 -f "$OUT" @"$LIST"
EOS
  chmod +x "$CONDOR_EXEC"

  # submission file
  SUB="partial_merge.sub"; rm -f "$SUB"
  cat > "$SUB" <<EOT
universe   = vanilla
executable = $CONDOR_EXEC
output     = $CONDOR_OUTDIR/merge.\$(Cluster).\$(Process).out
error      = $CONDOR_ERRDIR/merge.\$(Cluster).\$(Process).err
log        = $CONDOR_LOGDIR/merge.\$(Cluster).\$(Process).log
request_memory = 2000MB
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
stream_output = True
stream_error  = True
EOT

  i=0
  for SL in "$TMP_SUBDIR"/sublist_*; do
    ((i++))
    [[ $i -gt $queueN ]] && break
    echo "arguments = $SL $OUTPUT_DIR/${PARTIAL_PREFIX}_${i}.root" >> "$SUB"
    echo "queue" >> "$SUB"
  done

  condor_submit "$SUB"
  echo "[INFO] Submitted $queueN Condor jobs."
  exit 0
fi

###############################################################################
# 2) ADDCHUNKS  ---------------------------------------------------------------
###############################################################################
# build list of existing partials
partials=( "$OUTPUT_DIR"/${PARTIAL_PREFIX}_*.root )
[[ ! -e "${partials[0]}" ]] && { echo "[ERROR] No partial outputs"; exit 1; }
LIST="$OUTPUT_DIR/partialList.txt"
printf "%s\n" "${partials[@]}" > "$LIST"
FINAL="$OUTPUT_DIR/$MERGED_FILE"

if [[ "$SUBMODE" == "condor" ]]; then
  # helper for final merge (identical nounset fix)
  cat > "$CONDOR_EXEC" <<'EOS'
#!/usr/bin/env bash
set -eo pipefail
set +u
export USER="$(id -un)"; export LOGNAME="$USER"; export HOME="/sphenix/u/$USER"
MYINSTALL="/sphenix/user/$USER/install"
source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh "$MYINSTALL"
set -u
[[ $# -eq 2 ]] || { echo "[ERROR] Usage: $0 <list> <out>"; exit 1; }
hadd -v 3 -f "$2" @"$1"
EOS
  chmod +x "$CONDOR_EXEC"

  SUB="final_merge.sub"; rm -f "$SUB"
  cat > "$SUB" <<EOT
universe   = vanilla
executable = $CONDOR_EXEC
output     = $CONDOR_OUTDIR/final.\$(Cluster).\$(Process).out
error      = $CONDOR_ERRDIR/final.\$(Cluster).\$(Process).err
log        = $CONDOR_LOGDIR/final.\$(Cluster).\$(Process).log
request_memory = 2000MB
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
arguments  = $LIST $FINAL
queue
EOT
  condor_submit "$SUB"
  echo "[INFO] Submitted one Condor job to create $FINAL"
  exit 0
fi

echo "[INFO] Running local hadd -> $FINAL"
hadd -v 3 -f "$FINAL" @"$LIST"
echo "[INFO] Created $(ls -lh "$FINAL")"
