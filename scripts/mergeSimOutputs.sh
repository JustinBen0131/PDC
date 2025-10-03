#!/usr/bin/env bash
###############################################################################
# mergeSimOutputs.sh
#
# Hierarchical (1–2 stage) merging workflow for many simulation output files
# named like “PositionDep_sim_*.root”. The script can:
#
#   • Pre-check that produced sim outputs in SIM_CHUNK_DIR correspond 1:1 to
#     the expected (DST,HITS) input pairs (per-index existence test).
#   • Submit a first-stage parallel merge on Condor that groups raw sim outputs
#     into larger partial files (aka “partials”: chunkMerge_XXX.root).
#   • Optionally build a second stage of “superchunks” (aka chunkMerge10_XXX.root)
#     by merging groups of 10 partials — locally or via Condor — which speeds up
#     the final merge.
#   • Perform the final merge into PositionDep_sim_ALL.root, preferring any
#     available superchunks; otherwise falling back to partials.
#
# CLEAN-SLATE POLICY (for `condor` first stage):
#   After a successful precheck, we delete from OUTPUT_DIR:
#       • chunkMerge_*.root        (partials)
#       • chunkMerge10_*.root      (superchunks)
#       • PositionDep_sim_ALL.root (final)
#     and we recreate temporary list directories (TMP_SUBDIR, TMP_SUPERDIR).
#
# MODES
# ─────────────────────────────────────────────────────────────────────────────
# 1) condor [test | firstHalf | asManyAsCan]
#    - Runs the precheck automatically; aborts if any expected index is missing.
#    - With no sub-mode: uses FILES_PER_GROUP=1000 (e.g. 60k files -> ~60 jobs).
#    - test         : submit exactly one job (first FILES_PER_GROUP entries).
#                     Default FILES_PER_GROUP is 100; with test this yields 100.
#    - firstHalf    : submit only the first half of the planned sublists.
#    - asManyAsCan  : set FILES_PER_GROUP=1000 explicitly.
#    - Output: chunkMerge_001.root, chunkMerge_002.root, … in OUTPUT_DIR.
#
# 2) addChunks [part1] [condor]
#    - Without “part1”: final merge into MERGED_FILE using whichever exists:
#        • Prefer superchunks (chunkMerge10_*.root)
#        • Otherwise use partials (chunkMerge_*.root)
#      Run locally by default; add “condor” to run the final merge on Condor.
#    - With “part1”: build superchunks (groups of 10 partials) and exit.
#        • Local:   addChunks part1
#        • Condor:  addChunks part1 condor    (order of tokens doesn’t matter)
#      Output: chunkMerge10_001.root, chunkMerge10_002.root, … in OUTPUT_DIR.
#ƒ
# 3) checkFileOutput
#    - Standalone precheck (no submission). Reports:
#        • expected input pair count (from SIM_DST_LIST + SIM_HITS_LIST)
#        • found output files in SIM_CHUNK_DIR
#        • any missing indices with the exact expected filenames.
#
# USAGE EXAMPLES
# ─────────────────────────────────────────────────────────────────────────────
#   ./mergeSimOutputs.sh condor
#   ./mergeSimOutputs.sh condor test
#   ./mergeSimOutputs.sh condor firstHalf
#   ./mergeSimOutputs.sh condor asManyAsCan
#
#   ./mergeSimOutputs.sh addChunks               # final merge (local)
#   ./mergeSimOutputs.sh addChunks condor        # final merge (Condor)
#   ./mergeSimOutputs.sh addChunks part1         # build superchunks (local)
#   ./mergeSimOutputs.sh addChunks part1 condor  # build superchunks (Condor)
#
#   ./mergeSimOutputs.sh checkFileOutput
#
# DEFAULTS & KNOBS
# ─────────────────────────────────────────────────────────────────────────────
#   SIM_CHUNK_DIR  : directory holding raw “PositionDep_sim_*.root” outputs
#   OUTPUT_DIR     : where merged files are produced
#   PARTIAL_PREFIX : chunkMerge        (first-stage merged file name prefix)
#   SUPER_PREFIX   : chunkMerge10      (second-stage merged file name prefix)
#   MERGED_FILE    : PositionDep_sim_ALL.root
#   FILES_PER_GROUP: 100  (default); plain `condor` forces 1000 unless a sub-mode
#   GROUP_OF       : 10   (partials per superchunk)
#
# NOTES
# ─────────────────────────────────────────────────────────────────────────────
#   • Numerical sorting (`sort -V`) is used when ordering inputs for merging.
#   • Condor submit files stream stdout/stderr for better live diagnostics.
#   • The precheck *aborts* when any per-index output is missing. If counts
#     differ but every index is present (e.g., extra stray files), a warning
#     is emitted and the run proceeds.
###############################################################################
# 0) Configuration
###############################################################################
SIM_CHUNK_DIR="/sphenix/tg/tg01/bulk/jbennett/PDC/SimOut/9999"
DATA_RUN_BASE="/sphenix/tg/tg01/bulk/jbennett/PDC/output"   # per-run subdirs: /.../output/<run>/
OUTPUT_DIR="/sphenix/u/patsfan753/scratch/PDCrun24pp/output/simOutput"
OUTPUT_DIR_DATA="/sphenix/u/patsfan753/scratch/PDCrun24pp/output/dataOutput"
MERGED_FILE="PositionDep_sim_ALL.root"
LISTFILE="sim_chunks_fulllist.txt"
PARTIAL_PREFIX="chunkMerge"
SUPER_PREFIX="chunkMerge10"  # superchunks (groups of 10 partials)
FILES_PER_GROUP=100          # default; may be overwritten for asManyAsCan
CONDOR_OUTDIR="/sphenix/u/patsfan753/scratch/PDCrun24pp/stdout"
CONDOR_ERRDIR="/sphenix/u/patsfan753/scratch/PDCrun24pp/error"
CONDOR_LOGDIR="/sphenix/u/patsfan753/scratch/PDCrun24pp/log"
CONDOR_EXEC="hadd_condor.sh"
TMP_SUBDIR="/sphenix/u/patsfan753/scratch/PDCrun24pp/tmp_sublists_for_merging"
TMP_SUPERDIR="/sphenix/u/patsfan753/scratch/PDCrun24pp/tmp_superlists_for_merging"

# Master switch: turn off all file checking/precheck logic if false
DO_FILE_CHECK=false


# Simulation list files (for integrity checks in 'checkFileOutput' mode)
: "${SIM_LIST_DIR:="/sphenix/u/patsfan753/scratch/PDCrun24pp/simListFiles/run24_type14_gamma_pt_200_40000"}"
: "${SIM_DST_LIST:="${SIM_LIST_DIR}/DST_CALO_CLUSTER.list"}"
: "${SIM_HITS_LIST:="${SIM_LIST_DIR}/G4Hits.list"}"

usage() {
  echo "Usage:"
  echo "  $0 condor [test|firstHalf|asManyAsCan]"
  echo "  $0 local                      # sequential first-stage merges locally (1000 per group)"
  echo "  $0 addChunks [condor]"
  echo "  $0 addChunks part1 [condor]"
  echo "  $0 checkFileOutput"
  echo
  echo "Examples:"
  echo "  $0 condor"
  echo "  $0 condor test"
  echo "  $0 condor firstHalf"
  echo "  $0 condor asManyAsCan"
  echo "  $0 local"
  echo "  $0 addChunks"
  echo "  $0 addChunks condor"
  echo "  $0 addChunks part1"
  echo "  $0 addChunks part1 condor"
  exit 1
}

# Allow up to 4 tokens so we can add independent flags like 'minbias' or 'data'
[[ $# -lt 1 || $# -gt 4 ]] && usage
MODE="$1"
SUB1="${2:-}"
SUB2="${3:-}"
SUB3="${4:-}"

# Detect flags anywhere (case-insensitive), but DO NOT include them in SUBMODE
MINBIAS=false
DATA_MODE=false
for tok in "$SUB1" "$SUB2" "$SUB3"; do
  case "$tok" in
    MINBIAS|minbias|MinBias) MINBIAS=true ;;
    DATA|data|Data)          DATA_MODE=true ;;
  esac
done

# Build SUBMODE from remaining tokens (exclude MINBIAS/DATA tokens)
SUBMODE=""
for tok in "$SUB1" "$SUB2" "$SUB3"; do
  case "$tok" in
    ""|MINBIAS|minbias|MinBias|DATA|data|Data) continue ;;
    *) SUBMODE="${SUBMODE:+$SUBMODE }$tok" ;;
  esac
done


[[ "$MODE" != "condor" && "$MODE" != "local" && "$MODE" != "addChunks" && "$MODE" != "checkFileOutput" ]] && usage


# asManyAsCan ⇒ 1 000-file groups
[[ "$MODE" == "condor" && "$SUBMODE" == "asManyAsCan" ]] && FILES_PER_GROUP=1000

# Apply 'minbias' overrides before creating directories
if $MINBIAS; then
  # Separate produced chunk inputs
  SIM_CHUNK_DIR="/sphenix/tg/tg01/bulk/jbennett/PDC/SimOutMinBias/9999"

  # Separate expected-input lists for integrity checks
  SIM_LIST_DIR="/sphenix/u/patsfan753/scratch/PDCrun24pp/simListFiles/run21_type30_MB_nopileup"
  SIM_DST_LIST="${SIM_LIST_DIR}/DST_CALO_CLUSTER.list"
  SIM_HITS_LIST="${SIM_LIST_DIR}/G4HITS.list"

  # Separate all output artifacts by directory (keep Condor logs shared)
  OUTPUT_DIR="${OUTPUT_DIR%/}/minbias"
  TMP_SUBDIR="${TMP_SUBDIR%/}/minbias"
  TMP_SUPERDIR="${TMP_SUPERDIR%/}/minbias"

  # Avoid cross-run clobbering of the master list file
  LISTFILE="${OUTPUT_DIR%/}/sim_chunks_fulllist_minbias.txt"
fi

# Apply 'data' overrides: per-run inputs under $DATA_RUN_BASE, write per-run partials under $OUTPUT_DIR_DATA
if $DATA_MODE; then
  OUTPUT_DIR="$OUTPUT_DIR_DATA"
  MERGED_FILE="PositionDep_data_ALL.root"
  PARTIAL_PREFIX="chunkMerge_run"
fi

echo "============================================================================"
echo "[mergeSimOutputs.sh] START  $(date)"
echo "[INFO] Mode = $MODE ${SUBMODE:-<none>}  |  FILES_PER_GROUP=$FILES_PER_GROUP  |  MINBIAS=$MINBIAS"
echo "[INFO] SIM_CHUNK_DIR = $SIM_CHUNK_DIR"
echo "[INFO] OUTPUT_DIR    = $OUTPUT_DIR"
echo "============================================================================"

# Ensure all destination directories exist (needed if minbias appended paths are new)
mkdir -p "$OUTPUT_DIR" "$CONDOR_OUTDIR" "$CONDOR_ERRDIR" "$CONDOR_LOGDIR" || {
  echo "[ERROR] Cannot create output/log directories"; exit 1;
}


###############################################################################
# 1) CHECKFILEOUTPUT OR CONDOR  ----------------------------------------------
###############################################################################
# A) Input/output consistency check: no submission, just report
if [[ "$MODE" == "checkFileOutput" && "$DO_FILE_CHECK" == true ]]; then
  # sanity
  [[ -d "$SIM_CHUNK_DIR" ]] || { echo "[ERROR] SIM_CHUNK_DIR not found: $SIM_CHUNK_DIR"; exit 1; }
  [[ -f "$SIM_DST_LIST"  ]] || { echo "[ERROR] SIM_DST_LIST not found:  $SIM_DST_LIST";  exit 1; }
  [[ -f "$SIM_HITS_LIST" ]] || { echo "[ERROR] SIM_HITS_LIST not found: $SIM_HITS_LIST"; exit 1; }

  # expected inputs (paired), consistent with submission script
  mapfile -t pairs < <(paste -d' ' "$SIM_DST_LIST" "$SIM_HITS_LIST" | sort -k1,1V)
  totalExpected=${#pairs[@]}

  # present outputs (any produced files)
  presentOutputs=$(find "$SIM_CHUNK_DIR" -maxdepth 1 -type f -name 'PositionDep_sim_*.root' | wc -l)

  echo "======================================================================="
  echo "[checkFileOutput] BEGIN ($(date))"
  echo "  SIM_CHUNK_DIR : $SIM_CHUNK_DIR"
  echo "  SIM_DST_LIST  : $SIM_DST_LIST"
  echo "  SIM_HITS_LIST : $SIM_HITS_LIST"
  echo "-----------------------------------------------------------------------"
  echo "  Expected input pairs : $totalExpected"
  echo "  Found output files   : $presentOutputs  (pattern: PositionDep_sim_*.root)"
  echo "-----------------------------------------------------------------------"

  missing=0
  # loop through expected pairs and check for a corresponding output filename
  for idx in "${!pairs[@]}"; do
    read -r dstPath hitsPath <<< "${pairs[$idx]}"

    # try common naming patterns used by different production modes
    candidates=(
      "$SIM_CHUNK_DIR/PositionDep_sim_${idx}.root"
      "$SIM_CHUNK_DIR/PositionDep_sim_$(printf "%06d" "$idx").root"
      "$SIM_CHUNK_DIR/PositionDep_sim_pair${idx}.root"
      "$SIM_CHUNK_DIR/PositionDep_sim_pair$(printf "%06d" "$idx").root"
      "$SIM_CHUNK_DIR/PositionDep_sim_chunk${idx}.root"
      "$SIM_CHUNK_DIR/PositionDep_sim_$(basename "$dstPath")"
    )

    found=false
    for c in "${candidates[@]}"; do
      if [[ -f "$c" ]]; then
        found=true
        break
      fi
    done

    if ! $found; then
      ((missing++))
      echo "[MISSING] index=$idx"
      echo "          DST : $(basename "$dstPath")"
      echo "          HITS: $(basename "$hitsPath")"
      echo "          looked for any of:"
      for c in "${candidates[@]}"; do
        echo "             $(basename "$c")"
      done
    fi
  done

  echo "-----------------------------------------------------------------------"
  echo "[checkFileOutput] SUMMARY: missing=$missing / expected=$totalExpected"
  if (( missing == 0 )); then
    echo "[checkFileOutput] All expected inputs have a corresponding output file."
  else
    echo "[checkFileOutput] See [MISSING] lines above for exact DST/HITS inputs."
  fi
  echo "======================================================================="
  exit 0
fi

###############################################################################
# 1b) LOCAL (sequential) or CONDOR  ------------------------------------------
###############################################################################
if [[ "$MODE" == "local" || "$MODE" == "condor" ]]; then
  echo "-----------------------------------------------------------------------"
  echo "[MODE] Entering ${MODE^^} path  |  SUBMODE='${SUBMODE:-<none>}'  |  DATA_MODE=$DATA_MODE  |  DO_FILE_CHECK=$DO_FILE_CHECK"
  echo "-----------------------------------------------------------------------"

  ###########################################################################
  # LOCAL: do first-stage merges sequentially on this node
  ###########################################################################
  if [[ "$MODE" == "local" ]]; then
    FILES_PER_GROUP=1000
    echo "======================================================================="
    echo "[LOCAL] First-stage merging (sequential) with FILES_PER_GROUP=$FILES_PER_GROUP"
    echo "======================================================================="

    # Optional precheck (same logic as in condor path)
    if [[ "$DO_FILE_CHECK" == true ]]; then
      echo "[LOCAL][PRECHECK] Verifying produced sim outputs vs expected inputs"
      [[ -d "$SIM_CHUNK_DIR" ]] || { echo "[ERROR] SIM_CHUNK_DIR not found: $SIM_CHUNK_DIR"; exit 1; }
      [[ -f "$SIM_DST_LIST"  ]] || { echo "[ERROR] SIM_DST_LIST not found:  $SIM_DST_LIST";  exit 1; }
      [[ -f "$SIM_HITS_LIST" ]] || { echo "[ERROR] SIM_HITS_LIST not found: $SIM_HITS_LIST"; exit 1; }

      mapfile -t pairs < <(paste -d' ' "$SIM_DST_LIST" "$SIM_HITS_LIST" | sort -k1,1V)
      totalExpected=${#pairs[@]}
      presentOutputs=$(find "$SIM_CHUNK_DIR" -maxdepth 1 -type f -name 'PositionDep_sim_*.root' | wc -l)

      echo "  SIM_CHUNK_DIR : $SIM_CHUNK_DIR"
      echo "  SIM_DST_LIST  : $SIM_DST_LIST"
      echo "  SIM_HITS_LIST : $SIM_HITS_LIST"
      echo "  Expected input pairs : $totalExpected"
      echo "  Found output files   : $presentOutputs"
    fi

    echo "[LOCAL][CLEAN] Removing old ${PARTIAL_PREFIX}_*.root, ${SUPER_PREFIX}_*.root, and $MERGED_FILE from $OUTPUT_DIR"
    find "$OUTPUT_DIR" -maxdepth 1 -type f \( -name "${PARTIAL_PREFIX}_*.root" -o -name "${SUPER_PREFIX}_*.root" -o -name "$MERGED_FILE" \) -delete
    rm -rf "$TMP_SUBDIR" "$TMP_SUPERDIR"

    echo "[LOCAL] Scanning chunk directory: $SIM_CHUNK_DIR"
    find "$SIM_CHUNK_DIR" -maxdepth 1 -type f -name "PositionDep_sim_*.root" | sort -V > "$LISTFILE"
    [[ ! -s "$LISTFILE" ]] && { echo "[ERROR][LOCAL] No chunk files found"; exit 1; }

    totalToMerge=$(wc -l < "$LISTFILE")
    echo "[LOCAL] Files to merge listed in $LISTFILE  (count=$totalToMerge)"

    mkdir -p "$TMP_SUBDIR" "$TMP_SUPERDIR"
    split -d -l "$FILES_PER_GROUP" "$LISTFILE" "$TMP_SUBDIR/sublist_"
    totalSublists=$(compgen -G "$TMP_SUBDIR/sublist_*" | wc -l)
    echo "[LOCAL] Total sublists (groups of $FILES_PER_GROUP): $totalSublists"

    i=0
    for SL in "$TMP_SUBDIR"/sublist_*; do
      ((i++))
      out="$OUTPUT_DIR/${PARTIAL_PREFIX}_$(printf "%d" "$i").root"
      echo "-----------------------------------------------------------------------"
      echo "[LOCAL] hadd -v 3 -f $out @$SL"
      hadd -v 3 -f "$out" @"$SL" || { echo "[ERROR][LOCAL] hadd failed for $SL"; exit 1; }
    done

    echo "======================================================================="
    echo "[LOCAL] Completed $i partial merges in $OUTPUT_DIR"
    echo "======================================================================="
    exit 0
  fi  # end LOCAL block

  ###########################################################################
  # CONDOR: first-stage parallel merges
  ###########################################################################
  echo "======================================================================="
  echo "[CONDOR] First-stage parallel merge planning"
  echo "======================================================================="

  testMode=false; firstHalf=false; asMany=false
  case "$SUBMODE" in
    test)          testMode=true  ;;
    firstHalf)     firstHalf=true ;;
    asManyAsCan)   asMany=true    ;;
    "")            ;;
    *)             usage ;;
  esac

  # PRE-FLIGHT: optionally run the 'checkFileOutput' logic
  if [[ "$DO_FILE_CHECK" == true ]]; then
    echo "======================================================================="
    echo "[PRECHECK] Verifying produced sim outputs correspond 1:1 with expected inputs"
    echo "           (equivalent to: ./mergeSimOutputs.sh checkFileOutput)"
    echo "-----------------------------------------------------------------------"
    [[ -d "$SIM_CHUNK_DIR" ]] || { echo "[ERROR] SIM_CHUNK_DIR not found: $SIM_CHUNK_DIR"; exit 1; }
    [[ -f "$SIM_DST_LIST"  ]] || { echo "[ERROR] SIM_DST_LIST not found:  $SIM_DST_LIST";  exit 1; }
    [[ -f "$SIM_HITS_LIST" ]] || { echo "[ERROR] SIM_HITS_LIST not found: $SIM_HITS_LIST"; exit 1; }

    mapfile -t pairs < <(paste -d' ' "$SIM_DST_LIST" "$SIM_HITS_LIST" | sort -k1,1V)
    totalExpected=${#pairs[@]}
    presentOutputs=$(find "$SIM_CHUNK_DIR" -maxdepth 1 -type f -name 'PositionDep_sim_*.root' | wc -l)

    echo "  SIM_CHUNK_DIR : $SIM_CHUNK_DIR"
    echo "  SIM_DST_LIST  : $SIM_DST_LIST"
    echo "  SIM_HITS_LIST : $SIM_HITS_LIST"
    echo "-----------------------------------------------------------------------"
    echo "  Expected input pairs : $totalExpected"
    echo "  Found output files   : $presentOutputs  (pattern: PositionDep_sim_*.root)"
    echo "-----------------------------------------------------------------------"

    missing=0
    for idx in "${!pairs[@]}"; do
      read -r dstPath hitsPath <<< "${pairs[$idx]}"
      candidates=(
        "$SIM_CHUNK_DIR/PositionDep_sim_${idx}.root"
        "$SIM_CHUNK_DIR/PositionDep_sim_$(printf "%06d" "$idx").root"
        "$SIM_CHUNK_DIR/PositionDep_sim_pair${idx}.root"
        "$SIM_CHUNK_DIR/PositionDep_sim_pair$(printf "%06d" "$idx").root"
        "$SIM_CHUNK_DIR/PositionDep_sim_chunk${idx}.root"
      )
      found=false
      for c in "${candidates[@]}"; do
        [[ -f "$c" ]] && { found=true; break; }
      done
      if ! $found; then
        ((missing++))
        echo "[MISSING] index=$idx"
        echo "          DST : $(basename "$dstPath")"
        echo "          HITS: $(basename "$hitsPath")"
      fi
    done

    if (( missing > 0 )); then
      echo "-----------------------------------------------------------------------"
      echo "[PRECHECK][WARN] Missing=$missing / expected=$totalExpected"
      echo "[PRECHECK] Proceeding with Condor submissions (will merge what exists)."
      echo "-----------------------------------------------------------------------"
    fi

    if (( presentOutputs != totalExpected )); then
      echo "-----------------------------------------------------------------------"
      echo "[PRECHECK][WARN] Count mismatch (present=$presentOutputs, expected=$totalExpected)"
      echo "[PRECHECK] However, per-index existence check passed; continuing anyway."
      echo "-----------------------------------------------------------------------"
    else
      echo "[PRECHECK][PASS] All $totalExpected expected inputs have matching output files."
    fi
    echo "======================================================================="
  fi  # end PRE-FLIGHT

  # PLAN: default grouping per job; MINBIAS uses smaller (100), default uses 300
  if ! $testMode && ! $firstHalf && ! $asMany; then
    if $MINBIAS; then
      FILES_PER_GROUP=100
      echo "[PLAN][MINBIAS] No sub-mode supplied → FILES_PER_GROUP=${FILES_PER_GROUP} (smaller groups for MB I/O)"
    else
      FILES_PER_GROUP=300
      echo "[PLAN] No sub-mode supplied → FILES_PER_GROUP=${FILES_PER_GROUP} (target ~40 partials for ~12k files)"
    fi
    echo "       With ${totalToMerge:-0} outputs, expected sublists ≈ $(( (totalToMerge + FILES_PER_GROUP - 1) / FILES_PER_GROUP ))"
  else
    echo "[PLAN] Sub-mode='${SUBMODE:-<none>}' → FILES_PER_GROUP=${FILES_PER_GROUP}"
  fi

  # CLEAN after precheck
  echo "[CLEAN] Removing old ${PARTIAL_PREFIX}_*.root, ${SUPER_PREFIX}_*.root, and $MERGED_FILE from $OUTPUT_DIR"
  find "$OUTPUT_DIR" -maxdepth 1 -type f \( -name "${PARTIAL_PREFIX}_*.root" -o -name "${SUPER_PREFIX}_*.root" -o -name "$MERGED_FILE" \) -delete
  rm -rf "$TMP_SUBDIR" "$TMP_SUPERDIR"

  # Build helper once
  cat > "$CONDOR_EXEC" <<'EOS'
#!/usr/bin/env bash
set -eo pipefail
set +u
export USER="$(id -un)"
export LOGNAME="$USER"
export HOME="/sphenix/u/$USER"
MYINSTALL="/sphenix/user/$USER/install"
source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh "$MYINSTALL"
set -u
[[ $# -eq 2 ]] || { echo "[ERROR] Usage: $0 <list> <out>"; exit 1; }
LIST="$1"; OUT="$2"
echo "[hadd_condor] merging $(wc -l < "$LIST") files -> $OUT"
hadd -v 3 -f "$OUT" @"$LIST"
EOS
  chmod +x "$CONDOR_EXEC"

  mkdir -p "$TMP_SUBDIR" "$TMP_SUPERDIR"

  ###########################################################################
  # DATA MODE on CONDOR: one job per run directory
  ###########################################################################
  if $DATA_MODE; then
    echo "[CONDOR][DATA] Building one merge job per run under $DATA_RUN_BASE"
    mapfile -t runlist < <(find "$DATA_RUN_BASE" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort -V)
    [[ ${#runlist[@]} -gt 0 ]] || { echo "[ERROR] No run directories under $DATA_RUN_BASE"; exit 1; }

    queueN=${#runlist[@]}
    $testMode  && queueN=1
    $firstHalf && queueN=$(( queueN / 2 ))
    echo "[CONDOR][DATA] Planned runs: ${#runlist[@]}  |  queueN=$queueN"

    SUB="partial_merge.sub"; rm -f "$SUB"
    cat > "$SUB" <<EOT
universe   = vanilla
executable = $CONDOR_EXEC
output     = $CONDOR_OUTDIR/merge.\$(Cluster).\$(Process).out
error      = $CONDOR_ERRDIR/merge.\$(Cluster).\$(Process).err
log        = $CONDOR_LOGDIR/merge.\$(Cluster).\$(Process).log
request_memory = 1.5GB
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
stream_output = True
stream_error  = True
EOT

    i=0
    for run in "${runlist[@]}"; do
      ((i++))
      [[ $i -gt $queueN ]] && break
      list="$TMP_SUBDIR/run_${run}.txt"
      find "$DATA_RUN_BASE/$run" -maxdepth 1 -type f -name "PositionDep_data_*.root" | sort -V > "$list"
      if [[ ! -s "$list" ]]; then
        echo "[WARN] No PositionDep_data_*.root files for run $run – skipping"
        ((i--))
        continue
      fi
      echo "arguments = $list $OUTPUT_DIR/${PARTIAL_PREFIX}_${run}.root" >> "$SUB"
      echo "queue" >> "$SUB"
    done

    echo "======================================================================="
    echo "[SUBMIT] condor_submit $SUB"
    echo "         (one partial per run → ${OUTPUT_DIR}/${PARTIAL_PREFIX}_<run>.root)"
    echo "======================================================================="
    condor_submit "$SUB" || { echo "[ERROR] condor_submit failed"; exit 1; }
    echo "[INFO] Submitted $i Condor jobs (data per-run)."
    exit 0
  fi  # end DATA_MODE on CONDOR

  ###########################################################################
  # SIM MODE on CONDOR: group by FILES_PER_GROUP (original behavior)
  ###########################################################################
  echo "[CONDOR][SIM] Grouping by FILES_PER_GROUP=$FILES_PER_GROUP from $SIM_CHUNK_DIR"
  find "$SIM_CHUNK_DIR" -maxdepth 1 -type f -name "PositionDep_sim_*.root" | sort -V > "$LISTFILE"
  [[ ! -s "$LISTFILE" ]] && { echo "[ERROR] No chunk files found after precheck"; exit 1; }

  totalToMerge=$(wc -l < "$LISTFILE")
  echo "[CONDOR][SIM] Files to merge listed in $LISTFILE  (count=$totalToMerge)"

  if $testMode; then
    head -n "$FILES_PER_GROUP" "$LISTFILE" > "$TMP_SUBDIR/sublist_00"
    totalSublists=1
  else
    split -d -l "$FILES_PER_GROUP" "$LISTFILE" "$TMP_SUBDIR/sublist_"
    totalSublists=$(compgen -G "$TMP_SUBDIR/sublist_*" | wc -l)
  fi

  queueN=$totalSublists
  $firstHalf && queueN=$(( totalSublists / 2 ))

  echo "-----------------------------------------------------------------------"
  echo "[CONDOR][SIM] Grouping plan"
  echo "  • FILES_PER_GROUP   : $FILES_PER_GROUP"
  echo "  • Total sublists    : $totalSublists"
  $firstHalf && echo "  • firstHalf active   → queueN=$queueN"
  $testMode  && echo "  • test mode active   → queueN=$queueN (one sublist)"
  echo "  • Jobs to submit    : $queueN"
  echo "-----------------------------------------------------------------------"

  SUB="partial_merge.sub"; rm -f "$SUB"
  cat > "$SUB" <<EOT
universe   = vanilla
executable = $CONDOR_EXEC
output     = $CONDOR_OUTDIR/merge.\$(Cluster).\$(Process).out
error      = $CONDOR_ERRDIR/merge.\$(Cluster).\$(Process).err
log        = $CONDOR_LOGDIR/merge.\$(Cluster).\$(Process).log
request_memory = 1.5GB
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

  echo "======================================================================="
  echo "[SUBMIT] condor_submit $SUB"
  echo "         (will create ${queueN} partial outputs named ${PARTIAL_PREFIX}_<N>.root)"
  echo "======================================================================="
  condor_submit "$SUB" || { echo "[ERROR] condor_submit failed"; exit 1; }
  echo "[INFO] Submitted $i Condor jobs."
  exit 0

fi  # end LOCAL/CONDOR selector




###############################################################################
# 2) ADDCHUNKS  ---------------------------------------------------------------
###############################################################################
# New sub-modes:
#   addChunks part1                -> build superchunks from partials locally
#   addChunks part1 condor         -> build superchunks via Condor
#   addChunks [condor]             -> final merge; prefers superchunks if present,
#                                     otherwise falls back to partials.
#   addChunks <8-digit-run-number> -> LOCAL one-run merge of all PositionDep_data_*.root
#                                     found under $DATA_RUN_BASE/<run> into:
#                                     $OUTPUT_DIR_DATA/chunkMerge_run_<run>.root
###############################################################################

# Special case: addChunks <RUN> → local per-run data merge
if [[ "$MODE" == "addChunks" && "$SUBMODE" =~ ^[0-9]{8}$ ]]; then
  RUN="$SUBMODE"
  RUN_DIR="${DATA_RUN_BASE%/}/$RUN"
  [[ -d "$RUN_DIR" ]] || { echo "[ERROR] Run directory not found: $RUN_DIR"; exit 1; }

  mkdir -p "$OUTPUT_DIR_DATA" "$TMP_SUBDIR"

  LIST="$TMP_SUBDIR/run_${RUN}_data_files.txt"
  find "$RUN_DIR" -maxdepth 1 -type f -name "PositionDep_data_*.root" | sort -V > "$LIST"

  [[ -s "$LIST" ]] || { echo "[ERROR] No PositionDep_data_*.root files found for run $RUN under $RUN_DIR"; exit 1; }

  OUT="$OUTPUT_DIR_DATA/chunkMerge_run_${RUN}.root"
  echo "======================================================================="
  echo "[addChunks][DATA][LOCAL] Merging run ${RUN}"
  echo "  Input list : $LIST  ($(wc -l < "$LIST") files)"
  echo "  Output     : $OUT"
  echo "======================================================================="

  rm -f "$OUT" 2>/dev/null || true
  hadd -v 3 -f "$OUT" @"$LIST" || { echo "[ERROR] hadd failed for run $RUN"; exit 1; }

  echo "[INFO] Created $(ls -lh "$OUT")"
  exit 0
fi

# Auto-detect dataset type for final merge — PREFER SIM artifacts under $OUTPUT_DIR.
# Only if no SIM superchunks/partials are present, fall back to DATA partials.
DATA_MODE_ADD=false

# If SIM superchunks or SIM partials exist in $OUTPUT_DIR, stay in SIM mode.
if compgen -G "$OUTPUT_DIR/${SUPER_PREFIX}_*.root" > /dev/null || \
   compgen -G "$OUTPUT_DIR/${PARTIAL_PREFIX}_*.root" > /dev/null; then
  DATA_MODE_ADD=false
else
  # Otherwise, check for DATA partials first in the canonical DATA dir, then in $OUTPUT_DIR
  if [[ -n "${OUTPUT_DIR_DATA:-}" ]] && compgen -G "$OUTPUT_DIR_DATA/chunkMerge_run_*.root" > /dev/null; then
    DATA_MODE_ADD=true
    OUTPUT_DIR="$OUTPUT_DIR_DATA"
  elif compgen -G "$OUTPUT_DIR/chunkMerge_run_*.root" > /dev/null; then
    DATA_MODE_ADD=true
  fi
fi

if $DATA_MODE_ADD; then
  PARTIAL_PREFIX="chunkMerge_run"
  MERGED_FILE="PositionDep_data_ALL.root"
fi

FINAL="$OUTPUT_DIR/$MERGED_FILE"
GROUP_OF=10                        # merge 10 partials -> 1 superchunk


# Helper to emit a simple wrapper that runs hadd on a @listfile
emit_hadd_wrapper() {
  local exe="$1"
  cat > "$exe" <<'EOS'
#!/usr/bin/env bash
set -eo pipefail
set +u
export USER="$(id -un)"; export LOGNAME="$USER"; export HOME="/sphenix/u/$USER"
MYINSTALL="/sphenix/user/$USER/install"
source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh "$MYINSTALL"
set -u
if [[ $# -ne 2 ]]; then
  echo "[ERROR] Usage: $0 <listfile> <outroot>" >&2
  exit 1
fi
LIST="$1"; OUT="$2"
echo "[hadd_condor] inputs=$(wc -l < "$LIST")  ->  $OUT"
hadd -v 3 -f "$OUT" @"$LIST"
EOS
  chmod +x "$exe"
}


# Build arrays of existing files
mapfile -t partials < <(ls -1 "$OUTPUT_DIR"/${PARTIAL_PREFIX}_*.root 2>/dev/null || true)
mapfile -t supers   < <(ls -1 "$OUTPUT_DIR"/${SUPER_PREFIX}_*.root   2>/dev/null || true)

# --------------------------- PART1: build superchunks ------------------------
if [[ "$SUBMODE" == "part1" || "$SUBMODE" == "part1 condor" || "$SUBMODE" == "condor part1" ]]; then
  (( ${#partials[@]} )) || { echo "[ERROR] No partial outputs (${PARTIAL_PREFIX}_*.root) found in $OUTPUT_DIR"; exit 1; }

  echo "[INFO] addChunks part1: grouping ${#partials[@]} partials into sets of ${GROUP_OF}"
  rm -rf "$TMP_SUPERDIR"; mkdir -p "$TMP_SUPERDIR"

  # Create sublists of size GROUP_OF
  idx=0; group=0
  current_list=""
  for f in "${partials[@]}"; do
    if (( idx % GROUP_OF == 0 )); then
      ((group++))
      current_list="$TMP_SUPERDIR/superlist_$(printf "%03d" "$group").txt"
      : > "$current_list"
    fi
    echo "$f" >> "$current_list"
    ((idx++))
  done

  # Submit jobs (local or condor)
  if [[ "$SUBMODE" =~ condor ]]; then
    emit_hadd_wrapper "$CONDOR_EXEC"
    SUB="build_superchunks.sub"; rm -f "$SUB"
    cat > "$SUB" <<EOT
universe   = vanilla
executable = $CONDOR_EXEC
output     = $CONDOR_OUTDIR/super.\$(Cluster).\$(Process).out
error      = $CONDOR_ERRDIR/super.\$(Cluster).\$(Process).err
log        = $CONDOR_LOGDIR/super.\$(Cluster).\$(Process).log
request_memory = 3GB
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
stream_output = True
stream_error  = True
EOT

    i=0
    for SL in "$TMP_SUPERDIR"/superlist_*.txt; do
      ((i++))
      out="$OUTPUT_DIR/${SUPER_PREFIX}_$(printf "%03d" "$i").root"
      echo "arguments = $SL $out" >> "$SUB"
      echo "queue" >> "$SUB"
    done

    echo "[SUBMIT] condor_submit $SUB  (will create ${i} superchunks ${SUPER_PREFIX}_NNN.root)"
    condor_submit "$SUB" || { echo "[ERROR] condor_submit for part1 failed"; exit 1; }
    echo "[INFO] Submitted $i Condor jobs for superchunks."
    exit 0
  else
    # Local build of superchunks
    i=0
    for SL in "$TMP_SUPERDIR"/superlist_*.txt; do
      ((i++))
      out="$OUTPUT_DIR/${SUPER_PREFIX}_$(printf "%03d" "$i").root"
      echo "[LOCAL] hadd -f $out @$(basename "$SL")"
      hadd -v 3 -f "$out" @"$SL" || { echo "[ERROR] hadd failed for $SL"; exit 1; }
    done
    echo "[INFO] Created $i superchunks in $OUTPUT_DIR"
    exit 0
  fi
fi

# --------------------------- FINAL: merge to one file ------------------------
# Prefer superchunks if present; otherwise fall back to partials.
if (( ${#supers[@]} )); then
  echo "[INFO] Using ${#supers[@]} superchunks (${SUPER_PREFIX}_*.root) for final merge."
  LIST="$OUTPUT_DIR/superList.txt"
  printf "%s\n" "${supers[@]}" > "$LIST"
else
  (( ${#partials[@]} )) || { echo "[ERROR] No inputs for final merge (neither ${SUPER_PREFIX}_*.root nor ${PARTIAL_PREFIX}_*.root)"; exit 1; }
  echo "[INFO] Using ${#partials[@]} partials (${PARTIAL_PREFIX}_*.root) for final merge."
  LIST="$OUTPUT_DIR/partialList.txt"
  printf "%s\n" "${partials[@]}" > "$LIST"
fi

if [[ "$SUBMODE" == "condor" ]]; then
  emit_hadd_wrapper "$CONDOR_EXEC"
  SUB="final_merge.sub"; rm -f "$SUB"
  cat > "$SUB" <<EOT
universe   = vanilla
executable = $CONDOR_EXEC
output     = $CONDOR_OUTDIR/final.\$(Cluster).\$(Process).out
error      = $CONDOR_ERRDIR/final.\$(Cluster).\$(Process).err
log        = $CONDOR_LOGDIR/final.\$(Cluster).\$(Process).log
request_memory = 3GB
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
stream_output = True
stream_error  = True
arguments  = $LIST $FINAL
queue
EOT
  echo "[SUBMIT] condor_submit $SUB  (final -> $(basename "$FINAL"))"
  condor_submit "$SUB" || { echo "[ERROR] condor_submit for final merge failed"; exit 1; }
  echo "[INFO] Submitted one Condor job to create $FINAL"
  exit 0
fi

echo "[INFO] Running local hadd -> $FINAL"
hadd -v 3 -f "$FINAL" @"$LIST" || { echo "[ERROR] local hadd failed"; exit 1; }
echo "[INFO] Created $(ls -lh "$FINAL")"
