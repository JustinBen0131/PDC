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
OUTPUT_DIR="/sphenix/u/patsfan753/scratch/PDCrun24pp/output/simOutput"
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

[[ $# -lt 1 || $# -gt 3 ]] && usage
MODE="$1"
SUB1="${2:-}"
SUB2="${3:-}"
# Preserve existing two-arg behavior; allow combined submodes for addChunks
SUBMODE="$SUB1"
[[ -n "$SUB2" ]] && SUBMODE="$SUB1 $SUB2"
[[ "$MODE" != "condor" && "$MODE" != "local" && "$MODE" != "addChunks" && "$MODE" != "checkFileOutput" ]] && usage


# asManyAsCan ⇒ 1 000-file groups
[[ "$MODE" == "condor" && "$SUBMODE" == "asManyAsCan" ]] && FILES_PER_GROUP=1000

echo "============================================================================"
echo "[mergeSimOutputs.sh] START  $(date)"
echo "[INFO] Mode = $MODE $SUBMODE  |  FILES_PER_GROUP=$FILES_PER_GROUP"
echo "============================================================================"
mkdir -p "$OUTPUT_DIR" || { echo "[ERROR] Cannot create $OUTPUT_DIR"; exit 1; }

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

  if [[ "$MODE" == "local" ]]; then
    # -------- LOCAL: do first-stage merges sequentially on this node --------
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

    # CLEAN old outputs (partials, superchunks, final)
    echo "[LOCAL][CLEAN] Removing old ${PARTIAL_PREFIX}_*.root, ${SUPER_PREFIX}_*.root, and $MERGED_FILE from $OUTPUT_DIR"
    find "$OUTPUT_DIR" -maxdepth 1 -type f \( -name "${PARTIAL_PREFIX}_*.root" -o -name "${SUPER_PREFIX}_*.root" -o -name "$MERGED_FILE" \) -delete
    rm -rf "$TMP_SUBDIR" "$TMP_SUPERDIR"

    # Build master list of chunk files
    find "$SIM_CHUNK_DIR" -maxdepth 1 -type f -name "PositionDep_sim_*.root" \
      | sort -V > "$LISTFILE"
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
  fi

  # ----------------------------- CONDOR path (unchanged) ---------------------
  testMode=false; firstHalf=false; asMany=false
  case "$SUBMODE" in
    test)          testMode=true  ;;
    firstHalf)     firstHalf=true ;;
    asManyAsCan)   asMany=true    ;;
    "")            ;;
    *) usage ;;
  esac

# -------------------------------------------------------------------------
# PRE-FLIGHT: optionally run the 'checkFileOutput' logic
# -------------------------------------------------------------------------
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

  # Deep check: verify each expected index can be matched to an output file
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
      if [[ -f "$c" ]]; then found=true; break; fi
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
fi

  # -------------------------------------------------------------------------
  # PLAN: For plain 'condor' with no sub-mode, use 1000 files per job (60k -> 60 jobs)
  # -------------------------------------------------------------------------
  if ! $testMode && ! $firstHalf && ! $asMany; then
    FILES_PER_GROUP=1000
    echo "[PLAN] No sub-mode supplied → forcing FILES_PER_GROUP=${FILES_PER_GROUP}"
    echo "       With ${presentOutputs} outputs, expected sublists ≈ $(( (presentOutputs + FILES_PER_GROUP - 1) / FILES_PER_GROUP ))"
  else
    # keep whatever FILES_PER_GROUP is (default 100 or 1000 for asManyAsCan)
    echo "[PLAN] Sub-mode='$SUBMODE' → FILES_PER_GROUP=${FILES_PER_GROUP}"
  fi

  # -------------------------------------------------------------------------
  # CLEAN old outputs (partials, superchunks, final) only after successful precheck
  # -------------------------------------------------------------------------
  echo "[CLEAN] Removing old ${PARTIAL_PREFIX}_*.root, ${SUPER_PREFIX}_*.root, and $MERGED_FILE from $OUTPUT_DIR"
  find "$OUTPUT_DIR" -maxdepth 1 -type f \( -name "${PARTIAL_PREFIX}_*.root" -o -name "${SUPER_PREFIX}_*.root" -o -name "$MERGED_FILE" \) -delete
  rm -rf "$TMP_SUBDIR" "$TMP_SUPERDIR"

find "$SIM_CHUNK_DIR" -maxdepth 1 -type f -name "PositionDep_sim_*.root" \
       | sort -V > "$LISTFILE"
  [[ ! -s "$LISTFILE" ]] && { echo "[ERROR] No chunk files found after precheck"; exit 1; }

  totalToMerge=$(wc -l < "$LISTFILE")
  echo "[INFO] Files to merge listed in $LISTFILE  (count=$totalToMerge)"

  mkdir -p "$TMP_SUBDIR" "$TMP_SUPERDIR"
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
  echo "[INFO] Grouping plan"
  echo "  • FILES_PER_GROUP   : $FILES_PER_GROUP"
  echo "  • Total sublists    : $totalSublists"
  $firstHalf && echo "  • firstHalf active   → queueN=$queueN"
  $testMode  && echo "  • test mode active   → queueN=$queueN (one sublist)"
  echo "  • Jobs to submit    : $queueN"
  echo "-----------------------------------------------------------------------"

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
  echo "[INFO] Submitted $queueN Condor jobs."
  exit 0
fi



###############################################################################
# 2) ADDCHUNKS  ---------------------------------------------------------------
###############################################################################
# New sub-modes:
#   addChunks part1          -> build superchunks from partials locally
#   addChunks part1 condor   -> build superchunks via Condor
# Plain:
#   addChunks [condor]       -> final merge; prefers superchunks if present,
#                               otherwise falls back to partials.
###############################################################################
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

