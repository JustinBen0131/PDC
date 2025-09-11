#!/usr/bin/env bash
###############################################################################
# PositionDependentCorrect_Condor_submit.sh
#
# Driver for Position-Dependent Correction (PDC) production on sPHENIX:
#   • Data-side Condor submissions (chunked by file lists)
#   • Simulation-side submissions (Condor or Local), with support for:
#       - round slicing (first/second/third/fourthRound)
#       - sample <N> subsets
#       - full “doAll” processing
#       - optional tracking of submitted pairs across rounds
#   • Golden-run segmentation helper for data (split into job-balanced segments)
#
# PATHS (key defaults)
#   MACRO_PATH            : /sphenix/u/patsfan753/scratch/PDCrun24pp/macros/Fun4All_PDC.C
#   DST_LIST_DIR          : /sphenix/u/patsfan753/scratch/PDCrun24pp/dst_list
#   OUTDIR_DATA           : /sphenix/tg/tg01/bulk/jbennett/PDC/output
#   OUTDIR_SIM            : /sphenix/tg/tg01/bulk/jbennett/PDC/SimOut
#   LOG_DIR               : /sphenix/user/${USER}/PositionDependentCorrect
#   CONDOR_LISTFILES_DIR  : ${LOG_DIR}/condorListFiles
#
# SPECIAL OUTPUT (LOCAL doAll for SIM ONLY)
#   • When you run:  ./PositionDependentCorrect_Condor_submit.sh simOnly local doAll
#     the script processes ALL sim pairs locally (single macro call) and writes:
#       /sphenix/u/patsfan753/scratch/PDCrun24pp/output/PositionDep_sim_ALL.root
#     This is equivalent to Condor “doAll” followed by the final addChunks merger.
#
# TOP-LEVEL MODES (first argument)
# -----------------------------------------------------------------------------
# 1)  (no args)           Submit BOTH data & sim on Condor (all events).
#       Example:  ./PositionDependentCorrect_Condor_submit.sh
#
# 2)  noSim               Data only on Condor (all events).
#       Example:  ./PositionDependentCorrect_Condor_submit.sh noSim
#
# 3)  dataOnly <submode>  Data-only helpers & wrappers:
#       submodes:
#         • local              Run locally (≤ LOCAL_NEVENTS) using your local helper.
#         • condorTest         Submit a very small test (first run only).
#         • condorFirstTen     Cap total submitted jobs at MAX_JOBS_DATA.
#         • condor             Full scan (legacy behavior).
#             - Optional:  round <N>   Submit only the Nth golden-run segment.
#         • splitGoldenRunList Build goldenRuns_segment*.txt from Final_RunNumbers.
#         • condorAll          Alias of plain “condor”.
#       Examples:
#         ./PositionDependentCorrect_Condor_submit.sh dataOnly condor
#         ./PositionDependentCorrect_Condor_submit.sh dataOnly condor round 3
#         ./PositionDependentCorrect_Condor_submit.sh dataOnly splitGoldenRunList
#
# 4)  condor [firstTen]   Shorthand for data-only Condor submit.
#       • firstTen         Cap total submitted jobs at MAX_JOBS_DATA.
#       Examples:
#         ./PositionDependentCorrect_Condor_submit.sh condor
#         ./PositionDependentCorrect_Condor_submit.sh condor firstTen
#
# 5)  condorTest          Shorthand for data-only “condorTest”.
#       Example:  ./PositionDependentCorrect_Condor_submit.sh condorTest
#
# 6)  simOnly <args…>     Simulation-side driver (passes remaining args to
#                         submit_sim_condor). Supported tokens (any order):
#       Execution mode (choose one):
#         • local              Run on this node.
#         • condor             Submit grouped jobs to Condor.
#       Round selection (choose at most one):
#         • firstRound         Use pairs [0 .. 14999], wipe outputs.
#         • secondRound        Use pairs [15000 .. 29999].
#         • thirdRound         Use pairs [30000 .. 44999].
#         • fourthRound        Use pairs [45000 .. 59999].
#         • testSubmit         Use exactly 5 pairs (quick smoke test).
#         • sample <N>         Use the first N pairs, wipe outputs.
#         • firstRoundCheck    Inventory only; no submission.
#       Global switches:
#         • doAll              Override selection and process ALL available pairs.
#
#       Behavior notes:
#         • simOnly condor doAll
#             - Groups all pairs into Condor jobs (GROUP_SIZE_SIM=5), writes
#               per-job ROOT files into OUTDIR_SIM. Use external merger if needed.
#         • simOnly local doAll
#             - Runs one macro over ALL pairs locally and writes a single final:
#               /sphenix/u/patsfan753/scratch/PDCrun24pp/output/PositionDep_sim_ALL.root
#
#       Examples:
#         ./PositionDependentCorrect_Condor_submit.sh simOnly condor firstRound
#         ./PositionDependentCorrect_Condor_submit.sh simOnly condor sample 200
#         ./PositionDependentCorrect_Condor_submit.sh simOnly condor doAll
#         ./PositionDependentCorrect_Condor_submit.sh simOnly local  doAll
#
# 7)  submitTestSimCondor  Partial sim submission helper (small subset).
#       Example:  ./PositionDependentCorrect_Condor_submit.sh submitTestSimCondor
#
# DATA-SIDE CONSTANTS
#   CHUNK_SIZE_DATA = 10     # files per data job
#   MAX_JOBS_DATA   = 10000  # cap used by “firstTen” mode
#
# SIM-SIDE CONSTANTS (submit_sim_condor)
#   GROUP_SIZE_SIM  = 5      # (DST,G4) pairs per Condor job
#   LOCAL_NEVENTS   = 10000  # used by dataOnly local helper
#
# GOLDEN-RUN SEGMENTER
#   dataOnly splitGoldenRunList
#     - Reads $GOLDEN_RUN_LIST (one run per line)
#     - Counts expected jobs per run (ceil(Nfiles / CHUNK_SIZE_DATA))
#     - Emits goldenRuns_segment*.txt so each segment’s total job count
#       stays ≤ MAX_JOBS_DATA; prints a detailed summary table.
#
# LOGGING / STAGING
#   Logs:   ${LOG_DIR}/{stdout,error,log}
#   Staging for list files:  ${CONDOR_LISTFILES_DIR}
#
# EXAMPLES (quick reference)
#   • BOTH (default):               ./PositionDependentCorrect_Condor_submit.sh
#   • DATA (full Condor):           ./PositionDependentCorrect_Condor_submit.sh noSim
#   • DATA (segment #3):            ./PositionDependentCorrect_Condor_submit.sh dataOnly condor round 3
#   • SIM  (all via Condor):        ./PositionDependentCorrect_Condor_submit.sh simOnly condor doAll
#   • SIM  (all locally → ALL):     ./PositionDependentCorrect_Condor_submit.sh simOnly local doAll
###############################################################################


#############################
# 0) CONFIGURATION
#############################
CLR_R='\033[0m'   # reset
CLR_B='\033[1m'   # bold
CLR_G='\033[1;32m'
CLR_Y='\033[1;33m'
CLR_C='\033[1;36m'
CLR_Rd='\033[1;31m'


# Path to your Fun4All macro:
MACRO_PATH="/sphenix/u/patsfan753/scratch/PDCrun24pp/macros/Fun4All_PDC.C"

# Where to send Condor output/logs
LOG_DIR="/sphenix/user/${USER}/PositionDependentCorrect"
CONDOR_LISTFILES_DIR="${LOG_DIR}/condorListFiles"

# Default output directories (if needed by your actual code)
OUTDIR_DATA="/sphenix/tg/tg01/bulk/jbennett/PDC/output"
OUTDIR_SIM="/sphenix/tg/tg01/bulk/jbennett/PDC/SimOut"

# Where data DST .list files are
DST_LIST_DIR="/sphenix/u/patsfan753/scratch/PDCrun24pp/dst_list"
# e.g. "dst_calo_run2pp-00047289.list"

# Simulation DST + G4Hits .list files
# Where sim DST + G4Hits lists were generated (the PhotonJet pT=5 sample).
: "${SIM_LIST_DIR:="/sphenix/u/patsfan753/scratch/PDCrun24pp/simListFiles/run24_type14_gamma_pt_200_40000"}"
: "${SIM_DST_LIST:="${SIM_LIST_DIR}/DST_CALO_CLUSTER.list"}"
: "${SIM_HITS_LIST:="${SIM_LIST_DIR}/G4Hits.list"}"

# How many files per Condor chunk
FILES_PER_CHUNK=10

# For local modes, limit to 10k events:
LOCAL_NEVENTS=10000

#############################
# 1) ARGUMENT PARSING
#############################
MODE="$1"  # "local", "noSim", "simOnly", "localSimTest", "submitTestSimCondor" or empty
WHAT="$2"  # "both", "noSim", or empty

usage()
{
  echo "Usage examples:"
  echo "  $0 dataOnly splitGoldenRunList      # create goldenRuns_segment*.txt"
  echo "  $0 dataOnly condor round <N>        # submit runs listed in segment N"
  echo "  $0 local both"
  echo "  $0 local noSim"
  echo "  $0 localSimTest"
  echo "  $0 noSim"
  echo "  $0 simOnly"
  echo "  $0 submitTestSimCondor"
  echo "  $0 dataOnly  local|condorTest|condorAll|condorFirstTen"
  echo "  $0                      (submits both data & sim by default)"
  exit 1
}


echo "===================================================================="
echo "[DEBUG] Starting PositionDependentCorrect_Condor_submit.sh"
echo "[DEBUG] MODE='$MODE'  WHAT='$WHAT'"
echo
echo "[DEBUG] Macro path       : $MACRO_PATH"
echo "[DEBUG] DST_LIST_DIR     : $DST_LIST_DIR"
echo "[DEBUG] SIM_LIST_DIR     : $SIM_LIST_DIR"
echo "[DEBUG] SIM_DST_LIST     : $SIM_DST_LIST"
echo "[DEBUG] SIM_HITS_LIST    : $SIM_HITS_LIST"
echo "[DEBUG] LOG_DIR          : $LOG_DIR"
echo "[DEBUG] LOCAL_NEVENTS    : $LOCAL_NEVENTS"
echo "===================================================================="
echo

###############################################################################
#   dataOnly splitGoldenRunList
#        • Reads   $GOLDEN_RUN_LIST   (one run number / line)
#        • Computes how many Condor jobs each run would generate
#          (#jobs  =  ceil(  Nfiles / CHUNK_SIZE_DATA  ))
#        • Writes successive segments
#                 goldenRuns_segment1.txt, goldenRuns_segment2.txt, …
#          so that       Σ #jobs   ≤  MAX_JOBS_DATA      in each file
#        • Output location =  $GOLDEN_SPLIT_DIR
###############################################################################

GOLDEN_SPLIT_DIR="/sphenix/u/patsfan753/scratch/PDCrun24pp"
GOLDEN_RUN_LIST="${GOLDEN_SPLIT_DIR}/Final_RunNumbers_After_All_Cuts.txt"
GOLDEN_SEGMENT_PREFIX="${GOLDEN_SPLIT_DIR}/goldenRuns_segment"

split_golden_run_list()(
    set -e  # <- subshell only

    [[ -f "$GOLDEN_RUN_LIST" ]] ||
        { echo -e "${CLR_Rd}[ERROR]${CLR_R} Golden‑run list not found → $GOLDEN_RUN_LIST"; return 1; }

    echo -e "${CLR_B}────────────────────────────────────────────────────────────────────────────"
    echo -e "          GOLDEN‑RUN SPLITTER – VERBOSE DIAGNOSTICS"
    echo -e "          Source file  : ${CLR_C}$(basename "$GOLDEN_RUN_LIST")${CLR_R}"
    echo -e "          Target files : ${CLR_C}$(basename "${GOLDEN_SEGMENT_PREFIX}")N.txt${CLR_R}"
    echo -e "          Parameters   : CHUNK_SIZE_DATA=${CLR_Y}${CHUNK_SIZE_DATA}${CLR_R},  "\
            "MAX_JOBS_DATA=${CLR_Y}${MAX_JOBS_DATA}${CLR_R}"
    echo -e "────────────────────────────────────────────────────────────────────────────${CLR_R}"

    # ───────────────— global counters ──────────────────────────────────────
    local segment_no=1
    local jobs_in_segment=0  runs_in_segment=0
    local tot_runs=0  tot_jobs=0  tot_files=0
    local current_file="${GOLDEN_SEGMENT_PREFIX}${segment_no}.txt"
    : > "$current_file"

    # keep stats per finished segment
    declare -a segJobs  segRuns  segFiles

    # ---------- helper: flush current segment and start a new one ----------
    flush_segment() {
        segJobs[segment_no]=$jobs_in_segment
        segRuns[segment_no]=$runs_in_segment
        segFiles[segment_no]=$files_in_segment

        printf "${CLR_G}[SEGMENT %2d]${CLR_R} %-25s  runs=%3d  jobs=%5d  files=%6d\n" \
               "$segment_no" "$(basename "$current_file")" \
               "$runs_in_segment" "$jobs_in_segment" "$files_in_segment"

        (( ++segment_no ))
        current_file="${GOLDEN_SEGMENT_PREFIX}${segment_no}.txt"
        : > "$current_file"
        jobs_in_segment=0
        runs_in_segment=0
        files_in_segment=0
    }

    # header for run‑level table
    printf "${CLR_B}%-10s│%10s│%10s│%10s│%9s│%10s${CLR_R}\n" \
           "Run" "Files" "Jobs(run)" "Seg#Now" "J(seg)" "J(seg)%"
    printf "──────────┼──────────┼──────────┼──────────┼─────────┼──────────\n"

    # ───────────────— MAIN LOOP OVER RUN NUMBERS ───────────────────────────
    local files_in_segment=0
    while IFS= read -r runRaw; do
        [[ -z "$runRaw" || "$runRaw" =~ ^# ]] && continue
        local run=$(printf "%08d" "$runRaw")
        local list_file="${DST_LIST_DIR}/dst_calo_run2pp-${run}.list"

        if [[ ! -f "$list_file" ]]; then
            echo -e "${CLR_Y}[WARN]${CLR_R} No DST list for run ${run} – skipped"
            continue
        fi

        local nFiles; nFiles=$(wc -l < "$list_file")
        (( nFiles )) || { echo -e "${CLR_Y}[WARN]${CLR_R} Empty list for run ${run} – skipped"; continue; }

        local nJobs=$(( (nFiles + CHUNK_SIZE_DATA - 1) / CHUNK_SIZE_DATA ))

        # overflow? → close segment first
        if (( jobs_in_segment + nJobs > MAX_JOBS_DATA )); then
            flush_segment
        fi

        # append run to current segment file
        echo "$run" >> "$current_file"

        # update counters
        (( jobs_in_segment += nJobs ))
        (( runs_in_segment += 1 ))
        (( files_in_segment += nFiles ))
        (( tot_jobs += nJobs ))
        (( tot_runs += 1 ))
        (( tot_files += nFiles ))

        # live line
        local pct=$(( 100*jobs_in_segment / MAX_JOBS_DATA ))
        printf "%-10s│%10d│%10d│%10d│%9d│%9s%%\n" \
               "$run" "$nFiles" "$nJobs" "$segment_no" "$jobs_in_segment" "$pct"
    done < "$GOLDEN_RUN_LIST"

    # flush last segment
    flush_segment

    # ───────────────— SUMMARY ──────────────────────────────────────────────
    echo -e "${CLR_B}\n────────────────────────  SUMMARY  ────────────────────────${CLR_R}"
    printf "Segments produced : %s\n" "${#segJobs[@]}"
    printf "Total runs        : %s\n" "$tot_runs"
    printf "Total DST files   : %s\n" "$tot_files"
    printf "Total Condor jobs : %s  (theoretical)\n" "$tot_jobs"
    echo "-----------------------------------------------------------"
    printf "%-10s | %-5s jobs | %-5s runs | %-6s files\n" "Segment" "≈" "≈" "≈"
    echo "-----------------------------------------------------------"
    for s in "${!segJobs[@]}"; do
        printf "segment%-3d | %7d | %7d | %9d\n" \
               "$s" "${segJobs[$s]}" "${segRuns[$s]}" "${segFiles[$s]}"
    done
    echo -e "${CLR_G}[OK]${CLR_R} Golden‑run list split complete → files in  ${CLR_C}$GOLDEN_SPLIT_DIR${CLR_R}"
)

##############################################################################
# 2)  DATA‑SIDE CONDOR SUBMISSION (full‑run, chunked)
##############################################################################

# ---------- user‑tunable constants ----------
CHUNK_SIZE_DATA=10          # N  – files per Condor job
MAX_JOBS_DATA=10000         # Q  – global cap for “condor firstTen”
# ---------------------------------------------------------------------------


submit_data_condor() {
    # Signature:
    #   submit_data_condor  <runMode>  <limitSwitch>  [<runListFile>]
    #
    #     runMode       :  condor | condorTest
    #     limitSwitch   :  "" | firstTen      (ignored for condorTest)
    #     runListFile   :  optional flat file with run numbers (one / line);
    #                     if supplied we *only* submit the runs listed there.

    local runMode="${1:-condor}"
    local limitSwitch="${2:-}"
    local runListFile="${3:-}"        # ← NEW
    local jobLimit=0

    # ----------------  fancy verbosity flag --------------------------------
    local VERBOSE=0
    [[ "$runMode" == "condorTest" || "$limitSwitch" == "firstTen" ]] && VERBOSE=1
    vecho() { (( VERBOSE )) && echo "$@"; }

    # ----------------  global‑cap logic (old behaviour unchanged) ----------
    if [[ "$runMode" == "condor" ]]; then
    if [[ "$limitSwitch" == "firstTen" ]]; then
        jobLimit=$MAX_JOBS_DATA
    elif [[ "$limitSwitch" =~ ^[0-9]+$ ]]; then
        jobLimit=$limitSwitch          # user-supplied numeric limit
    fi
fi

    mkdir -p "$LOG_DIR"/{stdout,error} "$CONDOR_LISTFILES_DIR"

    # ----------------  (1) build the array   listFiles[]  -------------------
    declare -a listFiles

    if [[ -n "$runListFile" ]]; then
        # ---------- ROUND N MODE -------------------------------------------
        [[ -f "$runListFile" ]] || {
            echo "[ERROR] Supplied run‑list file not found → $runListFile"; return 1; }
        vecho "[VERBOSE] Using external run‑list file → $(basename "$runListFile")"

        while IFS= read -r rn; do
            [[ -z "$rn" || "$rn" =~ ^# ]] && continue
            rn=$(printf "%08d" "$rn")      # normalise to 8 digits
            lf="${DST_LIST_DIR}/dst_calo_run2pp-${rn}.list"
            if [[ -f "$lf" ]]; then
                listFiles+=( "$lf" )
            else
                echo "[WARN] Missing DST list for run ${rn} (skipped)"
            fi
        done < "$runListFile"

    else
        # ---------- ORIGINAL DIRECTORY SCAN --------------------------------
        mapfile -t listFiles < <(ls -1 "${DST_LIST_DIR}"/dst_calo_run2pp-*.list 2>/dev/null)
    fi

    (( ${#listFiles[@]} )) || { echo "[ERROR] No run lists selected → aborting."; return 1; }

    # ----------------  (2) info banner -------------------------------------
    vecho "[VERBOSE] Total run‑lists selected      : ${#listFiles[@]}"
    vecho "[VERBOSE] CHUNK_SIZE_DATA               : $CHUNK_SIZE_DATA"
    (( jobLimit )) && vecho "[VERBOSE] Global job cap (Q)        : $jobLimit"

    # ----------------  (3) loop over runs ----------------------------------
    local submitted=0  runCounter=0

    for dlist in "${listFiles[@]}"; do
        ((++runCounter))
        local runBase runNum
        runBase="$(basename "$dlist")"
        runNum="${runBase#dst_calo_run2pp-}"
        runNum="${runNum%.*}"

        mapfile -t allfiles < "$dlist"
        local total=${#allfiles[@]}
        (( total )) || { echo "[WARN] $runBase is empty – skipping"; continue; }
        vecho "[VERBOSE] ---- run $runNum  (files=$total)"

        # --- split the file into   CHUNK_SIZE_DATA‑line   parts -------------
        rm -f "${CONDOR_LISTFILES_DIR}/run${runNum}_chunk_"* 2>/dev/null || true
        split -l "$CHUNK_SIZE_DATA" -d -a 3 "$dlist" "${CONDOR_LISTFILES_DIR}/run${runNum}_chunk_"

        mapfile -t chunks < <(ls "${CONDOR_LISTFILES_DIR}/run${runNum}_chunk_"* 2>/dev/null)
        local nChunks=${#chunks[@]}
        (( nChunks )) || { echo "[ERROR] split produced zero chunks for run $runNum"; return 2; }

        # --- global‑cap guard ----------------------------------------------
        if (( jobLimit && submitted + nChunks > jobLimit )); then
            echo "[INFO] Adding run $runNum would exceed cap ($jobLimit) – stop."
            break
        fi

        # --- create one .sub per run (atomic‑run rule) ----------------------
        local subFile="PositionDependentCorrect_data_${runNum}.sub"
        cat > "$subFile" <<EOL
universe      = vanilla
executable    = PositionDependentCorrect_Condor.sh
log           = $LOG_DIR/job.\$(Cluster).\$(Process).log
output        = $LOG_DIR/stdout/job.\$(Cluster).\$(Process).out
error         = $LOG_DIR/error/job.\$(Cluster).\$(Process).err
request_memory= 1000MB
# Args: <runNum> <listFile> <data|sim> <clusterID> <nEvents> <chunkIdx> <hitsList> <destBase>
EOL

        local chunkIdx=0
        for chunkFile in "${chunks[@]}"; do
            ((++chunkIdx))
            printf 'arguments = %s %s data $(Cluster) 0 %d "" %s\nqueue\n\n' \
                   "$runNum" "$chunkFile" "$chunkIdx" "$OUTDIR_DATA" >> "$subFile"
        done

        echo "[INFO] condor_submit → $subFile  (jobs=$nChunks)"
        condor_submit "$subFile"

        (( submitted += nChunks ))

        [[ "$runMode" == "condorTest" ]] && { vecho "[VERBOSE] condorTest stops after first run"; break; }
    done

    echo "[INFO] Grand‑total data jobs submitted : $submitted"
    (( jobLimit )) && echo "[INFO] Global‑cap mode active (cap=$jobLimit)"
}



###############################################################################
# 3) submit_sim_condor  –  drive Position‑Dependent‑Correction jobs
#                      • “first|second|third|fourthRound|sample|testSubmit”
#                      • optional integer for sample size
#                      • execution mode:  condor | local
#
#   simOnly            – accepted as an alias that simply forces execMode=condor
#
# CLEANNING POLICY
# ─────────────────────────────────────────────────────────────────────────────
#   round keyword   what is wiped *before* submission
#   --------------------------------------------------------------------------
#   firstRound      $OUTDIR_SIM  +  log/ stdout/ error/
#   sample          $OUTDIR_SIM  +  log/ stdout/ error/
#   secondRound…fourthRound   log/ stdout/ error/   (ROOT files kept)
#   testSubmit      nothing
###############################################################################
submit_sim_condor() {

  set -Eeuo pipefail
  trap 'echo "[FATAL] Unexpected error near line $LINENO" >&2' ERR

  # ─────────────── 0) CONSTANT PATHS ───────────────────────────────────────
  local LOG_BASE="/sphenix/u/patsfan753/scratch/PDCrun24pp"
  local LOG_DIR_LOG="${LOG_BASE}/log"
  local LOG_DIR_OUT="${LOG_BASE}/stdout"
  local LOG_DIR_ERR="${LOG_BASE}/error"
  local OUTDIR_SIM="/sphenix/tg/tg01/bulk/jbennett/PDC/SimOut"
  local CONDOR_LISTFILES_DIR="${LOG_BASE}/condorListFiles"

  # round‑tracking (pairs remaining to be submitted across rounds)
  local TRACK_DIR="${LOG_BASE}/round_tracking"
  local TRACK_FILE="${TRACK_DIR}/sim_pairs_remaining.txt"

  mkdir -p "$LOG_DIR_LOG" "$LOG_DIR_OUT" "$LOG_DIR_ERR" "$CONDOR_LISTFILES_DIR"

  clean_logs() {
      echo "[CLEAN] Removing and recreating Condor log directories under: $LOG_BASE"
      echo "[CLEAN] rm -rfv \"$LOG_DIR_OUT\" \"$LOG_DIR_ERR\" \"$LOG_DIR_LOG\""
      rm -rfv "$LOG_DIR_OUT" "$LOG_DIR_ERR" "$LOG_DIR_LOG" 2>/dev/null || true
      echo "[CLEAN] mkdir -p \"$LOG_DIR_OUT\" \"$LOG_DIR_ERR\" \"$LOG_DIR_LOG\""
      mkdir -p "$LOG_DIR_OUT" "$LOG_DIR_ERR" "$LOG_DIR_LOG"
  }

  local roundArg="firstRound"     # default slice keyword
  local execMode="condor"         # default execution mode
  local sampleN=""                # set only if user supplied a pure integer
  local checkOnly=false           # when true, only count inputs; no submission
  local doAll=false               # when true, submit ALL pairs in one go

  for tok in "$@"; do
      case "$tok" in
          local|condor)
              execMode="$tok"
              ;;
          simOnly)
              execMode="condor"   # alias
              ;;
          doAll)
              doAll=true          # submit ALL pairs in one submission
              ;;
          firstRound|secondRound|thirdRound|fourthRound|testSubmit|sample)
              roundArg="$tok"
              ;;
          firstRoundCheck)
              roundArg="firstRound"   # use firstRound slice
              checkOnly=true          # but do not submit
              ;;
          ''|*[!0-9]*)               # ignore non‑numeric tokens
              ;;
          *)
              sampleN="$tok"         # pure integer => sample size
              ;;
      esac
  done

  # numeric without “sample” ⇒ treat as   sample <N>
  if [[ -n "$sampleN" && "$roundArg" != "sample" ]]; then
      roundArg="sample"
  fi

  # ─────────────── 2) DETERMINE OFFSET & CHUNK SIZE ───────────────────────
  local chunkSize=15000
  local offset=0
  local wipe_outputs=false        # whether $OUTDIR_SIM is deleted

  case "$roundArg" in
      firstRound)  offset=0        ; wipe_outputs=true ;;
      secondRound) offset=15000 ;;
      thirdRound)  offset=30000 ;;
      fourthRound) offset=45000 ;;
      testSubmit)  chunkSize=5 ;;
      sample)
          [[ "$sampleN" =~ ^[1-9][0-9]*$ ]] \
              || { echo "[ERROR] sample requires a positive integer"; return 1; }
          chunkSize="$sampleN"
          wipe_outputs=true ;;
      *) echo "[ERROR] Unknown round keyword '$roundArg'"; return 1 ;;
  esac

  echo "######################################################################"
  echo "[INFO] execMode='$execMode'   round='$roundArg'"
  echo "[INFO] offset=$offset   chunkSize=$chunkSize"
  echo "######################################################################"

  # ─────────────── 3) CLEAN-UP PRIOR TO SUBMISSION ────────────────────────
  if [[ "$execMode" == "condor" ]]; then
      if [[ "$wipe_outputs" == true ]]; then
          echo "[STEP] Removing previous ROOT outputs from $OUTDIR_SIM"
          [[ -d "$OUTDIR_SIM" ]] && find "$OUTDIR_SIM" -mindepth 1 -delete
      fi

      echo "[STEP] Cleaning old Condor logs (stdout, error, log) and recreating"
      echo "[CLEAN] rm -rfv \"$LOG_DIR_OUT\" \"$LOG_DIR_ERR\" \"$LOG_DIR_LOG\""
      rm -rfv "$LOG_DIR_OUT" "$LOG_DIR_ERR" "$LOG_DIR_LOG" 2>/dev/null || true

      echo "[CLEAN] mkdir -p \"$LOG_DIR_OUT\" \"$LOG_DIR_ERR\" \"$LOG_DIR_LOG\""
      mkdir -p "$LOG_DIR_OUT" "$LOG_DIR_ERR" "$LOG_DIR_LOG"

      # also reset the listfiles staging dir used by submissions
      echo "[CLEAN] rm -rfv \"$CONDOR_LISTFILES_DIR\" && mkdir -p \"$CONDOR_LISTFILES_DIR\""
      rm -rfv "$CONDOR_LISTFILES_DIR" 2>/dev/null || true
      mkdir -p "$CONDOR_LISTFILES_DIR"

      # visible proof the directories exist and are empty
      echo "[CLEAN] Post-cleanup directory state:"
      ls -ld "$LOG_DIR_OUT" "$LOG_DIR_ERR" "$LOG_DIR_LOG" "$CONDOR_LISTFILES_DIR"
      for d in "$LOG_DIR_OUT" "$LOG_DIR_ERR" "$LOG_DIR_LOG"; do
          echo "[CLEAN] Contents of $d (should be empty):"
          ls -la "$d"
      done
  fi

  # ─────────────── 4) BUILD THE LIST OF (DST,G4) PAIRS ───────────────────
  [[ -f "$SIM_DST_LIST"  ]] || { echo "[ERROR] Missing DST list $SIM_DST_LIST";  return 1; }
  [[ -f "$SIM_HITS_LIST" ]] || { echo "[ERROR] Missing G4 list  $SIM_HITS_LIST"; return 1; }

  mapfile -t simPairs < <(paste -d' ' "$SIM_DST_LIST" "$SIM_HITS_LIST" | sort -k1,1V)

  local maxN=${#simPairs[@]}
  (( offset < maxN )) || { echo "[WARN] Offset beyond list size – nothing to do."; return 0; }

  if [[ "$doAll" == true ]]; then
      echo "[INFO] doAll requested → overriding slice to submit ALL ${maxN} pairs in one submission."
      offset=0
      chunkSize=$maxN
      wipe_outputs=true      # ensure first-round-style full reset semantics
  fi

  local endIndex=$(( offset + chunkSize ))
  if (( endIndex > maxN )); then
      endIndex=$maxN
  fi
  local nJobs=$(( endIndex - offset ))
  (( nJobs > 0 )) || { echo "[ERROR] Selected zero pairs"; return 1; }

  # Always show an inventory of available inputs; and optionally short-circuit
  local dstLines;  dstLines=$(wc -l < "$SIM_DST_LIST"  || echo 0)
  local hitsLines; hitsLines=$(wc -l < "$SIM_HITS_LIST" || echo 0)

  # Preview the number of Condor jobs: ceil(nPairs / groupSizePreview).
  # Keep this preview in sync with the grouping used later (GROUP_SIZE_SIM).
  # If you use a larger grouping in doAll, reflect it here.
  local groupSizePreview=5
  local plannedJobs=$(( (nJobs + groupSizePreview - 1) / groupSizePreview ))

  echo "=================================================================="
  echo "[INVENTORY] Simulation input inventory"
  echo "  DST list : $(basename "$SIM_DST_LIST")"
  echo "    └─ total lines = ${dstLines}"
  echo "  HITS list: $(basename "$SIM_HITS_LIST")"
  echo "    └─ total lines = ${hitsLines}"
  echo "  Pairs available (post-merge) = ${maxN}"
  echo "  Round selection indices      = [${offset} .. $((endIndex-1))]"
  echo "  Selected pairs this round    = ${nJobs}"
  echo "  Planned Condor jobs          = ${plannedJobs}  (group size = ${groupSizePreview})"
  echo "  Mode flags                   = doAll=${doAll}  FAST_SUBMIT=$([[ "$doAll" == true ]] && echo true || echo false)  execMode=${execMode}"
  echo "  Preview first/last 3 selected pairs:"
  for (( k=0; k<3 && (offset+k)<endIndex; ++k )); do
      echo "    • $(printf "%6d" $((offset+k))) : ${simPairs[$((offset+k))]}"
  done
  if (( endIndex-offset > 6 )); then
      echo "    …"
  fi
  for (( k=3; k>=1; --k )); do
      if (( endIndex-k >= offset )); then
        echo "    • $(printf "%6d" $((endIndex-k))) : ${simPairs[$((endIndex-k))]}"
      fi
  done
  echo "=================================================================="

  # ---- tracking file initialization (skip entirely in doAll fast mode) ----
  local FAST_SUBMIT=false
  [[ "$doAll" == true ]] && FAST_SUBMIT=true

  if ! $FAST_SUBMIT; then
        echo "------------------------------------------------------------------"
        if [[ "$roundArg" == "firstRound" ]]; then
            mkdir -p "$TRACK_DIR"
            echo "[TRACK] firstRound: removing any previous tracking file and recreating fresh one"
            rm -f "$TRACK_FILE" 2>/dev/null || true
            printf "%s\n" "${simPairs[@]}" > "$TRACK_FILE"
        else
            if [[ ! -f "$TRACK_FILE" ]]; then
                mkdir -p "$TRACK_DIR"
                printf "%s\n" "${simPairs[@]}" > "$TRACK_FILE"
                echo "[TRACK][WARN] Tracking file was missing; re-created full list:"
            else
                echo "[TRACK] Continuing with existing tracking file:"
            fi
        fi

        # Light snapshot only
        local track_before_lines=0
        [[ -f "$TRACK_FILE" ]] && track_before_lines=$(wc -l < "$TRACK_FILE" 2>/dev/null || echo 0)
        echo "  • Tracking path : $TRACK_FILE"
        echo "  • Lines (BEFORE): ${track_before_lines}"
        echo "------------------------------------------------------------------"
  fi

  if [[ "$checkOnly" == true ]]; then
      echo "[INFO] firstRoundCheck: count-only mode; no Condor submission."
      echo "[INFO] Exiting after inventory + tracking snapshot."
      return 0
  fi

  # ─────────────── 5) LOCAL MODE  (bundle X pairs and run once) ───────────
  if [[ "$execMode" == "local" ]]; then
      echo "[STEP] LOCAL mode – bundling $nJobs pair(s) into one macro run"

      # Build temp list files for this local bundle
      local stamp; stamp="$(date +%Y%m%d_%H%M%S)"
      local dstList="${CONDOR_LISTFILES_DIR}/local_bundle_DST_${offset}_${endIndex}_${stamp}.list"
      local hitsList="${CONDOR_LISTFILES_DIR}/local_bundle_HITS_${offset}_${endIndex}_${stamp}.list"
      rm -f "$dstList" "$hitsList" 2>/dev/null || true

      local idx=0
      for (( i=offset; i<endIndex; i++ )); do
          idx=$(( idx + 1 ))
          IFS=' ' read -r dstFile hitFile <<< "${simPairs[$i]}"

          [[ -f "$dstFile" ]] || { echo "[WARN] Missing DST  $dstFile – skip"; continue; }
          [[ -f "$hitFile" ]] || { echo "[WARN] Missing HITS $hitFile – skip"; continue; }

          echo "$dstFile"  >> "$dstList"
          echo "$hitFile"  >> "$hitsList"
      done

      # Quick sanity
      local nd=$(wc -l < "$dstList"  || echo 0)
      local nh=$(wc -l < "$hitsList" || echo 0)
      if (( nd == 0 || nh == 0 || nd != nh )); then
          echo "[ERROR] Local bundle lists invalid (DST=$nd, HITS=$nh) – abort local run"
          return 1
      fi

      # If user requested "doAll" in local mode, write the final ALL file directly
      local finalOutDir="/sphenix/u/patsfan753/scratch/PDCrun24pp/output"
      local outRoot
      if $FAST_SUBMIT; then
          mkdir -p "$finalOutDir"
          outRoot="${finalOutDir}/PositionDep_sim_ALL.root"
          rm -f "$outRoot" 2>/dev/null || true
      else
          outRoot="${OUTDIR_SIM}/PositionDep_sim_bundle_${offset}_${endIndex}_${stamp}.root"
      fi

      echo "------------------------------------------------------------------"
      echo "[RUN]  single macro call over ${nd} pairs"
      echo "   DST-LIST : $(basename "$dstList")"
      echo "   HITS-LIST: $(basename "$hitsList")"
      echo "   OUT      : $(basename "$outRoot")"
      echo "------------------------------------------------------------------"

      # nevents=0 → process all events across all files in the lists
      root -b -q -l "${MACRO_PATH}(0, \"${dstList}\", \"${hitsList}\", \"${outRoot}\")"

      if $FAST_SUBMIT; then
          echo "[INFO] LOCAL doAll complete → $outRoot"
      else
          echo "[INFO] LOCAL bundle complete → $outRoot"
      fi
      return 0
  fi

  # ─────────────── 6) CONDOR SUBMISSION (group N pairs per job) ───────────
  local submitFile="PositionDependentCorrect_sim.sub"
  cat > "$submitFile" <<EOL
universe   = vanilla
executable = PositionDependentCorrect_Condor.sh
log        = ${LOG_DIR_LOG}/job.\$(Cluster).\$(Process).log
output     = ${LOG_DIR_OUT}/job.\$(Cluster).\$(Process).out
error      = ${LOG_DIR_ERR}/job.\$(Cluster).\$(Process).err
request_memory = 1500MB
# Args: <runNum=9999> <listFileData> <sim> <Cluster> <nEvents=0> <processID> [<listFileHits>]
EOL

  # How many (DST,G4) pairs per Condor job:
  local GROUP_SIZE_SIM=5

  # Verbose grouping plan
  echo "[PLAN] Building grouped lists for submission:"
  echo "       • GROUP_SIZE_SIM=${GROUP_SIZE_SIM}"
  echo "       • Pair range     = [${offset} .. $((endIndex-1))]"
  echo "       • Total pairs    = ${nJobs}"
  echo "       • Expected jobs  = $(( (nJobs + GROUP_SIZE_SIM - 1) / GROUP_SIZE_SIM ))"

  local groupIdx=0
  local i=$offset
  while (( i < endIndex )); do
        local gStart=$i
        local gEnd=$(( i + GROUP_SIZE_SIM ))
        if (( gEnd > endIndex )); then
            gEnd=$endIndex
        fi
        local pairsInGroup=$(( gEnd - gStart ))

        local dstGroup="${CONDOR_LISTFILES_DIR}/sim_group_DST_${gStart}_${gEnd}.list"
        local hitGroup="${CONDOR_LISTFILES_DIR}/sim_group_HITS_${gStart}_${gEnd}.list"
        rm -f "$dstGroup" "$hitGroup" 2>/dev/null || true

        # Fill the group list files
        local j
        for (( j=gStart; j<gEnd; j++ )); do
            IFS=' ' read -r dstFile hitFile <<< "${simPairs[$j]}"
            echo "$dstFile" >> "$dstGroup"
            echo "$hitFile" >> "$hitGroup"
        done

        # Log every group (comment out to throttle)
        echo "[GROUP $(printf "%6d" "${groupIdx}")] indices=${gStart}..$((gEnd-1))  pairs=${pairsInGroup}"
        echo "    dstList: ${dstGroup}"
        echo "    hitList: ${hitGroup}"

        # One job per group: pass the two list files to the wrapper
        printf 'arguments = 9999 %s sim $(Cluster) 0 %d %s\nqueue\n\n' \
               "$dstGroup" "$groupIdx" "$hitGroup" >> "$submitFile"

        (( ++groupIdx ))
        i=$gEnd
  done
  echo "[PLAN] Finished grouping: total jobs queued in submit file = ${groupIdx}"

  echo "[STEP] condor_submit → $submitFile  (jobs=$groupIdx)"
  if condor_submit "$submitFile"; then
      echo "=================================================================="
      echo "[INFO] condor_submit succeeded."
      echo "=================================================================="

      # ── round-tracking: remove submitted pairs from tracking file (skip in FAST_SUBMIT) ──
      if ! $FAST_SUBMIT; then
            if [[ -f "$TRACK_FILE" ]]; then
                local stamp; stamp="$(date +%Y%m%d_%H%M%S)"
                local submitted_list="${TRACK_DIR}/submitted_${roundArg}_${stamp}.txt"
                local tmp_remove="${TRACK_DIR}/to_remove.tmp"
                local tmp_remaining="${TRACK_DIR}/remaining.tmp"

                rm -f "$tmp_remove" "$tmp_remaining" 2>/dev/null || true

                echo "------------------------------------------------------------------"
                echo "[TRACK] Building submitted-this-round list"
                echo "  • Range indices: [${offset} .. $((endIndex-1))]"
                echo "  • Output file  : $submitted_list"
                echo "------------------------------------------------------------------"

                for (( i=offset; i<endIndex; i++ )); do
                    printf "%s\n" "${simPairs[$i]}" >> "$tmp_remove"
                    printf "%s\n" "${simPairs[$i]}" >> "$submitted_list"
                done

                grep -F -x -v -f "$tmp_remove" "$TRACK_FILE" > "$tmp_remaining" || true
                mv -f "$tmp_remaining" "$TRACK_FILE"
                rm -f "$tmp_remove" 2>/dev/null || true

                local after_lines=0
                after_lines=$(wc -l < "$TRACK_FILE" 2>/dev/null || echo 0)

                local total_pairs="$maxN"
                local processed_pairs=$(( total_pairs - after_lines ))
                local pct_done=0
                if [[ "$total_pairs" -gt 0 ]]; then
                    pct_done=$(( 100 * processed_pairs / total_pairs ))
                fi
                printf "[SUMMARY] Round: %-12s | submitted: %6d | total: %6d | remaining: %6d | done: %3d%%\n" \
                       "$roundArg" "$nJobs" "$total_pairs" "$after_lines" "$pct_done"
            else
                echo "[TRACK][WARN] Tracking file not found; skipping removal step."
            fi
      fi
  else
      echo "[ERROR] condor_submit failed"
      return 1
  fi
}


#############################
# 4) MAIN LOGIC
#############################
case "$MODE" in

  # -----------------------------------------------------------------------
  #  LEGACY SHORT‑CUTS  (kept for backward‑compatibility)
  # -----------------------------------------------------------------------
  noSim)         # data only, unlimited jobs
    echo "===================================================================="
    echo "[INFO] Condor => data only (ALL events)."
    echo "===================================================================="
    submit_data_condor          # defaults → runMode=condor limitSwitch=""
    exit 0
    ;;

  simOnly)       # simulation path is unchanged
    echo "===================================================================="
    echo "[INFO] Condor => sim only. round='$WHAT'"
    echo "===================================================================="
    submit_sim_condor "${@:2}"          # forward all args after MODE (e.g., condor firstRound doAll)
    exit 0
    ;;

  submitTestSimCondor)
    echo "===================================================================="
    echo "[INFO] Condor => TEST sim only (partial)."
    echo "===================================================================="
    submit_test_sim_condor
    exit 0
    ;;

  # -----------------------------------------------------------------------
  #  NEW ALIASES / WRAPPERS  (data‑side convenience)
  # -----------------------------------------------------------------------
  condor)          # syntax:  condor [firstTen]
    submit_data_condor "condor" "$WHAT"
    exit 0
    ;;

  condorTest)
    submit_data_condor "condorTest" ""
    exit 0
    ;;

  dataOnly)
      case "$WHAT" in
            local)
                echo "===================================================================="
                echo "[INFO] dataOnly local – running ≤ $LOCAL_NEVENTS events"
                echo "===================================================================="
                local_run_data ;;

            condorTest)
                submit_data_condor "condorTest" "" ;;

            condorFirstTen)
                submit_data_condor "condor" "firstTen" ;;

            condor)
                # -----------------------------------------------
                #   • dataOnly condor           → full scan (old)
                #   • dataOnly condor round N   → submit segment N
                # -----------------------------------------------
                if [[ "$3" == "round" && "$4" =~ ^[0-9]+$ ]]; then
                    segFile="${GOLDEN_SEGMENT_PREFIX}${4}.txt"
                    submit_data_condor "condor" "" "$segFile"
                else
                    submit_data_condor "condor" ""          # legacy path
                fi ;;

            splitGoldenRunList)
                split_golden_run_list ;;

            ""|condorAll)
                submit_data_condor "condor" "" ;;

            *)
                echo "[ERROR] Unknown dataOnly sub‑mode “$WHAT”"; usage ;;
        esac
        exit 0 ;;


  # -----------------------------------------------------------------------
  #  DEFAULT: submit both data & sim
  # -----------------------------------------------------------------------
  "")
    echo "===================================================================="
    echo "[INFO] Condor => BOTH data & sim (ALL events)."
    echo "===================================================================="
    submit_data_condor
    submit_sim_condor
    exit 0
    ;;

  *)
    echo "[ERROR] Unrecognized usage: $MODE"
    usage
    ;;
esac

echo "[DEBUG] End of script."
