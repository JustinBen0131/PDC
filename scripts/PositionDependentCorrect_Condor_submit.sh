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
#   OUTDIR_SIM            : /sphenix/tvig/tg01/bulk/jbennett/PDC/SimOut
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
LOG_DIR="/sphenix/u/patsfan753/scratch/PDCrun24pp"
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
CHUNK_SIZE_DATA=2          # N  – files per Condor job
MAX_JOBS_DATA=10000         # Q  – global cap for “condor firstTen”
# ---------------------------------------------------------------------------


submit_data_condor() {
    # Signature:
    #   submit_data_condor  <runMode>  <limitSwitch>  [<runListFile>]
    #
    #     runMode       :  condor | condorTest
    #     limitSwitch   :  "" | firstTen      (ignored for condorTest)
    #     runListFile   :  optional flat file with run numbers (one / line)

    local runMode="${1:-condor}"
    local limitSwitch="${2:-}"
    local runListFile="${3:-}"
    local jobLimit=0

    # ---- helpers -----------------------------------------------------------
    vecho() { (( VERBOSE )) && echo -e "$@"; }
    die() { echo -e "[ERROR] $*" >&2; echo 1; }
    check_dir_writable() {
      local d="$1"
      [[ -d "$d" ]] || die "Missing directory: $d"
      [[ -w "$d" ]] || die "Directory not writable: $d"
    }

    # ---- verbosity & tools sanity -----------------------------------------
    local VERBOSE=0
    [[ "$runMode" == "condorTest" || "$limitSwitch" == "firstTen" ]] && VERBOSE=1

    command -v condor_submit >/dev/null 2>&1 || return $(die "condor_submit not in PATH")
    command -v split         >/dev/null 2>&1 || return $(die "split not in PATH")

    # ---- global cap logic (unchanged semantics) ---------------------------
    if [[ "$runMode" == "condor" ]]; then
      if   [[ "$limitSwitch" == "firstTen"  ]]; then jobLimit=$MAX_JOBS_DATA
      elif [[ "$limitSwitch" =~ ^[0-9]+$   ]]; then jobLimit=$limitSwitch
      fi
    fi

    # ---- directories & environment ----------------------------------------
    mkdir -p "$LOG_DIR"/{stdout,error,log} "$CONDOR_LISTFILES_DIR" || return $(die "mkdir failed for logs/list dir")
    check_dir_writable "$LOG_DIR"
    check_dir_writable "$LOG_DIR/stdout"
    check_dir_writable "$LOG_DIR/error"
    check_dir_writable "$LOG_DIR/log"
    check_dir_writable "$CONDOR_LISTFILES_DIR"
    mkdir -p "$OUTDIR_DATA" || return $(die "mkdir failed for OUTDIR_DATA=$OUTDIR_DATA")
    vecho "[VERBOSE] Log root         : $LOG_DIR"
    vecho "[VERBOSE] Staging (lists)  : $CONDOR_LISTFILES_DIR"
    vecho "[VERBOSE] OUTDIR_DATA      : $OUTDIR_DATA"

    # ---- build listFiles[] -------------------------------------------------
    declare -a listFiles
    if [[ -n "$runListFile" ]]; then
      [[ -f "$runListFile" ]] || return $(die "Run-list file not found → $runListFile")
      vecho "[VERBOSE] Using external run-list file → $(basename "$runListFile")"
      while IFS= read -r rn; do
            [[ -z "$rn" || "$rn" =~ ^# ]] && continue
            rn=$(printf "%08d" "$((10#$rn))")   # force base-10; keep 8-digit zero-pad
            lf="${DST_LIST_DIR}/dst_calo_run2pp-${rn}.list"
            if [[ -f "$lf" ]]; then
                listFiles+=( "$lf" )
            else
                echo "[WARN] Missing DST list for run ${rn} (skipped)"
            fi
      done < "$runListFile"
    else
      mapfile -t listFiles < <(ls -1 "${DST_LIST_DIR}"/dst_calo_run2pp-*.list 2>/dev/null)
    fi
    (( ${#listFiles[@]} )) || return $(die "No run lists selected under $DST_LIST_DIR")

    vecho "[VERBOSE] Total run-lists selected : ${#listFiles[@]}"
    vecho "[VERBOSE] CHUNK_SIZE_DATA          : $CHUNK_SIZE_DATA"
    (( jobLimit )) && vecho "[VERBOSE] Global job cap (Q)     : $jobLimit"

    # ---- loop over runs ----------------------------------------------------
    local submitted=0 runCounter=0

    for dlist in "${listFiles[@]}"; do
      ((++runCounter))
      local runBase runNum
      runBase="$(basename "$dlist")"
      runNum="${runBase#dst_calo_run2pp-}"
      runNum="${runNum%.*}"

      mapfile -t allfiles < "$dlist"
      local total=${#allfiles[@]}
      if (( total == 0 )); then echo "[WARN] $runBase is empty – skipping"; continue; fi

      echo "────────────────────────────────────────────────────────────────────"
      echo "[RUN] ${runNum}  | files=${total}  | list=${runBase}"

      # clean old fragments for this run in staging
      rm -f "${CONDOR_LISTFILES_DIR}/run${runNum}_chunk_"* 2>/dev/null || true

      # split into chunks of CHUNK_SIZE_DATA
      split -l "$CHUNK_SIZE_DATA" -d -a 3 "$dlist" "${CONDOR_LISTFILES_DIR}/run${runNum}_chunk_" \
        || return $(die "split failed for run $runNum")
   
      # ensure chunk files are recognized as list files by the macro
      for f in "${CONDOR_LISTFILES_DIR}/run${runNum}_chunk_"*; do
        mv -f "$f" "${f}.list"
      done

      mapfile -t chunks < <(ls -1 "${CONDOR_LISTFILES_DIR}/run${runNum}_chunk_"*.list 2>/dev/null)
      local nChunks=${#chunks[@]}
      (( nChunks )) || return $(die "split produced zero chunks for run $runNum")


      # preview first two chunks (lines)
      for p in "${chunks[@]:0:2}"; do
        local lc; lc=$(wc -l < "$p" 2>/dev/null || echo 0)
        vecho "[VERBOSE]   chunk $(basename "$p") : ${lc} line(s)"
      done

      # cap guard
      if (( jobLimit && submitted + nChunks > jobLimit )); then
        echo "[INFO] Adding run $runNum would exceed cap ($jobLimit). Stopping submission loop."
        break
      fi

      # ---- write submit file (atomic) -------------------------------------
      local subFile="PositionDependentCorrect_data_${runNum}.sub"
      cat > "$subFile" <<EOL
universe      = vanilla
executable    = /sphenix/u/patsfan753/scratch/PDCrun24pp/PositionDependentCorrect_Condor.sh
initialdir    = /sphenix/u/patsfan753/scratch/PDCrun24pp
getenv        = True
log           = /sphenix/u/patsfan753/scratch/PDCrun24pp/log/job.\$(Cluster).\$(Process).log
output        = /sphenix/u/patsfan753/scratch/PDCrun24pp/stdout/job.\$(Cluster).\$(Process).out
error         = /sphenix/u/patsfan753/scratch/PDCrun24pp/error/job.\$(Cluster).\$(Process).err
request_memory= 1000MB
should_transfer_files   = NO
stream_output           = True
stream_error            = True
# Minimal mode hints for the macro:
environment   = PDC_MODE=DATA;PDC_RUNNUMBER=24
# Args: <runNum> <listFile> <data|sim> <clusterID> <nEvents> <chunkIdx> <hitsList> <destBase>
EOL

      # queue one proc per chunk; pass NONE for hits-list to avoid quote issues
      local chunkIdx=0
      for chunkFile in "${chunks[@]}"; do
        ((++chunkIdx))
        # Validate chunk file lines to catch empty-split edge cases early
        local ln; ln=$(wc -l < "$chunkFile" 2>/dev/null || echo 0)
        if (( ln == 0 )); then
          echo "[WARN] Empty chunk file: $chunkFile (skipping)"
          continue
        fi
        printf 'arguments = %s %s data $(Cluster) 0 %d NONE %s\nqueue\n\n' \
               "$runNum" "$chunkFile" "$chunkIdx" "$OUTDIR_DATA" >> "$subFile"
      done

      # show a short preview of the submit file for transparency
      vecho "[VERBOSE] Submit file head:"
      vecho "$(sed -n '1,25p' "$subFile")"
      vecho "[VERBOSE] Submit file tail:"
      vecho "$(tail -n 5 "$subFile")"

      echo "[INFO] condor_submit → $subFile  (jobs queued from file = $nChunks)"
      if condor_submit "$subFile"; then
        (( submitted += nChunks ))
      else
        echo "[ERROR] condor_submit failed for $subFile"
        echo "[HINT] Check syntax above; verify log/stdout/error directories exist and are writable."
        # show the exact first two 'arguments' lines to spot quoting issues
        echo "[DEBUG] First two 'arguments' lines:"
        grep -m2 '^arguments =' "$subFile" || true
        return 1
      fi

      if [[ "$runMode" == "condorTest" ]]; then
        vecho "[VERBOSE] condorTest: submitted first run only; stopping."
        break
      fi
    done

    echo "===================================================================="
    echo "[INFO] Grand-total data jobs submitted : $submitted"
    (( jobLimit )) && echo "[INFO] Global-cap mode active (cap=$jobLimit)"
    return 0
}



local_run_data() {
  # Usage:
  #   local_run_data              -> run on the FIRST file of the FIRST run-list
  #   local_run_data 1            -> same as above (first file)
  #   local_run_data <N>          -> run on the first N files of the FIRST run-list
  #   local_run_data /path/run.list [N]  -> run on the first N files of the given run-list
  #
  # Notes:
  #   - Writes a temporary one-line (or N-line) list under ${CONDOR_LISTFILES_DIR}
  #   - Invokes the same wrapper (PositionDependentCorrect_Condor.sh) the Condor path uses
  #   - Respects LOCAL_NEVENTS for the event cap

  local arg1="${1:-1}"               # default: 1 file
  local chosenList=""
  local nFiles=1

  # If arg1 looks like a list path, use it; otherwise treat as a file count
  if [[ -f "$arg1" ]]; then
    chosenList="$arg1"
    nFiles=1
    # Optional second argument can override number of files
    if [[ -n "$2" && "$2" =~ ^[0-9]+$ ]]; then
      nFiles="$2"
    fi
  else
    # arg1 is a count (or empty)
    if [[ "$arg1" =~ ^[0-9]+$ ]]; then
      nFiles="$arg1"
    fi
    # pick the FIRST run-list in DST_LIST_DIR
    chosenList=$(ls -1 "${DST_LIST_DIR}"/dst_calo_run2pp-*.list 2>/dev/null | head -n 1)
  fi

  if [[ -z "$chosenList" || ! -f "$chosenList" ]]; then
    echo "[ERROR] No data run-list found under DST_LIST_DIR=${DST_LIST_DIR}"
    return 1
  fi

  # Derive run number from the chosen list
  local runBase runNum
  runBase="$(basename "$chosenList")"
  runNum="${runBase#dst_calo_run2pp-}"
  runNum="${runNum%.*}"

  # Build a temporary list with the first N files of the chosen run-list
  mkdir -p "$CONDOR_LISTFILES_DIR"
  local stamp; stamp="$(date +%Y%m%d_%H%M%S)"
  local tmpList="${CONDOR_LISTFILES_DIR}/local_single_${runNum}_${nFiles}files_${stamp}.list"

  # Guard against empty/short lists; head -n handles N>len gracefully (just trims to available)
  head -n "$nFiles" "$chosenList" > "$tmpList"

  local lines
  lines=$(wc -l < "$tmpList" 2>/dev/null || echo 0)
  if (( lines == 0 )); then
    echo "[ERROR] Temporary local list is empty → $tmpList"
    return 1
  fi

  echo "===================================================================="
  echo "[INFO] dataOnly local – running ≤ ${LOCAL_NEVENTS} events"
  echo "       run       : ${runNum}"
  echo "       run-list  : ${chosenList}"
  echo "       tmp list  : ${tmpList}  (files=${lines})"
  echo "       output to : ${OUTDIR_DATA}"
  echo "===================================================================="

  mkdir -p "$OUTDIR_DATA"

  # Use PDC_RUNNUMBER=24 for Run-2 pp data (macro also defaults to 24 unless it detects run21/minbias)
  PDC_MODE=DATA PDC_RUNNUMBER=24 \
  bash "$(dirname "$0")/PositionDependentCorrect_Condor.sh" \
         "$runNum" "$tmpList" data LOCAL "$LOCAL_NEVENTS" 0 "" "$OUTDIR_DATA"
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

  # NEW: targeted resubmission of specific groups (no cleanup)
  # Accept either a path to a list file OR a numeric ClusterId after 'doAllResSubmit'
  local resubmitMode=false
  local resubmitFile=""
  local resubmitClusterId=""
  local expectResubmitArg=false

  # NEW: MINBIAS selector (switch input lists + output base when present)
  local useMinBias=false

  for tok in "$@"; do
      if $expectResubmitArg; then
          # One token allowed after doAllResSubmit: either a ClusterId (digits) or a path
          if [[ "$tok" =~ ^[0-9]+$ ]]; then
              resubmitClusterId="$tok"
          else
              resubmitFile="$tok"
          fi
          expectResubmitArg=false
          continue
      fi
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
          doAllResSubmit)
              resubmitMode=true
              expectResubmitArg=true
              ;;
          MINBIAS|minbias|MinBias)
              useMinBias=true     # activate MinBias inputs/outputs
              ;;
          firstRound|secondRound|thirdRound|fourthRound|testSubmit|sample)
              roundArg="$tok"
              ;;
          firstRoundCheck)
              roundArg="firstRound"   # use firstRound slice
              checkOnly=true          # but do not submit
              ;;
          ''|*[!0-9]*)               # ignore non‑numeric tokens here
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

  if $useMinBias; then
      SIM_LIST_DIR="/sphenix/u/patsfan753/scratch/PDCrun24pp/simListFiles/run21_type30_MB_nopileup"
      SIM_DST_LIST="${SIM_LIST_DIR}/DST_CALO_CLUSTER.list"
      SIM_HITS_LIST="${SIM_LIST_DIR}/G4Hits.list"
      if [[ ! -f "$SIM_HITS_LIST" ]]; then
          alt="${SIM_LIST_DIR}/G4HITS.list"
          [[ -f "$alt" ]] && SIM_HITS_LIST="$alt"
      fi
      OUTDIR_SIM="/sphenix/tg/tg01/bulk/jbennett/PDC/SimOutMinBias"
      TRACK_FILE="${TRACK_DIR}/sim_pairs_remaining_minbias.txt"
      CONDOR_LISTFILES_DIR="${CONDOR_LISTFILES_DIR}_minbias"
      echo "[MODE] MINBIAS active → SIM_LIST_DIR=${SIM_LIST_DIR}"
      echo "[MODE] MINBIAS active → OUTDIR_SIM=${OUTDIR_SIM}"
      echo "[MODE] MINBIAS active → TRACK_FILE=${TRACK_FILE}"
      echo "[MODE] MINBIAS active → CONDOR_LISTFILES_DIR=${CONDOR_LISTFILES_DIR}"
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
  if [[ "$execMode" == "condor" && "$resubmitMode" != true ]]; then
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
  # NEW BRANCH: targeted resubmission without any cleanup or slicing
  if [[ "$resubmitMode" == true ]]; then
      # Modes supported:
      #   doAllResSubmit </path/to/list.txt>
      #   doAllResSubmit <ClusterId>
      #
      # In the ClusterId case, auto‑discover the sim_group_DST list paths from condor_q,
      # verify strict consistency vs. the number of jobs, ensure the files exist,
      # persist to stuckJobs_*.txt, condor_rm the whole cluster, then resubmit.

      local OWNER="patsfan753"

      # Pick dataset-specific default output filename for the auto-generated list
      local stuckFile=""
      if $useMinBias; then
          stuckFile="${LOG_BASE}/stuckJobs_MB.txt"
      else
          stuckFile="${LOG_BASE}/stuckJobs_singlePhoton.txt"
      fi

      # If a ClusterId was supplied (and no file was given), discover list paths from condor_q
      if [[ -n "${resubmitClusterId:-}" && -z "${resubmitFile:-}" ]]; then
          # Basic guardrails
          if ! command -v condor_q >/dev/null 2>&1; then
              echo "[RESUBMIT][ERROR] condor_q not found in PATH"
              return 1
          fi
          if ! command -v condor_rm >/dev/null 2>&1; then
              echo "[RESUBMIT][ERROR] condor_rm not found in PATH"
              return 1
          fi
          if [[ ! "${resubmitClusterId}" =~ ^[0-9]+$ ]]; then
              echo "[RESUBMIT][ERROR] Invalid ClusterId: ${resubmitClusterId}"
              return 1
          fi

          echo "══════════════════════════════════════════════════════════════════════════"
          echo "[RESUBMIT][AUTO] Discovering list paths from condor queue"
          echo "  • Owner       : ${OWNER}"
          echo "  • ClusterId   : ${resubmitClusterId}"
          echo "  • Dataset     : $([[ $useMinBias == true ]] && echo 'MINBIAS' || echo 'single‑photon')"
          echo "  • Output list : ${stuckFile}"

          # Build the absolute-path regex for the right listfiles directory
          local listdir="condorListFiles"
          $useMinBias && listdir="condorListFiles_minbias"
          local rx="/sphenix/u/${OWNER}/scratch/PDCrun24pp/${listdir}/sim_group_DST_[^[:space:]]+[.]list"


          # Count how many procs are currently in the cluster (ground truth)
          local qcount
          qcount=$(condor_q -nobatch -const "Owner==\"${OWNER}\" && ClusterId==${resubmitClusterId}" -af:tq ProcId 2>/dev/null | wc -l | awk '{print $1}')
          qcount=${qcount:-0}

          if (( qcount <= 0 )); then
              echo "[RESUBMIT][ERROR] ClusterId=${resubmitClusterId} has 0 jobs in condor_q for Owner=${OWNER}."
              echo "                    Nothing to harvest; aborting."
              return 1
          fi
          echo "[RESUBMIT][AUTO] condor_q reports ${qcount} job(s) in the cluster."

          # Write discovered paths atomically
          local tmpfile
          if ! tmpfile="$(mktemp -p "${LOG_BASE}" .resubmit_auto_XXXXXX)"; then
              echo "[RESUBMIT][ERROR] mktemp failed (LOG_BASE=${LOG_BASE})"
              return 1
          fi

          # Pull Arguments/Args for the cluster, extract absolute sim_group_DST paths, uniq-sort
          if ! condor_q -nobatch -const "Owner==\"${OWNER}\" && ClusterId==${resubmitClusterId}" \
                        -af:tq Arguments Args 2>/dev/null \
             | awk -v rx="$rx" -F'\t' '
                {
                  a1=$1; a2=$2
                  gsub(/^"+|"+$/, "", a1); gsub(/^"+|"+$/, "", a2)
                  args = (a1!="" && a1!="(null)" && a1!="undefined") ? a1 : a2
                  if (args=="") next
                  while (match(args, rx)) {
                    print substr(args, RSTART, RLENGTH)
                    args = substr(args, RSTART+RLENGTH)
                  }
                }
              ' | sort -u > "$tmpfile"
          then
              echo "[RESUBMIT][ERROR] condor_q pipeline failed while discovering stuck jobs."
              rm -f "$tmpfile" >/dev/null 2>&1 || true
              return 1
          fi

          # Confirm we actually found paths
          local nuniq
          nuniq="$(wc -l < "$tmpfile" 2>/dev/null | awk '{print ($1==""?0:$1)}')"
          echo "[RESUBMIT][AUTO] Unique sim_group_DST list paths harvested: ${nuniq}"

          if (( nuniq != qcount )); then
              echo "[RESUBMIT][ERROR] Mismatch: extracted unique list paths (${nuniq}) != condor_q jobs (${qcount})."
              echo "                    Check that each queued job has exactly one sim_group_DST_* argument."
              echo "                    Stopping to avoid partial resubmission."
              echo "                    First few extracted paths (up to 10):"
              head -n 10 "$tmpfile" | sed 's/^/                      • /'
              rm -f "$tmpfile" >/dev/null 2>&1 || true
              return 1
          fi

          # Verify every DST list file exists; also verify companion HITS list exists
          local missing_dst=0 missing_hits=0 good=0
          local tmpverified
          if ! tmpverified="$(mktemp -p "${LOG_BASE}" .resubmit_verified_XXXXXX)"; then
            echo "[RESUBMIT][ERROR] mktemp failed (LOG_BASE=${LOG_BASE})"
            rm -f "$tmpfile" >/dev/null 2>&1 || true
            return 1
          fi

          while IFS= read -r p; do
              [[ -z "$p" ]] && continue
              if [[ ! -f "$p" ]]; then
                  echo "[RESUBMIT][ERROR] Missing DST list file on disk: $p"
                  (( ++missing_dst ))
                  continue
              fi
              local hp="${p/sim_group_DST_/sim_group_HITS_}"
              if [[ ! -f "$hp" ]]; then
                  echo "[RESUBMIT][ERROR] Missing HITS list file on disk: $hp"
                  (( ++missing_hits ))
                  continue
              fi
              echo "$p" >> "$tmpverified"
              (( ++good ))
          done < "$tmpfile"
          rm -f "$tmpfile" >/dev/null 2>&1 || true

          if (( missing_dst > 0 || missing_hits > 0 )); then
              echo "[RESUBMIT][ERROR] File verification failed:"
              echo "                    missing DST:  ${missing_dst}"
              echo "                    missing HITS: ${missing_hits}"
              rm -f "$tmpverified" >/dev/null 2>&1 || true
              return 1
          fi
          if (( good != qcount )); then
              echo "[RESUBMIT][ERROR] Verified good paths (${good}) != condor_q jobs (${qcount}). Abort."
              rm -f "$tmpverified" >/dev/null 2>&1 || true
              return 1
          fi

          # Publish verified list atomically
          mv -f "$tmpverified" "$stuckFile"
          sync
          echo "[RESUBMIT][AUTO] Verified ${good}/${qcount} paths; wrote to: $stuckFile"
          resubmitFile="$stuckFile"

          # Remove the entire cluster before resubmitting
          echo "[RESUBMIT][CLEANUP] Removing existing Condor jobs for ClusterId=${resubmitClusterId}"
          if ! condor_rm -const "Owner==\"${OWNER}\" && ClusterId==${resubmitClusterId}"; then
              echo "[RESUBMIT][WARN] condor_rm reported a non-zero status; will still verify removal."
          fi

          # Wait until the cluster disappears from the queue (bounded wait)
          local max_wait=60
          local waited=0
          while :; do
              local remain
              remain=$(condor_q -nobatch -const "Owner==\"${OWNER}\" && ClusterId==${resubmitClusterId}" -af:tq ProcId 2>/dev/null | wc -l | awk '{print $1}')
              remain=${remain:-0}
              if (( remain == 0 )); then
                  echo "[RESUBMIT][CLEANUP] condor_q shows 0 remaining jobs for ClusterId=${resubmitClusterId}."
                  break
              fi
              (( ++waited ))
              if (( waited > max_wait )); then
                  echo "[RESUBMIT][ERROR] Jobs still present after condor_rm (remain=${remain}) – aborting resubmission."
                  return 1
              fi
              echo "[RESUBMIT][CLEANUP] Waiting for queue to drain… remaining=${remain}"
              sleep 1
          done

      elif [[ -z "${resubmitFile:-}" ]]; then
          # Backward-compatible default: user didn’t pass a file nor a cluster → expect an existing file
          resubmitFile="$stuckFile"
          echo "[RESUBMIT][INFO] No ClusterId or file given; expecting an existing file:"
          echo "                 $resubmitFile"
      fi

      echo "[RESUBMIT] Using list file: ${resubmitFile}"

      # Validate provided or auto-generated list file
      [[ -f "$resubmitFile" ]] || { echo "[ERROR] doAllResSubmit: file not found → $resubmitFile"; return 1; }

      # Support BOTH Condor and Local execution for doAllResSubmit.
      if [[ "$execMode" == "local" ]]; then
          echo "[RESUBMIT][LOCAL] Sequentially running pairs listed in: $resubmitFile"
          mkdir -p "$OUTDIR_SIM"

          local idx=0
          while IFS= read -r line; do
              [[ -z "$line" || "$line" =~ ^# ]] && continue
              # Each line is the DST list path; derive HITS list path automatically
              local dstPath="$line"
              local hitPath="${dstPath/sim_group_DST_/sim_group_HITS_}"

              if [[ ! -f "$dstPath" ]]; then echo "[WARN] Missing DST list:  $dstPath (skip)"; continue; fi
              if [[ ! -f "$hitPath" ]]; then echo "[WARN] Missing HITS list: $hitPath (skip)"; continue; fi

              # Derive a stable output filename from the group indices, e.g. 945_950
              local base; base="$(basename "$dstPath")"                         # sim_group_DST_945_950.list
              local stampPart="${base#sim_group_DST_}"                          # 945_950.list
              stampPart="${stampPart%.list}"                                    # 945_950
              local outRoot="${OUTDIR_SIM}/PositionDep_sim_resubmit_${stampPart}.root"

              echo "------------------------------------------------------------------"
              echo "[RUN][LOCAL] idx=${idx}  DST=$(basename "$dstPath")"
              echo "             HITS=$(basename "$hitPath")"
              echo "             OUT =$(basename "$outRoot")"
              echo "------------------------------------------------------------------"

              rm -f "$outRoot" 2>/dev/null || true

              # Set runnumber for macro: 21 for MinBias, else 24
              if $useMinBias; then
                  export PDC_RUNNUMBER=21
              else
                  export PDC_RUNNUMBER=24
              fi
              echo "[RESUBMIT][LOCAL] PDC_RUNNUMBER=${PDC_RUNNUMBER}"

              # nevents=0 → process all events across all files in the lists
              PDC_RUNNUMBER=${PDC_RUNNUMBER} root -b -q -l "${MACRO_PATH}(0, \"${dstPath}\", \"${hitPath}\", \"${outRoot}\")"

              (( ++idx ))
          done < "$resubmitFile"

          (( idx > 0 )) || { echo "[ERROR] doAllResSubmit: no valid lines in $resubmitFile"; return 1; }
          echo "[RESUBMIT][LOCAL] Completed ${idx} local run(s)."
          return 0
      fi

      # Default: Condor resubmission (one job per line)
      # Ensure destination exists and propagate OUTDIR_SIM to the job
      mkdir -p "$OUTDIR_SIM"

      local submitFile="PositionDependentCorrect_sim_resubmit.sub"
      cat > "$submitFile" <<EOL
universe   = vanilla
executable = PositionDependentCorrect_Condor.sh
log        = ${LOG_DIR_LOG}/job.\$(Cluster).\$(Process).log
output     = ${LOG_DIR_OUT}/job.\$(Cluster).\$(Process).out
error      = ${LOG_DIR_ERR}/job.\$(Cluster).\$(Process).err
request_memory = 1500MB
EOL

      # Global environment for all queued procs in this submit file:
      if $useMinBias; then
            echo "environment = OUTDIR_SIM=/sphenix/tg/tg01/bulk/jbennett/PDC/SimOutMinBias;PDC_RUNNUMBER=21" >> "$submitFile"
      else
            echo "environment = PDC_RUNNUMBER=24" >> "$submitFile"
      fi

      echo "# Args: <runNum=9999> <listFileData> <sim> <Cluster> <nEvents=0> <processID> [<listFileHits>]" >> "$submitFile"

      local idx=0
      while IFS= read -r line; do
          [[ -z "$line" || "$line" =~ ^# ]] && continue
          local dstPath="$line"
          local hitPath="${dstPath/sim_group_DST_/sim_group_HITS_}"

          if [[ ! -f "$dstPath" ]]; then echo "[WARN] Missing DST list:  $dstPath (skip)"; continue; fi
          if [[ ! -f "$hitPath" ]]; then echo "[WARN] Missing HITS list: $hitPath (skip)"; continue; fi

          # Use the function's OUTDIR_SIM (MinBias path if MINBIAS is active)
          local destBasePath="$OUTDIR_SIM"
          printf 'arguments = 9999 %s sim $(Cluster) 0 %d %s %s\nqueue\n\n' \
                   "$dstPath" "$idx" "$hitPath" "$destBasePath" >> "$submitFile"
          (( ++idx ))
      done < "$resubmitFile"

      (( idx > 0 )) || { echo "[ERROR] doAllResSubmit: no valid lines in $resubmitFile"; return 1; }

      echo "══════════════════════════════════════════════════════════════════════════"
      echo "[STEP] condor_submit (resubmit) → $submitFile"
      echo "       Jobs queued from list    → $idx"
      echo "       List source              → $resubmitFile"
      echo "══════════════════════════════════════════════════════════════════════════"
      condor_submit "$submitFile"
      return 0
  fi

  # ORIGINAL PATH (unchanged below)
  [[ -f "$SIM_DST_LIST"  ]] || { echo "[ERROR] Missing DST list $SIM_DST_LIST";  return 1; }
  [[ -f "$SIM_HITS_LIST" ]] || { echo "[ERROR] Missing G4 list  $SIM_HITS_LIST"; return 1; }

  # Build robust (DST, G4) pairing by joining on the common filename stem:
  #   DST_CALO_CLUSTER_XXXX.root  <->  G4Hits_XXXX.root  (case-insensitive)
  dst_idx=$(mktemp) ; hits_idx=$(mktemp)

    awk '{
      full=$0;
      base=$0; sub(/.*\//,"",base); sub(/\.[Rr][Oo][Oo][Tt]$/,"",base);
      key=tolower(base);
      sub(/^dst_calo_cluster_/,"",key);
      print key, full
    }' "$SIM_DST_LIST" | sort -k1,1 > "$dst_idx"

    awk '{
      full=$0;
      base=$0; sub(/.*\//,"",base); sub(/\.[Rr][Oo][Oo][Tt]$/,"",base);
      key=tolower(base);
      sub(/^g4hits_/,"",key);  # matches both G4Hits_ and G4HITS_
      print key, full
    }' "$SIM_HITS_LIST" | sort -k1,1 > "$hits_idx"

  pairs_file=$(mktemp)
  join -j1 -o 1.2,2.2 "$dst_idx" "$hits_idx" > "$pairs_file" || true
  rm -f "$dst_idx" "$hits_idx"

  mapfile -t simPairs < "$pairs_file"
  rm -f "$pairs_file"

  if (( ${#simPairs[@]} == 0 )); then
      echo "[ERROR] No (DST,G4) pairs could be matched. Check lists under $SIM_LIST_DIR" >&2
      return 1
  fi

  echo "[PAIRING] Matched ${#simPairs[@]} (DST,G4) pair(s)."

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

      # Ensure staging dir exists & is writable
      if [[ ! -d "$CONDOR_LISTFILES_DIR" ]]; then
        echo "[DEBUG] Creating CONDOR_LISTFILES_DIR: $CONDOR_LISTFILES_DIR"
        mkdir -p "$CONDOR_LISTFILES_DIR" || { echo "[ERROR] mkdir failed: $CONDOR_LISTFILES_DIR"; return 1; }
      fi
      if [[ ! -w "$CONDOR_LISTFILES_DIR" ]]; then
        echo "[ERROR] Not writable: $CONDOR_LISTFILES_DIR"; return 1
      fi

      # Build temp list files for this local bundle
      local stamp; stamp="$(date +%Y%m%d_%H%M%S)"
      local tag=""
      if $useMinBias; then tag="_MinBias"; fi
      local dstList="${CONDOR_LISTFILES_DIR}/local_bundle_DST${tag}_${offset}_${endIndex}_${stamp}.list"
      local hitsList="${CONDOR_LISTFILES_DIR}/local_bundle_HITS${tag}_${offset}_${endIndex}_${stamp}.list"

      rm -f "$dstList" "$hitsList" 2>/dev/null || true

      echo "[DEBUG] Writing list files:"
      echo "       DST  → $dstList"
      echo "       HITS → $hitsList"

      local idx=0
      for (( i=offset; i<endIndex; i++ )); do
          idx=$(( idx + 1 ))
          IFS=' ' read -r dstFile hitFile <<< "${simPairs[$i]}"

          # Avoid slow Lustre metadata checks; just enqueue the paths.
          # The ROOT macro will report any missing/corrupt files.
          if [[ -z "$dstFile" || -z "$hitFile" ]]; then
              echo "[ERROR] Malformed pair at index $i: '${simPairs[$i]}'"
              return 1
          fi

          printf '%s\n' "$dstFile"  >> "$dstList"  || { echo "[ERROR] write failed: $dstList"; return 1; }
          printf '%s\n' "$hitFile"  >> "$hitsList" || { echo "[ERROR] write failed: $hitsList"; return 1; }
      done

      sync

      # Quick sanity with visible counts + heads
      local nd nh
      nd=$(wc -l < "$dstList"  2>/dev/null || echo 0)
      nh=$(wc -l < "$hitsList" 2>/dev/null || echo 0)
      echo "[DEBUG] Sanity counts: DST=$nd  HITS=$nh"
      if (( nd == 0 || nh == 0 || nd != nh )); then
          echo "[ERROR] Local bundle lists invalid (DST=$nd, HITS=$nh) – abort local run"
          echo "[DEBUG] Head of $dstList:";  head -n 3 "$dstList" || true
          echo "[DEBUG] Head of $hitsList:"; head -n 3 "$hitsList" || true
          return 1
      fi

      # Choose output location — ALWAYS under OUTDIR_SIM (respects MINBIAS)
      mkdir -p "$OUTDIR_SIM" || { echo "[ERROR] mkdir failed: $OUTDIR_SIM"; return 1; }

      # If this round implies a wipe (firstRound/sample or doAll), clean the target tree
      if [[ "$wipe_outputs" == true ]]; then
          echo "[CLEAN][LOCAL] Removing previous ROOT outputs from $OUTDIR_SIM"
          find "$OUTDIR_SIM" -mindepth 1 -delete 2>/dev/null || true
      fi

      local outRoot
      if $FAST_SUBMIT; then
          # Single-file ALL output lives in the MINBIAS-specific OUTDIR_SIM
          if $useMinBias; then
              outRoot="${OUTDIR_SIM}/PositionDep_sim_ALL_MinBias.root"
          else
              outRoot="${OUTDIR_SIM}/PositionDep_sim_ALL.root"
          fi
      else
          outRoot="${OUTDIR_SIM}/PositionDep_sim_bundle_${offset}_${endIndex}_${stamp}.root"
      fi
      rm -f "$outRoot" 2>/dev/null || true

      echo "------------------------------------------------------------------"
      echo "[RUN][LOCAL]  pairs=${nd}"
      echo "  ROOT macro : $MACRO_PATH"
      echo "  DST list   : $dstList"
      echo "  HITS list  : $hitsList"
      echo "  OUT file   : $outRoot"
      echo "------------------------------------------------------------------"

      # Set runnumber for macro: 21 for MinBias, else 24
      if $useMinBias; then
          export PDC_RUNNUMBER=21
      else
          export PDC_RUNNUMBER=24
      fi
      echo "[RUN][LOCAL] PDC_RUNNUMBER=${PDC_RUNNUMBER}"

      # Run through the same wrapper the Condor path uses, so env + invocation are identical
      local WRAPPER_SCRIPT="$(dirname "$0")/PositionDependentCorrect_Condor.sh"
      if [[ ! -x "$WRAPPER_SCRIPT" ]]; then
          echo "[ERROR] Wrapper not found or not executable: $WRAPPER_SCRIPT"
          return 1
      fi

      # Arguments: <runNum> <listFileData> <sim> <Cluster> <nEvents> <processID> [<listFileHits>] <destBase>
      # Use a sentinel run number for local (e.g., 9999). nEvents=0 → process all events.
      PDC_RUNNUMBER=${PDC_RUNNUMBER} bash "$WRAPPER_SCRIPT" 9999 "${dstList}" sim LOCAL 0 0 "${hitsList}" "${OUTDIR_SIM}"
      local rc=$?
      if (( rc != 0 )); then
            echo "[ERROR] PositionDependentCorrect_Condor.sh returned non-zero exit code: $rc"
            return $rc
      fi


      # Report the actual output created by the wrapper:
      # PositionDependentCorrect_Condor.sh writes to: ${OUTDIR_SIM}/9999/PositionDep_sim_<basename(dstList .list)>.root
      local produced="${OUTDIR_SIM}/9999/PositionDep_sim_$(basename "${dstList%.*}").root"
      if $FAST_SUBMIT; then
            echo "[INFO] LOCAL doAll complete → ${produced}"
      else
            echo "[INFO] LOCAL bundle complete → ${produced}"
      fi
      return 0
  fi


  # ─────────────── 6) CONDOR SUBMISSION (group N pairs per job) ───────────
  local submitFile="PositionDependentCorrect_sim.sub"

  local ENV_LINE="environment = "
  if $useMinBias; then
        ENV_LINE+="OUTDIR_SIM=/sphenix/tg/tg01/bulk/jbennett/PDC/SimOutMinBias;PDC_RUNNUMBER=21"
  else
        ENV_LINE+="PDC_RUNNUMBER=24"
  fi


  cat > "$submitFile" <<EOL
universe   = vanilla
executable = PositionDependentCorrect_Condor.sh
log        = ${LOG_DIR_LOG}/job.\$(Cluster).\$(Process).log
output     = ${LOG_DIR_OUT}/job.\$(Cluster).\$(Process).out
error      = ${LOG_DIR_ERR}/job.\$(Cluster).\$(Process).err
request_memory = 1500MB
${ENV_LINE}
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

        # Use the function's OUTDIR_SIM (MinBias path if MINBIAS is active)
        local destBasePath="$OUTDIR_SIM"
        printf 'arguments = 9999 %s sim $(Cluster) 0 %d %s %s\nqueue\n\n' \
               "$dstGroup" "$groupIdx" "$hitGroup" "$destBasePath" >> "$submitFile"

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
                # Forward an optional count or list path in $3 (e.g. '1' → first file)
                local_run_data "$3" ;;

            condorTest)
                submit_data_condor "condorTest" "" ;;

            condorFirstTen)
                submit_data_condor "condor" "firstTen" ;;

             condor)
                # -----------------------------------------------------------------
                #   • dataOnly condor                 → full scan (legacy)
                #   • dataOnly condor round N         → submit entire segment N
                #   • dataOnly condor startFrom RUN   → submit from RUN → end of its segment
                #     (NEVER falls back to full scan; errors out if not found)
                # -----------------------------------------------------------------

                if [[ "$3" == "startFrom" && "$4" =~ ^[0-9]+$ ]]; then
                    # Normalize to 8 digits to match goldenRuns_segment*.txt format
                    startRun=$(printf "%08d" "$((10#$4))")
                    mkdir -p "$CONDOR_LISTFILES_DIR"

                    # Collect segment files
                    shopt -s nullglob
                    segs=( "${GOLDEN_SEGMENT_PREFIX}"*.txt )
                    shopt -u nullglob
                    if (( ${#segs[@]} == 0 )); then
                        echo "[ERROR] startFrom: no segment files found at ${GOLDEN_SEGMENT_PREFIX}*.txt"
                        exit 1
                    fi

                    # Locate the segment containing startRun (exact line match)
                    segMatch=""
                    for s in "${segs[@]}"; do
                        if grep -qx "$startRun" "$s"; then
                            segMatch="$s"
                            break
                        fi
                    done

                    if [[ -z "$segMatch" ]]; then
                        echo "[ERROR] startFrom: run ${startRun} not found in any segment file."
                        echo "[HINT] Your segments are:"
                        for s in "${segs[@]}"; do
                            echo "  - $(basename "$s"): first/last 5 lines →"
                            awk 'NR<=5{print "    ",$0} END{print "    ..."; for(i=NR-4;i<=NR;i++) if(i>5&&i>0) print "    ",$i}' "$s"
                        done
                        exit 1
                    fi

                    # Build tail list from startRun → EOF of the matched segment (strip comments/blanks)
                    tmpRunList="${CONDOR_LISTFILES_DIR}/startFrom_${startRun}_$(date +%s).txt"
                    awk -v s="$startRun" '
                        BEGIN { emit=0 }
                        /^[[:space:]]*#/ || /^[[:space:]]*$/ { next }
                        $0==s { emit=1 }
                        emit { print $0 }
                    ' "$segMatch" > "$tmpRunList"

                    lines=$(wc -l < "$tmpRunList" 2>/dev/null || echo 0)
                    if (( lines == 0 )); then
                        echo "[ERROR] startFrom: built list is empty; segment=$(basename "$segMatch"), start=${startRun}"
                        exit 1
                    fi

                    echo "[INFO] startFrom: segment=$(basename "$segMatch"), first=${startRun}, lines=${lines}"
                    submit_data_condor "condor" "" "$tmpRunList"
                    exit 0

                elif [[ "$3" == "round" && "$4" =~ ^[0-9]+$ ]]; then
                    segFile="${GOLDEN_SEGMENT_PREFIX}${4}.txt"
                    echo "[INFO] round: submitting entire segment file $(basename "$segFile")"

                    # One-time cleanup ONLY for round 1 (data output + condor log paths)
                    if [[ "$4" -eq 1 ]]; then
                      echo "[CLEAN] Round 1 detected → clearing previous DATA outputs and Condor logs"
                      # DATA output path: remove everything inside (keep the top directory)
                      if [[ -d "$OUTDIR_DATA" ]]; then
                        echo "        OUTDIR_DATA → $OUTDIR_DATA"
                        find "$OUTDIR_DATA" -mindepth 1 -print -delete 2>/dev/null || true
                      else
                        echo "        [WARN] OUTDIR_DATA does not exist: $OUTDIR_DATA"
                      fi
                      # Condor paths: remove files only (keep directories)
                      for d in "$LOG_DIR/stdout" "$LOG_DIR/log" "$LOG_DIR/error"; do
                        if [[ -d "$d" ]]; then
                          echo "        CONDOR DIR  → $d"
                          find "$d" -type f -print -delete 2>/dev/null || true
                        else
                          echo "        [WARN] Condor dir not found: $d"
                        fi
                      done
                    fi

                    submit_data_condor "condor" "" "$segFile"
                    exit 0

                else
                    echo "[INFO] legacy: submitting full scan (no round/startFrom provided)"
                    submit_data_condor "condor" ""
                    exit 0
                fi
                ;;



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
