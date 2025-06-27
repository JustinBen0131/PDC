#!/usr/bin/env bash
###############################################################################
# PositionDependentCorrect_Condor_submit.sh
#
# Final reference script combining:
#   - local mode (both/noSim/localSimTest) with up to 10k events
#   - Condor submission mode (data only / sim only / both)
#   - "localSimTest" picks just the first line from each sim list
#   - merges for local modes so the macro stops after 10k events
#   - chunking for condor modes to run all events
#   - "submitTestSimCondor" => partial Condor sim submission using only 4 lines
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
    [[ "$runMode" == "condor" && "$limitSwitch" == "firstTen" ]] && jobLimit=$MAX_JOBS_DATA

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
  mkdir -p "$LOG_DIR_LOG" "$LOG_DIR_OUT" "$LOG_DIR_ERR" "$CONDOR_LISTFILES_DIR"

  clean_logs() {
      rm -vf "$LOG_DIR_LOG"/job.* "$LOG_DIR_OUT"/job.* "$LOG_DIR_ERR"/job.* \
         2>/dev/null || true
  }

  # ─────────────── 1) FLEXIBLE CLI PARSING ────────────────────────────────
  local roundArg="firstRound"     # default slice keyword
  local execMode="condor"         # default execution mode
  local sampleN=""                # set only if user supplied a pure integer

  for tok in "$@"; do
      case "$tok" in
          local|condor) execMode="$tok" ;;
          simOnly)      execMode="condor" ;;        # alias
          firstRound|secondRound|thirdRound|fourthRound|testSubmit|sample)
                       roundArg="$tok" ;;
          ''|*[!0-9]*) ;;                           # ignore non‑numeric tokens
          *)           sampleN="$tok" ;;            # pure integer
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

  # ─────────────── 3) CLEAN‑UP PRIOR TO SUBMISSION ────────────────────────
  if [[ "$execMode" == "condor" ]]; then
      if $wipe_outputs; then
          echo "[STEP] Removing previous ROOT outputs from $OUTDIR_SIM"
          [[ -d "$OUTDIR_SIM" ]] && find "$OUTDIR_SIM" -mindepth 1 -delete
      fi
      echo "[STEP] Cleaning old Condor logs"
      clean_logs
  fi

  # ─────────────── 4) BUILD THE LIST OF (DST,G4) PAIRS ───────────────────
  [[ -f "$SIM_DST_LIST"  ]] || { echo "[ERROR] Missing DST list $SIM_DST_LIST";  return 1; }
  [[ -f "$SIM_HITS_LIST" ]] || { echo "[ERROR] Missing G4 list  $SIM_HITS_LIST"; return 1; }

  mapfile -t simPairs < <(paste -d' ' "$SIM_DST_LIST" "$SIM_HITS_LIST" | sort -k1,1V)

  local maxN=${#simPairs[@]}
  (( offset < maxN )) || { echo "[WARN] Offset beyond list size – nothing to do."; return 0; }

  local endIndex=$(( offset + chunkSize ))
  (( endIndex > maxN )) && endIndex=$maxN
  local nJobs=$(( endIndex - offset ))
  (( nJobs > 0 )) || { echo "[ERROR] Selected zero pairs"; return 1; }

  # ─────────────── 5) LOCAL MODE  (single‑core loop) ──────────────────────
  if [[ "$execMode" == "local" ]]; then
      echo "[STEP] LOCAL mode – processing $nJobs pair(s)"
      local idx=0
      for (( i=offset; i<endIndex; i++ )); do
          idx=$(( idx + 1 ))
          IFS=' ' read -r dstFile hitFile <<< "${simPairs[$i]}"

          [[ -f "$dstFile" ]] || { echo "[WARN] Missing DST  $dstFile – skip"; continue; }
          [[ -f "$hitFile" ]] || { echo "[WARN] Missing HITS $hitFile – skip"; continue; }

          local outRoot="${OUTDIR_SIM}/PositionDep_sim_pair$(printf "%06d" "$i").root"
          echo "------------------------------------------------------------------"
          echo "[RUN]  $idx / $nJobs"
          echo "   DST : $(basename "$dstFile")"
          echo "   HITS: $(basename "$hitFile")"
          echo "   OUT : $(basename "$outRoot")"
          echo "------------------------------------------------------------------"

          root -b -q -l "${MACRO_PATH}(0, \"${dstFile}\", \"${hitFile}\", \"${outRoot}\")"
      done
      echo "[INFO] LOCAL slice complete → $OUTDIR_SIM"
      return 0
  fi

  # ─────────────── 6) CONDOR SUBMISSION ──────────────────────────────────
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

  for (( i=offset; i<endIndex; i++ )); do
      IFS=' ' read -r dstFile hitFile <<< "${simPairs[$i]}"
      printf 'arguments = 9999 %s sim $(Cluster) 0 %d %s\nqueue\n\n' \
             "$dstFile" "$i" "$hitFile" >> "$submitFile"
  done

  echo "[STEP] condor_submit → $submitFile  (jobs=$nJobs)"
  condor_submit "$submitFile" &&
      echo "[INFO] condor_submit succeeded." ||
      { echo "[ERROR] condor_submit failed"; return 1; }
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
    submit_sim_condor "$WHAT" "$3"      # pass optional count
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
