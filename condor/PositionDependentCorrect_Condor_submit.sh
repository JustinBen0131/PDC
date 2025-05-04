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
###############################################################################

#############################
# 0) CONFIGURATION
#############################

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
# ---------------------------------------------------------------------
# By default, we point to gamma_dst_calo_cluster.list & gamma_g4hits.list.
# If you want others, set these environment variables externally or modify below.
# ---------------------------------------------------------------------
: "${SIM_LIST_DIR:="/sphenix/u/bseidlitz/work/analysis/EMCal_pi0_Calib_2023/macros/condor/listFiles"}"
: "${SIM_DST_LIST:="${SIM_LIST_DIR}/gamma_dst_calo_cluster.list"}"
: "${SIM_HITS_LIST:="${SIM_LIST_DIR}/gamma_g4hits.list"}"

# How many files per Condor chunk
FILES_PER_CHUNK=10

# For local modes, limit to 10k events:
LOCAL_NEVENTS=10000

#############################
# 1) ARGUMENT PARSING
#############################
MODE="$1"  # "local", "noSim", "simOnly", "localSimTest", or empty
WHAT="$2"  # "both", "noSim", or empty

usage()
{
  echo "Usage examples:"
  echo "  $0 local both"
  echo "  $0 local noSim"
  echo "  $0 localSimTest"
  echo "  $0 noSim"
  echo "  $0 simOnly"
  echo "  $0            (submits both data & sim by default)"
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

#############################
# 2) LOCAL MODE FUNCTIONS
#############################

##
# Merge all data lists from DST_LIST_DIR => /tmp/local_data_merged.list
##
merge_files_for_local_data() {
  echo "-----------------------------------------------------"
  echo "[DEBUG] In merge_files_for_local_data()"
  echo "[DEBUG] Checking DST_LIST_DIR contents: $DST_LIST_DIR"
  ls -l "$DST_LIST_DIR"
  echo "-----------------------------------------------------"

  local merged="/tmp/local_data_merged.list"
  echo "[DEBUG] Creating/overwriting $merged"
  rm -f "$merged"
  touch "$merged"

  local -a listFiles
  # gather all that match 'dst_calo_run2pp-*.list'
  listFiles=( $(ls -1 "${DST_LIST_DIR}"/dst_calo_run2pp-*.list 2>/dev/null) )
  echo "[DEBUG] Found ${#listFiles[@]} data .list files"

  if [ ${#listFiles[@]} -eq 0 ]; then
    echo "[ERROR] No data .list files found matching dst_calo_run2pp-* in $DST_LIST_DIR"
    return 1
  fi

  echo "[INFO] Merging all data .list files => $merged"
  for oneList in "${listFiles[@]}"; do
    echo "   -> Adding from: $oneList"
    cat "$oneList" >> "$merged"
  done

  echo "[DEBUG] Combined data line count => $(wc -l < "$merged")"
  echo "$merged"
}

##
# Merge all SIM DST lines => /tmp/local_sim_dst_merged.list
# Merge all SIM hits lines => /tmp/local_sim_hits_merged.list
##
merge_files_for_local_sim() {
  echo "-----------------------------------------------------"
  echo "[DEBUG] In merge_files_for_local_sim()"
  echo "[DEBUG] Checking SIM_LIST_DIR: $SIM_LIST_DIR"
  ls -l "$SIM_LIST_DIR"
  echo "-----------------------------------------------------"

  local mergedDst="/tmp/local_sim_dst_merged.list"
  local mergedHits="/tmp/local_sim_hits_merged.list"

  echo "[DEBUG] Creating/overwriting $mergedDst and $mergedHits"
  rm -f "$mergedDst" "$mergedHits"
  touch "$mergedDst" "$mergedHits"

  # check if $SIM_DST_LIST exists
  if [ ! -f "$SIM_DST_LIST" ]; then
    echo "[ERROR] SIM_DST_LIST not found: $SIM_DST_LIST"
    return 1
  fi

  # check if $SIM_HITS_LIST exists
  if [ ! -f "$SIM_HITS_LIST" ]; then
    echo "[ERROR] SIM_HITS_LIST not found: $SIM_HITS_LIST"
    return 1
  fi

  echo "[INFO] Merging sim DST => $mergedDst  (from $SIM_DST_LIST)"
  cat "$SIM_DST_LIST" >> "$mergedDst"

  echo "[INFO] Merging sim G4Hits => $mergedHits (from $SIM_HITS_LIST)"
  cat "$SIM_HITS_LIST" >> "$mergedHits"

  echo "[DEBUG] Final line count => DST: $(wc -l < "$mergedDst"), HITS: $(wc -l < "$mergedHits")"
  # Return "dstFile hitsFile"
  echo "$mergedDst $mergedHits"
}

##
# Local run for data => merges everything, runs 10k events
##
local_run_data() {
  echo "######################################################################"
  echo "[INFO] local_run_data => up to $LOCAL_NEVENTS events"

  local dataMerged
  dataMerged="$(merge_files_for_local_data)"
  if [ -z "$dataMerged" ] || [ ! -s "$dataMerged" ]; then
    echo "[ERROR] Merged data .list is empty => cannot proceed"
    exit 1
  fi

  local emptyFile="/tmp/local_empty_hits.list"
  rm -f "$emptyFile"
  touch "$emptyFile"
  echo "[DEBUG] Using empty hits file => $emptyFile"

  echo "[INFO] Now calling root for data => will stop after $LOCAL_NEVENTS events"
  root -b -q -l "${MACRO_PATH}(${LOCAL_NEVENTS}, \"${dataMerged}\", \"${emptyFile}\")"
  local rc=$?
  if [ $rc -ne 0 ]; then
    echo "[ERROR] local_run_data ended with code=$rc"
    exit $rc
  fi
  echo "[INFO] local_run_data completed OK."
}

##
# Local run for simulation => merges everything, runs 10k events
##
local_run_sim() {
  echo "######################################################################"
  echo "[INFO] local_run_sim => up to $LOCAL_NEVENTS events"

  local mergedSim
  mergedSim="$(merge_files_for_local_sim)"
  if [ -z "$mergedSim" ]; then
    echo "[ERROR] Merging sim DST/HITS returned empty => cannot proceed"
    exit 1
  fi

  local simDstFile simHitsFile
  simDstFile="$(echo "$mergedSim" | awk '{print $1}')"
  simHitsFile="$(echo "$mergedSim" | awk '{print $2}')"

  if [ ! -s "$simDstFile" ] || [ ! -s "$simHitsFile" ]; then
    echo "[ERROR] Merged sim DST/HITS list is empty => cannot proceed"
    exit 1
  fi

  echo "[INFO] Now calling root for sim => will stop after $LOCAL_NEVENTS events"
  root -b -q -l "${MACRO_PATH}(${LOCAL_NEVENTS}, \"${simDstFile}\", \"${simHitsFile}\")"
  local rc=$?
  if [ $rc -ne 0 ]; then
    echo "[ERROR] local_run_sim ended with code=$rc"
    exit $rc
  fi
  echo "[INFO] local_run_sim completed OK."
}

##
# localSimTest => read only the *first line* from each sim list
#                 process up to 10k events
##
local_sim_test_first_pair() {
  echo "######################################################################"
  echo "[INFO] localSimTest => *only* the first line from each sim list"

  # check if sim lists exist
  if [ ! -f "$SIM_DST_LIST" ]; then
    echo "[ERROR] SIM_DST_LIST not found => $SIM_DST_LIST"
    exit 1
  fi
  if [ ! -f "$SIM_HITS_LIST" ]; then
    echo "[ERROR] SIM_HITS_LIST not found => $SIM_HITS_LIST"
    exit 1
  fi

  mapfile -t simDstAll < "$SIM_DST_LIST"
  mapfile -t simHitsAll < "$SIM_HITS_LIST"

  local totalDst=${#simDstAll[@]}
  local totalHits=${#simHitsAll[@]}
  local maxN=$(( totalDst < totalHits ? totalDst : totalHits ))

  if [ "$maxN" -eq 0 ]; then
    echo "[ERROR] No lines found in either sim DST list or hits list => cannot proceed"
    exit 1
  fi

  # just pick the first line from each
  local dstFile="${simDstAll[0]}"
  local hitFile="${simHitsAll[0]}"
  echo "[INFO] Using first pair:"
  echo "  DST => $dstFile"
  echo "  HIT => $hitFile"

  # write them to /tmp
  local tmpDst="/tmp/localSimTest_dst_firstpair.list"
  local tmpHits="/tmp/localSimTest_hits_firstpair.list"
  echo "$dstFile"  > "$tmpDst"
  echo "$hitFile" > "$tmpHits"

  echo "[INFO] Will run up to $LOCAL_NEVENTS events"
  root -b -q -l "${MACRO_PATH}(${LOCAL_NEVENTS}, \"${tmpDst}\", \"${tmpHits}\")"
  local rc=$?
  if [ $rc -ne 0 ]; then
    echo "[ERROR] localSimTest ended with code=$rc"
    exit $rc
  fi
  echo "[INFO] localSimTest (first pair) completed OK."
}

#############################
# 3) CONDOR SUBMISSION FOR ALL EVENTS
#############################

submit_data_condor() {
  echo "######################################################################"
  echo "[INFO] Submitting Condor jobs => data"
  mkdir -p "$LOG_DIR"/{stdout,error} "$CONDOR_LISTFILES_DIR"

  # find data .list files => e.g. "dst_calo_run2pp-00047289.list"
  local listFiles
  listFiles=( $(ls -1 "${DST_LIST_DIR}"/dst_calo_run2pp-*.list 2>/dev/null) )
  echo "[DEBUG] Found ${#listFiles[@]} data run-lists in $DST_LIST_DIR"

  if [ ${#listFiles[@]} -eq 0 ]; then
    echo "[ERROR] No data .list found => cannot submit"
    return
  fi

  for dlist in "${listFiles[@]}"; do
    local runBase
    runBase="$(basename "$dlist")"
    local runNum="${runBase#dst_calo_run2pp-}"   # e.g. 00047289.list
    runNum="${runNum%.*}"                       # e.g. 00047289

    mapfile -t allfiles < "$dlist"
    local total=${#allfiles[@]}
    if [ "$total" -eq 0 ]; then
      echo "[WARNING] $dlist is empty => skipping"
      continue
    fi

    echo "[INFO] Data run=$runNum => $total lines"

    local chunkIndex=0
    local i=0
    local chunkListOfLists="${CONDOR_LISTFILES_DIR}/data_run_${runNum}_chunks.txt"
    > "$chunkListOfLists"

    while [ $i -lt $total ]; do
      chunkIndex=$(( chunkIndex + 1 ))
      local chunkFile="${CONDOR_LISTFILES_DIR}/data_run_${runNum}_chunk${chunkIndex}.list"
      > "$chunkFile"

      for (( cc=0; cc<FILES_PER_CHUNK && i<total; cc++ )); do
        echo "${allfiles[$i]}" >> "$chunkFile"
        i=$(( i + 1 ))
      done

      echo "$chunkFile" >> "$chunkListOfLists"
      echo "[DEBUG] Created chunk => $chunkFile"
    done

    # each job => process all events => nevents=0
    cat > PositionDependentCorrect_data_${runNum}.sub <<EOL
universe                = vanilla
executable              = PositionDependentCorrect_Condor.sh
# Args: <runNum> <listFileData> <"data"|"sim"> <clusterID> <nEvents> [<listFileHits>]
arguments               = ${runNum} \$(filename) data \$(Cluster) 0
log                     = ${LOG_DIR}/job.\$(Cluster).\$(Process).log
output                  = ${LOG_DIR}/stdout/job.\$(Cluster).\$(Process).out
error                   = ${LOG_DIR}/error/job.\$(Cluster).\$(Process).err
request_memory          = 1000MB
queue filename from ${chunkListOfLists}
EOL

    echo "[INFO] condor_submit => PositionDependentCorrect_data_${runNum}.sub"
    condor_submit PositionDependentCorrect_data_${runNum}.sub
  done
}

submit_sim_condor() {
  echo "######################################################################"
  echo "[INFO] Submitting Condor jobs => simulation"
  mkdir -p "$LOG_DIR"/{stdout,error} "$CONDOR_LISTFILES_DIR"

  if [ ! -f "$SIM_DST_LIST" ]; then
    echo "[ERROR] SIM_DST_LIST not found => $SIM_DST_LIST"
    return
  fi
  if [ ! -f "$SIM_HITS_LIST" ]; then
    echo "[ERROR] SIM_HITS_LIST not found => $SIM_HITS_LIST"
    return
  fi

  # read lines
  mapfile -t simDstAll < "$SIM_DST_LIST"
  mapfile -t simHitsAll < "$SIM_HITS_LIST"
  local totalDst=${#simDstAll[@]}
  local totalHits=${#simHitsAll[@]}
  local maxN=$(( totalDst < totalHits ? totalDst : totalHits ))

  echo "[INFO] Found $totalDst DST lines, $totalHits G4Hit lines => up to $maxN pairs"
  if [ "$maxN" -eq 0 ]; then
    echo "[WARNING] zero pairs => no sim submission"
    return
  fi

  local chunkIndex=0
  local i=0
  local chunkListOfLists="${CONDOR_LISTFILES_DIR}/sim_chunks.txt"
  > "$chunkListOfLists"

  # chunk them in sets of $FILES_PER_CHUNK
  while [ $i -lt $maxN ]; do
    chunkIndex=$(( chunkIndex + 1 ))

    local chunkDst="${CONDOR_LISTFILES_DIR}/sim_dst_chunk${chunkIndex}.list"
    local chunkHits="${CONDOR_LISTFILES_DIR}/sim_hits_chunk${chunkIndex}.list"
    > "$chunkDst"
    > "$chunkHits"

    for (( cc=0; cc<FILES_PER_CHUNK && i<maxN; cc++ )); do
      echo "${simDstAll[$i]}"   >> "$chunkDst"
      echo "${simHitsAll[$i]}" >> "$chunkHits"
      i=$(( i + 1 ))
    done

    echo "$chunkDst $chunkHits" >> "$chunkListOfLists"
    echo "[DEBUG] Created chunk pair => $chunkDst + $chunkHits"
  done

  # each job => all events => nevents=0
  cat > PositionDependentCorrect_sim.sub <<EOL
universe                = vanilla
executable              = PositionDependentCorrect_Condor.sh
# Args: <runNum=9999> <listFileData> <"sim"> <clusterID> <nEvents=0> [<listFileHits>]
arguments               = 9999 \$(filenameA) sim \$(Cluster) 0 \$(filenameB)
log                     = ${LOG_DIR}/job.\$(Cluster).\$(Process).log
output                  = ${LOG_DIR}/stdout/job.\$(Cluster).\$(Process).out
error                   = ${LOG_DIR}/error/job.\$(Cluster).\$(Process).err
request_memory          = 1000MB
queue filenameA filenameB from ${chunkListOfLists}
EOL

  echo "[INFO] condor_submit => PositionDependentCorrect_sim.sub"
  condor_submit PositionDependentCorrect_sim.sub
}

#############################
# 4) MAIN LOGIC
#############################

case "$MODE" in
  local)
    # local => must specify 'both' or 'noSim'
    if [ "$WHAT" == "both" ]; then
      echo "===================================================================="
      echo "[INFO] local => data + sim, each up to $LOCAL_NEVENTS events"
      echo "===================================================================="
      local_run_data
      local_run_sim
      exit 0

    elif [ "$WHAT" == "noSim" ]; then
      echo "===================================================================="
      echo "[INFO] local => data only, up to $LOCAL_NEVENTS events"
      echo "===================================================================="
      local_run_data
      exit 0

    else
      echo "[ERROR] For 'local' mode, specify 'both' or 'noSim'"
      usage
    fi
    ;;

  localSimTest)
    echo "===================================================================="
    echo "[INFO] localSimTest => first pair only, up to $LOCAL_NEVENTS events"
    echo "===================================================================="
    local_sim_test_first_pair
    exit 0
    ;;

  noSim)
    echo "===================================================================="
    echo "[INFO] Condor => data only (ALL events)."
    echo "===================================================================="
    submit_data_condor
    exit 0
    ;;

  simOnly)
    echo "===================================================================="
    echo "[INFO] Condor => sim only (ALL events)."
    echo "===================================================================="
    submit_sim_condor
    exit 0
    ;;

  "")
    # no arguments => condor => data + sim
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
