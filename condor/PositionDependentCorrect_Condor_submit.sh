#!/usr/bin/env bash
###############################################################################
# PositionDependentCorrect_Condor_submit.sh
#
# Usage:
#   1) ./PositionDependentCorrect_Condor_submit.sh local both
#      => Run locally on both data & sim (just one file from each).
#
#   2) ./PositionDependentCorrect_Condor_submit.sh local noSim
#      => Run locally on data only (just one file).
#
#   3) ./PositionDependentCorrect_Condor_submit.sh noSim
#      => Submit Condor jobs for data only.
#
#   4) ./PositionDependentCorrect_Condor_submit.sh simOnly
#      => Submit Condor jobs for simulation only (gamma lists).
#
#   5) ./PositionDependentCorrect_Condor_submit.sh
#      => Default: submit Condor jobs for both data & sim in parallel.
#
#   6) ./PositionDependentCorrect_Condor_submit.sh localSimTest
#      => **NEW**: Run a local test on the first pair of simulation files,
#                  *without* submitting anything to Condor.
###############################################################################

#############################
# 0) PATHS, MACRO NAMES
#############################

# Where the macro is located:
MACRO_PATH="/sphenix/u/patsfan753/scratch/PDCrun24pp/macros/Fun4All_PDC.C"

# Output directories
OUTDIR_DATA="/sphenix/tg/tg01/bulk/jbennett/PDC/output"
OUTDIR_SIM="/sphenix/tg/tg01/bulk/jbennett/PDC/SimOut"

# Where data DST lists are found:
DST_LIST_DIR="/sphenix/u/patsfan753/scratch/PDCrun24pp/dst_list"
# Typically named like "dst_calo_run2pp-00047289.list"

# Where sim DST + G4Hits lists are found:
# ------------------------------------------------------------------
# Instead of always using "mb_pythia_...", we can override with single‐photon:
# For example:
#   export SIM_DST_LIST="listFiles/gamma_dst_calo_cluster.list"
#   export SIM_HITS_LIST="listFiles/gamma_g4hits.list"
# if not defined, fallback to defaults
# ------------------------------------------------------------------
: "${SIM_LIST_DIR:="/sphenix/u/patsfan753/scratch/PDCrun24pp/simLists"}"
: "${SIM_DST_LIST:="${SIM_LIST_DIR}/mb_pythia_dst_calo_cluster.list"}"
: "${SIM_HITS_LIST:="${SIM_LIST_DIR}/mb_pythia_g4hits.list"}"

# Condor logs (stdout/err/log)
LOG_DIR="/sphenix/user/${USER}/PositionDependentCorrect"
CONDOR_LISTFILES_DIR="${LOG_DIR}/condorListFiles"

# how many files per chunk for condor
FILES_PER_CHUNK=10

#############################
# 1) PARSE ARGUMENTS
#############################
MODE="$1"     # "local" / "noSim" / "simOnly" / "localSimTest" / or empty
WHAT="$2"     # "both" / "noSim" / or empty

#############################
# 2) UTILITY: LOCAL TESTS
#############################

local_test_data() {
  echo "[INFO] Running local test on data..."
  local runList
  runList="$(ls -1 "${DST_LIST_DIR}"/dst_calo_run2pp-*.list 2>/dev/null | head -n 1)"
  if [ -z "$runList" ]; then
    echo "[ERROR] No data list found in ${DST_LIST_DIR}."
    exit 1
  fi
  local firstRoot
  firstRoot="$(head -n 1 "$runList")"
  echo "[INFO] Using data file: $firstRoot"
  echo "$firstRoot" > /tmp/inputdata.txt
  echo "" > /tmp/inputdatahits.txt

  # run macro locally, example with "nevents=10"
  root -b -q -l "${MACRO_PATH}(10, \"/tmp/inputdata.txt\", \"/tmp/inputdatahits.txt\")"
}

local_test_sim() {
  echo "[INFO] Running local test on sim..."
  if [ ! -f "${SIM_DST_LIST}" ]; then
    echo "[ERROR] SIM_DST_LIST not found: ${SIM_DST_LIST}"
    exit 1
  fi
  if [ ! -f "${SIM_HITS_LIST}" ]; then
    echo "[ERROR] SIM_HITS_LIST not found: ${SIM_HITS_LIST}"
    exit 1
  fi
  local firstSimRoot
  local firstSimHit
  firstSimRoot="$(head -n 1 "${SIM_DST_LIST}")"
  firstSimHit="$(head -n 1 "${SIM_HITS_LIST}")"

  echo "[INFO] Using sim files: $firstSimRoot , $firstSimHit"
  echo "$firstSimRoot" > /tmp/inputdata.txt
  echo "$firstSimHit"  > /tmp/inputdatahits.txt

  root -b -q -l "${MACRO_PATH}(10, \"/tmp/inputdata.txt\", \"/tmp/inputdatahits.txt\")"
}

#############################
# 3) UTILITY: SUBMIT CONDOR (DATA)
#############################
submit_data_condor() {
  echo "[INFO] Submitting Condor for data runs..."
  mkdir -p "$LOG_DIR"/{stdout,error}
  mkdir -p "$CONDOR_LISTFILES_DIR"

  local listFiles
  # example: "dst_calo_run2pp-00047289.list"
  listFiles=($(ls -1 "${DST_LIST_DIR}"/dst_calo_run2pp-*.list 2>/dev/null))
  if [ ${#listFiles[@]} -eq 0 ]; then
    echo "[ERROR] No data .list files found in $DST_LIST_DIR"
    return
  fi

  for dlist in "${listFiles[@]}"; do
    local runBase
    runBase="$(basename "$dlist")"
    local runNum="${runBase#dst_calo_run2pp-}"
    runNum="${runNum%.*}"  # remove .list extension

    echo "[INFO] Found data run list: $dlist -> runNumber=$runNum"

    # read lines
    mapfile -t allfiles < "$dlist"
    local total=${#allfiles[@]}
    if [ $total -eq 0 ]; then
      echo "[WARNING] $dlist is empty. skipping."
      continue
    fi

    local chunkIndex=0
    local i=0
    local chunkListOfLists="${CONDOR_LISTFILES_DIR}/data_run_${runNum}_chunks.txt"
    > "$chunkListOfLists"

    while [ $i -lt $total ]; do
      chunkIndex=$(( chunkIndex + 1 ))
      local chunkFile="${CONDOR_LISTFILES_DIR}/data_run_${runNum}_chunk${chunkIndex}.list"
      > "$chunkFile"
      for((cc=0; cc<FILES_PER_CHUNK && i<total; cc++)); do
        echo "${allfiles[$i]}" >> "$chunkFile"
        i=$(( i + 1 ))
      done
      echo "$chunkFile" >> "$chunkListOfLists"
    done

    # create the .sub file
    cat > PositionDependentCorrect_data_${runNum}.sub <<EOL
universe                = vanilla
executable              = PositionDependentCorrect_Condor.sh
arguments               = ${runNum} \$(filename) data \$(Cluster)
log                     = ${LOG_DIR}/job.\$(Cluster).\$(Process).log
output                  = ${LOG_DIR}/stdout/job.\$(Cluster).\$(Process).out
error                   = ${LOG_DIR}/error/job.\$(Cluster).\$(Process).err
request_memory          = 1000MB
queue filename from ${chunkListOfLists}
EOL

    condor_submit PositionDependentCorrect_data_${runNum}.sub
  done
}

#############################
# 4) UTILITY: SUBMIT CONDOR (SIM)
#############################
submit_sim_condor() {
  echo "[INFO] Submitting Condor for SIM..."

  mkdir -p "$LOG_DIR"/{stdout,error}
  mkdir -p "$CONDOR_LISTFILES_DIR"

  # Make sure these are set for single gamma or whichever sample you want:
  echo "[INFO] Using SIM_DST_LIST=${SIM_DST_LIST}"
  echo "[INFO] Using SIM_HITS_LIST=${SIM_HITS_LIST}"

  if [ ! -f "$SIM_DST_LIST" ]; then
    echo "[ERROR] Sim DST list missing: $SIM_DST_LIST"
    return
  fi
  if [ ! -f "$SIM_HITS_LIST" ]; then
    echo "[ERROR] Sim Hits list missing: $SIM_HITS_LIST"
    return
  fi

  # chunk the sim DST
  mapfile -t simDstAll < "$SIM_DST_LIST"
  local totalDst=${#simDstAll[@]}
  if [ $totalDst -eq 0 ]; then
    echo "[WARNING] $SIM_DST_LIST is empty"
    return
  fi

  # chunk the sim hits
  mapfile -t simHitsAll < "$SIM_HITS_LIST"
  local totalHits=${#simHitsAll[@]}
  if [ $totalHits -eq 0 ]; then
    echo "[WARNING] $SIM_HITS_LIST is empty"
    return
  fi

  # assume they have same length or at least min(totalDst, totalHits)
  local maxN=$(( totalDst < totalHits ? totalDst : totalHits ))
  echo "[INFO] Found $totalDst DST lines, $totalHits G4Hit lines => pairing up to $maxN"

  local chunkIndex=0
  local i=0
  local chunkListOfLists="${CONDOR_LISTFILES_DIR}/sim_chunks.txt"
  > "$chunkListOfLists"

  while [ $i -lt $maxN ]; do
    chunkIndex=$(( chunkIndex + 1 ))
    local chunkDstFile="${CONDOR_LISTFILES_DIR}/sim_dst_chunk${chunkIndex}.list"
    local chunkHitsFile="${CONDOR_LISTFILES_DIR}/sim_hits_chunk${chunkIndex}.list"
    > "$chunkDstFile"
    > "$chunkHitsFile"

    for((cc=0; cc<FILES_PER_CHUNK && i<maxN; cc++)); do
      echo "${simDstAll[$i]}" >> "$chunkDstFile"
      echo "${simHitsAll[$i]}" >> "$chunkHitsFile"
      i=$(( i + 1 ))
    done

    # store the chunk file pair in chunkListOfLists
    echo "$chunkDstFile $chunkHitsFile" >> "$chunkListOfLists"
  done

  cat > PositionDependentCorrect_sim.sub <<EOL
universe                = vanilla
executable              = PositionDependentCorrect_Condor.sh
arguments               = 9999 \$(filenameA) sim \$(Cluster) \$(filenameB)
log                     = ${LOG_DIR}/job.\$(Cluster).\$(Process).log
output                  = ${LOG_DIR}/stdout/job.\$(Cluster).\$(Process).out
error                   = ${LOG_DIR}/error/job.\$(Cluster).\$(Process).err
request_memory          = 1000MB
queue filenameA filenameB from ${chunkListOfLists}
EOL

  condor_submit PositionDependentCorrect_sim.sub
}

#############################
# 5) MAIN LOGIC
#############################

if [ "$MODE" == "local" ]; then
  # local usage
  if [ "$WHAT" == "both" ]; then
    echo "[INFO] Local test => data + sim"
    local_test_data
    local_test_sim
    exit 0
  elif [ "$WHAT" == "noSim" ]; then
    echo "[INFO] Local test => data only"
    local_test_data
    exit 0
  else
    echo "[ERROR] Usage for local is: ./PositionDependentCorrect_Condor_submit.sh local [both|noSim]"
    exit 1
  fi

elif [ "$MODE" == "localSimTest" ]; then
  # NEW: local test for SIM ONLY (no Condor submission)
  echo "[INFO] Local simulation-only test => first pair of sim files"
  local_test_sim
  exit 0

elif [ "$MODE" == "noSim" ]; then
  # condor data only
  echo "[INFO] Submitting condor => data only"
  submit_data_condor
  exit 0

elif [ "$MODE" == "simOnly" ]; then
  # condor sim only
  echo "[INFO] Submitting condor => SIM only"
  submit_sim_condor
  exit 0

elif [ -z "$MODE" ]; then
  # no arguments => data + sim
  echo "[INFO] Submitting condor => BOTH data & sim"
  submit_data_condor
  submit_sim_condor
  exit 0

else
  echo "[ERROR] Unrecognized usage."
  echo "Try one of:"
  echo "  ./PositionDependentCorrect_Condor_submit.sh local both"
  echo "  ./PositionDependentCorrect_Condor_submit.sh local noSim"
  echo "  ./PositionDependentCorrect_Condor_submit.sh noSim"
  echo "  ./PositionDependentCorrect_Condor_submit.sh simOnly"
  echo "  ./PositionDependentCorrect_Condor_submit.sh   (submits both by default)"
  echo "  ./PositionDependentCorrect_Condor_submit.sh localSimTest  (new mode)"
  exit 1
fi
