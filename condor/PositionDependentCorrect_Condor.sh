#!/usr/bin/env bash
###############################################################################
# PositionDependentCorrect_Condor.sh
#
# This is the per‐job script Condor runs. It:
#  1) Sets up the environment
#  2) Reads arguments: runNumber, chunkFile, dataOrSim, clusterID, [optional chunkFile2 for sim].
#  3) Copies chunkFile(s) into local node or uses them directly.
#  4) Runs the Fun4All macro with appropriate input lists.
#  5) Places output in the designated directory.
###############################################################################

# 1) Setup environment
echo "[INFO] PositionDependentCorrect_Condor.sh: Starting..."
export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}
# Example local installation
MYINSTALL="/sphenix/user/${USER}/install"
source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh "$MYINSTALL"

# Parse arguments
runNumber="$1"
chunkFile1="$2"
dataOrSim="$3"
clusterID="$4"
nEvents="$5"         # e.g. "0"
chunkIndex="$6"      # condor's Process ID => chunk number
chunkFile2="$7"      # G4Hits chunk file


echo "[INFO]  runNumber   = $runNumber"
echo "[INFO]  chunkFile1  = $chunkFile1"
echo "[INFO]  dataOrSim   = $dataOrSim"
echo "[INFO]  clusterID   = $clusterID"
echo "[INFO]  nEvents     = $nEvents"
echo "[INFO]  chunkIndex  = $chunkIndex"
echo "[INFO]  chunkFile2  = $chunkFile2"

# 2) Define macro path (same as in the submit script)
MACRO_PATH="/sphenix/u/patsfan753/scratch/PDCrun24pp/macros/Fun4ALL_PDC.C"

# 3) Decide output directory
OUTDIR_DATA="/sphenix/tg/tg01/bulk/jbennett/PDC/output"
OUTDIR_SIM="/sphenix/tg/tg01/bulk/jbennett/PDC/SimOut"

if [ "$dataOrSim" == "data" ]; then
  outDir="$OUTDIR_DATA"
elif [ "$dataOrSim" == "sim" ]; then
  outDir="$OUTDIR_SIM"
else
  echo "[WARNING] dataOrSim not recognized => defaulting to data output."
  outDir="$OUTDIR_DATA"
fi

mkdir -p "$outDir"

# 4) Move chunkFile(s) local
if [ ! -f "$chunkFile1" ]; then
  echo "[ERROR] chunkFile1 missing: $chunkFile1"
  exit 1
fi
cp "$chunkFile1" inputdata.txt

# For sim, we also have chunkFile2
if [ "$dataOrSim" == "sim" ]; then
  if [ ! -f "$chunkFile2" ]; then
    echo "[ERROR] chunkFile2 missing for sim: $chunkFile2"
    exit 1
  fi
  cp "$chunkFile2" inputdatahits.txt
else
  # data => no second file
  echo "" > inputdatahits.txt
fi

# 5) Construct an output name
# ------------------------------------
firstRoot=$(head -n 1 inputdata.txt)
fileBaseName=$(basename "$firstRoot")
fileTag="${fileBaseName%.*}"

if [ "$dataOrSim" = "data" ]; then
  outFile="${OUTDIR_DATA}/PositionDep_data_chunk${chunkIndex}.root"
else
  outFile="${OUTDIR_SIM}/PositionDep_sim_chunk${chunkIndex}.root"
fi
mkdir -p "$(dirname "$outFile")"
echo "[INFO] Output file will be: $outFile"
echo "[INFO] Running macro with input data/hits..."

root -b -l -q "${MACRO_PATH}+(${nEvents}, \"inputdata.txt\", \"inputdatahits.txt\", \"${outFile}\")"
