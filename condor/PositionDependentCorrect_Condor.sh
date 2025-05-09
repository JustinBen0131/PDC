#!/usr/bin/env bash
###############################################################################
# PositionDependentCorrect_Condor.sh
#
# This is the per‐job script Condor runs. It:
#  1) Sets up the environment (SPhenix).
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
chunkFile2="$5"  # only used if dataOrSim=sim

echo "[INFO]  runNumber  = $runNumber"
echo "[INFO]  chunkFile1 = $chunkFile1"
echo "[INFO]  dataOrSim  = $dataOrSim"
echo "[INFO]  clusterID  = $clusterID"
echo "[INFO]  chunkFile2 = $chunkFile2"

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
  outFile="${OUTDIR_DATA}/PositionDep_${fileTag}.root"
else
  outFile="${OUTDIR_SIM}/PositionDep_${fileTag}.root"
fi
mkdir -p "$(dirname "$outFile")"
echo "[INFO] Output file will be: $outFile"
echo "[INFO] Running macro with input data/hits..."

root -b -l -q \
  "${MACRO_PATH}(0, \"inputdata.txt\", \"inputdatahits.txt\", \"${outFile}\")"


# 6) Move result. The macro might produce an OUTHIST*.root.
#    If the macro directly writes to e.g. Fun4All_EMCal_sp something,
#    you can rename it here to outFile or do in the macro itself.

# For example, if you know it produces: "OUTHIST_iter_..."
# Move or rename to $outFile

foundHist=$(ls OUTHIST_iter_*.root 2>/dev/null | head -n 1)
if [ -f "$foundHist" ]; then
  mv "$foundHist" "$outFile"
  echo "[INFO] Moved histogram file => $outFile"
else
  echo "[WARNING] No OUTHIST_iter_*.root found. Possibly the macro named differently?"
fi

echo "[INFO] PositionDependentCorrect_Condor.sh finished!"
