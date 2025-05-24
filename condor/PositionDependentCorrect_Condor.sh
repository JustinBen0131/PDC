#!/usr/bin/env bash
###############################################################################
# PositionDependentCorrect_Condor.sh
#
# Reduced-verbosity version focusing on why the macro might not be read properly.
# Rewritten to directly call macros/Fun4All_PDC.C without a MACRO_PATH variable.
###############################################################################

echo "====================================================================="
echo "[INFO] PositionDependentCorrect_Condor.sh: Starting..."
echo "[INFO] Host: $(hostname -f), Working dir: $(pwd)"
echo "====================================================================="

#############################
# 1) Setup environment
#############################
export USER="$(id -u -n)"
export LOGNAME="${USER}"
export HOME="/sphenix/u/${LOGNAME}"

MYINSTALL="/sphenix/user/${USER}/install"
source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh "$MYINSTALL"

#############################
# 2) Parse arguments
#############################
runNumber="$1"
chunkFile1="$2"
dataOrSim="$3"
clusterID="$4"
nEvents="$5"
chunkIndex="$6"   # offset + loop index
chunkFile2="$7"
echo "---------------------------------------------------------------------"
echo "[INFO] runNumber   = $runNumber"
echo "[INFO] chunkFile1  = $chunkFile1"
echo "[INFO] dataOrSim   = $dataOrSim"
echo "[INFO] clusterID   = $clusterID"
echo "[INFO] nEvents     = $nEvents"
echo "[INFO] chunkIndex  = $chunkIndex"
echo "[INFO] chunkFile2  = $chunkFile2"
echo "---------------------------------------------------------------------"

#############################
# 3) Check that macros/Fun4All_PDC.C exists
#############################
if [ ! -f "macros/Fun4All_PDC.C" ]; then
  echo "[ERROR] Cannot find macros/Fun4All_PDC.C in $(pwd) or the script’s directory."
  exit 1
fi

#############################
# 4) Decide output directory
#############################
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

##############################
## 5) Copy chunk files locally
##############################
#echo "[DEBUG] Checking chunkFile1 => $chunkFile1"
#if [ ! -f "$chunkFile1" ]; then
#  echo "[ERROR] chunkFile1 missing: $chunkFile1"
#  exit 1
#fi
#cp "$chunkFile1" inputdata.txt

#if [ "$dataOrSim" == "sim" ]; then
#  if [ ! -f "$chunkFile2" ]; then
#    echo "[ERROR] chunkFile2 missing for sim: $chunkFile2"
#    exit 1
#  fi
#  cp "$chunkFile2" inputdatahits.txt
#else
#  echo "" > inputdatahits.txt
#fi

#echo "[DEBUG] inputdata.txt (first few lines):"
#head -n 5 inputdata.txt
#echo "[DEBUG] inputdatahits.txt (first few lines):"
#head -n 5 inputdatahits.txt

#############################
# 6) Construct output filename
#############################
fileBaseName="$(basename "$chunkFile1")"
fileTag="${fileBaseName%.*}"

if [ "$dataOrSim" = "data" ]; then
  outFile="${outDir}/PositionDep_data_chunk${chunkIndex}.root"
else
  outFile="${outDir}/PositionDep_sim_chunk${chunkIndex}.root"
fi
mkdir -p "$(dirname "$outFile")"
echo "[INFO] Final output file: $outFile"

#############################
# 7) Run the macro directly
#############################
echo "[INFO] Invoking ROOT macro with command:"
echo "root -b -l -q \"macros/Fun4All_PDC.C(${nEvents}, \\\"${chunkFile1}\\\", \\\"${chunkFile2}\\\", \\\"${outFile}\\\")\""

root -b -l -q "macros/Fun4All_PDC.C(${nEvents}, \"${chunkFile1}\", \"${chunkFile2}\", \"${outFile}\")"
rc=$?

echo "---------------------------------------------------------------------"
echo "[INFO] ROOT return code: $rc"
if [ $rc -ne 0 ]; then
  echo "[ERROR] ROOT macro failed with return code $rc"
  exit $rc
fi


echo "[INFO] PositionDependentCorrect_Condor.sh finished OK."
exit 0
