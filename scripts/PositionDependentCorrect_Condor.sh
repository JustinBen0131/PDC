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
destBase="$8"

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
# Priority:
#   (1) caller supplied destBase → use it
#   (2) built‑in defaults      → OUTDIR_DATA / OUTDIR_SIM
if [[ -n "$destBase" ]]; then
  outDir="${destBase}/${runNumber}"
else
  if [[ "$dataOrSim" == "sim" ]]; then
    outDir="/sphenix/tg/tg01/bulk/jbennett/PDC/SimOut/${runNumber}"
  else
    outDir="/sphenix/tg/tg01/bulk/jbennett/PDC/output/${runNumber}"
  fi
fi
mkdir -p "$outDir"
echo "[INFO] Output directory = $outDir"


#############################
# 5) Construct output filename
#############################
fileBaseName=$(basename "$chunkFile1")     # DST_CALO_run2pp_…‑00000.root
fileTag=${fileBaseName%.*}

if [[ "$dataOrSim" == "sim" ]]; then
  outFile="${outDir}/PositionDep_sim_${fileTag}.root"
else
  outFile="${outDir}/PositionDep_data_${fileTag}.root"
fi
echo "[INFO] Final output file: $outFile"

#############################
# 6) Run the macro directly
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
