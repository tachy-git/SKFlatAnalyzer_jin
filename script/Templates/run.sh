#!/bin/bash
SECTION=`printf $1`
WORKDIR=`pwd`
Trail=0

export SKFlat_WD=[SKFlat_WD]
export SKFlat_LIB_PATH=$SKFlat_WD/lib
export SKFlatV="Run2UltraLegacy_v3"
export DATA_DIR="$SKFlat_WD/data/$SKFlatV"
export MYBIN="$SKFlat_WD/bin"
export PYTHONDIR="$SKFlat_WD/python"
export PATH=${MYBIN}:${PYTHONDIR}:${PATH}
export PYTHONPATH="${PYTHONPATH}:${PYTHONDIR}"
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$SKFlat_WD/DataFormats/include/:$SKFlat_WD/AnalyzerTools/include/:$SKFlat_WD/Analyzers/include/
export LHAPDFDIR=$SKFlat_WD/external/lhapdf
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SKFlat_LIB_PATH:$LHAPDFDIR/lib

source $SKFlat_WD/bin/BashColorSets.sh

#### make sure use C locale
export LC_ALL=C

#### Don't make root history
export ROOT_HIST=0

#### User conda env for root
source /opt/conda/bin/activate
conda activate pyg

#### modifying LD_LIBRARY_PATH to use libraries in baseRunDir
export LD_LIBRARY_PATH=$(echo $LD_LIBRARY_PATH|sed 's@'$SKFlat_WD'/lib@[masterJobDir]/lib@')

while [[ "$Trial" -lt 3 ]]; do
    echo "#### running ####"
    echo "root -l -b -q [baseRunDir]/run_${SECTION}.C"
    root -l -b -q [baseRunDir]/run_${SECTION}.C 2> err.log
    EXITCODE=$?
    if [ "$EXITCODE" -eq 5 ]; then
        echo "IO error occured.. running again in 300 seconds.."
        Trial=$((Trial+=1))
        sleep 300
    else
        break
    fi
done

if [ "$EXITCODE" -ne 0 ]; then
    echo "ERROR errno=$EXITCODE" >> err.log
fi

cat err.log >&2
exit $EXITCODE
