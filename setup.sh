RELEASE="`cat /etc/redhat-release`"
if [[ $HOSTNAME == *"tamsa"* ]]; then
  echo "@@@@ Working in tamsa"
  export SKFlat_WD="/data6/Users/$USER/SKFlatAnalyzer"
  export SKFlatRunlogDir="/data6/Users/$USER/SKFlatRunlog"
  export SKFlatOutputDir="/data6/Users/$USER/SKFlatOutput"

  # root configuration
  # Singlarity image
  if [[ $RELEASE == *"Fedora"* ]]; then
    source /opt/conda/bin/activate
    conda activate torch
  # using from host
  else
    # temporarily use ROOT and python from LCG environment
    source /cvmfs/sft.cern.ch/lcg/views/LCG_102cuda/x86_64-centos7-gcc8-opt/setup.sh
  fi
else
  echo "@@@@ Working in local"
  export SKFlat_WD="/home/$USER/workspace/SKFlatAnalyzer"
  export SKFlatRunlogDir="/home/$USER/workspace/SKFlatRunlog"
  export SKFlatOutputDir="/home/$USER/workspace/SKFlatOutput"
  # root configuration
  source ~/.conda-activate
  conda activate torch
fi

export SKFlat_LIB_PATH=$SKFlat_WD/lib/
mkdir -p $SKFlat_LIB_PATH
mkdir -p $SKFlat_WD/tar

export SKFlatV="Run2UltraLegacy_v3"
mkdir -p $SKFlat_WD/data/$SKFlatV
export DATA_DIR=$SKFlat_WD/data/$SKFlatV

alias skout="cd $SKFlatOutputDir/$SKFlatV"

export MYBIN=$SKFlat_WD/bin/
export PYTHONDIR=$SKFlat_WD/python/
export PATH=${MYBIN}:${PYTHONDIR}:${PATH}

export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$SKFlat_WD/DataFormats/include/:$SKFlat_WD/AnalyzerTools/include/:$SKFlat_WD/Analyzers/include/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SKFlat_LIB_PATH

source $SKFlat_WD/bin/BashColorSets.sh

## submodules ##
#source bin/CheckSubmodules.sh

if [ "$1" = "-q" ];then
    return
fi

## Todo list ##
python python/PrintToDoLists.py
source $SKFlat_WD/tmp/ToDoLists.sh
rm $SKFlat_WD/tmp/ToDoLists.sh

CurrentGitBranch=`git branch | grep \* | cut -d ' ' -f2`
printf "> Current SKFlatAnalyzer branch : "${BRed}$CurrentGitBranch${Color_Off}"\n"
echo "-----------------------------------------------------------------"
## Log Dir ##
echo "* Your Log Directory Usage (ctrl+c to skip)"
du -sh $SKFlatRunlogDir
