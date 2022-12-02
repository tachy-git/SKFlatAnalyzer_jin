export SKFlat_WD=`pwd`
export SKFlat_LIB_PATH=$SKFlat_WD/lib/
mkdir -p $SKFlat_LIB_PATH
mkdir -p $SKFlat_WD/tar

export SKFlatV="Run2UltraLegacy_v3"
mkdir -p $SKFlat_WD/data/$SKFlatV
export DATA_DIR=$SKFlat_WD/data/$SKFlatV

#### use cvmfs for root ####
RELEASE="`cat /etc/redhat-release`"
if [[ $RELEASE == *"Fedora"* ]]; then
    # working in local
    if [[ $HOSTNAME == *"fedora"* ]]; then
        source ~/.conda-activate
        conda activate torch
    # Using singularity image, source conda environment
    else
        source /opt/conda/bin/activate
        conda activate torch
    fi
else
    # Using from host, temporarily use ROOT from lcg
    source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.24.06/x86_64-centos7-gcc48-opt/bin/thisroot.sh
fi

if [[ $HOSTNAME == *"ui"*".sdfarm.kr"* ]]; then
  echo "@@@@ Working on KISTI"
  export SKFlatRunlogDir="/cms/ldap_home/$USER/SKFlatRunlog"
  export SKFlatOutputDir="/cms/ldap_home/$USER/SKFlatOutput"
elif [[ $HOSTNAME == *"tamsa1"* ]]; then
  echo "@@@@ Working on tamsa1"
  export SKFlatRunlogDir="/data6/Users/$USER/SKFlatRunlog"
  export SKFlatOutputDir="/data6/Users/$USER/SKFlatOutput"
elif [[ $HOSTNAME == *"tamsa2"* ]]; then
  echo "@@@@ Working on tamsa2"
  export SKFlatRunlogDir="/data6/Users/$USER/SKFlatRunlog"
  export SKFlatOutputDir="/data6/Users/$USER/SKFlatOutput"
elif [[ $HOSTNAME == *"knu"* ]]; then
  echo "@@@@ Working on KNU"
  export SKFlatRunlogDir="/u/user/$USER/scratch/SKFlatRunlog"
  export SKFlatOutputDir="/u/user/$USER/scratch/SKFlatOutput"
else
  echo "@@@@ Working in local"
  export SKFlatRunlogDir="/home/$USER/workspace/SKFlatRunlog"
  export SKFlatOutputDir="/home/$USER/workspace/SKFlatOutput"
fi

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
