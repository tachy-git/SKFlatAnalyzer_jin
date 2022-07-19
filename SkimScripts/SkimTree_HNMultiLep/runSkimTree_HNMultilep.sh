analyzer=SkimTree_HNMultiLep
rundir=runSkims
mcpath=${SKFlat_WD}/SkimScripts/${analyzer}/mc_lists_multilep/
datapath=${SKFlat_WD}/SkimScripts/${analyzer}/data_lists_multilep/
njobs=50
njobs_data=100
nmax=250
skim=' '
declare  -a era_list=("2016postVFP" "2016preVFP" "2017" )

if [[ $1 == "MC" ]]; then
    for i in "${era_list[@]}"
    do
        #SKFlat.py -a $analyzer  -l $mcpath/MC_${i}.txt  -n ${njobs}  --nmax ${nmax} -e ${i} &
	SKFlat.py -a $analyzer  -l $mcpath/MC.txt  -n ${njobs}  --nmax ${nmax} -e ${i} &
    done

    python $SKFlat_WD/script/BadFileChecker/runTAMSA_Skim.py -s $analyzer
fi

if [[ $1 == "DATA" ]]; then
    for i in "${era_list[@]}"
    do
        SKFlat.py -a $analyzer  -l $datapath/DATA_${i}.txt  -n ${njobs_data}  --nmax ${nmax}   -e ${i} &
        SKFlat.py -a $analyzer  -l $datapath/DATA_l_${i}.txt   -n ${njobs_data}  --nmax ${nmax}    -e ${i} &
    done

    SKFlat.py -a $analyzer  -l $datapath/DATA_2018.txt  -n ${njobs_data}  --nmax ${nmax}   -e 2018 &
    SKFlat.py -a $analyzer  -l $datapath/DATA_l_2018.txt   -n ${njobs_data}  --nmax ${nmax}    -e 2018 

    python $SKFlat_WD/script/BadFileChecker/runTAMSA_SkimData.py -s $analyzer
fi



