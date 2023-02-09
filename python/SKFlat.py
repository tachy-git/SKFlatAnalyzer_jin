#!/usr/bin/env python3
import os, sys
import shutil
import time
import argparse
import datetime
from CheckJobStatus import *
from GetXSECTable import *
from TimeTools import *
import random
import subprocess

## Arguments
parser = argparse.ArgumentParser(description='SKFlat Command')
parser.add_argument('-a', dest='Analyzer', default="")
parser.add_argument('-i', dest='InputSample', default="")
parser.add_argument('-p', dest='DataPeriod', default="ALL")
parser.add_argument('-l', dest='InputSampleList', default="")
parser.add_argument('-n', dest='NJobs', default=1, type=int)
parser.add_argument('-o', dest='Outputdir', default="")
parser.add_argument('-q', dest='Queue', default="fastq")
parser.add_argument('-e', dest='Era', default="2017",help="2016preVFP(2016a), 2016postVFP(2016b), 2017, 2018")
parser.add_argument('-y', dest='Year', default="",help="deprecated. use -e")
parser.add_argument('--skim', dest='Skim', default="", help="ex) SkimTree_Dilepton")
parser.add_argument('--python', action='store_true')
parser.add_argument('--no_exec', action='store_true')
parser.add_argument('--FastSim', action='store_true')
parser.add_argument('--userflags', dest='Userflags', default="")
parser.add_argument('--tagoutput', dest='TagOutput', default="")
parser.add_argument('--nmax', dest='NMax', default=0, type=int, help="maximum running jobs")
parser.add_argument('--reduction', dest='Reduction', default=1, type=float)
parser.add_argument('--memory', dest='Memory', default=0, type=float)
parser.add_argument('--batchname',dest='BatchName', default="")
args = parser.parse_args()

if args.Year:
    print("-y is depreciated. Using -e (Era) instead") 
    args.Era = args.Year

if args.Era == "2016a": args.Era = "2016preVFP"
if args.Era == "2016b": args.Era = "2016postVFP"

## make userflags as a list
Userflags = []
if args.Userflags != "":
    Userflags = (args.Userflags).split(',')

## Add absolute path for outputdir
if args.Outputdir and (not os.path.isabs(args.Outputdir)):
    args.Outputdir = f"{os.getcwd()}/{args.Outputdir}"

## TimeStamp
# 1) directory / file name style
timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H%M%S")
# 2) log style
jobStartTime = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
thisTime = ""

## Environment Variables
USER = os.environ["USER"]
if os.path.exists(os.environ['SKFlat_WD']+'/python/UserInfo_'+USER+'.py'):
    exec('from UserInfo_'+USER+' import *')
else:
    print("No UserInfo file")
    exit(1)
SKFlatLogEmail = UserInfo['SKFlatLogEmail']
SKFlatLogWebDir = UserInfo['SKFlatLogWebDir']
logEvery = UserInfo['LogEvery']

SKFlat_WD = os.environ['SKFlat_WD']
SKFlatV = os.environ['SKFlatV']
SAMPLE_DATA_DIR = f"{SKFlat_WD}/data/{SKFlatV}/{args.Era}/Sample"
SKFlatRunlogDir = os.environ['SKFlatRunlogDir']
SKFlatOutputDir = os.environ['SKFlatOutputDir']
SKFlat_LIB_PATH = os.environ['SKFlat_LIB_PATH']
UID = str(os.getuid())
HOSTNAME = os.environ['HOSTNAME']

## Check joblog email
sendLogToEmail = False if SKFlatLogEmail  == "" else True
sendLogToWeb   = False if SKFlatLogWebDir == "" else True

## Check hostname
isTAMSA1 = ("tamsa1" in HOSTNAME)
isTAMSA2 = ("tamsa2" in HOSTNAME)
isTAMSA = isTAMSA1 or isTAMSA2
if isTAMSA:
    if isTAMSA1: HOSTNAME = "TAMSA1"
    if isTAMSA2: HOSTNAME = "TAMSA2"
    sampleHOSTNAME = "SNU"
else:
    print("Working in local...")
    sampleHOSTNAME = "SNU"

## Are you skimming trees?
isSkimTree = "SkimTree" in args.Analyzer
if isSkimTree:
    if args.NMax == 0: args.NMax = 100  # Preventing from too heavy IO
    if args.NJobs == 1: args.NJobs = 0  # NJobs = 0 means NJobs -> NFiles
    if not isTAMSA:
        print("Skimming is only possible in SNU")
        exit(1)

## Make sample list
inputDataSampleList = ["DoubleMuon", "DoubleEG", "SingleMuon", "SingleElectron",
                       "SinglePhoton", "MuonEG", "EGamma"]
availableDataPeriods = []
if args.Era == "2016preVFP":
    availableDataPeriods = ["B_ver2","C","D","E","F"]
elif args.Era == "2016postVFP":
    availableDataPeriods = ["F","G","H"]
elif args.Era == "2017":
    availableDataPeriods = ["B","C","D","E","F"]
elif args.Era == "2018":
    availableDataPeriods = ["A", "B","C","D"]
else:
    print(f"[SKFlat.py] Wrong Era: {args.Era}")
    exit(1)
    
inputSampleList = []
stringForHash = ""
## When using txt file for input (i.e., -l option)
if args.InputSampleList:
    lines = open(args.InputSampleList)
    for line in lines:
        if line[0] == "#":
            continue
        line = line.strip("\n")
        inputSampleList.append(line)
        stringForHash += line
else:
    if args.InputSample in inputDataSampleList:
        if args.DataPeriod == "ALL":
            for period in availableDataPeriods:
                inputSampleList.append(f"{args.InputSample}:{period}")
                stringForHash += f"{args.InputSample}:{period}"
        elif args.DataPeriod in availableDataPeriods:
            inputSampleList.append(f"{args.InputSample}:{args.DataPeriod}")
            stringForHash += f"{args.InputSample}:{args.DataPeriod}"
        else:
            pass
    else:
        inputSampleList.append(args.InputSample)
        stringForHash += args.InputSample
fileRangesForEachSample = []

## add flags to hash
for flag in Userflags:
    stringForHash += flag

## Get random number for webdir
random.seed(hash(stringForHash+timestamp+args.Era))
randomNumber = int(random.random()*1000000)
webdirname = f"{timestamp}_{randomNumber}"
webdirpathbase = f"{SKFlatRunlogDir}/www/SKFlatAnalyzerJobLogs/{webdirname}"
while os.path.isdir(webdirpathbase):
    webdirpathbase += "_"

## Skim string
skimString = f"{args.Skim}_" if args.Skim else ""

## Difine master job directory
masterJobDir = f"{SKFlatRunlogDir}/{timestamp}__{randomNumber}__{args.Analyzer}__Era{args.Era}"
if args.Skim:
    masterJobDir = f"{masterJobDir}__{args.Skim}"
for flag in Userflags:
    masterJobDir += f"__{flag}"
masterJobDir += f"__{HOSTNAME}"

## copy library
os.makedirs(masterJobDir)
shutil.copytree(SKFlat_LIB_PATH, f"{masterJobDir}/lib")

## Loop over samples
# mask for each sample
isDoneForEachSample = []
isPostJobDoneForEachSample = []
baseDirForEachSample = []
xsecForEachSample = []
for inputSample in inputSampleList:
    njobs = args.NJobs
    isDoneForEachSample.append(False)
    isPostJobDoneForEachSample.append(False)
    
    # gloval variables
    isDATA = False
    dataPeriod = ""
    if ":" in inputSample:
        isDATA = True
        tmp = inputSample.split(":")
        inputSample, dataPeriod = tmp[0], tmp[1]
    
    # prepare output
    baseRunDir = f"{masterJobDir}/{inputSample}"
    if isDATA:
        baseRunDir = f"{baseRunDir}_period{dataPeriod}"
    os.makedirs(f"{baseRunDir}/output")
    
    # create webdir
    # cf) baseRunDir = $SKFlatRunlogDir/2019_02_26_222038__GetEffLumi__Era2016__KISTI/WW_pythia/
    thisWebDir = f"{webdirpathbase}/{baseRunDir.replace(SKFlatRunlogDir, '').replace(HOSTNAME+'/', HOSTNAME+'__')}"
    os.makedirs(thisWebDir)
    
    # get sample path
    samplePath = f"{SAMPLE_DATA_DIR}/For{sampleHOSTNAME}/{skimString}{inputSample}.txt"
    if isDATA:
        samplePath = f"{SAMPLE_DATA_DIR}/For{sampleHOSTNAME}/{skimString}{inputSample}_{dataPeriod}.txt"
    shutil.copy(samplePath, f"{baseRunDir}/input_filelist.txt")
    
    totalFiles = os.popen("sed 's/#.*//' "+samplePath+"|grep '.root'").readlines()
    nTotalFiles = len(totalFiles)
    if njobs > nTotalFiles or njobs == 0:
        njobs = nTotalFiles
    nfileperjob = int(nTotalFiles/njobs)
    remainder = nTotalFiles - njobs*nfileperjob
    
    fileRanges = []
    tmpEndLargerJob = 0
    nfileCheckSum = 0
    ## First remainder jobs will have (nfileperjob+1) files per job
    for it_job in range(remainder):
        fileRanges.append(range(it_job*(nfileperjob+1), (it_job+1)*(nfileperjob+1)))
        tmpEndLargerJob = (it_job+1)*(nfileperjob+1)
        nfileCheckSum += len(range(it_job*(nfileperjob+1), (it_job+1)*(nfileperjob+1)))
    
    ## Remaining njobs - remainder jobs will have (nfileperjob) files per job
    for it_job in range(njobs - remainder):
        fileRanges.append(range(tmpEndLargerJob+(it_job*nfileperjob), tmpEndLargerJob+((it_job+1)*(nfileperjob))))
        nfileCheckSum += len(range(tmpEndLargerJob+(it_job*nfileperjob), tmpEndLargerJob+((it_job+1)*(nfileperjob))))
    
    submitOutput = open(f"{baseRunDir}/SubmitOutput.log", "w")
    submitOutput.write(f"<SKFlat> nTotlaFiles = {nTotalFiles}\n")
    submitOutput.write(f"<SKFlat> njobs = {njobs}\n")
    submitOutput.write(f"<SKFlat> --> # of files per jobs = {nfileperjob}\n")
    if remainder >= njobs:
        submitOutput.write(f"<SKFlat> remainder = {remainder}\n")
        submitOutput.write(f"<SKFlat> while, njobs = {njobs}\n")
        submitOutput.write(f"<SKFlat> --> exit\n")
        print(f"Wrong matched job numbers, see {baseRunDir}/SubmitOutput.log for more detail")
        exit(1)
    submitOutput.write(f"nfileCheckSum = {nfileCheckSum}\n")
    submitOutput.write(f"nTotalFiles = {nTotalFiles}\n")
    assert nfileCheckSum == nTotalFiles
    submitOutput.close()
    fileRangesForEachSample.append(fileRanges)
    
    ## Get xsec ans sumW
    thisDasname = ""
    thisXsec = -1
    thisSumsign = -1
    thisSumW = -1
    if not isDATA and args.Analyzer != "GetEffLumi":
        sampleFilePath = f"{SAMPLE_DATA_DIR}/CommonSampleInfo/{inputSample}.txt"
        if not os.path.exists(sampleFilePath):
            print(f"No {sampleFilePath}")
            exit(1)
        for line in open(sampleFilePath).readlines():
            if line[0] == "#":
                continue
            words = line.split()
            assert inputSample == words[0]
            thisDasName, thisXsec, thisSumsign, thisSumW = words[1], words[2], words[4], words[5]
            break
    xsecForEachSample.append(thisXsec)
    
    ## Write run script
    commandsFileName = f"{args.Analyzer}_{args.Era}_{inputSample}"
    if isDATA:
        commandsFileName += f"_{dataPeriod}"
    for flag in Userflags:
        commandsFileName += f"__{flag}"
    runCommands = open(f"{baseRunDir}/{commandsFileName}.sh", "w")
    runCommands.write("#!/bin/bash\n")
    runCommands.write("SECTION=`printf $1`\n")
    runCommands.write("WORKDIR=`pwd`\n")
    runCommands.write("Trial=0\n")
    runCommands.write("\n")
    runCommands.write(f'export SKFlat_WD="{SKFlat_WD}"\n')
    runCommands.write("export SKFlat_LIB_PATH=$SKFlat_WD/lib\n")
    runCommands.write('export SKFlatV="Run2UltraLegacy_v3"\n')
    runCommands.write('export DATA_DIR="$SKFlat_WD/data/$SKFlatV"\n')
    runCommands.write('export MYBIN="$SKFlat_WD/bin"\n')
    runCommands.write('export PYTHONDIR="$SKFlat_WD/python"\n')
    runCommands.write("export PATH=${MYBIN}:${PYTHONDIR}:${PATH}\n")
    runCommands.write('export PYTHONPATH="${PYTHONPATH}:${PYTHONDIR}"\n')
    runCommands.write("export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$SKFlat_WD/DataFormats/include/:$SKFlat_WD/AnalyzerTools/include/:$SKFlat_WD/Analyzers/include/\n")
    runCommands.write("export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SKFlat_LIB_PATH\n")
    runCommands.write("\n")
    runCommands.write("source $SKFlat_WD/bin/BashColorSets.sh\n")
    runCommands.write("\n")
    runCommands.write("#### make sure use C locale\n")
    runCommands.write("export LC_ALL=C\n")
    runCommands.write("\n")
    runCommands.write("#### Don't make root history\n")
    runCommands.write("export ROOT_HIST=0\n")
    runCommands.write("\n")
    runCommands.write("#### User conda env for root\n")
    runCommands.write("source /opt/conda/bin/activate\n")
    runCommands.write("conda activate torch\n")
    runCommands.write("\n")    
    runCommands.write("#### modifying LD_LIBRARY_PATH to use libraries in baseRunDir\n")
    runCommands.write(f"export LD_LIBRARY_PATH=$(echo $LD_LIBRARY_PATH|sed 's@'$SKFlat_WD'/lib@{masterJobDir}/lib@')\n")
    runCommands.write("\n")
    runCommands.write('while [ "$Trial" -lt 3 ]; do\n')
    runCommands.write('  echo "#### running ####"\n')
    if args.python:
        runCommands.write(f'  echo "python {baseRunDir}'+'/run_${SECTION}.py"\n')
        runCommands.write(f'  echo Processing python {baseRunDir}'+'/run_${SECTION}.py...\n')
        runCommands.write(f"  python {baseRunDir}"+"/run_${SECTION}.py 2> err.log\n")
    else:
        runCommands.write(f'  echo "root -l -b -q {baseRunDir}'+'/run_${SECTION}.C"\n')
        runCommands.write(f"  root -l -b -q {baseRunDir}"+"/run_${SECTION}.C 2> err.log\n")
    runCommands.write("  EXITCODE=$?\n")
    runCommands.write('  if [ "$EXITCODE" -eq 5 ]; then\n')
    runCommands.write('    echo "IO error occured.. running again in 300 seconds.."\n')
    runCommands.write("    Trial=$((Trial+=1))\n")
    runCommands.write("    sleep 300\n")
    runCommands.write("  else\n")
    runCommands.write("    break\n")
    runCommands.write("  fi\n")
    runCommands.write("done\n")
    runCommands.write("\n")
    runCommands.write('if [ "$EXITCODE" -ne 0 ]; then\n')
    runCommands.write('  echo "ERROR errno=$EXITCODE" >> err.log\n')
    runCommands.write("fi\n")
    runCommands.write("")
    runCommands.write("cat err.log >&2\n")
    runCommands.write("exit $EXITCODE\n")
    runCommands.close()
    
    concurencyLimit = ""
    if args.NMax:
        concurrencyLimits = f"concurrencyLimits = n{args.NMax}.{os.genenv('USER')}"
    request_memory = ""
    if args.Memory:
        request_memory = f"request_memory = {args.Memory}"
    submitCommands = open(f"{baseRunDir}/submit.jds", "w")
    submitCommands.write(f"universe                = vanilla\n")
    submitCommands.write(f"executable              = {commandsFileName}.sh\n")
    submitCommands.write(f"jobbatchname            = {commandsFileName}\n")
    submitCommands.write(f'+singularityimage       = "/data6/Users/choij/Singularity/SKFlat"\n')
    submitCommands.write(f'+singularitybind        = "/cvmfs, /cms, /share"\n')
    submitCommands.write(f"requirements            = HasSingularity\n")
    submitCommands.write(f"arguments               = $(Process)\n")
    submitCommands.write(f"log                     = condor.log\n")
    # submitCommands.write(f"getenv                  = True\n")
    submitCommands.write(f"should_transfer_files   = Yes\n")
    submitCommands.write(f"when_to_transfer_output = ON_EXIT\n")
    submitCommands.write(f"output                  = job_$(Process).log\n")
    submitCommands.write(f"error                   = job_$(Process).err\n")
    submitCommands.write(f'transfer_output_remaps  = "hists.root = output/hists_$(Process).root"\n')
    submitCommands.write(f"{concurencyLimit}\n")
    submitCommands.write(f"{request_memory}\n")
    submitCommands.write(f"queue {njobs}\n")
    submitCommands.close()
    
    checkTotalNFile = 0
    for it_job in range(len(fileRanges)):
        checkToalNFile = checkTotalNFile + len(fileRanges[it_job])
        thisJobDir = f"{baseRunDir}/job_{it_job}"
        libDir = f"{masterJobDir}/lib".replace('///', '/').replace('//', '/')
        runFunctionName = f"run_{it_job}"
        runCfileFullPath = f"{baseRunDir}/run_{it_job}.C"
        runPyfileFullPath = f"{baseRunDir}/run_{it_job}.py"
        lhapdfpath = "/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/lhapdf/6.2.3/lib/libLHAPDF.so"
        if args.python:
            # first copy analyzer
            shutil.copy(f"PyAnalyzers/{args.Analyzer}.py", runPyfileFullPath)
            out = open(runPyfileFullPath, "a")
            out.write("\n\n\n")
            out.write('if __name__ == "__main__":\n')
            out.write(f"    m = {args.Analyzer}()\n")
            out.write(f'    m.SetTreeName("recoTree/SKFlat")\n')
            if isDATA:
                out.write(f'    m.IsDATA = True\n')
                out.write(f'    m.DataStream = "{inputSample}"\n')
            else:
                out.write(f'    m.IsDATA = False\n')
                out.write(f'    m.MCSample = "{inputSample}"\n')
                out.write(f'    m.xsec = {thisXsec}\n')
                out.write(f'    m.sumSign = {thisSumsign}\n')
                out.write(f'    m.sumW = {thisSumW}\n')
                if args.FastSim:
                    out.write(f"    m.IsFastSim = True\n")
                else:
                    out.write(f"    m.IsFastSim = False\n")
            out.write(f'    m.SetEra("{args.Era}")\n')
            if Userflags:
                out.write(f'    m.Userflags = std.vector[TString]()\n')
                for flag in Userflags: out.write(f'    m.Userflags.emplace_back("{flag}")\n')
            for it_file in fileRanges[it_job]:
                thisFileName = totalFiles[it_file].strip("\n")
                out.write(f'    if not m.AddFile("{thisFileName}"): exit(1)\n')
            if isSkimTree:
                tmpFileName = totalFiles[ fileRanges[it_job][0] ].strip("\n")
                skimOutDir = f"/gv0/DATA/SKFlat/{SKFlatV}/{args.Era}"
                if args.OutputDir:
                    skimOutDir = f"{args.Outputdir}/{SKFlatV}/{args.Era}"
                skimOutFileName = ""
                if isDATA:
                    skimOutDir = f"{skimOutDir}/DATA_{args.Analyzer}/{inputSample}/period{dataPeriod}"
                    skimOutFileName = f"SKFlatNtuple_{args.Era}_DATA_{it_job}.root"
                else:
                    skimOutDir = f"{skimOutDir}/MC_{args.Analyzer}/{thisDasName}"
                    skimOutFileName = f"SKFlatNtuple_{args.Era}_MC_{it_job}.root"
                skimOutDir += f"{skimOutDir}/{timestamp}/"
                os.makedirs(skimOutDir)
                out.write(f'    m.SetOutfilePath("{skimOutDir+skimOutFileName}")\n')
            else:
                out.write(f'    m.SetOutfilePath("hists.root")\n')
            if args.Reduction > 1:
                out.write(f"    m.MaxEvent = int(m.fChain.GetEntries()/{args.Reduction})\n")
            out.write(f"    m.Init()\n")
            out.write(f"    m.initializePyAnalyzer()\n")
            out.write(f"    m.initializeAnalyzerTools()\n")
            out.write(f"    m.SwitchToTempDir()\n")
            out.write(f"    m.Loop()\n")
            out.write(f"    m.WriteHist()\n")
            out.close()
        else:
            out = open(runCfileFullPath, "w")
            out.write(f"R__LOAD_LIBRARY({lhapdfpath})\n")
            out.write(f"void {runFunctionName}()" + "{\n")
            out.write(f"    {args.Analyzer} m;\n")
            out.write(f'    m.SetTreeName("recoTree/SKFlat");\n')
            if isDATA:
                out.write(f"    m.IsDATA = true;\n")
                out.write(f'    m.DataStream = "{inputSample}";\n')
            else:
                out.write(f'    m.MCSample = "{inputSample}";\n')
                out.write(f"    m.IsDATA = false;\n")
                out.write(f"    m.xsec = {thisXsec};\n")
                out.write(f"    m.sumSign = {thisSumsign};\n")
                out.write(f"    m.sumW = {thisSumW};\n")
                if args.FastSim:
                    out.write(f"    m.IsFastSim = true;\n")
                else:
                    out.write(f"    m.IsFastSim = false;\n")
            out.write(f'    m.SetEra("{args.Era}");\n')
            if Userflags:
                for flag in Userflags: out.write(f'    m.Userflags.emplace_back("{flag}");\n')
            for it_file in fileRanges[it_job]:
                thisFileName = totalFiles[it_file].strip("\n")
                out.write(f'    if(!m.AddFile("{thisFileName}")) exit(EIO);\n')
            if isSkimTree:
                tmpFileName = totalFiles[ fileRanges[it_job][0] ].strip("\n")
                skimOutDir = f"/gv0/DATA/SKFlat/{SKFlatV}/{args.Era}"
                if args.Outputdir:
                    skimOutDir = f"{args.Outputdir}/{SKFlatV}/{args.Era}"
                skimOutFileName = ""
                if isDATA:
                    skimOutDir = f"{skimOutDir}/DATA_{args.Analyzer}/{inputSample}/period{dataPeriod}"
                    skimOutFileName = f"SKFlatNtuple_{args.Era}_DATA_{it_job}.root"
                else:
                    skimOutDir = f"{skimOutDir}/MC_{args.Analyzer}/{thisDasName}"
                    skimOutFileName = f"SKFlatNtuple_{args.Era}_MC_{it_job}.root"
                skimOutDir += f"{skimOutDir}/{timestamp}/"
                os.makedirs(skimOutDir)
                out.write(f'    m.SetOutfilePath("{skimOutDir+skimOutFileName}");\n')
            else:
                out.write(f'    m.SetOutfilePath("hists.root");\n')
            if args.Reduction > 1:
                out.write(f"    m.MaxEvent = m.fChain->GetEntries()/{args.Reduction};\n")
            out.write(f"    m.Init();\n")
            out.write(f"    m.initializeAnalyzer();\n")
            out.write(f"    m.initializeAnalyzerTools();\n")
            out.write(f"    m.SwitchToTempDir();\n")
            out.write(f"    m.Loop();\n")
            out.write(f"    m.WriteHist();\n")
            out.write("}")
            out.close()
        
    # submit the jobs
    if not args.no_exec:
        cwd = os.getcwd()
        os.chdir(baseRunDir)
        condorOptions = ""
        if args.BatchName:
            condorOptions = f" -batch-name {args.BatchName}"
        os.system(f"condor_submit submit.jds {condorOptions}")
        os.chdir(cwd)

if args.no_exec:
    print(f"- RunDir = {baseRunDir}")
    exit()
    
## Set output directory
## if args.Outputdir is not set, go to default setting
finalOutputPath = args.Outputdir
if not args.Outputdir:
    finalOutputPath = f"{SKFlatOutputDir}/{SKFlatV}/{args.Analyzer}/{args.Era}/"
    for flag in Userflags:
        finalOutputPath += f"{flag}__"
    if isDATA:
        finalOutputPath += "/DATA"
    if isSkimTree:
        finalOutputPath = f"/gv0/DATA/SKFlat/{SKFlatV}/{args.Era}"
if not os.path.exists(finalOutputPath):
    os.makedirs(finalOutputPath)

print(f"#################################################")
print(f"Submission Finished")
print(f"- JobID = {randomNumber}")
print(f"- Analyzer = {args.Analyzer}")
print(f"- Skim = {args.Skim}")
print(f"- inputSampleList = {inputSampleList}")
print(f"- njobs = {njobs}")
print(f"- Era = {args.Era}")
print(f"- UserFlags = {Userflags}")
print(f"- RunDir = {baseRunDir}")
print(f"- output will be send to {finalOutputPath}")
print(f"#################################################")

##########################
## Submittion all done. ##
## Now monitor job      ##
##########################

## Loop over samples again
isAllSampleDone = False
gotError = False
errorLog = ""

try:
    while not isAllSampleDone:
        if gotError: break
        isAllSampleDone = True
        for it_sample in range(len(inputSampleList)):
            inputSample = inputSampleList[it_sample]
            isDone = isDoneForEachSample[it_sample]
            isPostJobDone = isPostJobDoneForEachSample[it_sample]
            
            if isPostJobDone: continue
            else: isAllSampleDone = False
            
            ## Global Variables
            isDATA = False
            dataPeriod = ""
            if ":" in inputSample:
                isDATA = True
                tmp = inputSample.split(":")
                inputSample, dataPeriod = tmp[0], tmp[1]
            
            ## Prepare output
            ## This should be copied from above
            baseRunDir = f"{masterJobDir}/{inputSample}"
            if isDATA:
                baseRunDir += f"_period{dataPeriod}"
            thisWebDir = f"{webdirpathbase}/{baseRunDir.replace(SKFlatRunlogDir,'').replace(HOSTNAME+'/',HOSTNAME+'__')}"
            
            if not isDone:
                ## This sample was not done in the previous monitoring
                ## Monitor again this time
                thisSampleDone = True
                
                ## Write job status until it's done
                statusLog = open(f"{baseRunDir}/JobStatus.log", "w")
                statusLog.write(f"Job submitted at {jobStartTime}\n")
                statusLog.write(f"JobNumber\t| Status\n")
                
                toStatusLog = []
                nEventRan = 0
                finished = []
                eventDone = 0
                eventTotal = 0
                
                totalEventRunTime = 0
                maxTimeLeft = 0
                maxEventRunTime = 0
                
                fileRanges = fileRangesForEachSample[it_sample]
                
                for it_job in range(len(fileRanges)):
                    thisJobDir = baseRunDir
                    thisStatus = CheckJobStatus(thisJobDir, args.Analyzer, it_job, HOSTNAME)
                    
                    if "ERROR" in thisStatus:
                        gotError = True
                        statusLog.write("#### ERROR OCCURED ####\n") 
                        statusLog.write(f"{thisStatus}\n") 
                        errorLog = thisStatus
                        break
                
                    if "FINISHED" not in thisStatus:
                        thisSampleDone = False
                        
                    outLog = ""
                    if "FINISHED" in thisStatus:
                        finished.append("Finished")
                        eventInfo = thisStatus.split()[1].split(":")
                        thisEventDone = int(eventInfo[2])
                        thisEventTotal = int(eventInfo[2])
                        
                        eventDone += thisEventDone
                        eventTotal += thisEventTotal
                        
                        #### start
                        lineEventRunTime = f"{thisStatus.split()[2]} {thisStatus.split()[3]}"
                        thisJobStartTime = GetDatetimeFromMyFormat(lineEventRunTime)
                        #### end
                        lineEventEndTime = f"{thisStatus.split()[4]} {thisStatus.split()[5]}"
                        thisJobEndTime = GetDatetimeFromMyFormat(lineEventEndTime)
                        
                        thisDiff = thisJobEndTime - thisJobStartTime
                        thisEventRunTime = 86400*thisDiff.days + thisDiff.seconds
                        
                        thisTimePerEvent = float(thisEventRunTime)/float(thisEventDone)
                        thisTimeLeft = (thisEventTotal-thisEventDone)*thisTimePerEvent
                        
                        totalEventRunTime += thisEventRunTime
                        maxTimeLeft = max(maxTimeLeft, thisTimeLeft)
                    elif "RUNNING" in thisStatus:
                        outLog = f"{it_job}\t| {this_status.split()[1]} %"
                        eventInfo = thisStatus.split()[2].spllit(":")
                        
                        thisEventDone, thisEventTotal = int(eventInfo[1]), int(eventInfo[2])
                        eventDone += thisEventDone
                        eventTotal += thisEventTotal
                        
                        lineEventRunTime = f"{thisStatus.split()[3]} {thisStatus.split()[4]}"
                        thisJobStartTime = GetDatetimeFromMyFormat(lineEventRunTime)
                        thisDiff = datetime.datetime.now() - thisJobStartTime
                        thisEventRunTime = 86400*thisDiff.days + thisDiff.seconds
                        
                        if thisEventDone == 0: thisEventDone = 1
                        thisTimePerEvent = float(thisEventRunTime) / float(thisEventDone)
                        thisTimeLeft = (thisEventTotal-thisEventDone)*thisTimePerEvent
                        
                        totalEventRunTime += thisEventRunTime
                        maxTimeLeft = max(maxTimeLeft, thisTimeLeft)
                        maxEventRunTime = max(maxEventRunTime, thisEventRunTime)
                        
                        outLog += f"({thisTimeLeft:.1}+ s ran, and {thisEventRunTime:.1} s left)"
                        toStatusLog.append(outlog)
                        nEventRan += 1
                    else:
                        outLog = f"{it_job}\t| {thisStatus}"
                        toStatusLog.append(outLog)
                        
                    ##---- END it_job loop
                    
                if gotError:
                    ## When error occured, change both isDone/postJob flag to True
                    isDoneForEachSample[it_sample] = True
                    isPostJobDoneForEachSample[it_sample] = True
                    break

                for line in toStatusLog: 
                    statusLog.write(f"{line}\n")
                statusLog.write("\n=====================================\n")
                statusLog.write(f"HOSTNAME = {HOSTNAME}\n")
                statusLog.write(f"queue = {args.Queue}\n")
                statusLog.write(f"{len(fileRanges)} job submitted\n")
                statusLog.write(f"{nEventRan} jobs are running\n")
                statusLog.write(f"{len(finished)} jobs are finished\n")
                thisTime = datetime.datetime.now()
                statusLog.write(f"XSEC = {xsecForEachSample[it_sample]}\n")
                statusLog.write(f"eventDone = {eventDone}\n")
                statusLog.write(f"eventTotal = {eventTotal}\n")
                statusLog.write(f"eventLeft = {eventTotal-eventDone}\n")
                statusLog.write(f"totalEventRunTime = {totalEventRunTime}\n")
                statusLog.write(f"maxTimeLeft = {maxTimeLeft}\n")
                statusLog.write(f"maxEventRunTime = {maxEventRunTime}\n")
                
                timePerEvent = 1
                if eventDone:       # exist finished events
                    timePerEvent = float(totalEventRunTime)/float(eventDone)
                statusLog.write(f"timePerEvent = {timePerEvent}\n")
                estTime = thisTime+datetime.timedelta(0, maxTimeLeft)
                statusLog.write(f"Estimated finishing time: {estTime.strftime('%Y-%m-%d %H:%M:%S')}\n")
                statusLog.write(f"Last checked at {thisTime.strftime('%Y-%m-%d %H:%M:%S')}\n")
                statusLog.close()
                
                ## copy statusLog to webdir
                shutil.copy(f"{baseRunDir}/JobStatus.log", thisWebDir)
                
                ## This time, it is found to be finished
                ## change the flag
                if thisSampleDone:
                    isDoneForEachSample[it_sample] = True
                
                ##---- END if finished
            else:
                ## Job was finished in the previous monitoring
                ## check if postjob is also finished
                if not isPostJobDone:
                    # copy output and change the postjob flag
                    # if skim, no need to hadd
                    if isSkimTree:
                        isPostJobDoneForEachSample[it_sample] = True
                        continue
                    
                    outputName = f"{args.Analyzer}_{skimString}{inputSample}"
                    if isDATA:
                        outputName += f"_{dataPeriod}"
                    
                    if args.TagOutput:
                        outputName += f"_{args.TagOutput}"
                        
                    if not gotError:
                        cwd = os.getcwd()
                        os.chdir(baseRunDir)
                        
                        ## if number of job is 1, we can just move the file
                        nFiles = len(fileRangesForEachSample[it_sample])
                        if nFiles == 1:
                            os.system('echo "nFiles = 1, so skipping hadd and just move the file" >> JobStatus.log')
                            os.system('ls -1 output/*.root >> JobStatus.log')
                            os.system(f'mv output/hists_0.root {outputName}.root')
                        else:
                            while True:
                                nhadd=int(os.popen("pgrep -x hadd -u $USER |wc -l").read().strip())
                                if nhadd<4: break
                                os.system(f'echo "Too many hadd currently (nhadd={nhadd}). Sleep 60s" >> JobStatus.log')
                                time.sleep(60)
                            print(f"hadd target {outputName}.root")
                            os.system(f'hadd -f {outputName}.root output/*.root >> JobStatus.log')
                            os.system('rm output/*.root')

                        ## Final Outputpath
                        shutil.move(f"{outputName}.root", finalOutputPath)
                        os.chdir(cwd)
                    isPostJobDoneForEachSample[it_sample] = True
        
        if sendLogToWeb:
            shutil.copytree(webdirpathbase, SKFlatLogWebDir)
        time.sleep(20)

except KeyboardInterrupt:
    print("interrupted!")
    
## Send Email now
if sendLogToEmail:
    from SendEmail import *
    jobFinishedEmail = f"""#### Job Info ####
    HOST = {HOSTNAME}
    JobID = {randomNumber}
    Analyzer = {args.Analyzer}
    Era = {args.Era}
    Skim = {args.Skim}
    # of jobs = {njobs}
    inputSample = {inputSampleList}
    {GetXSECTable(inputSampleList, xsecForEachSample)}
    Output sent to: {finalOutputPath}
    """
    emailTitle = f"[{HOSTNAME}] Summary of JobID {randomNumber}"
    if gotError:
        jobFinishedEmail = f"#### ERROR OCCURED ####\n" + jobFinishedEmail
        jobFinishedEmail = errorLog+"\n--------------------------------------\n"+jobFinishedEmail
        emailTitle = f"[ERROR] Summary of JobID {randomNumber}"
    SendEmail(USER, SKFlatLogEmail, emailTitle, jobFinishedEmail)
