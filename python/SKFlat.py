#!/usr/bin/env python
import os
import logging 
import shutil
import argparse
import datetime 
import random

from Submission.helper import SampleListHandler, SampleProcessor

## Arguments
parser = argparse.ArgumentParser(description='SKFlat Command')
parser.add_argument("-a", dest="Analyzer", required=True, help="Analyzer name, defined in Analyzers/")
parser.add_argument("-i", dest="InputSample", default="", help="Input sample name")
parser.add_argument("-p", dest="DataPeriod", default="ALL", help="data period")
parser.add_argument("-l", dest="InputSampleList", default="", help="Input sample list")
parser.add_argument("-n", dest="njobs", default=1, type=int, help="no. of jobs to be splitted")
parser.add_argument("-o", dest="OutputDir", default="")
parser.add_argument("-e", dest="Era", default="2017", help="2016preVFP / 2016postVFP / 2017 / 2018")
parser.add_argument('--skim', dest='Skim', default="", help="ex) SkimTree_Dilepton")
parser.add_argument("--python", action="store_true", default=False, help="Run python analyzers")
parser.add_argument("--no_exec", action="store_true", default=False, help="Prepare running directory")
parser.add_argument("--fastsim", action="store_true", default=False, help="FastSim")
parser.add_argument("--debug", action="store_true", default=False, help="Debug mode")
parser.add_argument("--userflags", dest="Userflags", default="", help="User pre-defined flags")
parser.add_argument("--tagoutput", dest="TagOutput", default="")
parser.add_argument("--nmax", dest="nmax", default=0, type=int, help="maximum running jobs")
parser.add_argument("--reduction", dest="Reduction", default=1, type=int)
parser.add_argument("--memory", dest="Memory", default=2048, type=int)
parser.add_argument("--batchname", dest="BatchName", default="")
parser.add_argument("--logevery", dest="logEvery", default=1000, type=int)
args = parser.parse_args()

logging.getLogger().setLevel(logging.INFO)
if args.debug:
    logging.getLogger().setLevel(logging.DEBUG)

if args.Era == "2016a": args.Era = "2016preVFP"
if args.Era == "2016b": args.Era = "2016postVFP"

## make userflags as a list
Userflags = []
if args.Userflags != "":
    Userflags = args.Userflags.replace(' ', '').split(',')

## Add absolute path for outputdir
if args.OutputDir and (not os.path.isabs(args.OutputDir)):
    args.OutputDir = f"{os.getcwd()}/{args.OutputDir}"

## TimeStamp
# 1) directory / file name style
timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H%M%S")
# 2) log style
jobStartTime = datetime.datetime.now().strftime("%Y-%m-%D %H:%M:%S")
thisTime = None

## Environment Variables
USER = os.environ["USER"]
SKFlat_WD = os.environ['SKFlat_WD']
SKFlatV = os.environ['SKFlatV']
SAMPLE_DATA_DIR = f"{SKFlat_WD}/data/{SKFlatV}/{args.Era}/Sample"
SKFlatRunlogDir = os.environ['SKFlatRunlogDir']
SKFlatOutputDir = os.environ['SKFlatOutputDir']
SKFlat_LIB_PATH = os.environ['SKFlat_LIB_PATH']
UID = str(os.getuid())
HOSTNAME = os.environ['HOSTNAME']

# read userinfo
path_userinfo = f"{SKFlat_WD}/python/UserInfo_{USER}.py"
try:
    with open(path_userinfo, 'r') as f:
        exec(f.read())
except:
    logging.error(f"UserInfo file not found: {path_userinfo}")
    raise FileNotFoundError
SKFlatLogEmail = UserInfo["SKFlatLogEmail"]
SKFlatLogWebDir = UserInfo["SKFlatLogWebDir"]

## check job log destination
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
    logging.info("Working in local...")
    sampleHOSTNAME = "SNU"
    
## Are you skimming trees?
isSkimTree = "SkimTree" in args.Analyzer
if isSkimTree:
    if args.nmax == 0: args.nmax = 100  # Preventing from too heavy IO
    if args.njobs == 1: args.njobs = 0  # NJobs = 0 means NJobs -> NFiles
    if not isTAMSA:
        logging.error("Skimming is only possible in SNU")
        exit(1)

## Make sample list
sampleListHandler = SampleListHandler(args.Era, args.DataPeriod)
inputSampleList = []
stringForHash = ""
if args.InputSample:
    inputSampleList, stringForHash = sampleListHandler.generateSampleListFromInputSample(args.InputSample)
elif args.InputSampleList:
    inputSampleList, stringForHash = sampleListHandler.generateSampleListFromInputSampleList(args.InputSampleList)
else:
    raise ValueError("No input sample")

## add flags to hash
for flag in Userflags:
    stringForHash += flag
    
## Skim string
skimString = f"{args.Skim}_" if args.Skim else ""

## Define master job directory
randomNumber = int(random.random()*1000000)
masterJobDir = f"{SKFlatRunlogDir}/{timestamp}__{randomNumber}__{args.Analyzer}__Era{args.Era}"
if args.Skim: masterJobDir += f"__{args.Skim}"
for flag in Userflags:  masterJobDir += f"__{flag}"
masterJobDir += f"__{HOSTNAME}"

## End of environment settings

## Copy libraries
## copy library
try:
    os.makedirs(masterJobDir)
    shutil.copytree(SKFlat_LIB_PATH, f"{masterJobDir}/lib")
except Exception as e:
    raise e

## Loop over samples
# mask for each sample
processorForEachSample = {}
isDATA = False

for inputSample in inputSampleList:
    args_dict = {
        "Analyzer": args.Analyzer,
        "Era": args.Era,
        "njobs": args.njobs,
        "nmax": args.nmax,
        "userflags": args.Userflags,
        "masterJobDir": masterJobDir,
        "Memory": args.Memory,
        "SKFlat_WD": SKFlat_WD,
        "SKFlatV": SKFlatV,
        "SAMPLE_DATA_DIR": SAMPLE_DATA_DIR,
        "skimString": skimString,
        "timestamp": timestamp,
        "Reduction": args.Reduction,
        "OutputDir": args.OutputDir,
        "lhapdfpath": "/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/lhapdf/6.2.3/lib/libLHAPDF.so",
    }
    processor = SampleProcessor(inputSample, **args_dict)
    if processor.isDATA: isDATA = True  # at least one sample is data
    processor.prepareRunDirectory()
    processor.generateSubmissionScripts(args.python, isSkimTree)
    
    if args.no_exec:
        logging.info(f"RunDir = {processor.baseRunDir}")
        exit()
    else:
        processor.submitJobs(args.BatchName)
    
    processorForEachSample[inputSample] = processor    

## Set output directory
## if args.OutputDir is not given, use default output directory
finalOutputPath = args.OutputDir
if not args.OutputDir:
    finalOutputPath = f"{SKFlatOutputDir}/{SKFlatV}/{args.Analyzer}/{args.Era}/"
    for flag in Userflags:
        finalOutputPath += f"{flag}__"
    if isDATA:
        finalOutputPath += "/DATA"
    if isSkimTree:
        finalOutputPath = f"/gv0/DATA/SKFlat/{SKFlatV}/{args.Era}"
os.makedirs(finalOutputPath, exist_ok=True)

print(f"#################################################")
print(f"Submission Finished")
print(f"JobID = {randomNumber}")
print(f"Analyzer = {args.Analyzer}")
print(f"Skim = {args.Skim}")
print(f"inputSampleList = {inputSampleList}")
print(f"njobs = {args.njobs}")
print(f"Era = {args.Era}")
print(f"UserFlags = {Userflags}")
print(f"RunDir = {processorForEachSample[inputSampleList[0]].baseRunDir}")
print(f"output will be send to {finalOutputPath}")
print(f"#################################################")

## Submission done. 
