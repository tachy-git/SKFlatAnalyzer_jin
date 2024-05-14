#!/usr/bin/env python
import os
import argparse
import logging
import shutil
import time
import datetime
import random
import importlib.util

from Submission import SampleListHandler, SampleProcessor
from Monitoring import CondorJobHandler

ENVs = {}

def updateTimeStamp():
    global ENVs
    ENVs["timestamp"] = datetime.datetime.now().strftime("%Y-%m-%d_%H%M%S")
    ENVs["jobstarttime"] = datetime.datetime.now().strftime("%Y-%m-%D %H:%M:%S") 
    
## Parse arguments
## fastsim, tagoutput, queue are depreciated in this version
def parseArguments():
    parser = argparse.ArgumentParser(description='SKFlat Command')
    parser.add_argument("--analyzer", "-a", required=True, type=str, help="name of analyzer")
    parser.add_argument("--input_sample", "-i", default="", type=str, help="input sample name")
    parser.add_argument("--input_sample_list", "-l", default="", type=str, help="input sample list")
    parser.add_argument("--era", "-e", required=True, type=str, help="2016preVFP / 2016postVFP / 2017 / 2018")
    parser.add_argument("--data_period", "-p", default="ALL", type=str, help="data period")
    parser.add_argument("--njobs", "-n", default=1, type=int, help="no. of jobs to be splitted")
    parser.add_argument("--userflags", default="", type=str, help="User pre-defined flags, splitted by ,")
    parser.add_argument("--reduction", default=1, type=int, help="factor to reduce total number of events")
    parser.add_argument("--output_dir", "-o", default="", type=str, help="output directory")
    parser.add_argument("--skim", default="", type=str, help="skim string, ex. SkimTree_Dilepton")
    parser.add_argument("--memory", default=2048, type=int, help="Memory in MB")
    parser.add_argument("--batchname", default="", type=str, help="Batch name")
    parser.add_argument("--log_every", default=1000, type=int, help="logging for ever N events")
    parser.add_argument("--nmax", default=0, type=int, help="maximum running jobs")
    parser.add_argument("--python", action="store_true", default=False, help="Run python analyzers")
    parser.add_argument("--no_exec", action="store_true", default=False, help="Prepare running directory")
    parser.add_argument("--debug", action="store_true", default=False, help="Debug mode")
    args = parser.parse_args()
    return args

def setLoggingLevel(args):
    logging.getLogger().setLevel(logging.INFO)
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)

def processEra(args):
    if args.era == "2016a": args.era = "2016preVFP"
    if args.era == "2016b": args.era = "2016postVFP"
    
## Make args.userflags as list
def processUserflags(args):
    if args.userflags != "":
        args.userflags = args.userflags.replace(' ', '').split(",")
    else:
        args.userflags = []

def processOutputDir(args):
    if args.output_dir and (not os.path.isabs(args.output_dir)):
        args.output_dir = f"{os.getcwd()}/{args.output_dir}"

def processENVs(args):
    global ENVs
    ENVs = {
        "USER": os.environ["USER"],
        "UID": str(os.getuid()),
        "HOSTNAME": os.environ["HOSTNAME"],
        "SKFlat_WD": os.environ["SKFlat_WD"],
        "SKFlatV": os.environ["SKFlatV"],
        "SAMPLE_DATA_DIR": f"{os.environ['SKFlat_WD']}/data/{os.environ['SKFlatV']}/{args.era}/Sample",
        "SKFlatRunlogDir": os.environ["SKFlatRunlogDir"],
        "SKFlatOutputDir": os.environ["SKFlatOutputDir"],
        "SKFlat_LIB_PATH": os.environ["SKFlat_LIB_PATH"],
    }

def processUserInfo():
    global ENVs
    path_userinfo = f"{ENVs['SKFlat_WD']}/python/UserInfo_{ENVs['USER']}.py"
    try:
        spec = importlib.util.spec_from_file_location("UserInfo", path_userinfo)
        user_info = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(user_info)
        ENVs["SKFlatLogEmail"] = user_info.UserInfo["SKFlatLogEmail"]
        ENVs["SKFlatLogWebDir"] = user_info.UserInfo["SKFlatLogWebDir"]
    except:
        raise FileNotFoundError(f"UserInfo file not found: {path_userinfo}")
        
def processIOForSkimming(args):
    if "SkimTree" in args.analyzer:
        if args.nmax == 0: args.nmax = 100  # preventing from too heavy IO
        if args.njobs == 1: args.njobs = 0  # njobs =0 means njobs -> nfiles
        if not "TAMSA" in ENVs['HOSTNAME'].upper():
            raise ValueError("Skimming is only possible in SNU")

def generateSampleList(args):
    handler = SampleListHandler(args.era, args.data_period)
    assert not (args.input_sample and args.input_sample_list)
    if args.input_sample:
        input_sample_list, _ = handler.generateSampleListFromInputSample(args.input_sample)
    elif args.input_sample_list:
        input_sample_list, _ = handler.generateSampleListFromInputSampleList(args.input_sample_list)
    else:
        raise ValueError("No input sample")
    logging.debug(f"input_sample_list: {input_sample_list}")
    
    return input_sample_list

def updateMasterJobDir(args):
    global ENVs
    job_id = int(random.random()*1000000)
    master_job_dir = f"{ENVs['SKFlatRunlogDir']}/{ENVs['timestamp']}__{job_id}__{args.analyzer}__Era{args.era}"
    if args.skim: master_job_dir += f"__{args.skim}"
    for flag in args.userflags:  
        master_job_dir += f"__{flag}"
    master_job_dir += f"__{ENVs['HOSTNAME']}"
    
    ENVs["MasterJobDir"] = master_job_dir
    ENVs["JobID"] = job_id
    logging.debug(f"MasterJobDir: {master_job_dir}")
    
def mkdirFinalOutputPath(args, includeDataSample):
    final_output_path = args.output_dir
    if not args.output_dir:
        final_output_path = f"{ENVs['SKFlatOutputDir']}/{ENVs['SKFlatV']}/{args.analyzer}/{args.era}"
    flags = ""
    for flag in args.userflags:
        flags += f"{flag}__"
    final_output_path = os.path.join(final_output_path, flags)
    if includeDataSample:
        final_output_path = os.path.join(final_output_path, "DATA")
    if "SkimTree" in args.analyzer:
        final_output_path = f"/gv0/DATA/SKFlat/{ENVs['SKFlatV']}/{args.era}"
    os.makedirs(final_output_path, exist_ok=True)
    return final_output_path
    
def preprocess():
    args = parseArguments()
    setLoggingLevel(args)
    processEra(args)
    processUserflags(args)
    processOutputDir(args)
    processENVs(args)       # Update global variables ENVs
    processUserInfo()       # Update global variables ENVs
    processIOForSkimming(args)
    logging.debug(f"ENVs: {ENVs}")
    logging.debug(f"args: {args}")

    return args

if __name__ == "__main__":
    ## preprocess args and ENVs
    args = preprocess()
    updateTimeStamp()  
    
    ## prepare master job directory
    updateMasterJobDir(args)
    master_job_dir, job_id = ENVs["MasterJobDir"], ENVs["JobID"]
    try:
        os.makedirs(master_job_dir)
        shutil.copytree(ENVs["SKFlat_LIB_PATH"], f"{master_job_dir}/lib")
    except Exception as e:
        logging.error(f"Error in making master job directory: {e}")
        raise e

    ## now submission 
    input_sample_list = generateSampleList(args)
    processor_holder = {}
    includeDataSample = False
    includeMCSample = True
    for sample in input_sample_list:
        processor = SampleProcessor(sample, args, ENVs)
        processor.prepareRunDirectory()
        processor.generateSubmissionScripts()
        
        if args.no_exec:
            logging.info(f"RunDir = {processor.baseRunDir}")
            exit()
        else:
            processor.submitJobs()
        
        if processor.isDATA:
            includeDataSample = True
            includeMCSample = False
        processor_holder[sample] = processor
    
    final_output_path = mkdirFinalOutputPath(args, not includeMCSample)
    print(f"#################################################")
    print(f"Submission Finished")
    print(f"JobID = {job_id}")
    print(f"Analyzer = {args.analyzer}")
    print(f"Skim = {args.skim}")
    print(f"input samples = {input_sample_list}")
    print(f"njobs = {args.njobs}")
    print(f"Era = {args.era}")
    print(f"UserFlags = {args.userflags}")
    print(f"RunDir = {processor_holder[input_sample_list[0]].baseRunDir}")
    print(f"output will be send to {final_output_path}")
    print(f"#################################################")
    
    ## Start monitoring    
    ## For each job, check status -> RUNNING / FINISHED / ERROR
    ## If RUNNING, just update the status log
    ## If ERROR, stop monitoring and send error email
    ## If FINISHED, hadd the output files and send email
    ## Loop until all jobs are finished
    isAllSampleDone = False
    while not isAllSampleDone:
        isAllSampleDone = True
        for sample, processor in processor_holder.items():
            if processor.isDone:
                continue
            
            isAllSampleDone = False
            final_output_path = mkdirFinalOutputPath(args, processor.isDATA)
            handler = CondorJobHandler(processor, ENVs["SKFlatLogEmail"])
            handler.monitorJobStatus()  # Update status flags and JobStatus.log 
            handler.postProcess(final_output_path)   

            if processor.isError:
                logging.error(f"Error in {processor.sampleName}")
                exit(1)
        time.sleep(10) 
        if isAllSampleDone:
            logging.info("All jobs are done")
            break
    
    ## No need to check postprocess if skimming
    if "SkimTree" in args.analyzer:
        exit()
    
    ## Now monitor the postprocesses, waiting for all samples to be hadd and send to final output path
    isAllPostJobDone = False
    while not isAllPostJobDone:
        isAllPostJobDone = True
        for sample, processor in processor_holder.items():
            if processor.isPostJobDone:
                continue
            if processor.isError:
                continue
            
            isAllPostJobDone = False
            handler = CondorJobHandler(processor, ENVs["SKFlatLogEmail"])
            handler.monitorPostProcess()
            
            if processor.isError:
                logging.error(f"Error in {processor.sampleName}")
                exit(1)
        time.sleep(10)
        
        if isAllPostJobDone:
            logging.info("All post processes are done")
            exit()
                
