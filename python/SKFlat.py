#!/usr/bin/env python
import os
import logging 
import shutil
import argparse
import datetime 
import random
import importlib.util

from Submission import SampleListHandler, SampleProcessor

def get_timestamp():
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H%M%S")
    thistime = datetime.datetime.now().strftime("%Y-%m-%D %H:%M:%S") 
    return timestamp, thistime

## Parse arguments
## fastsim and tagoutput are depreciated in this version
def parse_arguments():
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

def set_logging_level(args):
    logging.getLogger().setLevel(logging.INFO)
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    
def process_era(args):
    if args.era == "2016a": args.era = "2016preVFP"
    if args.era == "2016b": args.era = "2016postVFP"

## Make args.userflags as list
def process_userflags(args):
    userflags = []
    if args.userflags != "":
        userflags = args.userflags.replace(' ', '').split(',')
    args.userflags = userflags

## Add absolute path for outputdir
def process_output_dir(args):
    if args.output_dir and (not os.path.isabs(args.output_dir)):
        args.output_dir = f"{os.getcwd()}/{args.output_dir}"

def get_environment_variables(args):
    env_vars = {
        "USER": os.environ["USER"],
        "UID": str(os.getuid()),
        "SKFlat_WD": os.environ["SKFlat_WD"],
        "SKFlatV": os.environ["SKFlatV"],
        "SAMPLE_DATA_DIR": f"{os.environ['SKFlat_WD']}/data/{os.environ['SKFlatV']}/{args.era}/Sample",
        "SKFlatRunlogDir": os.environ["SKFlatRunlogDir"],
        "SKFlatOutputDir": os.environ["SKFlatOutputDir"],
        "SKFlat_LIB_PATH": os.environ["SKFlat_LIB_PATH"],
        "HOSTNAME": os.environ["HOSTNAME"]
    }
    return env_vars

def process_user_info(env_vars):
    user = env_vars["USER"]
    SKFlat_WD = env_vars["SKFlat_WD"]
    path_userinfo = f"{SKFlat_WD}/python/UserInfo_{user}.py"
    try:
        spec = importlib.util.spec_from_file_location("UserInfo", path_userinfo)
        user_info = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(user_info)
        env_vars["SKFlatLogEmail"] = user_info.UserInfo["SKFlatLogEmail"]
        env_vars["SKFlatLogWebDir"] = user_info.UserInfo["SKFlatLogWebDir"]
    except FileNotFoundError:
        logging.error(f"UserInfo file not found: {path_userinfo}")
        env_vars["SKFlatLogEmail"] = ""
        env_vars["SKFlatLogWebDir"] = ""

def determine_hostname(env_vars):
    hostname = env_vars["HOSTNAME"]
    is_tamsa1 = ("tamsa1" in hostname)
    is_tamsa2 = ("tamsa2" in hostname)
    is_tamsa = is_tamsa1 or is_tamsa2
    if is_tamsa:
        if is_tamsa1: hostname = "TAMSA1"
        if is_tamsa2: hostname = "TAMSA2"
        sample_hostname = "ForSNU"
    else:
        logging.info("Working in local...")
        sample_hostname = "ForSNU"
    return hostname, sample_hostname 

def process_jobs_for_skimming(args, hostname):
    if args.skim:
        if args.nmax == 0: args.nmax = 100  # preventing from too heavy IO
        if args.njobs == 1: args.njobs = 0  # njobs =0 means njobs -> nfiles
        if not "TAMSA" in hostname:
            raise ValueError("Skimming is only possible in SNU")

def generate_sample_list_and_hash(args):
    handler = SampleListHandler(args.era, args.data_period)
    if args.input_sample:
        input_sample_list, string_for_hash = handler.generateSampleListFromInputSample(args.input_sample)
    elif args.input_sample_list:
        input_sample_list, string_for_hash = handler.generateSampleListFromInputSampleList(args.input_sample_list)
    else:
        raise ValueError("No input sample")
    
    ## add flags to hash
    for flag in args.userflags:
        string_for_hash += flag
    
    return input_sample_list, string_for_hash

def generate_master_job_dir(args, env_vars, timestamp):
    random_number = int(random.random()*1000000)
    master_job_dir = f"{env_vars['SKFlatRunlogDir']}/{timestamp}__{random_number}__{args.analyzer}__Era{args.era}"
    if args.skim: master_job_dir += f"__{args.skim}"
    for flag in args.userflags:  
        master_job_dir += f"__{flag}"
    master_job_dir += f"__{env_vars['HOSTNAME']}"
    return master_job_dir, random_number

def submission():
    args = parse_arguments()
    set_logging_level(args)
    process_era(args)
    process_userflags(args)
    process_output_dir(args)
    env_vars = get_environment_variables(args)
    process_user_info(env_vars)
    
    hostname, sample_hostname = determine_hostname(env_vars)    
    process_jobs_for_skimming(args, hostname)
    timestamp, job_start_time = get_timestamp()
    
    input_sample_list, string_for_hash = generate_sample_list_and_hash(args)
    master_job_dir, job_id = generate_master_job_dir(args, env_vars, timestamp)
    ## skim string
    ##skim_string = f"{args.skim}_" if args.skim else ""
    
    ## End environment setting
    
    ## Prepare the master_job_dir
    try:
        os.makedirs(master_job_dir)
        shutil.copytree(env_vars["SKFlat_LIB_PATH"], f"{master_job_dir}/lib")
    except Exception as e:
        logging.error(f"Error in making master job directory: {master_job_dir}")
        raise e

    ## processor holder
    processorFor = {}
    ## loop over samples
    isDATA = False
    for input_sample in input_sample_list:
        args_dict = {
            "analyzer": args.analyzer,
            "era": args.era,
            "njobs": args.njobs,
            "nmax": args.nmax,
            "userflags": args.userflags,
            "masterJobDir": master_job_dir,
            "memory": args.memory,
            "skimString": args.skim,
            "SKFlat_WD": env_vars["SKFlat_WD"],
            "SKFlatV": env_vars["SKFlatV"],
            "SAMPLE_DATA_DIR": env_vars["SAMPLE_DATA_DIR"],
            "timestamp": timestamp,
            "reduction": args.reduction,
            "output_dir": args.output_dir,
            "lhapdfpath": f"{env_vars['SKFlat_WD']}/external/lhapdf/lib/libLHAPDF.so",
        }
        processor = SampleProcessor(input_sample, **args_dict)
        processor.prepareRunDirectory()
        processor.generateSubmissionScripts(args.python, bool(args.skim))

        if args.no_exec:
            logging.info(f"RunDir = {processor.baseRunDir}")
            exit()
        else:
            processor.submitJobs(args.batchName)

        if processor.isDATA:
            isDATA = True
        processorFor[input_sample] = processor
        
    final_output_path = args.output_dir
    if not args.output_dir:
        final_output_path = f"{env_vars['SKFlatOutputDir']}/{env_vars['SKFlatV']}/{args.analyzer}/{args.era}/"
    for flag in args.userflags:
        final_output_path += f"{flag}__"
    if isDATA:
        final_output_path += "/DATA"
    if args.skim:
        final_output_path = f"/gv0/DATA/SKFlat/{env_vars['SKFlatV']}/{args.era}"
    os.makedirs(final_output_path, exist_ok=True)

    print(f"#################################################")
    print(f"Submission Finished")
    print(f"JobID = {job_id}")
    print(f"Analyzer = {args.analyzer}")
    print(f"Skim = {args.skim}")
    print(f"input samples = {input_sample_list}")
    print(f"njobs = {args.njobs}")
    print(f"Era = {args.era}")
    print(f"UserFlags = {args.userflags}")
    print(f"RunDir = {processorFor[input_sample_list[0]].baseRunDir}")
    print(f"output will be send to {final_output_path}")
    print(f"#################################################")

    return processorFor

if __name__ == "__main__":
    processorFor = submission()
    print(processorFor)
    
## check job log destination
#sedLogToEmail = False if SKFlatLogEmail  == "" else True
#sendLogToWeb   = False if SKFlatLogWebDir == "" else True