#!/usr/bin/env python

### Include ROOT for TH1 GetEffectiveLumi integration 
import os,ROOT



def Exit(timelist, exit_code):
    
    Initial_time=0.
    for x in timelist:
        if Initial_time == 0:
            Initial_time=x[1]
        print str(x[0]) + '  ' + str(x[1]-Initial_time) +' s'
        Initial_time=x[1]

    if exit_code:
        exit()

### include time for printing of processing tiume
import thread,time

#----- start time for verbose output
vtimes=[]
start_time = time.time()
vtimes.append(["start_time",time.time()])

''' import user functions/classes '''
from os import listdir
from CheckLocalSKFlat import *
from SampleFileCore   import *


os.system("rm *~;rm -rf tmp*")

##---------------- SET ARG INPUT
import argparse
parser = argparse.ArgumentParser(description='option')
parser.add_argument('-e', dest='era', default="NULL")
parser.add_argument('-s', dest='skim', default="NULL")
parser.add_argument('--Full', action='store_true')
parser.add_argument('--Signal', action='store_true')
parser.add_argument('--UserInput', action='store_true')
parser.add_argument('--VERBOSE', action='store_true')
args = parser.parse_args()
input_era=args.era
input_skim=args.skim
update_script='script/UpdateSampleForSNULatest.sh'


user = os.getenv("USER")

if not "SKFlat_WD"  in os.environ:
    print 'Please setup area.....'
    exit()


print ("UpdateTamsaSKFlatSampleInput [config] ")
if args.era == "Null":

    print "Running All eras"
else:
    print "Running " +  args.era



GIT_BRANCH="Run2UltraLegacy-jalmond"
GIT_BRANCH_SIGNAL="HNL_Run2UltraLegacy-shared"

os.system("git br -vv | grep '*' > tempGIT.txt")
git_branch_infoline = open("tempGIT.txt").readlines()[0]
os.system('rm tempGIT.txt')
print git_branch_infoline.split()[1] 
if args.Signal:
    if not  git_branch_infoline.split()[1] == GIT_BRANCH_SIGNAL:
        print ("ERROR: Running SIGNAL option in wrong branch... exit ")
        print ("git co " + GIT_BRANCH_SIGNAL)
        exit()
else:
    if not  git_branch_infoline.split()[1] == GIT_BRANCH:
        print ("ERROR: Running in wrong branch... exit ")
        print ("git co " + GIT_BRANCH)

        exit()

if args.UserInput:
    update_script=os.getenv("SKFlat_WD")+'/bin/UpdateSampleForSNU.sh'

''' set up global vars '''
var_skflat_wd=os.getenv("SKFlat_WD")
var_sk_out=os.getenv("SKFlatOutputDir")

FLAT_Version= ReadConfig("SKFlatVersion")
path_dir_sk_datadir=var_skflat_wd+"/data/"+FLAT_Version+"/"
path_sklat_dir="/gv0/DATA/SKFlat/"+FLAT_Version+"/"

##------------------- set ERAS
''' set eras to run'''
arr_eras         = ReadConfig("Eras")
miniaod_versions = ReadConfig("MINIAODversions")

if not input_era  == "NULL":
    arr_eras=[input_era]

print (" ")*50
print ("------------------------------------------------------------------------------------------")
if args.Full:
    print ('@@ running UpdateTamsaSKFlatSampleInput.py [Full update/check] ..... this may take longer')
else:
    print ('@@ running UpdateTamsaSKFlatSampleInput.py [Partial update/check]')
print ("------------------------------------------------------------------------------------------")
''' loop over eras '''

skim_list = []

##------------------- SET UP SKIMS THAT CODE CAN RUN
if  input_skim == "HNL":
    skim_list   = ["SkimTree_HNFake","SkimTree_HNMultiLep"]
    skim_string = " "

    for x in skim_list:
        skim_string = skim_string + x + " "
    print ("@@ running UpdateTamsaSKFlatSampleInput.py skims = %s" % (skim_string))

elif input_skim=="NULL":
    print "@@ running UpdateTamsaSKFlatSampleInput.py running no skims"
else:
    print "@@ running UpdateTamsaSKFlatSampleInput.py invalid skim"
    exit()

if tmpDir(os.getenv("PWD")):
    os.system('rm -r tmp_*')

for era in arr_eras:

    print '-'*70
    print 'UpdateTamsaSKFlatSampleInput: Skipping the following DSN:'
    print '-'*70
            
    skip_dsn_list=open(os.getenv("SKFlat_WD")+"/script/SKFlatSampleList/config/DSN_Skip.txt","r")

    for line in skip_dsn_list:
        era_skip = line.split()[0]

        if era_skip == "201*":
            print ' -- ' + line.split()[1]
        elif era_skip == era:
            print ' -- ' + line.split()[1]

    skip_dsn_list.close()
    print '-'*70


    start_era_time = time.time()
    vtimes.append(["start_era_time",time.time()])
    

    #if args.VERBOSE:
    print 'UpdateTamsaSKFlatSampleInput: Processing time = ' + str(round(start_era_time - start_time,2)) 

    workdir=make_work_dir(args.VERBOSE)

    #if args.VERBOSE:
    print "UpdateTamsaSKFlatSampleInput:array_from_googledoc [%s] " % (era)



  ##--- Get list of samples on tamsa server (with samples remmoved if in SKip file)                                                        
    _skdir_list=get_skflat_on_tamsa(path_sklat_dir,era)


    up_start_era_time = time.time()
    vtimes.append(["get_list_time",time.time()])


    if check_processed_status(_skdir_list,era):
        exit ()

    print 'UpdateTamsaSKFlatSampleInput [post check_processed_status] : Processing time = ' + str(round(up_start_era_time - start_time,2)) + ' s'


    array_from_googledoc        = get_array_from_googledoc(era,miniaod_versions,False)
    array_from_googledoc_sig    = get_array_from_googledoc(era,miniaod_versions,True)

    if args.Signal:
        array_from_googledoc = array_from_googledoc_sig
    

    up2_start_era_time = time.time()
    vtimes.append(["get_array_time",time.time()])


    print 'UpdateTamsaSKFlatSampleInput [post get_array_from_googledoc : Processing time = ' + str(round(up2_start_era_time - start_time,2)) + ' s'


    AssignTimeStamp(path_sklat_dir+era+'/MC/',array_from_googledoc)

    MissingSamples(get_array_from_googledoc(era,miniaod_versions,False),array_from_googledoc_sig,path_sklat_dir+era+'/MC/',era)

    up3_start_era_time = time.time()
    vtimes.append(["Missing_file_time",time.time()])
    
    
    print 'UpdateTamsaSKFlatSampleInput [post MissingSamples : Processing time = ' + str(round(up3_start_era_time - start_time,2)) + ' s'


    if args.VERBOSE:
        for row in array_from_googledoc:
            print row[0] + ' '*(40-len(row[0])) + row[1] + ' '*(70-len(row[1]) ) + row[2]                                                                                                  
    os.system('cp -r  ' +var_skflat_wd+"/data/"+FLAT_Version+"/"+era+"/Sample/ForSNU/ " + workdir +'/' )

    if args.Full:
        print 'running bash ' + update_script
        os.system('bash ' + update_script + ' ' + path_sklat_dir+era+'/MC/')

    ##--- Get list of samples on tamsa server (with samples removed if in SKip file)


    era_time_v1 = time.time()
    vtimes.append(["Pre_sample_setup_time",time.time()])
    

    CopyCommonSampleFiles(era,var_skflat_wd+"/data/"+FLAT_Version+"/"+era+"/Sample/CommonSampleInfo/",var_sk_out+"/"+FLAT_Version+"/SampleCheck/"+era+"/", path_dir_sk_datadir+era+"/Sample/SampleSummary_MC.txt",_skdir_list)

    vtimes.append(["post_sample_setup_time",time.time()])

    UpdateFiles = update_mc_samplelist_from_googledoc(array_from_googledoc,era,_skdir_list,path_dir_sk_datadir+era+"/Sample/SampleSummary_MC.txt" ,path_sklat_dir, workdir)
    
    vtimes.append(["post_update_file_time",time.time()])


    check_file_diff(array_from_googledoc,era,UpdateFiles,path_sklat_dir,var_skflat_wd+"/data/"+FLAT_Version+"/"+era+"/Sample/ForSNU/",var_skflat_wd,var_sk_out+"/"+FLAT_Version+"/SampleCheck/"+era+"/", workdir, args.Full)


    vtimes.append(["post_filediff_time",time.time()])

    if len(UpdateFiles) >0:
        print "@@ running UpdateTamsaSKFlatSampleInput.py UpdateFiles:"
        _isample=0
        for x in UpdateFiles:
            print 'Update['+str(_isample)+']: ',x
            _isample=_isample+1


    if  len(UpdateFiles)==0:
        if not args.Full:
            print 'No updates ...'
            print '-'*75
            print ''
            
            PrintSampleInfo(array_from_googledoc)            

            CheckMultipleTimeStamps(path_sklat_dir+era+'/MC/')
            
            if CheckDuplicates(array_from_googledoc):
                print '-----------------------'*3
            update_summarymc_file(_era)

            continue
            
    #''' update input sample lists in ForSNU/ for new files added and in case some files removed'''

    up_start_us_time = time.time()

    for f in UpdateFiles:
        print 'UpdateTamsaSKFlatSampleInput: command [bash ' + update_script + ' '+path_sklat_dir+era+'/MC/'+f + ']'
        os.system('bash ' + update_script + ' '+path_sklat_dir+era+'/MC/'+f)
    up_end_us_time = time.time()

    vtimes.append(["post_update_UpdateSampleForSNU_time",time.time()])

    if args.VERBOSE:
        print 'UpdateTamsaSKFlatSampleInput:time ('+update_script+') '+ era +' = ' + str(round(up_end_us_time - up_start_us_time,2))


    up_start_sdiff_time = time.time()

    check_samplefile_diff(var_skflat_wd+"/data/"+FLAT_Version+"/"+era+"/Sample/ForSNU/", workdir + "/ForSNU/" , var_sk_out+"/"+FLAT_Version+"/SampleCheck/"+era+"/")

    vtimes.append(["post_check_samplefile_diff_time",time.time()])


    up_end_sdiff_time = time.time()

    if args.VERBOSE:
        print 'UpdateTamsaSKFlatSampleInput:time (check_samplefile_diff) '+ era +' = ' + str(round(up_end_sdiff_time - up_start_sdiff_time,2))
    

    #''' make make_common_sampleinfo files for new samples'''

    up_start_comm_time = time.time()


    make_common_sampleinfo(array_from_googledoc,era,_skdir_list,var_skflat_wd+"/data/"+FLAT_Version+"/"+era+"/Sample/CommonSampleInfo/",var_sk_out+"/"+FLAT_Version+"/SampleCheck/")
    up_end_comm_time = time.time()

    if args.VERBOSE:
        print 'UpdateTamsaSKFlatSampleInput:time (make_common_sampleinfo) '+ era +' = ' + str(round(up_end_comm_time - up_start_comm_time,2))
    
    vtimes.append(["post_make_common_sampleinfo_time",time.time()])

    up_start_eff_time = time.time()

    #''' check if sample is new and run eff lumi, and fill summary and input files with nevent info'''
    updated_list = get_effective_lumi(array_from_googledoc,era,var_sk_out+"/"+FLAT_Version+"/", var_sk_out+"/"+FLAT_Version+"/SampleCheck/"+era+"/", path_dir_sk_datadir, _skdir_list,path_dir_sk_datadir+era+"/Sample/SampleSummary_MC.txt",skim_list, workdir,args.Full)


    vtimes.append(["post_get_effective_lumi_diff_time",time.time()])

    up_end_eff_time = time.time()


    if args.VERBOSE:
        print 'UpdateTamsaSKFlatSampleInput:time (get_effective_lumi)  '+ era +' = ' + str(round(up_end_eff_time - up_start_eff_time,2))

    print 'UpdateTamsaSKFlatSampleInput:time (total)  '+ era +' = ' + str(round(up_end_eff_time - start_time,2))


    for skim in skim_list:
        if os.path.exists(path_sklat_dir+era+'/MC_'+skim+'/'):
            for x in updated_list:
                if os.path.exists(path_sklat_dir+era+'/MC_'+skim+'/'+x):
                    print 'UpdateTamsaSKFlatSampleInput: command[python '+var_skflat_wd+'/script/BadFileChecker/runTAMSA_Skim.py -e ' + era + ' -s  ' + skim + ' -d ' + x + ']'
                    os.system('python '+var_skflat_wd+'/script/BadFileChecker/runTAMSA_Skim.py -e ' + era + ' -s  ' + skim + ' -d ' + x)
            os.system('bash '+update_script+' '+path_sklat_dir+era+'/MC_'+skim+'/')
    os.system('rm -rf ' + workdir)

    PrintSampleInfo(array_from_googledoc)
    
    vtimes.append(["post_summary_time",time.time()])


    CheckMultipleTimeStamps(path_sklat_dir+era+'/MC/')

    if CheckDuplicates(array_from_googledoc):
        print '-----------------------'*3


    vtimes.append(["finished_ime",time.time()])

    update_summarymc_file(era)
            
Exit(vtimes,True)
