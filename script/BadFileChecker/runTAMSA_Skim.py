import os,csv
from CheckBadFile import *
from optparse import OptionParser

from urllib import urlopen


def get_url_from_era(_era,is_signal):

    """ depending on era loop over dir in tamsa and add line to Samples_MC.txt """                                                                                                                    
                                       
    var_format_val=10
    var_url_2016pre   = 'https://docs.google.com/spreadsheets/d/e/2PACX-1vSVnaxYKoksltw7L896wgWPAHuW2lnJfrFYWxIIIUHVZXJ-1IzmFd2BV2806B_RKwIgKP9rcHxTw242/pub?gid=1589532672&single=true&output=csv'
    var_url_2016post  = 'https://docs.google.com/spreadsheets/d/e/2PACX-1vSVnaxYKoksltw7L896wgWPAHuW2lnJfrFYWxIIIUHVZXJ-1IzmFd2BV2806B_RKwIgKP9rcHxTw242/pub?gid=1502984176&single=true&output=csv'
    var_url_2017      = 'https://docs.google.com/spreadsheets/d/e/2PACX-1vSVnaxYKoksltw7L896wgWPAHuW2lnJfrFYWxIIIUHVZXJ-1IzmFd2BV2806B_RKwIgKP9rcHxTw242/pub?gid=604962654&single=true&output=csv'
    var_url_2018      = 'https://docs.google.com/spreadsheets/d/e/2PACX-1vSVnaxYKoksltw7L896wgWPAHuW2lnJfrFYWxIIIUHVZXJ-1IzmFd2BV2806B_RKwIgKP9rcHxTw242/pub?gid=225171833&single=true&output=csv'

    if is_signal:
        var_url_2016pre   = 'https://docs.google.com/spreadsheets/d/e/2PACX-1vSVnaxYKoksltw7L896wgWPAHuW2lnJfrFYWxIIIUHVZXJ-1IzmFd2BV2806B_RKwIgKP9rcHxTw242/pub?gid=1966662212&single=true&output=csv'
        var_url_2016post = 'https://docs.google.com/spreadsheets/d/e/2PACX-1vSVnaxYKoksltw7L896wgWPAHuW2lnJfrFYWxIIIUHVZXJ-1IzmFd2BV2806B_RKwIgKP9rcHxTw242/pub?gid=277494749&single=true&output=csv'
        var_url_2017 = 'https://docs.google.com/spreadsheets/d/e/2PACX-1vSVnaxYKoksltw7L896wgWPAHuW2lnJfrFYWxIIIUHVZXJ-1IzmFd2BV2806B_RKwIgKP9rcHxTw242/pub?gid=1117199949&single=true&output=csv'
        var_url_2018 = 'https://docs.google.com/spreadsheets/d/e/2PACX-1vSVnaxYKoksltw7L896wgWPAHuW2lnJfrFYWxIIIUHVZXJ-1IzmFd2BV2806B_RKwIgKP9rcHxTw242/pub?gid=461082258&single=true&output=csv'

    if _era == "2016preVFP":
        return var_url_2016pre
    if _era == "2016postVFP":
        return var_url_2016post
    if _era == "2017":
        return var_url_2017
    if _era == "2018":
        return var_url_2018
    else: 
        print _era
    return ""


def get_dsn_from_alias(_era, _alias, _isSignal):

    _url= get_url_from_era(_era, _isSignal)
    var_index_var=1

    _value="NULL"
    cr = csv.reader(urlopen(_url).readlines())
    for row in cr:
        if row[0] == "END":
            break
        if not row[1] == "DS":
            if row[0] == _alias:
                _value=row[var_index_var]


    if _value == "NULL":
        _url_sig= get_url_from_era(_era,True)
        cr2 = csv.reader(urlopen(_url_sig).readlines())
        for row in cr2:
            if row[0] == "END":
                break
            if not row[1] == "DS":
                if row[0] == _alias:
                    _value=row[var_index_var]



    return _value



import argparse   
parser = argparse.ArgumentParser(description='option')
parser.add_argument('-e', dest='era', default="NULL")
parser.add_argument('-d', dest='mcdir', default="NULL")
parser.add_argument('-l', dest='mclist', default="NULL")
parser.add_argument('-s', dest='skimname', default="NULL")
parser.add_argument('--RunData', action='store_true')

args = parser.parse_args()

era=args.era
mcdir=args.mcdir
skimname = args.skimname
mclist=args.mclist

if skimname== "NULL":
  skimname="SkimTree_HNMultiLep"

Production="Run2UltraLegacy_v3"
tag="SKFlat_Run2UltraLegacy_v3"
Years = ["2016preVFP", "2016postVFP", "2017", "2018"]

if not era == "NULL":
  Years= [era]

  
from os import listdir
from os.path import isfile,isdir, join
 
out = open('summaryMC_'+skimname+'.sh','w')
out.write('#!/bin/sh \n')

for year in Years:
  
    _path='MC_'+skimname
    
    PD = "/gv0/DATA/SKFlat/"+Production+"/"+year+"/"+_path+"/"


    #os.system('bash '+os.getenv("SKFlat_WD")+'/bin/UpdateSampleForSNU.sh  ' + PD)
  
    DatasetList = [f for f in listdir(PD) if isdir(join(PD,f))]

    if not mcdir == "NULL":
        DatasetList=[mcdir]

        
    if not mclist == "NULL":
        DatasetList=[]
        dlist = open(mclist,"r")
        for dline in dlist:
            dline=dline.split()[0]
            dsn=get_dsn_from_alias(year, dline.split()[0], False)
            if not dsn == "NULL":
                DatasetList.append(dsn)
        dlist.close()

    for datasets in DatasetList:

        sample_list=[]
        sample_path= PD + datasets  + "/"
        print sample_path

        path_ext = [f for f in listdir(sample_path) if isdir(join(sample_path,f))]
          
        for tag_path in path_ext:
            sample_tag_path= PD + datasets  + "/" + tag_path+ "/";
        
            sample_list.append(sample_tag_path)
      
      
        for sample in sample_list:
            Files = [f for f in listdir(sample) if isfile(join(sample,f))]
      
            nBad=False
            to_update=True
            for File in Files:
                filepath_access= sample+File
                if CheckBadFile(filepath_access) or not CheckHasEntries(filepath_access):
                    print '@@ ' + filepath_access
                    print '@@ BAD FILE'
                    print '@@ --> remove'
                    os.system("rm " + filepath_access)
                    nBad=True
                    if to_update:
                        to_update=False
                        os.system('bash  '  + os.getenv("SKFlat_WD") + '/bin/UpdateSampleForSNU.sh  ' + sample_path)
                        
            if not nBad:
                print '-- NO BAD FILES in ' + sample

out.close()
