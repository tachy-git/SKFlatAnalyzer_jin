import os
from CheckBadFile import *
from optparse import OptionParser


import argparse   
parser = argparse.ArgumentParser(description='option')
parser.add_argument('-e', dest='era', default="NULL")
parser.add_argument('-d', dest='mcdir', default="NULL")
parser.add_argument('-s', dest='skimname', default="NULL")

args = parser.parse_args()

era=args.era
mcdir=args.mcdir
skimname = args.skimname
if skimname== "NULL":
  skimname="SkimTree_HNMultiLep"

Production="Run2UltraLegacy_v3"
tag="SKFlat_Run2UltraLegacy_v3"
Years = ["2016preVFP", "2016postVFP", "2017", "2018"]

if not era == "NULL":
  Years= [era]

  
from os import listdir
from os.path import isfile,isdir, join



out = open('summaryData_'+skimname+'.sh','w')
out.write('#!/bin/sh \n')

for year in Years:
  


  
  _path="DATA_"+skimname
    
  PD = "/gv0/DATA/SKFlat/"+Production+"/"+year+"/"+_path+"/"

  #os.system('bash '+os.getenv("SKFlat_WD")+'/bin/UpdateSampleForSNU.sh  ' + PD)

  DatasetList = [f for f in listdir(PD) if isdir(join(PD,f))]

  if not mcdir == "NULL":
    DatasetList=[mcdir]

  
  for datasets in DatasetList:

    sample_list=[]
    sample_path= PD + datasets  + "/"
    print sample_path

    path_pext = [f for f in listdir(sample_path) if isdir(join(sample_path,f))]
        
    for period in path_pext:
      
      period_path = PD + datasets  + "/" + period + "/"
      
      path_ext = [f for f in listdir(period_path) if isdir(join(period_path,f))]
    
      
      for tag_path in path_ext:
        sample_tag_path= period_path + tag_path+ "/"
      
        sample_list.append(sample_tag_path)
        
      
    for sample in sample_list:
      Files = [f for f in listdir(sample) if isfile(join(sample,f))]
      
      nBad=False
      for File in Files:
        
        filepath_access= sample+File
        if CheckBadFile(filepath_access) or not CheckHasEntries(filepath_access):
          print '@@ ' + filepath_access
          print '@@ BAD FILE'
          print '@@ --> remove'
          os.system("rm " + filepath_access)
          os.system('bash  '  + os.getenv("SKFlat_WD") + '/bin/UpdateSampleForSNU.sh  ' + PD + datasets)
          
          nBad=True
      if not nBad:
        print '-- NO BAD FILES in ' + sample
out.close()
