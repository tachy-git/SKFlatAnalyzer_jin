import os
from CheckBadFile import *
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-e", "--e", dest="era", default="NULL",help="tag")
parser.add_option("-d", "--d", dest="mcdir", default="NULL",help="tag")
(options, args) = parser.parse_args()
era=options.era
mcdir=options.mcdir


Production="Run2UltraLegacy_v3"
tag="SKFlat_Run2UltraLegacy_v3"
Years = ["2016preVFP", "2016postVFP", "2017", "2018"]

if not era == "NULL":
  Years= [era]

  
from os import listdir
from os.path import isfile,isdir, join



out = open('summary.sh','w')
out.write('#!/bin/sh \n')

for year in Years:

  #PDs = open('samples_'+Year+'_MC.txt').readlines()

  
  MCpath="MC"
  #PD = PD.splitlines()[0]
  PD = "/gv0/DATA/SKFlat/"+Production+"/"+year+"/"+MCpath+"/"

  DatasetList = [f for f in listdir(PD) if isdir(join(PD,f))]

  if not mcdir == "NULL":
    DatasetList=[mcdir]
  
  for datasets in DatasetList:

    #if not "WWW" in datasets:
     # continue

    sample_list=[]
    sample_path= PD + datasets  + "/"
    print sample_path

    path_ext = [f for f in listdir(sample_path) if isdir(join(sample_path,f))]
          
    for tag_path in path_ext:
      sample_tag_path= PD + datasets  + "/" + tag_path+ "/"
      
      path_ext2 = [f for f in listdir(sample_tag_path) if isdir(join(sample_tag_path,f))]
      for full_path in path_ext2:
        sample_list.append(sample_tag_path+full_path+"/")
      
      
    for sample in sample_list:
      print "samples = " + sample
      Files = [f for f in listdir(sample) if isfile(join(sample,f))]
      
      nBad=False
      for File in Files:
        #print sample + " ---> " + File
        filepath_access= sample+File
        if CheckBadFile(filepath_access) or not CheckHasEntries(filepath_access):
          print '@@ ' + filepath_access
          print '@@ BAD FILE'
          print '@@ --> remove'
          out.write("rm " + filepath_access)
          nBad=True
      if not nBad:
        print '-- NO BAD FILES in ' + sample
out.close()
