#!/usr/bin/env python

import os,ROOT

''' import reclared functions/classes '''
from os import listdir

from optparse import OptionParser

parser = OptionParser()
#parser.add_option("-a", "--a", dest="aliaslist", default="NULL",help="tag")

(options, args) = parser.parse_args()
#alias_list=options.aliaslist


from UpdateSampleFiles import *
from CheckSKFlat import *


''' set up global vars '''
var_skflat_wd=os.getenv("SKFlat_WD")
var_sk_out=os.getenv("SKFlatOutputDir")
path_dir_sk_datadir=var_skflat_wd+"/data/Run2UltraLegacy_v3/"
path_sklat_dir="/gv0/DATA/SKFlat/Run2UltraLegacy_v3/"

''' set eras to run'''
arr_eras= ["2016postVFP", "2016preVFP", "2017","2018"]

workdir=make_work_dir()

print " "*50
print "------------------------------------------------------------------------------------------"
print '@@ running UpdateTamsaSKFlatSampleInput.py....'
print "------------------------------------------------------------------------------------------"
''' loop over eras '''

for era in arr_eras:
    
    arr_alias = get_alias_perera(era,get_skflat_on_tamsa(path_sklat_dir,era))
    get_missing_alias_perera(era,get_skflat_on_tamsa(path_sklat_dir,era), "SkimTree_HNMultiLep")


os.rmdir(workdir)
