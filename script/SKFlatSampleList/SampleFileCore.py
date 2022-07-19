#!/usr/bin/env python

import os, sys, ROOT, csv

from urllib import urlopen


def IsData(alias):

    if "DoubleEG" in alias:
        return True
        
    if "DoubleMuon" in alias:
        return True

    if "EGamma" in alias:
        return True

    if "MuonEG" in alias:
        return True


    if "SingleElectron" in alias:
        return True


    if "SingleMuon" in alias:
        return True

    return False


def now():
    from datetime import datetime
    return str(datetime.now().month) + str(datetime.now().day)+ str(datetime.now().hour)+str(datetime.now().minute)


def make_work_dir(verbose):

    _dirname='tmp_'+now()
    if not os.path.exists(_dirname):
        os.mkdir(_dirname)
    if verbose:
        print "work directory = " + _dirname
    return _dirname

def print_message(_type, _outstring):

    var_format_val=10
    var_printformat="%%%%%%%%%%%%%%"
    
    if _type==1:
        print " "
        print var_printformat*var_format_val
        print "@@ "+_outstring
        print var_printformat*var_format_val
    
    else:
        print _outstring

    return 


def CheckSetup():

    if not os.getenv("SKFLAT_WD") == "":
        return True
    else:
        print ('set up SKFlatAnalyzer area before running........')
        return False
def ReadConfig(conf_var):

    arr_conf=[]
    CONFIG_File='config/config.txt'
    _conf = open(CONFIG_File,"r")
    for line in _conf:
        sline = line
        sline = sline.replace(',','')
        sline = sline.replace('=','')
        if 'arr_eras' in sline and conf_var == "Eras":
            for x in range(1,len(sline.split())):
                arr_conf.append(sline.split()[x])
        if 'miniaod_versions' in sline and conf_var == "MINIAODversions":
            for x in range(1,len(sline.split())):
                arr_conf.append(sline.split()[x])

        if 'SKFlatVersion' in sline and conf_var == "SKFlatVersion":
            return sline.split()[1]
            
                
    _conf.close()

    return arr_conf

def nevents(_filename):

    _localfile = open(_filename, "r")
    nevent_tree1=0
    
    for f in _localfile:
        #os.system("ls " + f.split()[0])
        _file = ROOT.TFile((f.split())[0])
        _tree  = _file.Get("recoTree/SKFlat")
        if _tree:
            nevent_tree1 = nevent_tree1 + _tree.GetEntries()
        _file.Close()
    _localfile.close()

    return nevent_tree1
    print "recoTree/SKFlat  has a total event count of " + str(nevent_tree1)



def get_url_from_era(_era,is_signal):

    """ depending on era loop over dir in tamsa and add line to Samples_MC.txt """                                                                                              

    var_format_val=10

    maod_versions=["MiniAODv2"]
    urls=[]
    for maod_version in maod_versions:
    
        GD_File='config/GoogleDocLink'+maod_version+'.txt'
        if is_signal:
            GD_File='config/GoogleDocLinkSignal'+maod_version+'.txt'

        GD_File_READ = open (GD_File,"r")

        for line in GD_File_READ:
            if len(line.split()) ==2 :
                if line.split()[0] == _era:
                    urls.append(line.split()[1])
            
        GD_File_READ.close()
    return urls
        



def get_url_from_era_def(_era,is_signal,maod_version):

    """ depending on era loop over dir in tamsa and add line to Samples_MC.txt """                                                                                                                                                           
    var_format_val=10

    GD_File='config/GoogleDocLink'+maod_version+'.txt'
    if is_signal:
        GD_File='config/GoogleDocLinkSignal'+maod_version+'.txt'

    GD_File_READ = open (GD_File,"r")

    for line in GD_File_READ:
        if len(line.split()) ==2 :
            if line.split()[0] == _era:
                GD_File_READ.close()
                return line.split()[1]
    print ('Error in assigning GD page from era')
    return '-11111'


def get_array_from_googledoc(_era, miniaod_versions,_is_signal):
    
    _full_gd_array = []

    for mver in miniaod_versions:
        _array_v=get_array_from_googledoc_def(_era, _is_signal, mver)
        
        for x in _array_v:
            _full_gd_array.append(x)

    return _full_gd_array

def get_array_from_googledoc_def(_era, _is_signal,v_maod):
    
    _url = get_url_from_era_def(_era, _is_signal,v_maod)

    arr_doc=[]
    arr_docN=[]

    if _url == "NULL":
        return arr_doc
    print 'Reading ' + _url
    cr = csv.reader(urlopen(_url).readlines())

    for row in cr:
        if row[0] == "END":
            break
        if not row[1] == "DS" and not row[1] == "" and row[2] == "Y":
            arr_doc.append(row )
        if not row[1] == "DS" and not row[1] == "" and row[2] == "N":
            arr_docN.append(row )
            
    print ("Number of samples with process status Y(es) in GoogleDoc = " + str(len(arr_doc)))
    print ' '*20
    print ("Number of samples with process status N(o)  in GoogleDoc = " + str(len(arr_docN)))
    for x in arr_docN:
        print x[0]
    print ' '*20

    return arr_doc

def n_directories(_dirname):
    
    return len(os.walk(_dirname).next()[1])



def find_googledoc_var_from_alias(_era,_varname , _alias,_urls,_urls_sig):

    for x in range(0,len(_urls)):
        _var=find_googledoc_var_from_alias_def(_era,_varname , _alias,_urls[x],_urls_sig[x])
        
        if not _var == "NULL":
            return _var

    print "UpdateSampleFiles:find_googledoc_var_from_alias [ERROR]"
    exit()

def find_googledoc_var_from_alias_def(_era,_varname , _alias,_url,_url_sig):

    var_index_var=1

    if _varname == "xsec":
        var_index_var=4

    _value="NULL"
    cr = csv.reader(urlopen(_url).readlines())
    for row in cr:
        if row[0] == "END":
            break
        if not row[1] == "DS" and not row[1] == "":
            if row[0] == _alias:
                _value=row[var_index_var]

    if _value == "NULL":
        cr2 = csv.reader(urlopen(_url_sig).readlines())
        for row in cr2:
            if row[0] == "END":
                break
            if not row[1] == "DS" and not row[1] == ""  and row[2] == "Y":
                if row[0] == _alias:
                    _value=row[var_index_var]

    return _value



def run_skim_from_googledoc(_era,_alias , _skim,_urls):

    for x in _urls:
        
        if not run_skim_from_googledoc_def(_era,_alias , _skim,x) == "NULL":
            return run_skim_from_googledoc_def(_era,_alias , _skim,x)
        

def run_skim_from_googledoc_def(_era,_alias , _skim,_url):
    
    skim_value="NULL"
    cr = csv.reader(urlopen(_url).readlines())
    skim_row=0
    for row in cr:
        if row[0] == "END":
            break
        if row[1] == "DS" :
            if row[10].split()[0]    == _skim:
                skim_row=10
            elif row[11].split()[0] == _skim:
                skim_row=11
            elif row[12].split()[0] == _skim:
                skim_row=12
            else:
                return "N"
        elif row[0].split()[0] == _alias and  not row[1] == ""  and row[2] == "Y":
            skim_value=row[skim_row].split()[0]
                
    return skim_value

def CheckMultipleTimeStamps(sk_path):

    from os import listdir
    from os.path import isfile,isdir, join
    arr_tamsadir  = [f for f in listdir(sk_path) if isdir(join(sk_path,f))]

    print sk_path
    
    for _dir in arr_tamsadir:
        
        sk_path_mc=sk_path+_dir
        arr_tamsadir_ts  = [f for f in listdir(sk_path_mc) if isdir(join(sk_path_mc,f))]
        
        if len(arr_tamsadir_ts) > 1:
            list_ts=''
            for x in arr_tamsadir_ts:
                list_ts=list_ts+x +', '
            print 'Multiple timestamps for ' + sk_path_mc + '  [ ' + list_ts+' ] '

def GetTimeStamp(sk_path,_dsn):

    from os import listdir
    from os.path import isfile,isdir, join
    arr_tamsadir  = [f for f in listdir(sk_path) if isdir(join(sk_path,f))]

    for _dir in arr_tamsadir:
        if not _dsn == _dir:
            continue
        
        sk_path_mc=sk_path+_dir
        arr_tamsadir_ts  = [f for f in listdir(sk_path_mc) if isdir(join(sk_path_mc,f))]
        
        if len(arr_tamsadir_ts) == 1:
            if arr_tamsadir_ts[0] == "SKFlat_Run2UltraLegacy_v3":
                sk_path_mc=sk_path+_dir+'/SKFlat_Run2UltraLegacy_v3/'
                arr_tamsadir_ts  = [f for f in listdir(sk_path_mc) if isdir(join(sk_path_mc,f))]

        
        if len(arr_tamsadir_ts) == 1:
            return arr_tamsadir_ts[0]

        elif len(arr_tamsadir_ts) > 1:
            return arr_tamsadir_ts[0]


def MissingSamples(array_gd,array_gd_sig,sk_path,_era):

    
    missing = []
    wrong_assign=[]
    dsn_list=[]
    for row in array_gd:
        if row[0] == "END":
            break

        if not row[1] == "":
            dsn_list.append(row[1])

        if not row[1] == "" and row[2] == "Y":
            if not os.path.exists(sk_path+"/"+row[1]):
                wrong_assign.append(row[0])
        else:
            if not os.path.exists(sk_path+"/"+row[1]):
                missing.append(row[0])

    for row in array_gd_sig:
        if row[0] == "END":
            break

        if not row[1] == "":
            dsn_list.append(row[1])

        if not row[1] == "" and row[2] == "Y":
            if not os.path.exists(sk_path+"/"+row[1]):
                wrong_assign.append(row[0])
        else:
            if not os.path.exists(sk_path+"/"+row[1]):
                missing.append(row[0])


    if len(wrong_assign)>0:
        print 'Samples that are set Y for processed on GD but missing in tamsa include'
    for x in wrong_assign:
        print x
    if len(missing) > 0:
        print 'Samples that are filled on GD but missing in tamsa include'
    for x in missing:
        print x


    from os import listdir
    from os.path import isfile,isdir, join
    arr_tamsadir  = [f for f in listdir(sk_path) if isdir(join(sk_path,f))]
    missing_dsn=[]
    for _dir in arr_tamsadir:
        _dir_tmp=_dir
        _dir_tmp.replace(sk_path+'/','')
        if not _dir_tmp in dsn_list:
            
            skip_dsn_list=open(os.getenv("SKFlat_WD")+"/script/SKFlatSampleList/config/DSN_Skip.txt","r")
            skip=False
            for line in skip_dsn_list:
                if "END" in line:
                    break
                if line.split()[1] == _dir_tmp and line.split()[0] == _era  :
                    skip=True
            skip_dsn_list.close()
            if not skip:
                print 'Sample on tamsa missing in googledoc ' + _dir_tmp
                missing_dsn.append(_dir_tmp)

    if len(missing_dsn)> 0:
        print ("Some samples not on google doc. Add following to config/DSN_Skip.txt OR add to googledoc")
            

    for x in missing_dsn:
        #skip_dsn_list=open(os.getenv("SKFlat_WD")+"/script/SKFlatSampleList/config/DSN_Skip.txt","r")
        #skip=False
        #for line in skip_dsn_list:
        #    if line.split()[1] == x:
        #        skip=True
        #skip_dsn_list.close()
        #if not skip:
        print 'echo "'+ _era+ ' ' +  x + '"  >> config/DSN_Skip.txt'
    if len(missing_dsn)> 0:
        var_userinput= raw_input ('If running on SM MC and these are signals you can ignore.... Should we continue [y/n]:')                                                              
        if var_userinput == "y":                                                                                                                        
            return 
        exit()


def AssignTimeStamp(sk_path,array_gd):
    
    array_gd_ts=[]
    for row in array_gd:
        if row[0] == "END":
            array_gd_ts.append(row)
            break
        if not row[1] == "":
            dsn=row[1]
            ts = GetTimeStamp(sk_path,dsn)
            row.append(ts)
            array_gd_ts.append(row)

        else:
            array_gd_ts.append(row)
                        
    return array_gd_ts
    
def CheckDuplicates(array_gd):
    
    alias_list=[]
    dsn_list=[]
    for row in array_gd:
        if row[0] == "END":
            break
        if not row[1] == "":
            if row[0] in alias_list:
                print 'Alias ['+row[0]+'] is used twice OR wrong set as processed twice, fix'
                return True

            if row[1] in dsn_list:
                print 'DATASET ['+row[1]+'] is used twice OR wrong set as processed twice, fix'
                return True
                                    
            alias_list.append(row[0])
            dsn_list.append(row[1])
    return False

def PrintSampleInfo(array_gd):

    print 'Summary for SKFlat TAMSA samples'
    nwidth=254
    print '-'*nwidth
    print 'Alias'+' '*(30) + ' | '  + 'MAOD Version' + ' '*(3)  + ' | ' + 'xsec' + ' '*(11)+ ' | '  + 'DAS Sample [1]=MINIAODSIM' + ' '*(142)+ ' | ' + 'TimeStamp   |' 
    print '-'*nwidth


    for row in array_gd:
        if row[0] == "END":
            break
        if row[1] == "":
            print '-'*nwidth
        if len(row) < 14:
            row.append('-')
            
        if not row[1] == "DS" and not row[1] == "" and row[2] == "Y":
            dasname=row[7]
            dasname=dasname.replace('MINIAODSIM','[1]')
            print row[0] + ' '*(35 - len(row[0])) + ' | ' + row[3] + ' '*(15 - len(row[3]))  + ' | ' + row[4] + ' '*(15 - len(row[4]))+ ' | ' + dasname + ' '*(165 - len(dasname))+ ' | ' + row[13]
    print '-'*nwidth

def find_googledoc_from_dsn(array_gd,_era,_varname , _dsn):

    _value="NULL"

    for row in array_gd:
        if row[0] == "END":
            break
        if not row[1] == "DS" and not row[1] == "" and row[2] == "Y":
            if row[1] == _dsn:
                _value=row[0]
                
    return _value

def find_googledoc_var_from_dsn(array_gd,_era,_varname , _dsn):
    
    var_index_var=0

    if _varname == "xsec":
        var_index_var=4

    _value="NULL"

    for row in array_gd:
        if row[0] == "END":
            break
        if not row[1] == "DS" and not row[1] == "" and row[2] == "Y":
            if row[1] == _dsn:
                _value=row[var_index_var]
            elif row[1] == _dsn+'/':
                _value=row[var_index_var]
                print('Fix '/' in DSN in google spreadsheet for ',dsn,' in era',era)

    return _value


def check_processed_status(_dsnlist,_era):
    
    print_message(1,"check_dsn_list_from_googledoc ["+_era+"]")

    _urls=get_url_from_era(_era,False)

    arr_googledoc_dsn=[]
    arr_notgoogledoc_dsn=[]
    for x in _urls:

        cr = csv.reader(urlopen(x).readlines())
        for row in cr:
            if row[0] == "END":
                break
            if not row[1] == "DS" and not row[1] == "" and row[2] == "Y":
                arr_googledoc_dsn.append(row[1])
            if not row[1] == "DS" and not row[1] == "" and row[2] == "N":
                arr_notgoogledoc_dsn.append(row[1])

    Incorrect=False
    for x in arr_googledoc_dsn:
        if not x in _dsnlist:
            print "check_processed_status: ["+_era+"]   " + x + " googledoc status incorrect"        
            Incorrect=True
    for x in arr_notgoogledoc_dsn:
        if x in _dsnlist:
            print "check_processed_status: ["+_era+"]   " + x + " googledoc status incorrect"
            Incorrect=True

    
    return Incorrect



def make_common_sampleinfo(array_from_googledoc,_era, _dirlist, _outfiledir,data_skoutput):

    print_message(1,"make_commonsample_file ["+_era+"]")

    for DSN in _dirlist:
        
        array_gd = array_from_googledoc

        var_alias = find_googledoc_var_from_dsn(array_gd,_era, "alias", DSN)
        var_xsec  = find_googledoc_var_from_dsn(array_gd,_era, "xsec" , DSN)
        
        if var_alias == "NULL":
            print "skipping " + DSN
            continue

        if os.path.exists(_outfiledir+"/"+var_alias+".txt"):
            #print "make_common_sampleinfo: "+_outfiledir+"/"+var_alias+".txt already exists..."  
            continue

        print_message(2,"make_commonsample_file: processing  " + _outfiledir+"/"+var_alias+".txt")

        commonfile= open(_outfiledir+"/"+ var_alias+".txt","w")
        commonfile.write("# alias PD xsec nmc sumsign sumw \n")
        commonfile.write(var_alias + "\t" + DSN + "\t" + var_xsec +   "\t" + "1." + "\t" + "1."  + "\t" + "1.\n" )
        commonfile.close()
        print "make_common_sampleinfo: added file " + _outfiledir+"/"+ var_alias+".txt"
        
        # Check If old Eff lumi file is present and delete 

        if os.path.exists(data_skoutput+"/"+ var_alias +".txt"):
            print ('rm ' + data_skoutput+"/"+ var_alias +".txt")
            os.system('rm ' + data_skoutput+"/"+ var_alias +".txt")


def update_summarymc_file(_era):
    
    print_message(1,"update_summarymc_file ")
    currentdir = os.getenv("PWD")
    var_skflat_wd=os.getenv("SKFlat_WD")
    FLAT_Version= ReadConfig("SKFlatVersion")

    os.chdir(var_skflat_wd+"/data/"+FLAT_Version+"/"+_era+"/Sample/")
    os.system("python Summarize.py" )
    os.chdir(currentdir)
    return 




def check_ds_in_current_list(_dslist, _varname):

    b_to_update=True

    for l in _dslist:
        _var=l.split()[1]
        if _var == _varname:
            b_to_update=False
    if b_to_update:
        print "check_ds_in_current_list:DSN " + _varname +" not in GD"
        for x in _dslist:
            print x
        
    return b_to_update

def check_xsec_in_current_list( _dslist, _dsn, _varname):

    b_to_update=True

    for l in _dslist:
        _var_dsn=l.split()[1]
        _var_xsec=l.split()[2]
        if _dsn == _var_dsn:
            if _var_xsec == _varname:
                b_to_update=False

    if b_to_update: 
        print "check_xsec_in_current_list: xsec not correct"
    

    return b_to_update


def update_mc_samplelist_from_googledoc(array_from_googledoc,_era, _dirlist, _path_mc_outfile, _path_sklat_dir,work_dir):

    
    update_list=[]
  
    """ depending on era loop over dir in tamsa and add line to Samples_MC.txt"""
    print_message(1,"update_mc_samplelist_from_googledoc: [Start]")

    """ fill arr_copy_file with lines in local MC summary file"""
    arr_copy_file=[]
   
    ToUpdate=False
    r_mc_outfile=open(_path_mc_outfile,"r")
    for line in r_mc_outfile:
        arr_copy_file.append(line)
    r_mc_outfile.close()
    
    """ setup w_outfile to fill with updated list"""
    _tmp_path_mc_outfile="tmp_summary_"+_era+".txt"
    w_mc_outfile = open(_tmp_path_mc_outfile, "w")
    w_mc_outfile.write('# alias PD xsec nmc sumsign sumw \n')

    """ loop over dataset list on tamsa"""

    for dsn in _dirlist:
 
        #""" get alias xsec from googledoc"""

        array_gd = array_from_googledoc


        var_alias = find_googledoc_var_from_dsn(array_gd,_era,"alias", dsn)
        var_xsec  = find_googledoc_var_from_dsn(array_gd,_era,"xsec" , dsn)

        if var_alias == "NULL" :
            continue

        if n_directories(_path_sklat_dir+_era+"/MC/"+dsn) == 0:
            print "update_mc_samplelist_from_googledoc ["+_era+"] directory "+ dsn + " is missing." 
            continue
        if var_alias == "NULL":
            print "Dirname " + dsn + " has no match from google doc"
            #print (array_from_googledoc,array_from_googledoc_sig,_era, _dirlist, _path_mc_outfile, _path_sklat_dir,work_dir)
            print ('DSN= ',dsn)
            print( _path_mc_outfile, _path_sklat_dir,work_dir)  
            
            skip_dsn_list=open(os.getenv("SKFlat_WD")+"/script/SKFlatSampleList/config/DSN_Skip.txt","r")
            skip=False
            for line in skip_dsn_list:
                if line.split()[1] == dsn and line.split()[0] == _era:
                    skip=True
            skip_dsn_list.close()

            if not skip:
                print ('Add sample to condif/DSNSkip.txt OR add to spreadsheet')
                exit()
            else:
                continue

        if var_xsec == "NULL":
            print"xsec not found"
            exit()


        _nevents_no_w="1."
        _nevents_sign="1."
        _nevents_w="1."
        sample_exists=False
        for mcline in arr_copy_file:
            if var_alias == mcline.split()[0]:
                sample_exists=True
                w_mc_outfile.write(mcline)

                _nevents_no_w=mcline.split()[3]
                _nevents_sign=mcline.split()[4]
                _nevents_w=mcline.split()[5]
                mc_xsec = mcline.split()[2]
                if not mc_xsec.split()[0]  == var_xsec.split()[0]:
                    print "Change in xsec for " +var_alias + " -->"  + mc_xsec+ ": " + var_xsec.split()[0]+":"
                    ToUpdate=True
                    update_list.append(dsn)
                    
        line_to_file = var_alias + "\t"  + dsn + "\t"+ var_xsec + "\t"+_nevents_no_w+"\t"+_nevents_sign+"\t"+_nevents_w+" \n"
        if not sample_exists:
            w_mc_outfile.write(line_to_file)
            update_list.append(dsn)

            ToUpdate=True

    w_mc_outfile.close()



    os.system("bash script/diff_script_mcfile.sh "+_tmp_path_mc_outfile + " " +_path_mc_outfile+" " + work_dir)

    if  os.path.getsize(work_dir+"/tmp.txt") > 1:
        var_flag= "NULL"
        ToUpdate=True
                    
        #arr_answers = ['y','n']
        #while not  var_flag in arr_answers:
        print 'Replace ' + _path_mc_outfile + ' with updated file'
        #var_userinput= raw_input ('Replace ' + _path_mc_outfile + ' with updated file [y/n]:')
            #var_flag=var_userinput
            #print var_userinput
            #if var_userinput == "y":
        os.system('mv ' + _tmp_path_mc_outfile + ' '  + _path_mc_outfile )
        print "Updated...... " +_path_mc_outfile
            #else:
            #    print "Did not update " + _path_mc_outfile
    os.system("rm "+work_dir+"/tmp.txt")

    if len(update_list) == 0:
        print ('No new samples to process.....')
    print_message(1,"update_mc_samplelist_from_googledoc: [end]")

    return update_list
    

def check_bad_files(badfile_dir,_era, _dir, RunFull):
    
    if not RunFull:
        return 
    print_message(1,"check_bad_files ["+_era+"] ["+ _dir +"]")
    currentdir = os.getenv("PWD")
    os.chdir(badfile_dir)
    os.system('python  runTAMSA.py -e ' + _era + ' -d ' + _dir)
    os.system("chmod u+x  summary.sh")
    os.system('./summary.sh')
    os.system('rm summary.sh')
    os.chdir(currentdir)
    

def CopyCommonSampleFiles(_era,commonpath, geteff_path, _path_mc_outfile,gd_list_wo_skipped):


    currentdir = os.getenv("PWD")
    var_skflat_wd=os.getenv("SKFlat_WD")
    FLAT_Version= ReadConfig("SKFlatVersion")

    print ("python Summarize.py")
    os.chdir(var_skflat_wd+"/data/"+FLAT_Version+"/"+_era+"/Sample/")
    os.system("python Summarize.py" )
    os.chdir(currentdir)
    os.system('cat ' + _path_mc_outfile)

    print('Removing files in ' + geteff_path)

    from os import listdir
    from os.path import isfile,isdir, join

    if len([f for f in listdir(geteff_path) if isfile(join(geteff_path,f))]) > 0:
        os.system('rm ' + geteff_path + '/*.txt')

    alias_copy_file=[]
    _path_mc_outfile_tmp = "tmp.txt"
    w_mc_outfile=open(_path_mc_outfile_tmp,"w")
    r_mc_outfile=open(_path_mc_outfile,"r")
    for line in r_mc_outfile:
        if "#" in line:
            w_mc_outfile.write(line)
        else:
            
            alias=line.split()[0]
            if float(line.split()[3]) == 1. and float(line.split()[4]) == 1. and float(line.split()[5]) == 1. :
                print ("Removing sample " + alias+ " from " + _path_mc_outfile)
            else:
                alias_copy_file.append(alias)
                w_mc_outfile.write(line)
            
    r_mc_outfile.close()
    w_mc_outfile.close()
    os.system('mv ' + _path_mc_outfile_tmp + ' ' + _path_mc_outfile)
    print('Updated ' + _path_mc_outfile)
    print '-'*40
    os.system('cat ' + _path_mc_outfile)
    print '-'*40

    forsnupath=commonpath
    forsnupath=forsnupath.replace('CommonSampleInfo','ForSNU')


    common_path  = [f for f in listdir(commonpath) if isfile(join(commonpath,f))]


    get_list  = [f for f in listdir(geteff_path) if isfile(join(geteff_path,f))]
    deleted_alias=[]

    for x in common_path:
        
        if IsData(x):
            continue
        alias_from_file=x
        alias_from_file=alias_from_file.replace('.txt','')

        if not alias_from_file in alias_copy_file:
            print ('Alias ' + alias_from_file + ' is not in ' + _path_mc_outfile)
            print ("removing from " + commonpath)
            print ('rm ' + commonpath + '/'+ x + '\n')
            print('rm ' + forsnupath + '/'+ x + '\n')
            os.system('rm ' + commonpath + '/'+ x)
            os.system('rm ' + forsnupath + '/'+ x)
            deleted_alias.append(x)
            continue

        if not os.path.exists(forsnupath+"/"+x):
            print ('Alias ' + alias_from_file + ' is missing forsnu file')
            print ('rm ' + commonpath + '/'+ x + '\n')
            print('rm ' + forsnupath + '/'+ x + '\n')
            os.system('rm ' + commonpath + '/'+ x)
            os.system('rm ' + forsnupath + '/'+ x)
            deleted_alias.append(x)

            continue


        tmp_file=False
        open_common = open(commonpath+"/"+x, "r")
        for  line in open_common:
            code_seg="\t" + "1." + "\t" + "1."  + "\t" + "1." 
            if code_seg in line:
                tmp_file=True
        open_common.close()

        if tmp_file:
            print 'File ' + commonpath+"/"+x  + ' is set to 1. as default \n'
            os.system('cat ' + commonpath+"/"+x)
            #l_userinput= raw_input ('rm ' + commonpath + '/'+ x +  '[y/n] ')
            #if l_userinput == "y":
            print ('rm ' + commonpath + '/'+ x + '\n')
            print('rm ' + forsnupath + '/'+ x + '\n')
            os.system('rm ' + commonpath + '/'+ x)
            os.system('rm ' + forsnupath + '/'+ x)
            deleted_alias.append(x)
        else:
            
            if not x in get_list:
                print 'cp ' + commonpath + '/'+ x  + ' ' + geteff_path
                os.system('cp ' + commonpath + '/'+ x  + ' ' + geteff_path)


    forsnu_path  = [f for f in listdir(forsnupath) if isfile(join(forsnupath,f))]
    for x in forsnu_path:
        if "Skim" in x:
            continue
        if "skim" in x:
            continue

        if IsData(x):
            continue
        if not os.path.exists(commonpath + '/'+ x ):
            print ("missing | " +commonpath + '/'+ x )
            print ('rm ' + commonpath + '/'+ x + '\n')
            print('rm ' + forsnupath + '/'+ x + '\n')
            os.system('rm ' + commonpath + '/'+ x)
            os.system('rm ' + forsnupath + '/'+ x)
            deleted_alias.append(x)

            continue

        nfiles=0
        readfile = open(forsnupath + '/'+ x,"r")
        for line in readfile:
            if ".root" in line:
                nfiles=nfiles+1
        readfile.close()
        
        if nfiles == 0:
            print ("Empty file " + forsnupath + '/'+ x)
            os.system("cat " + forsnupath + '/'+ x)
            print ('rm ' + commonpath + '/'+ x + '\n')
            print('rm ' + forsnupath + '/'+ x + '\n')
            os.system('rm ' + commonpath + '/'+ x)
            os.system('rm ' + forsnupath + '/'+ x)
            deleted_alias.append(x)

            continue


            
    print ("Delete removed files from samplelist")
    for dfile in deleted_alias:
        if  os.path.exists(geteff_path+"/"+ dfile):
            os.system('rm ' + geteff_path+"/"+ dfile)

    alias_copy_file=[]
    _path_mc_outfile_tmp = "tmp2.txt"
    w_mc_outfile=open(_path_mc_outfile_tmp,"w")
    r_mc_outfile=open(_path_mc_outfile,"r")
    for line in r_mc_outfile:
        if "#" in line:
            w_mc_outfile.write(line)
        else:
            alias=line.split()[0]+'.txt'
            if not  alias in deleted_alias:
                w_mc_outfile.write(line)
            else:
                print ("Removing sample " + alias+ " from " + _path_mc_outfile)

    r_mc_outfile.close()
    w_mc_outfile.close()
    os.system('mv ' + _path_mc_outfile_tmp + ' ' + _path_mc_outfile)
    print(' '*40)
    print('*'*40)
    os.system('cat ' + _path_mc_outfile)
    print('*'*40)


    print (geteff_path)
    os.system('ls ' + geteff_path)
    get_list2  = [f for f in listdir(geteff_path) if isfile(join(geteff_path,f))]

    print '-'*40
    print ('List of samples already processed')
    for x in get_list2:
        print x
    print ''*40
    print '-'*40
    print ('List of samples to be processed')
    for x in deleted_alias:
        print x
    print '-'*40

    print ("python Summarize.py")
    os.chdir(var_skflat_wd+"/data/"+FLAT_Version+"/"+_era+"/Sample/")
    os.system("python Summarize.py" )
    os.chdir(currentdir)
    os.system('cat ' + _path_mc_outfile)
    

        
    return


def get_effective_lumi(array_from_googledoc,_era,_skoutput ,data_skoutput, _skdatadir,_dirlist,_summary_path,skim_list, _workdir, RunFull):

    """  fill file with all alias that has not been ran previously"""

    print_message(1,"get_effective_lumi ["+_era+"]")

    """ loop over dataset list on tamsa"""
    var_url = get_url_from_era(_era,False)

    arr_alias=[]
    arr_alias_torun=[]
    array_gd = array_from_googledoc

    for dsn in _dirlist:
        
        array_gd = array_from_googledoc


        var_alias = find_googledoc_var_from_dsn(array_gd,_era,"alias", dsn)

        if var_alias == "NULL" :
            print "Skipping NULL [get_effective_lumi] " +dsn
            continue

        if not os.path.exists(data_skoutput+"/"+ var_alias +".txt"):
            print "get_effective_lumi: adding " + var_alias + " to processing list"
            arr_alias_torun.append(var_alias)

         
        arr_alias.append(var_alias)
    if os.path.exists( _workdir+"/MC"+_era+".txt"):
        os.system("rm " + _workdir +  "/MC"+_era+".txt")

    w_list=open(_workdir+"/MC"+_era+".txt","w")

    for x in arr_alias_torun:
        w_list.write(x.split()[0]+"\n")
    w_list.close()
    
    return_list=[]

    if len(arr_alias_torun) > 0:
        currentdir = os.getenv("PWD")
        print "SKFlat.py -a GetEffLumi -l "+currentdir+"/"+_workdir+"/MC"+_era+".txt -n 50  --nmax 300  -e "+_era

        os.chdir(os.getenv("SKFlat_WD"))
        os.system("SKFlat.py -a GetEffLumi -l "+currentdir+"/"+_workdir+"/MC"+_era+".txt -n 50  --nmax 300  -e "+_era )
        for x in arr_alias_torun:
            print "SKFlat.py -a GetEffLumi -i "+x.split()[0] +" -n 50  --nmax 300  -e "+_era 
            #print ('cp ' + _skoutput + "/GetEffLumi/"+_era + "/GetEffLumi_"+ x.split()[0] +".root " + data_skoutput+"/")
            #os.system('cp ' + _skoutput + "/GetEffLumi/"+_era + "/GetEffLumi_"+ x.split()[0] +".root " + data_skoutput+"/")

        #GetEFf=False
        #while not GetEFf:
        #    l_userinput= raw_input ('Check if Eff lumi job is ok?  [y/ MCname]:')   
        #    if l_userinput == "y" :                                                                                                                                                             
        #        print "Good"
        #        GetEFf=True
        #    else:
        #        os.system("SKFlat.py -a GetEffLumi -i "+l_userinput+" -n 50  --nmax 300  -e "+_era )
                

        for skim in skim_list:
            new_list=open(currentdir+"/"+_workdir+"/MC"+_era+".txt","r")
            new_skimlist=open(currentdir+"/"+_workdir+"/MC_"+skim+"_"+_era+".txt","w")
            runSkim=False
            for l in new_list:
                l = l.split()[0]
                #allowed_inputs=['y','n']
                #l_userinput ='NULL'
                #while not l_userinput in allowed_inputs:
                #    l_userinput= raw_input ('Sample to update ['+l+']: make skim ' + skim + '  [y/n]:')

                
                if run_skim_from_googledoc(_era,l , skim,var_url) == "Y": #l_userinput == "y":
                    print 'Sample to update ['+l+']: make skim ' + skim
                    new_skimlist.write(l+'\n')
                    runSkim=True
                    return_list.append(find_googledoc_var_from_alias(_era, "dsn", l.split()[0],var_url))

            new_list.close()
            new_skimlist.close()
            if runSkim:
                os.system("SKFlat.py -a "+skim+" -l "+currentdir+"/"+_workdir+"/MC_"+skim+"_"+_era+".txt -n 100  --nmax 300  -e "+_era )
                
            
            os.system("rm "+currentdir+"/"+_workdir+"/MC_"+skim+"_"+_era+".txt")
        os.chdir(currentdir)
        
    else:
        print "get_effective_lumi: all samples proccessed previously"

    ''' delete job submittion file '''
    os.system("rm "+ _workdir+ "/MC"+_era+".txt")    

    
    ''' run over ds list at tamsa and fill common samplefile'''

    print ('Fill CommonSampleFiles')

    update_array=[]

    for dsn in _dirlist:
        ''' access alias and xsec fmor google doc'''

        array_gd = array_from_googledoc

        var_alias = find_googledoc_var_from_dsn(array_gd,_era,"alias", dsn)
        if var_alias == "NULL" :
            continue


        if not RunFull:
            if not var_alias in arr_alias_torun:
                print ("skipping " + var_alias + " since not running Full mode")
                continue
                
        else:
            print 'Filling for ' + dsn

        var_xsec  = find_googledoc_var_from_dsn(array_gd,_era,"xsec" , dsn)

        #''' get nevents from GetEffLumi job'''

        dirpath   = _skoutput + "/GetEffLumi/"+_era + "/"
        _file     = ROOT.TFile(dirpath + "/GetEffLumi_"+ var_alias + ".root")
        hist      = _file.Get("sumW")
        nevents_w = hist.Integral()
        signhist      = _file.Get("sumSign")
        nevents_sign = signhist.Integral()
        _file.Close()
        nevents_no_w=0
        
        orig_xsec=""
        orig_nevent_no_w=""
        orig_nevent_sign=""
        orig_nevent_w=""

        
        print "Reading  : " + dsn
        print "Reading  " + _skdatadir + _era+ "/Sample/CommonSampleInfo/"+var_alias+".txt"

        orig_common_list = open(_skdatadir + _era+ "/Sample/CommonSampleInfo/"+var_alias+".txt","r")

        for line in orig_common_list:
            if not "#" in line:
                if len(line.split()) < 1:
                    continue
                orig_xsec=line.split()[2]
                orig_nevent_no_w=line.split()[3]
                orig_nevent_sign=line.split()[4]
                orig_nevent_w=line.split()[5]
        orig_common_list.close()

        print ('Filled for original values')
        update_file=False

        if not orig_xsec == var_xsec:
            update_file=True
            print "CommonSampleInfo xsec updated for " + var_alias

            #if not orig_nevent_no_w==str(nevents_no_w):
            #update_file=True
            #print "CommonSampleInfo xsec updated for nevents_no_w " + str(nevents_no_w)
        if not float(orig_nevent_w)==nevents_w:
            update_file=True
            nevents_no_w=nevents(_skdatadir + _era+ "/Sample/ForSNU/"+var_alias + ".txt")
            print "CommonSampleInfo updated for nevents_w " + str(nevents_w)

        elif not float(orig_nevent_sign)==nevents_sign:
            update_file=True
            nevents_no_w=nevents(_skdatadir + _era+ "/Sample/ForSNU/"+var_alias + ".txt")
            print "CommonSampleInfo updated for nevents_sign " + str(nevents_sign)
        else: 
            nevents_no_w=orig_nevent_no_w

        if update_file:
            ''' make commonfile for alias'''
            common_list=open(_skdatadir + _era+ "/Sample/CommonSampleInfo/"+var_alias+".txt","w")
            common_list.write("# alias PD xsec nmc sumsign sumw\n")
            common_list.write( var_alias + "\t" + dsn + "\t" + var_xsec + "\t" + str(nevents_no_w) +"\t"+ str(nevents_sign) +"\t"+ str(nevents_w)+" \n")   


            common_list.close()
            os.system("git diff " + _skdatadir + _era+ "/Sample/CommonSampleInfo/"+var_alias+".txt")

            update_array.append([var_alias,var_xsec, nevents_no_w, nevents_sign , nevents_w])


    if len(update_array) > 0:
        update_summarymc_file(_era)

    return return_list

def find_timestamp(_path_to_samplelist , _timestamps):

    list_used = []
    for ts in _timestamps:
        filelist = open(_path_to_samplelist,"r")
        for rline in filelist:
            if ts in rline:
                list_used.append(ts)
                break
        filelist.close()
    return list_used

def CopyGetEff(def_dir, new_dir):
    
    from os import listdir
    from os.path import isfile,isdir, join
    arr_tmp_dir  = [f for f in listdir(def_dir) if isfile(join(def_dir,f))]
    for _dir in arr_tmp_dir:
        if not os.path.exists(new_dir+'/'+_dir):
            os.system('cp ' + def_dir + '/'+_dir + ' ' + new_dir+'/')
            print('cp ' + def_dir + '/'+_dir + ' ' + new_dir+'/')

def tmpDir(cdir):
    
    from os import listdir
    from os.path import isfile,isdir, join
    arr_tmp_dir  = [f for f in listdir(cdir) if isfile(join(cdir,f))]
    for _dir in arr_tmp_dir:
        if 'tmp' in _dir:
            return True
    return False

def check_samplefile_diff(orig_dir, updated_dir, geteff_filepath):
    
    from os import listdir
    from os.path import isfile,isdir, join
    arr_orig_dir  = [f for f in listdir(orig_dir) if isfile(join(orig_dir,f))]
    arr_updated_dir  = [f for f in listdir(updated_dir) if isfile(join(updated_dir,f))]

    import filecmp

    for x in arr_updated_dir:
        if not filecmp.cmp(orig_dir+'/'+x, updated_dir+'/'+x):
            alias_changed = x
            alias_changed = alias_changed.replace('.txt','')
            if os.path.exists(geteff_filepath + "/"+ alias_changed +".txt"):
                #var_userinput= raw_input ('Remove ' + geteff_filepath + '/GetEffLumi_'+ alias_changed +'.root since ForSNU file changed  [y/n]:')
                print 'Remove ' + geteff_filepath + '/'+ alias_changed +'.txt since ForSNU file changed'
                #if var_userinput == "y" :
                os.system("rm " + geteff_filepath + "/"+ alias_changed +".txt")
                print "###############################"*2
                print "removing " + geteff_filepath + "/"+ alias_changed +"."
                print "###############################"*2
                
                
    
def check_file_diff(array_from_googledoc,_era, _dsn,  _path_to_skflat_mc, _path_to_samplefiles,_var_skflat_wd,_skoutput, work_dir, RunFull):
    
    from os import listdir
    from os.path import isfile,isdir, join
    from subprocess import check_output
    print_message(1,"check_file_diff ["+_era+"] [start]")

    """ loop over dataset list on tamsa"""
    #var_url = get_url_from_era(_era,False)
    
    path_mc_era=_path_to_skflat_mc+_era+"/MC/"
    arr_tamsadir  = [f for f in listdir(path_mc_era) if isdir(join(path_mc_era,f))]
    arr_snufiles  = [f for f in listdir(_path_to_samplefiles) if isfile(join(_path_to_samplefiles,f))]

    hasBadFile=False

    for _dir in _dsn:
        
        array_gd = array_from_googledoc

        var_alias = find_googledoc_var_from_dsn(array_gd,_era,"alias", _dir)
        if var_alias == "NULL" :
            continue

        runBadFile=False
        
        if var_alias+".txt" in arr_snufiles:
            path_mc_era_dir=path_mc_era+_dir
            arr_timestamp  = [f for f in listdir(path_mc_era_dir) if isdir(join(path_mc_era_dir,f))]
            timestamp=[]
            if len(arr_timestamp) > 0:
                #timestamp=arr_timestamp[0]
                timestamp=find_timestamp(_path_to_samplefiles+"/"+var_alias+".txt" , arr_timestamp)

            for _ts in timestamp:
                source = path_mc_era_dir+"/"+ _ts
                target = _path_to_samplefiles + var_alias+".txt"
                os.system("chmod u+x  script/diff_script.sh")
                print ("bash script/diff_script.sh " + target + " " + source+" " + var_alias+ " " + work_dir)
                os.system("bash script/diff_script.sh " + target + " " + source+" " + var_alias+ " " + work_dir)


            if  os.path.getsize(work_dir+"/tmp"+var_alias+".txt") > 1:
                print "check_file_diff: "+ var_alias+".txt has changed size... will rerun BadFileChecker and Effective luminosity" 

                runBadFile=True
            os.system("rm "+work_dir+"/tmp"+var_alias+".txt")
        else:
            print "check_file_diff: "+var_alias+".txt is new, will check for broken or empty files..."
            runBadFile=True
            
        
        if runBadFile:
            hasBadFile=True
            if os.path.exists(_skoutput+"/"+ var_alias +".txt"):
                os.system("rm " + _skoutput+"/"+ var_alias +".txt")
                print "###############################"*2
                print "removing " + _skoutput+"/"+ var_alias +".txt"
                print "###############################"*2

            check_bad_files(_var_skflat_wd+"/script/BadFileChecker/",_era, _dir, RunFull)                                                                                                                                                               

    if not hasBadFile:
        print 'no bad/corrupt files found....'

    print_message(1,"check_file_diff ["+_era+"] [end]")
    return


def get_alias_perera(_era, dirlist):

    arr_alias=[]
    var_url = get_url_from_era(_era,False)

    print "List of MC in ["+_era+"]"
    for x in dirlist:
        var_url = get_url_from_era(_era,False)
        if find_googledoc_var_from_dsn(_era,"alias", x, var_url) == "NULL":
            var_url = get_url_from_era(_era,True)

        var_alias = find_googledoc_var_from_dsn(_era,"alias", x, var_url)
        if var_alias == "NULL" :
            print "Skipping " +var_alias
            continue

        arr_alias.append(var_alias)
    return arr_alias

def get_missing_alias_perera(_era, dirlist, skimname):

    print "get_missing_alias_perera: era = ["+_era+"]"
    arr_alias=[]
    var_url = get_url_from_era(_era,False)

    skim_dir="/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC_"+ skimname
    if not os.path.exists(skim_dir):
        exit()


        
    from os import listdir
    from os.path import isfile,isdir, join
    arr_skimsamples  = [f for f in listdir(skim_dir) if isdir(join(skim_dir,f))]

    for _alias in dirlist:
        if not _alias in arr_skimsamples:
            print _alias



def print_alias_perera(_era, dirlist):

    var_url = get_url_from_era(_era,False)

    print "List of MC in ["+_era+"]"
    for x in dirlist:
        var_url = get_url_from_era(_era,False)
        if find_googledoc_var_from_dsn(_era,"alias", x, var_url) == "NULL":
            var_url = get_url_from_era(_era,True)

        var_alias = find_googledoc_var_from_dsn(_era,"alias", x, var_url)

        print var_alias
