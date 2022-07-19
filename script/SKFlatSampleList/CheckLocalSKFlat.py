import os
from os import listdir
from os.path import isfile,isdir, join

def get_dsn_list(_mclist_path):
    
    arr_dsn=[]
    r_localfile = open(_mclist_path, "r")

    for line in _rlocalfile:
        line_blocks = line.split()
        if  not len(line_blocks) == 3:
            continue
            arr_dsn.append(line_blocks[1])

    _localfile.close()
    return arr_dsn





def get_skflat_on_tamsa(_path,_era):
    
    path_era=_path+_era+"/MC/"
    arr_skflatdir = [f for f in listdir(path_era) if isdir(join(path_era,f))]
    arr_skflatdir.sort()


    skip_dsn_list=open(os.getenv("SKFlat_WD")+"/script/SKFlatSampleList/config/DSN_Skip.txt","r")

    for line in skip_dsn_list:
        if len(line.split()) ==  2:
            if line.split()[0] == "201*" or line.split()[0] ==  _era:
                arr_skflatdir.remove(line.split()[1] )
    skip_dsn_list.close()
    
    return arr_skflatdir
    
        
                
    
