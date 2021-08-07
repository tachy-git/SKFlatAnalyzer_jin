import os,ROOT
histmaps=os.popen("find $SKFlat_WD/data/$SKFlatV -name 'histmap*.txt'").readlines()
for histmap in histmaps:
    histmap=histmap.strip()
    lines=open(histmap).readlines()
    for line in lines:
        line=line.strip()
        line=line.split("#",1)[0]
        if line=="": continue
        words=line.split()
        if "PileUp" in histmap:
            if len(words)<3:
                print histmap+": Wrong format"
                print words
                continue
            filename=words[2]
            histname=words[0]+"_"+words[1]
        else:
            if len(words)<5:
                print histmap+": Wrong format"
                print words
                continue
            filename=words[3]
            histname=words[4]
        f=ROOT.TFile(os.path.dirname(histmap)+"/"+filename)
        h=f.Get(histname)
        if not h:
            print ":".join([histmap,filename,histname])+": Not found"

    
