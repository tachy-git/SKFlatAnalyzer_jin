import os

os.system('ls -1 CommonSampleInfo/*.txt &> temp.txt');
lines = open('temp.txt').readlines()
os.system('rm temp.txt')

out_MC = open('SampleSummary_MC.txt','w')
out_Sig_Type1 = open('SampleSummary_Signal_Type1.txt','w')

FirstLine = '# alias PD xsec nmc sumsign sumw\n'
out_MC.write(FirstLine)
out_Sig_Type1.write(FirstLine)

for line in lines:

  line = line.strip('\n')
  alias = line.replace('CommonSampleInfo/','').replace('.txt','')

  infoline = open(line).readlines()[1]

  if "Type" in alias or "private" in alias:
    out_Sig_Type1.write(infoline)
  else:
    out_MC.write(infoline)

out_MC.close()
out_Sig_Type1.close()
