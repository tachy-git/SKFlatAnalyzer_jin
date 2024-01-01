from datetime import datetime, timedelta
import os

def PutColor(a,color):
  return '${'+color+'}'+a+'${Color_Off}'

USER = os.environ['USER']
if os.path.exists(os.environ['SKFlat_WD']+'/python/UserInfo_'+USER+'.py'):
  exec('from UserInfo_'+USER+' import *')
else:
  print("No UserInfo file")
  exit(1)

try:
    LogRemoval = int(UserInfo['LogRemoval'])
except:
    LogRemoval = -1

if (LogRemoval > 0) and (LogRemoval < 3):
    LogRemoval = 3 # safety cut in order to not remove currently running jobs after date change

SKFlat_WD=os.environ['SKFlat_WD']
TMPDIR = SKFlat_WD+'/tmp/'
os.system('mkdir -p '+TMPDIR)
out = open(TMPDIR+'/OldLogs.sh','w')

SKFlatRunlogDir = os.environ['SKFlatRunlogDir']
dirnames = os.listdir(SKFlatRunlogDir)
dirnamesToRemove = []

for dirname in dirnames:
    if ("__" in dirname):
        TimeStamp = dirname.split("__")[0].split("_")
        LogYear = int(TimeStamp[0])
        LogMonth = int(TimeStamp[1])
        LogDay = int(TimeStamp[2])

        LogCreated = datetime(LogYear, LogMonth, LogDay)
        TimeToday = datetime.now()

        TimeDelta = TimeToday - LogCreated
        if (LogRemoval > 0):
            if (TimeDelta > timedelta(days=LogRemoval)):
                dirnamesToRemove.append(dirname)

if (len(dirnamesToRemove) > 0):
    out.write('printf "\\n"\n')
    msgPrint = "Old log directories to remove in "+SKFlatRunlogDir
    out.write('printf "'+PutColor(msgPrint,'BYellow')+'"\n\n')
    out.write('printf "\\n\\n\\n"\n')
    for dirname in dirnamesToRemove:
        out.write('printf "'+PutColor(dirname,'Yellow')+'"\n')
        out.write('rm -r ' + SKFlatRunlogDir + '/' + dirname + '\n')
        out.write('printf "\\n"\n')
out.write('printf "\\n\\n"')
out.close()
