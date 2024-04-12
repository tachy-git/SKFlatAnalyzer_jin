import os
import logging
import time, datetime
from CheckJobStatus import CheckJobStatus
from TimeTools import GetDatetimeFromMyFormat
from SendEmail import SendEmail


class CondorJobHandler:
    def __init__(self, processor, log_email):
        self.processor = processor
        self.log_email = log_email
        self.to_status_log = []
        self.running = []
        self.finished = []
        self.evt_done = 0
        self.evt_total = 0
        self.total_evt_runtime = 0
        self.max_evt_runtime = 0
        self.max_time_left = 0
        self.err_log = []
        self.jobbatchname = ""
    
    def reset(self):
        self.to_status_log = []
        self.running = []
        self.finished = []
        self.evt_done = 0
        self.evt_total = 0
        self.total_evt_runtime = 0
        self.max_evt_runtime = 0
        self.max_time_left = 0
        self.err_log = []
        self.jobbatchname = ""
    
    def monitorJobStatus(self):
        # self.isDone = False
        # is not done in the previous monitoring
        # monitor again
        status_log = open(f"{self.processor.baseRunDir}/JobStatus.log", "w")
        status_log.write(f"Job submitted at {self.processor.jobstarttime}\n")
        status_log.write(f"Job ID\t| Status\n")
        self.reset()
        
        ## loop over all jobs
        for iproc, _ in enumerate(self.processor.fileRanges):
            status = CheckJobStatus(self.processor.baseRunDir, 
                                    self.processor.analyzer, 
                                    iproc, 
                                    self.processor.hostname.upper())
            if "RUNNING" in status:
                self.running.append(iproc)
                out_log = self.RUNNING(iproc, status)
            elif "FINISHED" in status:
                self.finished.append(iproc)
                out_log = self.FINISHED(status)
            elif "ERROR" in status:
                status_log.write("#### ERROR OCCURED ####\n")
                out_log = self.ERROR(status)
            else:
                out_log = f"{iproc}\t| {status}"
            logging.debug(f"{self.processor.sampleName} - {status}: {out_log}")
            self.to_status_log.append(out_log)
        
        for line in self.to_status_log:
            status_log.write(f"{line}\n")
        
        status_log.write("\n=====================================\n")
        status_log.write(f"HOSTNAME = {self.processor.hostname}\n")
        status_log.write(f"{len(self.processor.fileRanges)} job submitted\n")
        status_log.write(f"{len(self.running)} jobs are running\n")
        status_log.write(f"{len(self.finished)} jobs are finished\n")
        status_log.write(f"XSEC = {self.processor.xsec}S\n")
        status_log.write(f"eventDone = {self.evt_done}\n")
        status_log.write(f"eventTotal = {self.evt_total}\n")
        status_log.write(f"eventLeft = {self.evt_total-self.evt_done}\n")
        status_log.write(f"totalEventRunTime = {self.total_evt_runtime}\n")
        status_log.write(f"maxTimeLeft = {self.max_time_left}\n")
        status_log.write(f"maxEventRunTime = {self.max_evt_runtime}\n")
        
        time_per_evt = 1
        if self.evt_done:    # exist finished events
            time_per_evt = float(self.total_evt_runtime)/self.evt_done
        status_log.write(f"Estimated time per event = {time_per_evt}\n")
        finish_time = datetime.datetime.now() + datetime.timedelta(seconds=self.max_time_left)
        status_log.write(f"Estimated finish time = {finish_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        status_log.write(f"Last checked at {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        status_log.close()
        
    def postProcess(self, final_output_path):
        if self.err_log:
            logging.error(f"Error in {self.processor.sampleName}")
            self.sendErrorMail()
            self.processor.isError = True
        elif len(self.processor.fileRanges) == len(self.finished):
            # no need to hadd files if skimming
            if self.processor.skim:
                logging.info(f"Skimmed output files are ready in {final_output_path}")
            else:    
                self.preparePostProcessing(final_output_path)
                cwd = os.getcwd()
                os.chdir(self.processor.baseRunDir)
                os.system(f"condor_submit postprocess.sub")
                logging.info(f"Postprocessing submitted for {self.jobbatchname}")
                os.chdir(cwd)
            self.processor.isDone = True
        else:
            # Job is still running
            # do noting in this stage
            pass
    
    def sendErrorMail(self):
        title = f"[ERROR] Summary of JobID {self.processor.job_id}"
        email = f"""
#### ERROR OCCURRED ####
{self.err_log[0]}
        
#### Job Info ####
HOST = {self.processor.hostname}
JobID = {self.processor.job_id}
Analyzer = {self.processor.analyzer}
Era = {self.processor.era}
Skim = {self.processor.skim}
# of jobs = {self.processor.njobs}
input sample = {self.processor.sampleName}
xsec = {self.processor.xsec}"""
        SendEmail(os.environ['USER'], self.log_email, title, email)
    
    def sendFinishMail(self):
        title = f"[{self.processor.hostname}] Summary of JobID {self.processor.job_id}"
        email = f"""
#### Job Info ####
HOST = {self.processor.hostname}
JobID = {self.processor.job_id}
Analyzer = {self.processor.analyzer}
Era = {self.processor.era}
Skim = {self.processor.skim}
# of jobs = {self.processor.njobs}
input sample = {self.processor.sampleName}
xsec = {self.processor.xsec}"""
        SendEmail(os.environ['USER'], self.log_email, title, email) 
    
    def preparePostProcessing(self, final_output_path):
        with open(f"script/Templates/PostProcess/condor.sub", "r") as f:
            template = f.read()
        with open(f"{self.processor.baseRunDir}/postprocess.sub", "w") as f:
            self.jobbatchname = f"{self.processor.era}_{self.processor.analyzer}_{self.processor.sampleName}"
            if self.processor.isDATA:
                self.jobbatchname += f"_{self.processor.dataPeriod}"
            template = template.replace("[SAMPLENAME]", self.jobbatchname)
            f.write(template)
        
        output_name = f"{self.processor.analyzer}_{self.processor.sampleName}"
        if self.processor.isDATA:
            output_name += f"_{self.processor.dataPeriod}"
        with open(f"script/Templates/PostProcess/hadd.sh", "r") as f:
            template = f.read()
        with open(f"{self.processor.baseRunDir}/hadd.sh", "w") as f:
            template = template.replace("[BASERUNDIR]", self.processor.baseRunDir)
            template = template.replace("[OUTPUTNAME]", output_name)
            template = template.replace("[FINALOUTPUTPATH]", final_output_path)
            f.write(template)
        # change permission to 755 using python module
        os.chmod(f"{self.processor.baseRunDir}/hadd.sh", 0o755)
    
    def monitorPostProcess(self):
        # check if hadd.log exists
        if os.path.exists(f"{self.processor.baseRunDir}/hadd.log"):
            logging.info(f"Postprocessing finished for {self.processor.sampleName}")
            time.sleep(1)
            while True:
                try:
                    with open(f"{self.processor.baseRunDir}/hadd.err", "r") as f:
                        for line in f.readlines():
                            if "WARNING" in line:
                                continue
                            else:
                                logging.error(f"{self.processor.sampleName} - {line}")
                                self.err_log.append(line)
                    break
                except:
                    continue
            
            if self.err_log:
                self.sendErrorMail()
                self.processor.isError = True
            else:    
                self.sendFinishMail()
                self.processor.isPostJobDone = True
        
    def RUNNING(self, iproc, status):
        out_log = f"{iproc}\t| {status.split()[1]} %"
        evt_info = status.split()[2].split(":")
        proc_evt_done, proc_evt_total = int(evt_info[1]), int(evt_info[2])
        self.evt_done += proc_evt_done
        self.evt_total += proc_evt_total
        
        line_evt_runtime = f"{status.split()[3]} {status.split()[4]}"
        proc_starttime = GetDatetimeFromMyFormat(line_evt_runtime) 
        diff = datetime.datetime.now() - proc_starttime
        proc_evt_runtime = diff.total_seconds()
        
        if proc_evt_done == 0: 
            proc_evt_done = 1
        proc_time_per_evt = float(proc_evt_runtime)/proc_evt_done
        proc_time_left = (proc_evt_total - proc_evt_done)*proc_time_per_evt
        
        self.total_evt_runtime += proc_evt_runtime
        self.max_time_left = max(self.max_time_left, proc_time_left)
        self.max_evt_runtime = max(self.max_evt_runtime, proc_evt_runtime)
        
        out_log += f"({proc_time_left:.1}s ran, and {proc_evt_runtime}s left)"
        return out_log
    
    def FINISHED(self, status):
        evt_info = status.split()[1].split(":")
        proc_evt_done, proc_evt_total = int(evt_info[2]), int(evt_info[2])
        self.evt_done += proc_evt_done
        self.evt_total += proc_evt_total
        
        line_evt_runtime = f"{status.split()[2]} {status.split()[3]}"
        proc_starttime = GetDatetimeFromMyFormat(line_evt_runtime)
        line_evt_endtime = f"{status.split()[4]} {status.split()[5]}"
        proc_endtime = GetDatetimeFromMyFormat(line_evt_endtime)
        diff = proc_endtime - proc_starttime
        proc_evt_runtime = diff.total_seconds()
        
        proc_time_per_evt = float(proc_evt_runtime)/proc_evt_done
        proc_time_left = (proc_evt_total - proc_evt_done)*proc_time_per_evt
        
        self.total_evt_runtime += proc_evt_runtime
        self.max_time_left = max(self.max_time_left, proc_time_left)
        self.isDone = True
        return ""
    
    def ERROR(self, status):
        self.err_log.append(status)
        return ""
