import os
import shutil
import logging
import importlib.util

class SampleListHandler:
    inputDataSampleList = ["DoubleMuon", "DoubleEG", "SingleMuon", "SingleElectron",
                           "SinglePhoton", "MuonEG", "EGamma"]
    availableDataPeriods = {
        "2016preVFP": ["B_ver2", "C", "D", "E", "F"],
        "2016postVFP": ["F", "G", "H"],
        "2017": ["B", "C", "D", "E", "F"],
        "2018": ["A", "B", "C", "D"]
    }
    def __init__(self, era, dataPeriod):
        self.era = era
        self.dataPeriod = dataPeriod
        
    def generateSampleListFromInputSample(self, inputSample):
        sampleList = []
        stringForHash = ""
        # data
        if inputSample in self.inputDataSampleList:
            if self.dataPeriod == "ALL":
                for period in self.availableDataPeriods[self.era]:
                    sampleList.append(f"{inputSample}:{period}")
                    stringForHash += f"{inputSample}:{period}"
            elif self.dataPeriod in self.availableDataPeriods[self.era]:
                sampleList.append(f"{inputSample}:{self.dataPeriod}")
                stringForHash += f"{inputSample}:{self.dataPeriod}"
            else:
                raise ValueError(f"Wrong data period: {self.dataPeriod}")
        # MC
        else:
            sampleList.append(inputSample)
            stringForHash += inputSample

        return (sampleList, stringForHash)
    
    def generateSampleListFromInputSampleList(self, inputSampleList):
        sampleList = []
        stringForHash = ""
        with open(inputSampleList, "r") as f:
            for line in f:
                if line.startswith("#"): continue
                inputSample = line.strip("\n")
                thisSampleList, thisStringForHash = self.generateSampleListFromInputSample(inputSample)
                sampleList += thisSampleList
                stringForHash += thisStringForHash
        return (sampleList, stringForHash)
            

class SampleProcessor:
    def __init__(self, sample, args, ENVs):
        self.job_id = ENVs["JobID"]
        self.hostname = ENVs["HOSTNAME"]
        self.masterJobDir = ENVs["MasterJobDir"]
        self.SKFlat_WD = ENVs["SKFlat_WD"]
        self.SKFlatV = ENVs["SKFlatV"]
        self.SAMPLE_DATA_DIR = ENVs["SAMPLE_DATA_DIR"]
        self.jobstarttime = ENVs["jobstarttime"]
        self.analyzer = args.analyzer
        self.era = args.era
        self.njobs = args.njobs
        self.nmax = args.nmax
        self.memory = args.memory
        self.userflags = args.userflags
        self.skim = args.skim
        self.reduction = args.reduction
        self.outputdir = args.output_dir
        self.timestamp = ENVs["timestamp"]
        self.python = args.python
        self.isDATA = ":" in sample
        self.lhapdfpath = f"{ENVs['SKFlat_WD']}/external/lhapdf/lib/libLHAPDF.so"
        self.sampleName = ""
        self.dataPeriod = ""

        if self.isDATA:
            self.sampleName, self.dataPeriod = sample.split(":")
            self.baseRunDir = f"{self.masterJobDir}/{self.sampleName}_period{self.dataPeriod}"
        else:
            self.sampleName = sample
            self.baseRunDir = f"{self.masterJobDir}/{self.sampleName}" 
        
        self.totalFiles = []            
        self.fileRanges = []
        self.dasName = ""
        self.xsec = -1
        self.sumSign = -1
        self.sumW = -1
        
        # for monitoring
        self.isDone = False
        self.isPostJobDone = False
        self.isError = False
    
    def prepareRunDirectory(self):
        os.makedirs(f"{self.baseRunDir}/output")
        ## get sample path
        prefix = ""
        suffix = ""
        if self.skim:
            prefix = f"{self.skim}_"
        if self.isDATA:
            suffix = f"_{self.dataPeriod}"
        samplePath = f"{self.SAMPLE_DATA_DIR}/ForSNU/{prefix}{self.sampleName}{suffix}.txt" 
        shutil.copy(samplePath, f"{self.baseRunDir}/input_filelist.txt")
        
        ## calculate how many files to be splitted into each jobs
        self.totalFiles = os.popen("sed 's/#.*//' "+samplePath+"|grep '.root'").readlines()
        nTotalFiles = len(self.totalFiles)
        if self.njobs > nTotalFiles or self.njobs == 0:
            self.njobs = nTotalFiles
        nfileperjob = int(nTotalFiles/self.njobs)
        remainder = nTotalFiles - self.njobs*nfileperjob

        tmpEndLargerJob = 0
        nfileCheckSum = 0
        # first remainder jobs will have (nfileperjob+1) files per job
        for ijob in range(remainder):
            self.fileRanges.append(range(ijob*(nfileperjob+1), (ijob+1)*(nfileperjob+1)))
            tmpEndLargerJob = (ijob+1)*(nfileperjob+1)
            nfileCheckSum += len(range(ijob*(nfileperjob+1), (ijob+1)*(nfileperjob+1)))
    
        ## Remaining njobs - remainder jobs will have (nfileperjob) files per job
        for ijob in range(self.njobs - remainder):
            self.fileRanges.append(range(tmpEndLargerJob+(ijob*nfileperjob), tmpEndLargerJob+((ijob+1)*(nfileperjob))))
            nfileCheckSum += len(range(tmpEndLargerJob+(ijob*nfileperjob), tmpEndLargerJob+((ijob+1)*(nfileperjob))))
        logging.debug(f"prepare run directory for {self.sampleName}")
        logging.debug(f"samplePath: {samplePath}")
        logging.debug(f"nTotalFiles: {nTotalFiles}")
        logging.debug(f"njobs: {self.njobs}")
        logging.debug(f"nfileperjob: {nfileperjob}")
        logging.debug(f"remainder: {remainder}")
        logging.debug(f"fileRanges: {self.fileRanges}")
        logging.debug(f"nfileCheckSum: {nfileCheckSum}")
        assert nfileCheckSum == nTotalFiles
        
        ## Get xsec and sumW
        if not self.isDATA:
            sampleInfoPath = f"{self.SAMPLE_DATA_DIR}/CommonSampleInfo/{self.sampleName}.txt"
            if not os.path.exists(sampleInfoPath):
                raise FileNotFoundError(f"Sample info file not found: {sampleInfoPath}")
            with open(sampleInfoPath, "r") as f:
                for line in f:
                    if line.startswith("#"): continue
                    words = line.split()
                    thisSample, self.dasName, self.xsec, _, self.sumSign, self.sumW = tuple(words)
                    assert self.sampleName == thisSample
                    break
            
    def generateSubmissionScripts(self):
        ## write run script
        commandsFileName = f"{self.analyzer}_{self.era}_{self.sampleName}"
        if self.isDATA:
            commandsFileName += f"_{self.dataPeriod}"
        if self.userflags:
            for flag in self.userflags:
                commandsFileName += f"__{flag}"
        if self.python:
            with open(f"script/Templates/Submission/run.python.sh", "r") as f:
                template = f.read()
        else:
            with open(f"script/Templates/Submission/run.sh", "r") as f:
                template = f.read()
        with open(f"{self.baseRunDir}/{commandsFileName}.sh", "w") as f:
            template = template.replace("[SKFlat_WD]", self.SKFlat_WD)
            template = template.replace("[baseRunDir]", self.baseRunDir)
            template = template.replace("[masterJobDir]", self.masterJobDir)
            f.write(template)
            
        ## write condor submission script
        concurrencyLimits = ""
        request_memory = ""
        if self.nmax:
            concurrencyLimits = f"concurrencyLimits = n{self.nmax}.{os.getenv('USER')}"
        with open(f"script/Templates/Submission/condor.sub", "r") as f:
            template = f.read()
        with open(f"{self.baseRunDir}/condor.sub", "w") as f:
            template = template.replace("[commandsFileName]", commandsFileName)
            template = template.replace("[concurrencyLimits]", concurrencyLimits)
            template = template.replace("[request_memory]", str(self.memory))
            template = template.replace("[njobs]", str(self.njobs))
            f.write(template)
            
        ## runners
        totalSampleCounter = 0
        for ijob in range(len(self.fileRanges)):
            totalSampleCounter += len(self.fileRanges[ijob])
            runScriptName = f"run_{ijob}"
            runScriptFullPath = f"{self.baseRunDir}/{runScriptName}.py" if self.python else f"{self.baseRunDir}/{runScriptName}.C"
            lhapdfpath = self.lhapdfpath
            if self.python:
                shutil.copy(f"{self.SKFlat_WD}/PyAnalyzers/{self.analyzer}.py", runScriptFullPath)
                # add main function
                out = open(runScriptFullPath, "a")
                out.write("\n\n\n")
                out.write('if __name__ == "__main__":\n')
                out.write(f'    gSystem.Load("{lhapdfpath}")\n')
                out.write(f"    m = {self.analyzer}()\n")
                out.write(f'    m.SetTreeName("recoTree/SKFlat")\n')
                if self.isDATA:
                    out.write(f'    m.IsDATA = True\n')
                    out.write(f'    m.DataStream = "{self.sampleName}"\n')
                else:
                    out.write(f'    m.IsDATA = False\n')
                    out.write(f'    m.MCSample = "{self.sampleName}"\n')
                    out.write(f'    m.xsec = {self.xsec}\n')
                    out.write(f'    m.sumSign = {self.sumSign}\n')
                    out.write(f'    m.sumW = {self.sumW}\n')
                    out.write(f"    m.IsFastSim = False\n")
                out.write(f'    m.SetEra("{self.era}")\n')
                if self.userflags:
                    out.write(f'    m.Userflags = vector[TString]()\n')
                    for flag in self.userflags: 
                        out.write(f'    m.Userflags.emplace_back("{flag}")\n')
                for ifile in self.fileRanges[ijob]:
                    thisFileName = self.totalFiles[ifile].strip("\n")
                    out.write(f'    if not m.AddFile("{thisFileName}"): exit(1)\n')
                if "SkimTree" in self.analyzer:
                    skimOutDir = f"/gv0/DATA/SKFlat/{self.SKFlatV}/{self.era}"
                    if self.outputdir:
                        skimOutDir = f"{self.outputdir}/{self.SKFlatV}/{self.era}"
                    skimOutFileName = ""
                    if self.isDATA:
                        skimOutDir = f"{skimOutDir}/DATA_{self.analyzer}/{self.sampleName}/period{self.dataPeriod}"
                        skimOutFileName = f"SKFlatNtuple_{self.era}_DATA_{ijob}.root"
                    else:
                        skimOutDir = f"{skimOutDir}/MC_{self.analyzer}/{self.dasName}"
                        skimOutFileName = f"SKFlatNtuple_{self.era}_MC_{ijob}.root"
                    skimOutDir = f"{skimOutDir}/{self.timestamp}"
                    os.makedirs(skimOutDir)
                    out.write(f'    m.SetOutfilePath("{skimOutDir}/{skimOutFileName}")\n')
                else:
                    out.write(f'    m.SetOutfilePath("hists.root")\n')
                if self.reduction > 1:
                    out.write(f"    m.MaxEvent = int(m.fChain.GetEntries()/{self.reduction})\n")
                out.write(f"    m.Init()\n")
                out.write(f"    m.initializePyAnalyzer()\n")
                out.write(f"    m.initializeAnalyzerTools()\n")
                out.write(f"    m.SwitchToTempDir()\n")
                out.write(f"    m.Loop()\n")
                out.write(f"    m.WriteHist()\n")
                #out.write(f"    del m\n")
                out.close()
            else:
                out = open(runScriptFullPath, "w")
                out.write(f"R__LOAD_LIBRARY({lhapdfpath})\n")
                out.write(f"void {runScriptName}()" + "{\n")
                out.write(f"    {self.analyzer} m;\n")
                out.write(f'    m.SetTreeName("recoTree/SKFlat");\n')
                if self.isDATA:
                    out.write(f"    m.IsDATA = true;\n")
                    out.write(f'    m.DataStream = "{self.sampleName}";\n')
                else:
                    out.write(f'    m.MCSample = "{self.sampleName}";\n')
                    out.write(f"    m.IsDATA = false;\n")
                    out.write(f"    m.xsec = {self.xsec};\n")
                    out.write(f"    m.sumSign = {self.sumSign};\n")
                    out.write(f"    m.sumW = {self.sumW};\n")
                    out.write(f"    m.IsFastSim = false;\n")
                out.write(f'    m.SetEra("{self.era}");\n')
                if self.userflags:
                    for flag in self.userflags: 
                        out.write(f'    m.Userflags.emplace_back("{flag}");\n')
                for ifile in self.fileRanges[ijob]:
                    thisFileName = self.totalFiles[ifile].strip("\n")
                    out.write(f'    if(!m.AddFile("{thisFileName}")) exit(EIO);\n')
                if "SkimTree" in self.analyzer:
                    skimOutDir = f"/gv0/DATA/SKFlat/{self.SKFlatV}/{self.era}"
                    if self.outputDir:
                        skimOutDir = f"{self.outputDir}/{self.SKFlatV}/{self.era}"
                    skimOutFileName = ""
                    if self.isDATA:
                        skimOutDir = f"{skimOutDir}/DATA_{self.analyzer}/{self.sampleName}/period{self.dataPeriod}"
                        skimOutFileName = f"SKFlatNtuple_{self.era}_DATA_{ijob}.root"
                    else:
                        skimOutDir = f"{skimOutDir}/MC_{self.analyzer}/{self.dasName}"
                        skimOutFileName = f"SKFlatNtuple_{self.era}_MC_{ijob}.root"
                    skimOutDir = f"{skimOutDir}/{self.timestamp}"
                    os.makedirs(skimOutDir, exist_ok=True)
                    out.write(f'    m.SetOutfilePath("{skimOutDir}/{skimOutFileName}");\n')
                else:
                    out.write(f'    m.SetOutfilePath("hists.root");\n')
                if self.reduction > 1:
                    out.write(f"    m.MaxEvent = m.fChain->GetEntries()/{self.reduction};\n")
                out.write(f"    m.Init();\n")
                out.write(f"    m.initializeAnalyzer();\n")
                out.write(f"    m.initializeAnalyzerTools();\n")
                out.write(f"    m.SwitchToTempDir();\n")
                out.write(f"    m.Loop();\n")
                out.write(f"    m.WriteHist();\n")
                out.write("}")
                out.close()
    
    def submitJobs(self, batchname=""):
        cwd = os.getcwd()
        os.chdir(self.baseRunDir)
        condorOptions = ""
        if batchname:
            condorOptions = f"-batch-name {batchname}"
        os.system(f"condor_submit {condorOptions} condor.sub")
        os.chdir(cwd)
