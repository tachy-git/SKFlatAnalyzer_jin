import os, shutil
import pandas as pd

csv = pd.read_csv("modelInfo.csv", sep=",\s", engine="python")
for idx in csv.index:
    sig, bkg = csv.loc[idx, 'signal'], csv.loc[idx, 'background']
    model = csv.loc[idx, 'model']
    optim = csv.loc[idx, 'optimizer']
    initLR = csv.loc[idx, 'initLR']
    scheduler = csv.loc[idx, 'scheduler']

    modelPath = f"/data6/Users/choij/ChargedHiggsAnalysis/triLepRegion/full/Combine__/{sig}_vs_{bkg}/models/{model}_{optim}_initLR-{str(initLR).replace('.','p')}_{scheduler}.pt"
    
    # check path
    if os.path.exists(modelPath):
        shutil.copy(modelPath, f"models/{sig}_vs_{bkg}.pt")
    else:
        print(modelPath)
