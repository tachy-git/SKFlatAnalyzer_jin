#!/bin/bash
BASERUNDIR=[BASERUNDIR]

cd $BASERUNDIR

NFILES=`ls -1 output | wc -l`

if [ "$NFILES" -eq 1 ]; then
    echo "$NFILES = 1, so skipping hadd and just move the file" >> JobStatus.log
    mv output/hist_0.root [OUTPUTNAME].root
else
    hadd -f [OUTPUTNAME].root output/hists_*.root >> JobStatus.log
    rm output/*.root
fi

# Move the final output
mv [OUTPUTNAME].root [FINALOUTPUTPATH]
echo "FINISHED hadd and moved the final output to [FINALOUTPUTPATH]"