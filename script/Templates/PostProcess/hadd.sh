#!/bin/bash
BASERUNDIR=[BASERUNDIR]

cd $BASERUNDIR

NFILES=`ls -1 output | wc -l`

if [ "$NFILES" -eq 1 ]; then
    echo "$NFILES = 1, so skipping hadd and just move the file" >> JobStatus.log
    mv output/hists_0.root [OUTPUTNAME].root
else
    hadd -f [OUTPUTNAME].root output/hists_*.root >> JobStatus.log
    rm output/*.root
    mv [OUTPUTNAME].root [FINALOUTPUTPATH]
fi

echo "FINISHED hadd and moved the final output to [FINALOUTPUTPATH]"
