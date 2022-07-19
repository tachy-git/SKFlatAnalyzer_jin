Usage:
--------------------------------------------
Run all eras:

python  UpdateTamsaSKFlatSampleInput.py

--------------------------------------------

Run all eras and check files are ok and not changed from previous run:

python  UpdateTamsaSKFlatSampleInput.py --Full

--------------------------------------------

-s option runs pre coded skims for samples
python  UpdateTamsaSKFlatSampleInput.py -s skimname 
--------------------------------------------

python  UpdateTamsaSKFlatSampleInput.py -e 'era' 

--------------------------------------------

usage: UpdateTamsaSKFlatSampleInput.py [-h] [-e ERA] [-s SKIM] [--Full]
                                       [--UserInput] [--VERBOSE]
option

optional arguments:
  -h, --help   show this help message and exit
  -e ERA         [ALLOW JUST ONE ERA TO BE RUN]
  -s SKIM	 [RUNS SKIM ON NEW SAMPLES IF SET IN GOOGLEDOC]
  --Full	 [THIS RUNS CHECK ON PREXISTING SAMPLES]
  --UserInput    [THIS MEANS USER MUST SPECIFY WHICH SAMPLE TO KEEP IN COMMAND LINE IF MORE TAHN ONE TIMESTAMP EXISTS]
  --VERBOSE      [FULL PRINT OUT DEBUG MODE]
