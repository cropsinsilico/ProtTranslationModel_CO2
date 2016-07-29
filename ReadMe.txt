Gene regulatory network model for photosynthesis driver:

----------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------
Linux:
For Biocluster only: Get out of head node: qsub -I -X
cd to the GrCM/ directory
----------------------------------------------------------------------------------------------------------
Method 1: Directly calling using the executable
Run:
Execute the program: src/GrMC.py ./Input/GrCM_input.txt ./Output/GrCM_output.txt &> ./Output/GrCM_log.txt
----------------------------------------------------------------------------------------------------------
Method 2: Using a bash script
cd to the GrCM/ directory
Run the script: 
. ../../raila/Plants_in_Silico/RMQ-support/interface/setup.sh
../../raila/Plants_in_Silico/RMQ-support/interface/PsiRun.py GrCM.yml > ./Output/GrCM_log.txt

----------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------
Windows:
To be added later
