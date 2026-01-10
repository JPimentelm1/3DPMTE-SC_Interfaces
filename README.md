# 3DPMTE-SC_Interfaces
## Abaqus-python, 3D PMTE-SC code, and user interface subroutine for linear elastic (brittle) interfaces.

1. PMTESC3D_DnBisectionLoadP_subprocess.py: Python script to be launched with Abaqus CAE or from the CMD; this file is tuned to run in a desktop machine with Windows OS. To launch from CMD prompt: abaqus cae noGUI=PMTESC3D_DnBisectionLoad_subprocess.py

3. PMTESCabaqusV0.py: script containing functions which are called by PMTESC3D_DnBisectionLoadP_subprocess.py after each simulation run for data post-processing and total energy (Pi+DeltaR) calculations for the various starting damage configurations.

4. As the Fortran file includes a URDFIL user subroutine for reading the energy records at specific points in the analysis procedure, the .inp file requires the following keyword for the efficient solver termination call inside the routine (the CALL XITT() command):

  ** OUTPUT REQUESTS
  
  ** 
  
  ** This block is for the URDFIL subroutine
  
  *ENERGY FILE, FREQUENCY=1
  
  **

