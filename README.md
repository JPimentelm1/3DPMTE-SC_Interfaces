# 3DPMTE-SC_Interfaces
## Abaqus-python, 3D PMTE-SC code, and user interface subroutine for linear elastic (brittle) interfaces.

## 1. PMTESTmain.py: python script to be launched from Abaqus CAE or from the CMD. If launched from CMD, type: abaqus cae noGUI=PMTESTmain.py

## 2. PMTESCabaqus.py: script containing functions that are called by the PMTESTmain.py script during the simulation for data post-processing and total energy (Pi+DeltaR) calculations for the various starting damage configurations.
