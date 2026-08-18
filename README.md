# 3DPMTE-SC_Interfaces
## Abaqus-python, 3D PMTE-SC code, with user interface subroutine for linear elastic (brittle) material interfaces.

In the context of Finite Fracture Mechanics, the classical formulation of the coupled criterion for predicting the failure load in a notched or un-notched composite component with, say, a constant thickness $$B$$, assuming the nucleation of a single mode-I crack emanating from the notch root or adhesive front, is given by the expression:

$$
\mathbf{
\begin{array}{c}\sigma(x)  \geq  \sigma_c,\\
B\int_0^{\Delta a}G(x)dx \geq B\int_0^{\Delta a} G_{Ic} dx = G_{Ic}B\Delta a,
\end{array}}
$$

where the system of two nonlinear equations *must be solved simultaneously for (i) the failure load and (ii) the associated crack nucleation length $$\Delta a$$.*

## Files:
1. PMTESC3D_DnBisectionLoadP_subprocess.py: Python script to be launched from Abaqus CAE or from the CMD; this file is tuned to run in a desktop machine with Windows OS. To launch from CMD prompt: abaqus cae noGUI=PMTESC3D_DnBisectionLoad_subprocess.py

3. PMTESCabaqusV0.py: script containing functions which are called by PMTESC3D_DnBisectionLoad_subprocess.py after each simulation run for data post-processing and total energy (Pi+DeltaR) calculations for the various starting damage configurations.

4. The Fortran file "UINTERLEBIM3D_kincxit.for" includes a URDFIL user subroutine for reading the energy records at specific points in the analysis procedure, the .inp file requires the following keyword for the efficient solver termination call inside the routine (the CALL XITT() command):

  ** OUTPUT REQUESTS
  
  ** 
  
  ** This block is for the URDFIL subroutine
  
  *ENERGY FILE, FREQUENCY=1
  **

5. **"UINTERLEBIM3D_kincxitV2.for"** is the latest and more efficient user subroutine for 3D interfaces, as the convergence termination criteria is given by the (absolute) change of the damage state (i.e., change of broken springs) between the current and previous increments in a given load step: no need to track energy variables with the URDFIL routine across iterations, thus eliminating the request for the abaqus .fil file and reading of energy records, making this routine appropriate for modeling brittle fracture of large scale industrial models.
