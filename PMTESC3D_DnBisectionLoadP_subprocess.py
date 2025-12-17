# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 15:12:59 2023

@authors: Mar Munoz Reja Moreno;
          Jose M. Pimentel: load bisection algorithm, starting damage dist.
        implementation (uinter), Python's subprocess module for abaqus job cluster
        execution, jobs parallelization approach.
"""
from __future__ import print_function
from os import chdir, path, mkdir, getcwd, remove
import os
from shutil import rmtree, move,copy
import shutil
from glob import glob
from subprocess import call
import subprocess
import fileinput
import sys
import time
import collections
import pickle
import numpy as np
import math as mth
import stat
# Redirect stdout to the command prompt
sys.stdout = sys.__stdout__
import PMTESCabaqusV0

from abaqus import *
from abaqusConstants import *
from odbAccess import *
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import __main__
import regionToolset
import regionToolset
import shutil
import threading


#%%funciones
def replacement(archivo, previousw, nextw):
    for line in fileinput.input(archivo, inplace=1):
        line = line.replace(previousw, nextw)
        sys.stdout.write(line)
def lineasdefinteraccion(archivo, cadenaprimera, nfilas):
    abrirarchivo = open(archivo, "rt")
    cadenascambiare=[]
    comienzo=0
    for line in abrirarchivo:
    	# replacing the string and write to output file
        if line.startswith(cadenaprimera):
            cadenascambiare.append(line)
            comienzo=1
        elif comienzo>0 and comienzo<nfilas:
            cadenascambiare.append(line)
            comienzo=comienzo+1
    abrirarchivo.close()
    return(cadenascambiare)
    
def borrar_archivos(list_delete):
    for files in list_delete:
        files_delete = glob.glob(files)
        for k in range(len(files_delete)):
            remove(files_delete[k])
            
def mover_archivos(list_mover,folder):
    for files in list_mover:
        files_move = glob.glob(files)
        for k in range(len(files_move)):
            move(files_move[k],folder)
            
def model_outputfile(actual_directory,sim_file,K,Dn,newload,TF,enerHtotalN1,sumaenerinterN1,lNeL2damageKN1,lNeL2damagetotal,energy_evol):
    chdir(actual_directory)
    # print(os.getcwd(), sim_file)
    outfile=open(actual_directory+'/'+sim_file,'a')
    outfile.write(str('%d\t %d\t %.8e\t %.8e\t %.8e\t %.8e\t %d\t'%(K,Dn,newload,TF,enerHtotalN1,sumaenerinterN1,lNeL2damageKN1)))
    outfile.write(str('%d\t'%(lNeL2damagetotal)))
    outfile.write(str('%d\t %.8e\t' %(energy_evol[-1][0],energy_evol[-1][5]))+'\n')
    
    outfile.close()

def KDn_outputfile(actual_directory,dirname,K,Dn,enerHtotalN,interf_energy,NeL2damagekn,tf,energy_evol):
    chdir(actual_directory)
    outfile=open(actual_directory+'/'+dirname+'/'+'KDn_recordfile.txt','a')
    for kinc in range(1,len(energy_evol)):
        outfile.write(str('%d\t %d\t %d\t %.8e\t %.8e\t %d\t'%(K,Dn,energy_evol[kinc][0],tf,interf_energy,NeL2damagekn)))
        outfile.write(str('%.8e\t %.8e\t' %(energy_evol[kinc][-2],energy_evol[kinc][-1])))
        outfile.write(str('%.8e\t %.8e\t %.8e\t %.8e\t' %(energy_evol[kinc][1],energy_evol[kinc][5],
                                                          energy_evol[kinc][6],energy_evol[kinc][7]))+'\n')
    outfile.close()
    
def add_urdfill_output(input_file_path):
    """
    Finds the '** OUTPUT REQUESTS' header in an Abaqus input file and
    inserts a block for the URDFIL subroutine output right after it.

    The function will insert the following lines:
    ** This block is for the URDFIL subroutine
    *ENERGY FILE, FREQUENCY=1

    Args:
        input_file_path: The full path to the Abaqus input file (.inp).
    """
    # The text block to insert, with newlines included
    urdfill_block = (
        "** This block is for the URDFIL subroutine\n"
        "*ENERGY FILE, FREQUENCY=1\n"
        "**\n"
        "*CONTACT FILE, FREQUENCY=1\n"
        "SDV\n"
        "**"
    )

    # Define the target header to search for
    target_header = "** OUTPUT REQUESTS"

    try:
        # Read all lines from the original file into a list
        with open(input_file_path, 'r') as file:
            lines = file.readlines()

        # Find the index for insertion
        insert_index = -1
        for i, line in enumerate(lines):
            # Make the search case-insensitive and ignore leading/trailing whitespace
            if line.strip().upper() == target_header:
                insert_index = i + 1  # We want to insert on the line AFTER the header
                break

        # If the header was found, insert the block and rewrite the file
        if insert_index != -1:
            lines.insert(insert_index, urdfill_block)
            
            # Overwrite the original file with the modified content
            with open(input_file_path, 'w') as file:
                file.writelines(lines)
            print('Success: URDFIL output block added')
        else:
            # If the loop finishes without finding the header
            "Warning: Header '{}' not found. File was not modified.".format(target_header)

    except FileNotFoundError:
        "Error: The file at '{}' could not be found.".format(input_file_path)
    except Exception as e:
        "An unexpected error occurred: {}".format(e)

def modify_amplitude_data(input_file_path, new_load):
    """
    Finds a specific amplitude header in an Abaqus input file and
    replaces the data on the immediately following line.

    Specifically, it searches for the line containing '*Amplitude, name=AMP-1'
    and replaces the next line with:
             0., {new_load}, 1., {new_load}

    Args:
        input_file_path: The full path to the Abaqus input file (.inp).
        new_load: The positive float value to insert into the amplitude data.
    """
    
    # The header line to search for
    target_header = "*AMPLITUDE, name=AMP-1"
    # The header string to search for, cleaned of spaces for a robust check
    target_header_clean = "*AMPLITUDE,NAME=AMP-1"
    
    # Ensure the new_load value is a positive float
    if not isinstance(new_load, (int, float)) or new_load < 0:
        print("Error: new_load must be a positive number.")
        return

    try:
        # Read all lines from the original file
        with open(input_file_path, 'r') as file:
            lines = file.readlines()

        # Find the index of the line to replace
        line_to_replace_index = -1
        for i, line in enumerate(lines):
            # Clean up the line for a robust, case-insensitive comparison
            cleaned_line = line.strip().upper().replace(" ", "")
            if cleaned_line == target_header_clean.replace(" ", ""):
                # The target is the line immediately after the header
                line_to_replace_index = i + 1
                break

        # If the header was found, replace the data line and rewrite the file
        if line_to_replace_index != -1 and line_to_replace_index < len(lines):
            # Format the new data line, including leading spaces and a newline character
            new_data_line = " 0., {:.3f}, 1., {:.3f}\n".format(new_load,new_load)
            lines[line_to_replace_index] = new_data_line

            # Overwrite the original file with the modified content
            with open(input_file_path, 'w') as file:
                file.writelines(lines)
            print("Success: Amplitude 'AMP-1' was updated.")
        else:
            # If the loop finishes without finding the header
            "Warning: Header '{}' not found. File was not modified.".format(target_header)

    except FileNotFoundError:
        "Error: The file at '{}' could not be found.".format(input_file_path)
    except Exception as e:
        "An unexpected error occurred: {}".format(e)

def modify_amplitude_dataV2(input_file_path, new_load, new_load_p1):
    """
    Finds the '*Amplitude, name=AMP-1' header in an Abaqus input file
    and replaces the data on the immediately following line. (Python 2.7 compatible)
    
    The search logic has been enhanced to find the header even if it's prefixed by
    other text or keywords on the same line.

    Args:
        input_file_path: The full path to the Abaqus input file (.inp).
        new_load: The float value for the amplitude at time 0.
        new_load_p1: The float value for the amplitude at time 1.
    """
    # The header string to search for, cleaned of spaces for a robust check
    target_header_clean = "*AMPLITUDE,NAME=AMP-1"

    # Input validation remains
    if not all(isinstance(v, (int, float)) for v in [new_load, new_load_p1]):
        print("Error: 'new_load' and 'new_load_p1' must be numeric values.")
        return

    try:
        with open(input_file_path, 'r') as file:
            lines = file.readlines()

        line_to_replace_index = -1
        # --- CORRECTION APPLIED HERE ---
        for i, line in enumerate(lines):
            # Clean up the line: uppercase, remove all spaces, remove leading/trailing whitespace
            cleaned_line_for_search = line.strip().upper().replace(" ", "")

            # Check if the target string is contained anywhere in the cleaned line
            if target_header_clean in cleaned_line_for_search:
                # The target is the line *after* the header
                line_to_replace_index = i + 1
                break
        # --------------------------------

        if line_to_replace_index != -1 and line_to_replace_index < len(lines):
            # Format the new data line using .format()
            new_data_line = " 0., {0}, 1., {1}\n".format(new_load, new_load_p1)
            
            lines[line_to_replace_index] = new_data_line

            with open(input_file_path, 'w') as file:
                file.writelines(lines)
            print("Success: Amplitude 'AMP-1' in '{0}' was updated.".format(input_file_path))
        else:
            print("Warning: Header '{0}' not found. File was not modified.".format(target_header_clean))

    except IOError:
        print("Error: The file at '{0}' could not be found.".format(input_file_path))
    except Exception as e:
        print("An unexpected error occurred: {0}".format(e))

# restart input file function definition:
def create_crack_restart_inp(file_name, n, displacement_mag):
    """
    Generates an Abaqus restart .inp file for crack propagation step 'n'.
    
    Args:
        file_name (str): Name of the new .inp file (without extension).
        n (int): The current loop iterator (step number). Starts at 2.
        displacement_mag (float): The target displacement value for this step 
        (e.g., 0.0395).
    """
    
    # Calculate the step to read from
    prev_step = n - 1
    
    # Create unique names for this iteration to avoid conflicts
    amp_name = "AMP-BIS-{}".format(n)
    step_name = "Step-{}".format(n)
    
    full_filename = file_name + '.inp'
    
    print ("Writing restart file: {}".format(full_filename))
    
    with open(full_filename, 'w') as f:
        # --- HEADER ---
        f.write("*HEADING\n")
        f.write("Crack Propagation Restart - Step {}, applied load={}\n".format(n, displacement_mag))
        
        # --- RESTART READ ---
        # The 'END STEP' is crucial here to reset KINC to 1
        f.write("*RESTART, READ, STEP={}, END STEP\n".format(prev_step))
        
        # --- AMPLITUDE DEFINITION ---
        # We define a NEW amplitude for this specific step.
        # Format: Time, Value, Time, Value
        # This creates a constant hold at 'displacement_mag' across the step (0.0 to 1.0)
        f.write("**\n")
        f.write("** New Amplitude for Step {}\n".format(n))
        f.write("*AMPLITUDE, NAME={}\n".format(amp_name))
        f.write("0.0, {}, 1.0, {}\n".format(displacement_mag, displacement_mag))
        
        # --- STEP DEFINITION ---
        f.write("**\n")
        f.write("** STEP: {}\n".format(step_name))
        f.write("**\n")
        f.write("*Step, name={}, nlgeom=NO, inc=25, unsymm=NO\n".format(step_name))
        f.write("*Static, direct\n")
        f.write("0.04, 1., \n")
        
        # --- BOUNDARY CONDITIONS ---
        f.write("**\n")
        f.write("** BOUNDARY CONDITIONS\n")
        f.write("**\n")
        # Apply Positive Displacement to TOPN
        f.write("** Name: Disp-BC-4 Type: Displacement/Rotation\n")
        f.write("*Boundary, amplitude={}\n".format(amp_name))
        f.write("TOPN, 2, 2, 1.\n")
        
        # Apply Negative Displacement to BOTN
        f.write("** Name: Disp-BC-5 Type: Displacement/Rotation\n")
        f.write("*Boundary, amplitude={}\n".format(amp_name))
        f.write("BOTN, 2, 2, -1.\n")
        
        # --- OUTPUT REQUESTS ---
        f.write("**\n")
        f.write("** OUTPUT REQUESTS\n")
        f.write("** This block is for the URDFIL subroutine\n")
        f.write("*ENERGY FILE, FREQUENCY=1\n")
        
        # Important: Write restart for the NEXT loop iteration
        f.write("*Restart, write, frequency=1\n")
        
        f.write("**\n")
        f.write("** FIELD OUTPUT: F-Output-1\n")
        f.write("**\n")
        f.write("*Output, field\n")
        f.write("*Node Output, nset=INITIALFRONT\n")
        f.write("COORD, \n")
        
        f.write("**\n")
        f.write("** FIELD OUTPUT: F-Output-2\n")
        f.write("**\n")
        f.write("*Contact Output, nset=\"TOP ARM-1\"\n")
        f.write("CDISP, CSTRESS, SDV\n")
        
        f.write("**\n")
        f.write("** FIELD OUTPUT: F-Output-3\n")
        f.write("**\n")
        f.write("*Element Output, elset=\"TOP ARM-1\", directions=NO\n")
        f.write("SDV, \n")
        
        f.write("**\n")
        f.write("** FIELD OUTPUT: F-Output-4\n")
        f.write("**\n")
        f.write("*Node Output, nset=\"TOP ARM-1\"\n")
        f.write("U, \n")
        
        f.write("**\n")
        f.write("** FIELD OUTPUT: F-Output-5\n")
        f.write("**\n")
        f.write("*Node Output, nset=TOPN\n")
        f.write("CF, TF, UT\n")
        
        f.write("**\n")
        f.write("** HISTORY OUTPUT: H-Output-1\n")
        f.write("**\n")
        f.write("*Output, history\n")
        f.write("*Energy Output\n")
        f.write("ALLIE, ALLSE, ALLWK\n")
        
        f.write("**\n")
        f.write("** HISTORY OUTPUT: H-Output-2\n")
        f.write("**\n")
        f.write("*Node Output, nset=TOPN\n")
        f.write("CF2, TF2, U2\n")
        
        f.write("*End Step\n")
    return full_filename           
#%% datos de entrada
call('cls',shell=True)
actual_directory = os.getcwd()
#cambiar a actual_directory la siguiente linea si quiero correrlo en el actual
working_directory=actual_directory
start = time.time()
list_delete = ('*.sta', '*.dat', '*.sim', '*.prt', '*.msg', '*.com', '*.jnl',\
               '*.mtx','*.pes','*.par','*.pmg','*.odb','*.ipm','*.pyc')
borrar_archivos(list_delete)

# force (1) or displacement (0) control:
control=0
# initial load increment:
incr=1
factorcarga=1.0

# initial load:
if control==1:
    cargainicial=float(1)
else:
    cargainicial=1.0E-3
# initial step:
k=1
# Bisection algorithm stop tolerance
bisectol=10.0E-3
# Crack advance bisection tolerance
bisectol2=float(2)*bisectol
# Initial bisection tolerance
eps0=1
# Bisection iteration counter
bisect_itercnt=1
# Bisection limit of iterations
numiter=50
# Bisection scheme: full (1) or crack onset only (0)
bis_option=0
#
Nmin = round(-mth.log10(bisectol)/mth.log10(2))

nameinp='DCBuinter_L237m0_Uctrl'
UINTER_lst=['UINTERLEBIM3D_kincxit.for','UINTERLEBIMAMA3D.for']
# nameUMAT=UINTER_lst[0]
nameSUBR=UINTER_lst[0]
amplname='AMP-1'
godb=1  #si quiero guardar todos los odb pongo 1, si no 0
#definicion de mi uinter para cambiar
cadenasUinter=['*Surface Interaction, name=IntProp-1, user, depvar=22, properties=9\n','\n']
kn = 300; kt = kn*1

cadenasUinterN=[kn,kt,kt,0.09375,0.0,0.0,8.]

# 11 starting damage configurations:
cadenasUinterN1=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(1)+", "+str(control)+'\n']))]
cadenasUinterN2=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(2)+", "+str(control)+'\n']))]
cadenasUinterN3=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(3)+", "+str(control)+'\n']))]
cadenasUinterN4=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(4)+", "+str(control)+'\n']))]
cadenasUinterN5=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(5)+", "+str(control)+'\n']))]
cadenasUinterN6=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(6)+", "+str(control)+'\n']))]
cadenasUinterN7=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(7)+", "+str(control)+'\n']))]
cadenasUinterN8=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(8)+", "+str(control)+'\n']))]
cadenasUinterN9=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(9)+", "+str(control)+'\n']))]
cadenasUinterN10=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(10)+", "+str(control)+'\n']))]
cadenasUinterN11=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(11)+", "+str(control)+'\n']))]

Uinterprops=[cadenasUinterN1,cadenasUinterN11]    
cadenaSDVcontact=['CDISP, CSTRESS, SDV\n']
# cadenanombreSDV='ASSEMBLY_VIGA_INF_SUPER_SUP/ASSEMBLY_VIGA_SUP_SUPER_INF'
cadenanombreSDV='ASSEMBLY_TOP_SURF/ASSEMBLY_BOTTOM_SURF'
Dirname='AdaptiveF_ODBs-'+nameinp+'_knn'+str(cadenasUinterN[0])+'_mu'+str(cadenasUinterN[-1])+'Ctrl'+str(control)+'-Bisopt'+str(bis_option)

##para cluster
#cpus=10
#memoria=30000
#unidades=MEGA_BYTES
#para pc
cpus=1
memoria=90
#unidades=PERCENTAGE

#%% corremos para sacar una carga cerca de la carga que rompe el primer PI por CT
chdir(working_directory)
if (path.exists(Dirname) == True):
    rmtree(Dirname)
if godb==1:
    mkdir(Dirname)
sim_file=Dirname+'/stepdata.txt'

nameodbcargainicial=nameinp+'-carga-inicial'
shutil.copy(nameinp+'.inp',nameodbcargainicial+'.inp')
modify_amplitude_data(nameodbcargainicial+'.inp', cargainicial)
#Buscamos uinter='*Surface Interaction'
cadenasCohesivo=lineasdefinteraccion(nameodbcargainicial+'.inp', '*Surface Interaction',3)  #buscamos la cadena a cambiar
#reemplazamos la primera linea en los dos archivos de inicio
replacement(nameodbcargainicial+'.inp', cadenasCohesivo[0], cadenasUinter[0])
#reemplazamos la segunda linea en los dos archivos de inicio
replacement(nameodbcargainicial+'.inp', cadenasCohesivo[1], cadenasUinter[1])
#reemplazamos la tercera linea en los dos archivos de inicio. Esta es distinta
replacement(nameodbcargainicial+'.inp', cadenasCohesivo[2], cadenasUinterN1[0][:])
#reemplazamos la linea que pide las salida de las superficies de contacto    
cadenaFOcontact=lineasdefinteraccion(nameodbcargainicial+'.inp', 'CDISP, CSTRESS,',1)  #buscamos la cadena a cambiar
#cadenaSDVcontact=['CDISP, CSTRESS, SDV\n']
replacement(nameodbcargainicial+'.inp', cadenaFOcontact[0], cadenaSDVcontact[0])
replacement(nameodbcargainicial+'.inp', cadenaFOcontact[0], cadenaSDVcontact[0])
signalFile = 'exit_signal.txt'
# --- Ensure the signal file does not exist before starting ---
# This is to avoid false positives from previous runs.
try:
    os.remove(signalFile)
except OSError:
    pass

myJob = mdb.JobFromInputFile(name=nameodbcargainicial, 
        inputFileName=nameodbcargainicial, type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None,
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, memory=memoria, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, userSubroutine=nameSUBR, 
        scratch=actual_directory, parallelizationMethodExplicit=DOMAIN, numDomains=cpus, 
        activateLoadBalancing=False, multiprocessingMode=DEFAULT, numCpus=cpus)

# myJob.submit(consistencyChecking=OFF)
# myJob.waitForCompletion()

# --- 1. SETUP ---
jobName = myJob.name
# The command to run Abaqus from the command line
# 'interactive' tells Abaqus to run the job and then exit.
# Define the absolute path to the executable
abaqus_executable = r'C:\SIMULIA\Abaqus\6.14-1\code\bin\abq6141.exe'

# Populate the template with all necessary variables
# Build a list of command arguments
command = [
    abaqus_executable,
    'job=' + jobName,
    'user=' + nameSUBR,
    'interactive'
]
# command = 'abaqus job={} user={} interactive'.format(jobName,nameSUBR)
try:
    os.remove(signalFile)
    # Also remove the Abaqus lock file if it exists
    if os.path.exists('{}.lck'.format(jobName)):
        os.remove('{}.lck'.format(jobName))
except OSError:
    pass # Ignore errors if files don't exist
    
# --- 2. LAUNCH ABAQUS (NON-BLOCKING) ---
# subprocess.Popen launches the command in a new process
# and the script immediately continues to the next line.
process = subprocess.Popen(command, cwd=actual_directory)

# --- 3. MONITOR THE JOB (POLLING LOOP) ---
print("--- Entering monitoring loop... ---")
# The 'process.poll()' method checks if the subprocess has terminated.
# It returns 'None' if it's still running.
while process.poll() is None:
    "Waiting for '{}' to finish... (checking again in N seconds)".format(jobName)
    # 'time.sleep()' makes the script pause, preventing it from
    # using 100% CPU while it waits.
    time.sleep(10)
    
"\n--- Abaqus job has finished ---"
"--- The Abaqus process returned exit code: {} ---".format(process.returncode)

# --- Check for the signal file after the job finishes ---
if os.path.exists(signalFile):
    with open(signalFile, 'r') as f:
        content = f.read().strip()
    print('Abaqus job terminated early by UINTER')
    print('Signal received: ',content)
    # Clean up the signal file
    os.remove(signalFile)
    
###obtener la nueva carga del odb
print(nameodbcargainicial, cargainicial, cadenanombreSDV)
newload = PMTESCabaqusV0.PMTESCcriTen_carga(nameodbcargainicial, cargainicial, cadenanombreSDV)
# 
newload=(newload*factorcarga)
print('Applied load: ',float(newload))
###borrar archivos
if godb==1:
    list_odbs=([nameodbcargainicial+'*.odb',nameodbcargainicial+'*.inp'])
    mover_archivos(list_odbs,Dirname)
#borramos el resto de archivos de ese paso                                          
list_delete = ([nameodbcargainicial+'.*'])
#borrar_archivos(list_delete)
###preparamos nuestro inp k=1
nameinpK=nameinp+str(1)
shutil.copy(nameinp+'.inp',nameinpK+'.inp')

###importo inp y impongo nueva carga
# print(newload)
mdb.ModelFromInputFile(inputFileName=nameinpK+'.inp', name=nameinpK)
mdb.models[nameinpK].amplitudes[amplname].setValues(timeSpan=STEP, smooth=SOLVER_DEFAULT, data=((0.0, newload), (1.0, newload)))
mdb.Job(name=nameinpK, model=nameinpK, description='', 
            type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, 
            memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, memory=memoria,
            explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
            modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine=nameSUBR, 
            scratch=actual_directory, resultsFormat=ODB, multiprocessingMode=DEFAULT, numDomains=cpus, numCpus=cpus, 
            numGPUs=0)

####escribo el inp de k=1 con la nueva carga
mdb.jobs[nameinpK].writeInput(consistencyChecking=OFF)
# Add .fil output requests for each input file:
add_urdfill_output(nameinpK+'.inp')
# Directly overwrite the *AMP-1 values in the .inp file:
modify_amplitude_data(nameinpK+'.inp', newload)

##

###inicializamos variables:
NeL2damageKm1=[]
datos_salida=[]
stepdata=[]
# sys.exit()
#%% LOAD BISECTION ALGORITHM STARTS HERE

while (k<=numiter):
    # nameinpK=nameinp+str(1)
    # shutil.copy(nameinp+'.inp',nameinpK+'.inp')
    #%% 
    #%
    
    subroutine_file = nameSUBR  # Update this path
    input_file = nameinpK+'.inp'  # Update this path
    common_file = 'my_common.inc'
    
    Joblist=[]
    Pmtesc_postp=[]
    for Dn in range(1,len(Uinterprops)+1):
        #creamos los inp nuevos
        nameinpDNk=nameinpK+'_N'+str(Dn)+'_K'+str(k)
        # nameinpN2k=nameinpK+'_N2_K'+str(k)
        copy(nameinpK+'.inp',nameinpDNk+'.inp')  
        # copy(nameinpK+'.inp',nameinpN2k+'.inp')  
        print(Uinterprops[Dn-1][:][0])
        #abrimos el inp y buscamos la frase con las pripiedades de la UMAT o Uinter
        #Buscamos uinter='*Surface Interaction'
        cadenasCohesivo=lineasdefinteraccion(nameinpK+'.inp', '*Surface Interaction',3) #buscamos la cadena a cambiar
        #reemplazamos la primera linea en los dos archivos de inicio
        replacement(nameinpDNk+'.inp', cadenasCohesivo[0], cadenasUinter[0])
        # replacement(nameinpN2k+'.inp', cadenasCohesivo[0], cadenasUinter[0])
        #reemplazamos la segunda linea en los dos archivos de inicio
        replacement(nameinpDNk+'.inp', cadenasCohesivo[1], cadenasUinter[1])
        # replacement(nameinpN2k+'.inp', cadenasCohesivo[1], cadenasUinter[1])
        #reemplazamos la tercera linea en los dos archivos de inicio. Esta es distinta
        replacement(nameinpDNk+'.inp', cadenasCohesivo[2], Uinterprops[Dn-1][:][0])
        # replacement(nameinpN2k+'.inp', cadenasCohesivo[2], cadenasUinterN2[0])
        #reemplazamos la linea que pide las salida de las superficies de contacto    
        cadenaFOcontact=lineasdefinteraccion(nameinpDNk+'.inp', 'CDISP, CSTRESS',1)  #buscamos la cadena a cambiar
        cadenaSDVcontact=['CDISP, CSTRESS, SDV\n']
        replacement(nameinpDNk+'.inp', cadenaFOcontact[0], cadenaSDVcontact[0])
        # replacement(nameinpN2k+'.inp', cadenaFOcontact[0], cadenaSDVcontact[0])
        
    
    Joblist=[]
    processlist=[]    
    # ***ABAQUS JOB PARALLELIZATION APPROACH:***
    for job in range(1,len(Uinterprops)+1):
        
        nameinpDNk=nameinpK+'_N'+str(job)+'_K'+str(k)
        # Create unique work directory for each job
        job_name = '{}'.format(nameinpDNk)
        # work_dir = os.path.join(base_dir, job_name)
        work_dir = os.path.join(actual_directory, os.path.basename(job_name + '_scratch'))
        
        if os.path.exists(work_dir):
            shutil.rmtree(work_dir)
        if not os.path.exists(work_dir):
            os.makedirs(work_dir)
        # Add .fil output requests for each input file:
        add_urdfill_output(nameinpDNk+'.inp')
        # Directly overwrite the *AMP-1 values in the .inp file:
        modify_amplitude_data(nameinpDNk+'.inp', newload)
        
        # Copy subroutine to job directory (CRITICAL STEP)
        job_sub_path = os.path.join(work_dir, os.path.basename(subroutine_file))
        job_commfile_path = os.path.join(work_dir, os.path.basename(common_file))
        shutil.copy(subroutine_file, job_sub_path)
        shutil.copy(common_file, job_commfile_path)
        # Copy .inp file into the job work directory
        shutil.copy(nameinpDNk+'.inp', work_dir)
        #
        
        myJobN = mdb.JobFromInputFile(name=nameinpDNk, 
            inputFileName=nameinpDNk, type=ANALYSIS, 
            atTime=None, waitMinutes=0, waitHours=0, queue=None,
            memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, memory=memoria, 
            explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, userSubroutine=nameSUBR, 
            scratch='', parallelizationMethodExplicit=DOMAIN, numDomains=cpus, 
            activateLoadBalancing=False, multiprocessingMode=DEFAULT, numCpus=cpus)
        
        
        myJobN.setValues(
        scratch=work_dir,
        userSubroutine=job_sub_path,  # Point to local subroutine,
        numCpus=1,          # MUST BE 1 FOR PARALLEL EXECUTION
        numDomains=1,        # Disable domain decomposition
        multiprocessingMode=DEFAULT)
        
        Joblist.append(myJobN)
        chdir(work_dir)
        # Create a copy of the current environment
        # job_env = os.environ.copy()
        # Set the ABAQUS_SCRATCH variable for this specific job
        # job_env["ABAQUS_SCRATCH"] = work_dir
        
        abaqus_executable = r'C:\SIMULIA\Abaqus\6.14-1\code\bin\abq6141.exe'

        # Populate the template with all necessary variables
        # Build a list of command arguments
        commandN = [
            abaqus_executable,
            'job=' + myJobN.name,
            'user=' + nameSUBR,
            'interactive']
            
        # commandN = 'abaqus job={} user={} interactive'.format(myJobN.name,nameSUBR)
        # myJobN.submit() #consistencyChecking=OFF
        # Add a print statement for verification during the run
        print(job_name,work_dir)
        # try:
            # os.remove(signalFile)
            # # Also remove the Abaqus lock file if it exists
            # if os.path.exists('{}.lck'.format(myJobN.name)):
                # os.remove('{}.lck'.format(myJobN.name))
        # except OSError:
            # pass # Ignore errors if files don't exist
            
        # --- 2. LAUNCH ABAQUS JOB(NON-BLOCKING) ---
        # subprocess.Popen launches the command in a new process
        # and the script immediately continues to the next line.
        process = subprocess.Popen(commandN, cwd=work_dir)
        processlist.append(process)
                            
        # return to the parent, main script directory:
        chdir(actual_directory)
    
    # =================================================================
    # WAITING LOOP: Wait for all jobs to complete
    # =================================================================
    print('\n--- All jobs submitted. Waiting for completion... ---')
    # --- LOOP 2: WAIT FOR JOBS TO FINISH ---
    # This loop will pause the script until every job is done.
    for i, proc in enumerate(processlist):
        # The .wait() method blocks execution until the process finishes.
        return_code = proc.wait()
        
        # Optional: Check if the job completed successfully
        if return_code == 0:
            "Job {0} completed successfully.".format(i + 1)
        else:
            "!!! ERROR: Job {0} failed with return code {1}.".format(i + 1, return_code)
    
    for job_obj in Joblist:
        job_name = '{}'.format(job_obj.name)
        # work_dir = os.path.join(base_dir, job_name)
        work_dir = os.path.join(actual_directory, os.path.basename(job_name + '_scratch'))
        chdir(work_dir)
        if k > 1:
            steps_data=np.load(actual_directory+'\sim_file.npy')
            # print(steps_data)
            kstepindx=-1
            km1stepindx=np.where(steps_data[:,0] == k-1)[0][-1]
            Fk = stepdata[kstepindx][1]
            Fkm1 = stepdata[km1stepindx][1]
            incre = Fk-Fkm1
            eps0 = abs(incre)/abs(Fk)
            NeL2damagek=stepdata[kstepindx][5]
            # NeL2damagekm1=stepdata[km1stepindx][5]
            if eps0<bisectol and NeL2damagek>0 and stepdata[-1][-1]<0:
                eps0file = open(work_dir+'\eps0_value.txt', 'w')
                eps0file.write(str(float(1)))
                eps0file.close()
                
                # List of file extensions to delete
                extensions_to_delete = ['*.sim', '*.com', '*.jnl',\
                                '*.mtx','*.pes','*.par','*.pmg','*.ipm','*.pyc']
                for ext in extensions_to_delete:
                    # glob.glob finds all files matching the pattern
                    files_to_remove = glob.glob(ext)
                    for f in files_to_remove:
                        try:
                            os.remove(f)
                            "  - Deleted: {}".format(f)
                        except OSError as e:
                            "  - Error deleting file {}: {}".format(f, e)
            else:
                eps0file = open(work_dir+'\eps0_value.txt', 'w')
                eps0file.write(str(float(0)))
                eps0file.close()
                # List of file extensions to delete
                extensions_to_delete = ['*.inp','*.dat', '*.sim', '*.com', '*.jnl',\
                                '*.mtx','*.pes','*.par','*.pmg','*.ipm','*.pyc','*.for','*.f','*.fil']
                
                for ext in extensions_to_delete:
                    # glob.glob finds all files matching the pattern
                    files_to_remove = glob.glob(ext)
                    for f in files_to_remove:
                        try:
                            os.remove(f)
                            "  - Deleted: {}".format(f)
                        except OSError as e:
                            "  - Error deleting file {}: {}".format(f, e)
        elif k==1:
            # List of file extensions to delete
            extensions_to_delete = ['*.inp','*.dat', '*.sim', '*.prt', '*.com', '*.jnl',\
                            '*.mtx','*.pes','*.par','*.pmg','*.ipm','*.pyc','*.res','*.mdl','*.stt','*.for','*.f','*.fil']
            
            for ext in extensions_to_delete:
                # glob.glob finds all files matching the pattern
                files_to_remove = glob.glob(ext)
                for f in files_to_remove:
                    try:
                        os.remove(f)
                        "  - Deleted: {}".format(f)
                    except OSError as e:
                        "  - Error deleting file {}: {}".format(f, e)
        # return to the parent, main script directory:
        chdir(actual_directory)
    
    # for Dn in range(1,len(Uinterprops)-4):
        # nameinpDNk=nameinpK+'_N'+str(Dn)+'_K'+str(k)
        # job_name = '{}'.format(nameinpDNk)
        # work_dir = os.path.join(actual_directory, os.path.basename(job_name + '_scratch'))
        # work_dirodb=actual_directory+'//'+job_name + '_scratch'
        # odb_file=work_dirodb+'//'+nameinpDNk
        # # ODB Post-processing; simulation data output:
        # enerHtotalN,sumaenerinterN,NeL2damageKN,NeL2damagetotalN,tf,energy_evol=PMTESCabaqus.PMTESCcritEneSubr(odb_file,control,cadenanombreSDV,work_dirodb)
        # model_outputfile(actual_directory,sim_file,k,Dn,newload,tf,enerHtotalN,sumaenerinterN,len(NeL2damageKN),len(NeL2damagetotalN),energy_evol)
        # Pmtesc_postp.append([enerHtotalN,sumaenerinterN,NeL2damageKN,NeL2damagetotalN,tf,nameinpDNk])
        # datos_salida.append([k,Dn,float(newload),float(tf),float(enerHtotalN),sumaenerinterN,len(NeL2damageKN),nameinpDNk,len(NeL2damagetotalN)])
        # print([k,Dn,tf,enerHtotalN,sumaenerinterN,len(NeL2damageKN),len(NeL2damagetotalN),nameinpDNk,energy_evol[-1][-1]]) 
        # # move the job directory to the parent script directory:
        # shutil.move(work_dir, Dirname)
        
    # Joblist = []
    # processlist = []
    # for job in range(7,len(Uinterprops)+1):
        # nameinpDNk=nameinpK+'_N'+str(job)+'_K'+str(k)
        # # Create unique work directory for each job
        # job_name = '{}'.format(nameinpDNk)
        # # work_dir = os.path.join(base_dir, job_name)
        # work_dir = os.path.join(actual_directory, os.path.basename(job_name + '_scratch'))
        
        # if os.path.exists(work_dir):
            # shutil.rmtree(work_dir)
        # if not os.path.exists(work_dir):
            # os.makedirs(work_dir)
        
        # # Copy subroutine to job directory (CRITICAL STEP)
        # job_sub_path = os.path.join(work_dir, os.path.basename(subroutine_file))
        # job_commfile_path = os.path.join(work_dir, os.path.basename(common_file))
        # shutil.copy(subroutine_file, job_sub_path)
        # shutil.copy(common_file, job_commfile_path)
        # # Copy .inp file into the job work directory
        # shutil.copy(nameinpDNk+'.inp', work_dir)
        # #
        
        # myJobN = mdb.JobFromInputFile(name=nameinpDNk, 
            # inputFileName=nameinpDNk, type=ANALYSIS, 
            # atTime=None, waitMinutes=0, waitHours=0, queue=None,
            # memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, memory=memoria, 
            # explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, userSubroutine=nameSUBR, 
            # scratch='', parallelizationMethodExplicit=DOMAIN, numDomains=cpus, 
            # activateLoadBalancing=False, multiprocessingMode=DEFAULT, numCpus=cpus)
        # Joblist.append(myJobN)
        
        # myJobN.setValues(
        # scratch=work_dir,
        # userSubroutine=job_sub_path,  # Point to local subroutine,
        # numCpus=1,          # MUST BE 1 FOR PARALLEL EXECUTION
        # numDomains=1,        # Disable domain decomposition
        # multiprocessingMode=DEFAULT)
        
        # chdir(work_dir)
        # # Create a copy of the current environment
        # # job_env = os.environ.copy()
        # # Set the ABAQUS_SCRATCH variable for this specific job
        # # job_env["ABAQUS_SCRATCH"] = work_dir
        
        # abaqus_executable = '/opt/abaqus/Commands/abq2022'

        # # Populate the template with all necessary variables
        # # Build a list of command arguments
        # commandN = [
            # abaqus_executable,
            # 'job=' + myJobN.name,
            # 'user=' + nameSUBR,
            # 'interactive']
            
        # # commandN = 'abaqus job={} user={} interactive'.format(myJobN.name,nameSUBR)
        # # myJobN.submit() #consistencyChecking=OFF
        # # Add a print statement for verification during the run
        # print(job_name,work_dir)
        # # try:
            # # os.remove(signalFile)
            # # # Also remove the Abaqus lock file if it exists
            # # if os.path.exists('{}.lck'.format(myJobN.name)):
                # # os.remove('{}.lck'.format(myJobN.name))
        # # except OSError:
            # # pass # Ignore errors if files don't exist
            
        # # --- 2. LAUNCH ABAQUS JOB(NON-BLOCKING) ---
        # # subprocess.Popen launches the command in a new process
        # # and the script immediately continues to the next line.
        # process = subprocess.Popen(commandN, cwd=work_dir)
        # processlist.append(process)
                            
        # # return to the parent, main script directory:
        # chdir(actual_directory)
    
    # # --- LOOP 2: WAIT FOR JOBS TO FINISH ---
    # # This loop will pause the script until every job is done.
    # for i, proc in enumerate(processlist):
        # # The .wait() method blocks execution until the process finishes.
        # return_code = proc.wait()
        
        # # Optional: Check if the job completed successfully
        # if return_code == 0:
            # print "Job {0} completed successfully.".format(i + 1)
        # else:
            # print "!!! ERROR: Job {0} failed with return code {1}.".format(i + 1, return_code)
    
    # for job_obj in Joblist:
        # # 
        # job_name = '{}'.format(job_obj.name)
        # # work_dir = os.path.join(base_dir, job_name)
        # work_dir = os.path.join(actual_directory, os.path.basename(job_name + '_scratch'))
        # chdir(work_dir)
        # # List of file extensions to delete
        # extensions_to_delete = ['*.inp','*.dat', '*.sim', '*.prt', '*.com', '*.jnl',\
                        # '*.mtx','*.pes','*.par','*.pmg','*.ipm','*.pyc','*.res','*.mdl','*.stt','*.for','*.f','*.fil']
        
        # for ext in extensions_to_delete:
            # # glob.glob finds all files matching the pattern
            # files_to_remove = glob.glob(ext)
            # for f in files_to_remove:
                # try:
                    # os.remove(f)
                    # "  - Deleted: {}".format(f)
                # except OSError as e:
                    # "  - Error deleting file {}: {}".format(f, e)
        # # return to the parent, main script directory:
        # chdir(actual_directory)  
        
    for Dn in range(1,len(Uinterprops)+1):
        nameinpDNk=nameinpK+'_N'+str(Dn)+'_K'+str(k)
        job_name = '{}'.format(nameinpDNk)
        work_dir = os.path.join(actual_directory, os.path.basename(job_name + '_scratch'))
        work_dirodb=actual_directory+'//'+job_name + '_scratch'
        odb_file=work_dirodb+'//'+nameinpDNk
        # ODB Post-processing; simulation data output:
        enerHtotalN,sumaenerinterN,NeL2damageKN,NeL2damagetotalN,tf,energy_evol=PMTESCabaqusV0.PMTESCcritEneSubr(odb_file,control,cadenanombreSDV,work_dirodb)
        KDn_outputfile(actual_directory,Dirname,k,Dn,enerHtotalN,sumaenerinterN,len(NeL2damageKN),tf,energy_evol)
        model_outputfile(actual_directory,sim_file,k,Dn,newload,tf,enerHtotalN,sumaenerinterN,len(NeL2damageKN),len(NeL2damagetotalN),energy_evol)
        Pmtesc_postp.append([enerHtotalN,sumaenerinterN,NeL2damageKN,NeL2damagetotalN,tf,nameinpDNk])
        datos_salida.append([k,Dn,float(newload),float(tf),float(enerHtotalN),sumaenerinterN,len(NeL2damageKN),nameinpDNk,len(NeL2damagetotalN)])
        # print([k,Dn,tf,enerHtotalN,sumaenerinterN,len(NeL2damageKN),len(NeL2damagetotalN),nameinpDNk,energy_evol[-1][5]]) 
        # move the job directory to the parent script directory:
        # shutil.move(work_dir, Dirname)
        
    #%% Algorithm flow control for the ODB file with minimum total energy PI+DeltaR:
    MindelPIpDeltaR=min(row[0] for row in Pmtesc_postp); print(MindelPIpDeltaR)
    min_row = min(Pmtesc_postp, key=lambda x:x[0])
    max_row = max(Pmtesc_postp, key=lambda x:x[0])
    indx = min(enumerate(Pmtesc_postp), key=lambda x:x[1][0])[0]
    # for row in range(len(Pmtesc_postp)):
    #     if len(Pmtesc_postp[row][2])>0.0 and Pmtesc_postp[row][0]<=min_row[0]:
    odbdef=Pmtesc_postp[indx][-1]
    NeL2damage=Pmtesc_postp[indx][2]
    NeL2damagetotal=Pmtesc_postp[indx][3]
    job_name = '{}'.format(odbdef)
        
        # elif (max_row>=min_row) and len(Pmtesc_postp[row][2])==0.0:
        #     # if len(Pmtesc_postp[row][2])==0.0:
        #     # 
        #     odbdef=Pmtesc_postp[row][-1]
        #     NeL2damage=Pmtesc_postp[row][2]
        #     NeL2damagetotal=Pmtesc_postp[row][3]							 
        #     job_name = '{}'.format(odbdef)
    print(odbdef, len(NeL2damage))
    work_dirodb=Dirname+'//'+job_name + '_scratch'
    odb_file=job_name +'_scratch'+'//'+odbdef
    enerHtotalN,sumaenerinterN,NeL2damageKN,NeL2damagetotalN,tf,energy_evol=PMTESCabaqusV0.PMTESCcritEneSubr(odb_file,control,cadenanombreSDV,work_dirodb)
    model_outputfile(actual_directory,sim_file,k,Dn,newload,tf,enerHtotalN,sumaenerinterN,len(NeL2damageKN),len(NeL2damagetotalN),energy_evol)
    datos_salida.append([float(newload),float(tf),float(enerHtotalN),sumaenerinterN,len(NeL2damageKN),odbdef,len(NeL2damagetotalN)])
    
    stepdata.append([k,newload,tf,enerHtotalN,sumaenerinterN,len(NeL2damageKN),len(NeL2damagetotalN),energy_evol[-1][5]])
    sim_filenp=np.array(stepdata)
    np.save('sim_file.npy', sim_filenp)
    
    if k==1:
        steps_data=np.load('sim_file.npy')
        print(steps_data)
        eps0=1.0; 
        if stepdata[0][5]>0:
            newload = 0.5*stepdata[0][1]
        else:
            newload= 2*stepdata[0][1]
        print(k,newload)
        #creamos los inp nuevos
        nameinpK=nameinp+str(k)
        #copy(nameinpN1k+'.inp',nameinpK+'.inp')
        mdb.ModelFromInputFile(inputFileName=nameinpK+'.inp', name=nameinpK)
        mdb.models[nameinpK].amplitudes[amplname].setValues(timeSpan=STEP, smooth=SOLVER_DEFAULT, data=((0.0, newload), (1.0, newload)))
        # job_name=nameinpK

        # write the .inp with newload
        mdb.jobs[nameinpK].writeInput(consistencyChecking=OFF)
        list_odbs=(['*_K'+str(k)+'.odb','*_K'+str(k)+'.inp',
                    '*_K'+str(k)+'.sta'])
        mover_archivos(list_odbs,Dirname)
        #borramos el resto de archivos de ese paso                                          
        list_delete = ('*.sta', '*.dat', '*.sim', '*.prt', '*.msg', '*.com', '*.jnl',\
                       '*.mtx','*.pes','*.par','*.pmg','*.odb','*.ipm','*.pyc','*.res','*.mdl','*.stt')
        borrar_archivos(list_delete)
        for Dn in range(1,len(Uinterprops)+1):
            nameinpDNk=nameinpK+'_N'+str(Dn)+'_K'+str(k)
            job_name = '{}'.format(nameinpDNk)
            work_dir = os.path.join(actual_directory, os.path.basename(job_name + '_scratch'))
            shutil.rmtree(work_dir)
    else:
        steps_data=np.load('sim_file.npy')
        print(steps_data)
        kstepindx=np.where(steps_data[:,0] == k)[0][-1]
        km1stepindx=np.where(steps_data[:,0] == k-1)[0][-1]
        Fk = stepdata[kstepindx][1]
        Fkm1 = stepdata[km1stepindx][1]
        incre = Fk-Fkm1
        eps0 = abs(incre)/abs(Fk)
        NeL2damagek=stepdata[kstepindx][5]
        NeL2damagekm1=stepdata[km1stepindx][5]
        signaldamk = 0
        matches = np.where(steps_data[:,0] == k - 1)[0]
        if matches.size > 0:
            km1_stepindx = matches[-1]
            print("Last matching index:", km1_stepindx)
        else:
            print("No match found.")
        
        if k==2 and abs(incre)!=0 and stepdata[-1][-1]<0 and stepdata[-1][5]!=0:
            if stepdata[-2][5]>0 and stepdata[-2][-1]<0:
                newload = 0.5*Fk
            elif stepdata[-2][5]==0:
                newload= Fkm1 + abs(0.5*incre)
            
            #creamos los inp nuevos
            nameinpK=nameinp+str(k)
            shutil.copy(nameinp+str(k-1)+'.inp',nameinpK+'.inp')
            mdb.ModelFromInputFile(inputFileName=nameinpK+'.inp', name=nameinpK)
            mdb.models[nameinpK].amplitudes[amplname].setValues(timeSpan=STEP, smooth=SOLVER_DEFAULT, data=((0.0, newload), (1.0, newload)))
            # Create the job from the model
            mdb.Job(name=nameinpK, model=nameinpK)
            # write the .inp with newload
            mdb.jobs[nameinpK].writeInput(consistencyChecking=OFF)
            list_odbs=(['*_K'+str(k)+'.odb','*_K'+str(k)+'.inp',
                    '*_K'+str(k)+'.sta'])
            mover_archivos(list_odbs,Dirname)
            #borramos el resto de archivos de ese paso                                          
            list_delete = ('*.sta', '*.dat', '*.sim', '*.prt', '*.msg', '*.com', '*.jnl',\
                           '*.mtx','*.pes','*.par','*.pmg','*.odb','*.ipm','*.pyc','*.res','*.mdl','*.stt')
            borrar_archivos(list_delete)
            for Dn in range(1,len(Uinterprops)+1):
                nameinpDNk=nameinpK+'_N'+str(Dn)+'_K'+str(k)
                job_name = '{}'.format(nameinpDNk)
                work_dir = os.path.join(actual_directory, os.path.basename(job_name + '_scratch'))
                shutil.rmtree(work_dir)
        elif k==2 and abs(incre)!=0 and stepdata[-1][5]==0:
            if stepdata[-2][5]==0:
                newload = 2*Fk
            elif stepdata[-2][5]>0 and stepdata[-2][-1]<0:
                newload= Fkm1 - 0.5*abs(incre)
            
            #creamos los inp nuevos
            nameinpK=nameinp+str(k)
            shutil.copy(nameinp+str(k-1)+'.inp',nameinpK+'.inp')
            mdb.ModelFromInputFile(inputFileName=nameinpK+'.inp', name=nameinpK)
            mdb.models[nameinpK].amplitudes[amplname].setValues(timeSpan=STEP, smooth=SOLVER_DEFAULT, data=((0.0, newload), (1.0, newload)))
            # Create the job from the model
            mdb.Job(name=nameinpK, model=nameinpK)
            # write the .inp with newload
            mdb.jobs[nameinpK].writeInput(consistencyChecking=OFF)
            list_odbs=(['*_K'+str(k)+'.odb','*_K'+str(k)+'.inp',
                    '*_K'+str(k)+'.sta'])
            mover_archivos(list_odbs,Dirname)
            #borramos el resto de archivos de ese paso                                          
            list_delete = ('*.sta', '*.dat', '*.sim', '*.prt', '*.msg', '*.com', '*.jnl',\
                           '*.mtx','*.pes','*.par','*.pmg','*.odb','*.ipm','*.pyc','*.res','*.mdl','*.stt')
            borrar_archivos(list_delete)
            for Dn in range(1,len(Uinterprops)+1):
                nameinpDNk=nameinp+str(k-1)+'_N'+str(Dn)+'_K'+str(k)
                job_name = '{}'.format(nameinpDNk)
                work_dir = job_name + '_scratch'
                shutil.rmtree(work_dir)
        elif k==2 and abs(incre)!=0 and stepdata[-1][-1]>0 and stepdata[-1][5]!=0:
            if stepdata[-2][5]==0:
                newload = 2*Fk
            elif stepdata[-2][5]!=0 and stepdata[-2][-1]<0:
                newload= Fkm1 - 0.5*abs(incre)
            elif stepdata[-2][5]!=0 and stepdata[-2][-1]>0:
                newload = 2*Fk
            
            #creamos los inp nuevos
            nameinpK=nameinp+str(k)
            shutil.copy(nameinp+str(k-1)+'.inp',nameinpK+'.inp')
            mdb.ModelFromInputFile(inputFileName=nameinpK+'.inp', name=nameinpK)
            mdb.models[nameinpK].amplitudes[amplname].setValues(timeSpan=STEP, smooth=SOLVER_DEFAULT, data=((0.0, newload), (1.0, newload)))
            # Create the job from the model
            mdb.Job(name=nameinpK, model=nameinpK)
            # write the .inp with newload
            mdb.jobs[nameinpK].writeInput(consistencyChecking=OFF)
            list_odbs=(['*_K'+str(k)+'.odb','*_K'+str(k)+'.inp',
                    '*_K'+str(k)+'.sta'])
            mover_archivos(list_odbs,Dirname)
            #borramos el resto de archivos de ese paso                                          
            list_delete = ('*.sta', '*.dat', '*.sim', '*.prt', '*.msg', '*.com', '*.jnl',\
                           '*.mtx','*.pes','*.par','*.pmg','*.odb','*.ipm','*.pyc','*.res','*.mdl','*.stt')
            borrar_archivos(list_delete)
            for Dn in range(1,len(Uinterprops)+1):
                nameinpDNk=nameinp+str(k-1)+'_N'+str(Dn)+'_K'+str(k)
                job_name = '{}'.format(nameinpDNk)
                work_dir = job_name + '_scratch'
                shutil.rmtree(work_dir)
        elif k>2 and NeL2damagek!=0 and stepdata[-1][-1]<0 and eps0<bisectol and k<=numiter:
            print(k, Fk, tf, eps0, NeL2damagek, job_name +'_scratch'+'//'+odbdef)
            work_dirodb=Dirname+'//'+job_name + '_scratch'
            odb_file=job_name +'_scratch'+'//'+odbdef
            work_dir = os.path.join(actual_directory, os.path.basename(job_name + '_scratch'))
            shutil.copy(nameinp+'.inp', work_dir)
            shutil.copy(nameSUBR, work_dir)
            # sys.exit()
            if control==0:
                # CRACK PROPAGATION AND ARREST STAGE: displacement control
                # -------------------------------------------------------------
                # -------------------------------
                # CRACK ADVANCES WITH BISECTION
                # -------------------------------
                if bis_option==1:
                # commands that copy the winner input file and create/launch restart
                # analysis files
                    
                    # name of the Abaqus worker script
                    worker_script = 'runrestart_job.bat' 
                    # define the new .inp and directory name
                    nameinp_k_mas_1=odbdef+'_crackprop'+'_m'+str(1)
                    work_dirkm = os.path.join(actual_directory, 
                                 os.path.basename(odbdef+'_crackprop'+'_m'+str(0) + '_scratch'))
                    # create a new crack propagation working directory
                    os.mkdir(work_dirkm)
                    # copy the worker script into the new work directory
                    shutil.copy(worker_script, work_dirkm)
                    
                    # files we need to export
                    inp_file = os.path.join(odbdef + '.inp')
                    mdl_file = os.path.join(work_dir+'\\'+odbdef + '.mdl')
                    prt_file = os.path.join(work_dir+'\\'+odbdef + '.prt')
                    res_file = os.path.join(work_dir+'\\'+odbdef + '.res')
                    odb_file = os.path.join(work_dir+'\\'+odbdef + '.odb')
                    stt_file = os.path.join(work_dir+'\\'+odbdef + '.stt')
                    export_fil = [mdl_file,prt_file,odb_file,nameSUBR,nameinp+'.inp',
                                  res_file,stt_file,common_file]
                    # export the files into the new directory
                    for f in export_fil:
                        shutil.copy(f, work_dirkm)
                    for Dn in range(1,len(Uinterprops)+1):
                        nameinpDNk=nameinp+str(k-1)+'_N'+str(Dn)+'_K'+str(k)
                        job_name = '{}'.format(nameinpDNk)
                        work_dir = job_name + '_scratch'
                        shutil.rmtree(work_dir)
                    """
                    # *************************************************************************
                    # NOTE: for any crack propagation step inside the loop, the load amplitude
                    # values must be the same across all increments of the step, i.e.,
                    # ((0,value),(1.0,value)), consistent with crack initiation steps,
                    # which will enable a most precise crack advance load prediction;
                    # should be invariant to the bisection algorithm scheme.
                    # *************************************************************************
                    """
                    # load vector preparation:
                    start = newload + float(0.05*newload)
                    end = 6.5*newload
                    num = 27
                    power = float(2.25) # no bias towards the start
                    # Create a normalized linear space (0 to 1)
                    if num == 1:
                        linear_norm = [start]
                    else:
                        linear_norm = [float(i) / (num - 1) for i in range(num)]
                    # Apply the power and scale
                    load_vector = [start + (end - start) * (norm_val ** power) for norm_val in linear_norm]
                    
                    load_list = []
                    load_list.append([newload, 0])
                    numiter_max = 25
                    newld_iter = load_vector[0] #starting load value
                    for indx,value in enumerate(load_vector):
                        # if indx>=4: sys.exit() # provisional termination
                        advance_iternum = 1
                        
                        if indx==0: 
                            n = 2
                            # newld_iter = value
                            inpfil_name = odbdef+'_m'+str(indx+1)+'_bisiter_'+str(advance_iternum)
                            inp_filename = create_crack_restart_inp(inpfil_name, n, newld_iter)
                            shutil.move(inp_filename, work_dirkm)
                        else:
                            n += 1
                        """
                        
                        START OF ADVANCE BISECTION ALGORITHM
                        
                        """
                        chdir(work_dirkm)
                        
                        epsb2 = float(1)
                        
                        while advance_iternum<numiter_max:
                            if indx==0 and advance_iternum==1: 
                                old_load=newload
                                # load_list.append(newld_iter)
                                odb_file = odbdef
                            if indx==0 and advance_iternum>1:
                                old_load=load_list[-1][0]
                                # load_list.append(newld_iter)
                            if indx>0 and advance_iternum==1:
                                newld_iter = load_list[-1][0] + 2.0*incre
                            
                            if indx==0: 
                                oldjob_name = odbdef
                            else: oldjob_name = odb_file
                            
                            # Define extensions generated by Abaqus
                            # extensions_to_clean = ['.odb', '.lck', '.com', '.dat', '.msg', '.sta', '.prt', '.res', '.mdl', '.stt']
                            # print ("Cleaning up old attempt for: {}".format(inpfil_name))
                            # for ext in extensions_to_clean:
                            #     # Construct the full file path
                            #     file_to_remove = os.path.join(work_dirkm, inpfil_name + ext)
                                
                            #     # Remove it if it exists
                            #     if os.path.exists(file_to_remove):
                            #         try:
                            #             os.remove(file_to_remove)
                            #         except OSError:
                            #             print ("Warning: Could not remove locked file: {}".format(file_to_remove))
                                        
                            print("START OF ADVANCE BISECTION ALGORITHM: base job {}".format(oldjob_name))
                                # oldjob_name is the last converged job after the 
                                # nucleation bisection algorithm 
                            # Define files to monitor
                            lock_file = os.path.join(inpfil_name + '.lck')
                            # ... (your pre-check to remove old .lck files) ...
                            # if os.path.exists(lock_file):
                            #     os.remove(lock_file)
                                    
                            # --- THIS IS THE MODIFIED PART ---
                            # Call the .bat file and pass the job_name as an argument
                            command = [
                                worker_script,
                                inpfil_name, # current job name
                                oldjob_name, # restart base analysis
                                work_dirkm # work directory argument
                            ]
                            
                            print("Submitting job via Worker .bat file...")
                            # We use .call() because the .bat file is fast
                            # and will exit as soon as the job is submitted.
                            subprocess.call(command, cwd=work_dirkm)
                            # --- Rest of your script ---
                            
                            # 1. Wait for the job to start
                            # print ("Waiting for job to start (for .lck file)...")
                            # time_to_start = 0
                            # ... (rest of your .lck polling logic) ...
                            # while not os.path.exists(lock_file):
                            #     time.sleep(1)
                                # ... (etc) ...
                            
                            # 2. Wait for the job to finish
                            # print ("Job is running. Polling {}.lck...".format(inpfil_name))
                            # while os.path.exists(lock_file):
                            #     time.sleep(10)
                            print ("--- Iteration {} complete. ---".format(inpfil_name))
                            chdir(work_dirkm)
                            enerHtotalN,sumaenerinterN,NeL2damageKN,NeL2damagetotalN,tf,energy_evol=PMTESCabaqusV0.PMTESCcritEneSubr(inpfil_name,control,cadenanombreSDV,work_dirkm)
                            if indx==0 and advance_iternum==1: epsb2=1.0
                            elif indx==0 and advance_iternum!=1: epsb2=abs(newld_iter-old_load)/abs(newld_iter)
                            elif indx>0:
                                epsb2=abs(newld_iter-old_load)/abs(newld_iter)

                            print(indx,advance_iternum,epsb2,enerHtotalN,len(NeL2damageKN),len(NeL2damagetotalN),newld_iter,old_load)
                            # model_outputfile(actual_directory,sim_file,k,Dn,newload,tf,enerHtotalN,sumaenerinterN,len(NeL2damageKN),len(NeL2damagetotalN),energy_evol)
                            
                            if len(NeL2damageKN)>0 and epsb2>bisectol2:
                                print('CRACK ADVANCES, WITH epsb2 {}> bisectol2 {}'.format(epsb2,bisectol2))
                                advance_iternum += 1
                                load_list.append([newld_iter, 1])
                                last_row = next((row for row in reversed(load_list) if row[1] == 0), None)
                                old_load=last_row[0]
                                chdir(actual_directory)
                                newld_iter=newld_iter - 0.5*abs(newld_iter-old_load)
                                inpfil_name = odbdef+'_m'+str(indx+1)+'_bisiter_'+str(advance_iternum)
                                inp_filename = create_crack_restart_inp(inpfil_name, n, newld_iter)
                                # odb_file = inpfil_name
                                shutil.move(inp_filename, work_dirkm)
                                
                                chdir(work_dirkm)
                                
                                # List of file extensions to delete
                                extensions_to_delete = ['*.dat', '*.sim', '*.com', '*.jnl'\
                                                '*.mtx','*.pes','*.par','*.pmg','*.ipm','*.pyc','*.fil']
                                for ext in extensions_to_delete:
                                    # glob.glob finds all files matching the pattern
                                    files_to_remove = glob.glob(ext)
                                    for f in files_to_remove:
                                        try:
                                            os.remove(f)
                                            "  - Deleted: {}".format(f)
                                        except OSError as e:
                                            "  - Error deleting file {}: {}".format(f, e)
                                chdir(actual_directory)
                            if len(NeL2damageKN)==0 and epsb2>bisectol2:
                                print('NO CRACK ADVANCE, WITH epsb2 {}> bisectol2 {}'.format(epsb2,bisectol2))
                                advance_iternum += 1
                                load_list.append([newld_iter, 0])
                                last_row = next((row for row in reversed(load_list) if row[1] == 1), None)
                                if last_row==None:
                                    last_row=next((row for row in reversed(load_list) if row[1] == 0), None)
                                old_load=last_row[0]
                                chdir(actual_directory)
                                newld_iter=old_load + 1.5*abs(newld_iter-old_load)
                                
                                inpfil_name = odbdef+'_m'+str(indx+1)+'_bisiter_'+str(advance_iternum)
                                inp_filename = create_crack_restart_inp(inpfil_name, n, newld_iter)
                                # odb_file = inpfil_name
                                shutil.move(inp_filename, work_dirkm)
                                
                                chdir(work_dirkm)
                                # List of file extensions to delete
                                extensions_to_delete = ['*.dat', '*.sim', '*.com', '*.jnl'\
                                                '*.mtx','*.pes','*.par','*.pmg','*.ipm','*.pyc','*.fil']
                                for ext in extensions_to_delete:
                                    # glob.glob finds all files matching the pattern
                                    files_to_remove = glob.glob(ext)
                                    for f in files_to_remove:
                                        try:
                                            os.remove(f)
                                            "  - Deleted: {}".format(f)
                                        except OSError as e:
                                            "  - Error deleting file {}: {}".format(f, e)
                                chdir(actual_directory)
                            if len(NeL2damageKN)==0 and epsb2<bisectol2:
                                print('NO CRACK ADVANCE, WITH epsb2 {}< bisectol2 {}'.format(epsb2,bisectol2))
                                advance_iternum += 1
                                load_list.append([newld_iter, 0])
                                last_row = next((row for row in reversed(load_list) if row[1] == 1), None)
                                if last_row==None:
                                    last_row=next((row for row in reversed(load_list) if row[1] == 0), None)
                                old_load=last_row[0]
                                chdir(actual_directory)
                                newld_iter=old_load + float(1.5)*abs(newld_iter-old_load)
                                
                                inpfil_name = odbdef+'_m'+str(indx+1)+'_bisiter_'+str(advance_iternum)
                                inp_filename = create_crack_restart_inp(inpfil_name, n, newld_iter)
                                # odb_file = inpfil_name
                                shutil.move(inp_filename, work_dirkm)
                                
                                chdir(work_dirkm)
                                # List of file extensions to delete
                                extensions_to_delete = ['*.sim', '*.com', '*.jnl'\
                                                '*.mtx','*.pes','*.par','*.pmg','*.ipm','*.pyc','*.fil']
                                for ext in extensions_to_delete:
                                    # glob.glob finds all files matching the pattern
                                    files_to_remove = glob.glob(ext)
                                    for f in files_to_remove:
                                        try:
                                            os.remove(f)
                                            "  - Deleted: {}".format(f)
                                        except OSError as e:
                                            "  - Error deleting file {}: {}".format(f, e)
                                chdir(actual_directory)
                            if len(NeL2damageKN)>0  and epsb2<bisectol2:
                                print('CRACK ADVANCES, WITH epsb2 {}< bisectol2 {}'.format(epsb2,bisectol2))
                                chdir(actual_directory)
                                ref_load = newld_iter
                                load_list.append([newld_iter, 1])
                                incre = abs(newld_iter-old_load)
                                newld_iter = newld_iter + 2.0*incre
                                last_iter = advance_iternum
                                inpfil_name = odbdef+'_m'+str(indx+2)+'_bisiter_'+str(1)
                                inp_filename = create_crack_restart_inp(inpfil_name, n+1, newld_iter)
                                if os.path.exists(work_dirkm+'\\'+inp_filename):
                                    print("Input file {} already exists in the work dir".format(inp_filename))
                                else:
                                    shutil.move(inp_filename, work_dirkm)
                                oldjob_base = odbdef+'_m'+str(indx+1)+'_bisiter_'+str(advance_iternum)
                                odb_file = odbdef+'_m'+str(indx+1)+'_bisiter_'+str(advance_iternum)
                                # ADVANCE ITERATION POSTPROCESS
                                enerHtotalN,sumaenerinterN,NeL2damageKN,NeL2damagetotalN,tf,energy_evol=PMTESCabaqusV0.PMTESCcritEneSubr(work_dirkm+'//'+odb_file,control,cadenanombreSDV,work_dirkm)
                                # model_outputfile(actual_directory,sim_file,k,11,val,tf,enerHtotalN,sumaenerinterN,len(NeL2damageKN),len(NeL2damagetotalN),energy_evol)
                                # datos_salida.append([float(val),float(tf),float(enerHtotalN),sumaenerinterN,len(NeL2damageKN),job_name,len(NeL2damagetotalN)])
                                # stepdata.append([k,val,tf,enerHtotalN,sumaenerinterN,len(NeL2damageKN),len(NeL2damagetotalN),energy_evol[-1][5]])
                                KDn_outputfile(actual_directory,Dirname,k,1,enerHtotalN,sumaenerinterN,len(NeL2damageKN),tf,energy_evol)
                                
                                chdir(work_dirkm)
                                # List of file extensions to delete
                                extensions_to_delete = ['*.sim', '*.com', '*.jnl','*.txt','*.png','*.dat'\
                                                '*.mtx','*.pes','*.par','*.pmg','*.ipm','*.pyc','*.fil']
                                for ext in extensions_to_delete:
                                    # glob.glob finds all files matching the pattern
                                    files_to_remove = glob.glob(ext)
                                    for f in files_to_remove:
                                        try:
                                            os.remove(f)
                                            "  - Deleted: {}".format(f)
                                        except OSError as e:
                                            "  - Error deleting file {}: {}".format(f, e)
                                # Remove unnecessary files for each iteration but the last one
                                for iterindx in range(1,advance_iternum):
                                    iter_filename=odbdef+'_m'+str(indx+1)+'_bisiter_'+str(iterindx)
                                    del_list = [iter_filename+'.stt',iter_filename+'.res',iter_filename+'.mdl',
                                                iter_filename+'.prt',iter_filename+'.msg',iter_filename+'.odb',iter_filename+'.inp']
                                    for file in del_list:
                                        os.remove(file)
                                # Change to the current working directory
                                chdir(actual_directory)
                                break
                            
                            # 
                            if len(NeL2damageKN)>=0 and advance_iternum==numiter_max:
                                chdir(work_dirkm)
                                # List of file extensions to delete
                                extensions_to_delete = ['*.dat', '*.sim', '*.com', '*.jnl'\
                                                '*.mtx','*.pes','*.par','*.pmg','*.ipm','*.pyc','*.fil']
                                for ext in extensions_to_delete:
                                    # glob.glob finds all files matching the pattern
                                    files_to_remove = glob.glob(ext)
                                    for f in files_to_remove:
                                        try:
                                            os.remove(f)
                                            "  - Deleted: {}".format(f)
                                        except OSError as e:
                                            "  - Error deleting file {}: {}".format(f, e)
                                chdir(actual_directory)
                                print("Simulation terminated after {} advance bisection iterations, EXIT".format(advance_iternum))
                                sys.exit()
                            
                            
                        chdir(actual_directory)
                        """
                        END OF ADVANCE BISECTION ALGORITHM
                        """
                    os.chdir(work_dirkm)
                    # List of file extensions to delete
                    extensions_to_delete = ['*.sim', '*.com', '*.jnl','*.stt','*.prt','*.mdl','*.res','*.msg'\
                                    '*.mtx','*.pes','*.par','*.pmg','*.ipm','*.pyc','*.for','*.sta']
                    for ext in extensions_to_delete:
                        # glob.glob finds all files matching the pattern
                        files_to_remove = glob.glob(ext)
                        for f in files_to_remove:
                            try:
                                os.remove(f)
                                "  - Deleted: {}".format(f)
                            except OSError as e:
                                "  - Error deleting file {}: {}".format(f, e)
                    sys.exit()
                else:
                    # CRACK ADVANCES WITH PRESCRIBED LOAD INCREMENTS
                    # This is the name of your NEW Abaqus worker script
                    worker_script = 'run_job.bat' 
                    # define the new .inp and directory name
                    nameinp_k_mas_1=odbdef+'_crackprop'+'_m'+str(1)
                    work_dirkm = os.path.join(actual_directory, 
                                 os.path.basename(odbdef+'_crackprop'+'_m'+str(0) + '_scratch'))
                    # create a new crack propagation working directory
                    os.mkdir(work_dirkm)
                    # copy the worker script into the new work directory
                    shutil.copy(worker_script, work_dirkm)
                    
                    # files we need to export
                    inp_file = os.path.join(odbdef + '.inp')
                    mdl_file = os.path.join(work_dir+'\\'+odbdef + '.mdl')
                    prt_file = os.path.join(work_dir+'\\'+odbdef + '.prt')
                    res_file = os.path.join(work_dir+'\\'+odbdef + '.res')
                    odb_file = os.path.join(work_dir+'\\'+odbdef + '.odb')
                    stt_file = os.path.join(work_dir+'\\'+odbdef + '.stt')
                    export_fil = [mdl_file,prt_file,odb_file,nameSUBR,nameinp+'.inp',
                                  res_file,stt_file,common_file]
                    
                    # export the files into the new directory
                    for f in export_fil:
                        shutil.copy(f, work_dirkm)
                    """
                    # HERE AFTERWARDS, notice the use of the mdb.models and the
                    # mdb.models[].InitialState procedures to create the model 
                    # and load predefined fields for the actual step increment.
                    """
                    # copy the original .inp file
                    copy(nameinp+'.inp', nameinp_k_mas_1+'.inp')
                    #cargamos el input original de ABAQUS
                    mdb.ModelFromInputFile(inputFileName=nameinp_k_mas_1+'.inp', name=nameinp_k_mas_1)
                    #imponemos el estado inicial del job ganador entre las n para que sea el inicio del paso anterior
                    instancesname=mdb.models[nameinpK].rootAssembly.instances.keys()
                    instanceslista=[]
                    
                    for inst in range(len(instancesname)):
                        instanceslista.append(mdb.models[nameinpK].rootAssembly.instances[instancesname[inst]])
                        
                    mdb.models[nameinpK].InitialState(updateReferenceConfiguration=OFF, fileName=odbdef, 
                              endStep=LAST_STEP, endIncrement=STEP_END, name='Predefined Field-1', 
                              createStepName='Initial', instances=tuple(instanceslista))
                    # cambiamos el nombre del step  
                    mdb.models[nameinpK].steps.changeKey(fromName='Step-1', toName='Step-'+str(1)+'_m'+str(1))
                    
                    cadenasUinter=['*Surface Interaction, name=IntProp-1, user, depvar=22, properties=9\n','\n']
                    kn = 300; kt = kn*1
                    cadenasUinterN=[kn,kt,kt,0.09375,0.0,0.0,8.]
                    # 11 starting damage configurations:
                    cadenasUinterN1=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(1)+", "+str(control)+'\n']))]
                    cadenasUinterN2=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(2)+", "+str(control)+'\n']))]
                    cadenasUinterN3=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(3)+", "+str(control)+'\n']))]
                    cadenasUinterN4=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(4)+", "+str(control)+'\n']))]
                    cadenasUinterN5=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(5)+", "+str(control)+'\n']))]
                    cadenasUinterN6=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(6)+", "+str(control)+'\n']))]
                    cadenasUinterN7=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(7)+", "+str(control)+'\n']))]
                    cadenasUinterN8=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(8)+", "+str(control)+'\n']))]
                    cadenasUinterN9=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(9)+", "+str(control)+'\n']))]
                    cadenasUinterN10=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(10)+", "+str(control)+'\n']))]
                    cadenasUinterN11=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(11)+", "+str(control)+'\n']))]
    
                    Uinterprops=[cadenasUinterN1,cadenasUinterN11]    
                    cadenaSDVcontact=['CDISP, CSTRESS, SDV\n']
                    cadenanombreSDV='ASSEMBLY_TOP_SURF/ASSEMBLY_BOTTOM_SURF'
                    
                    #Buscamos uinter='*Surface Interaction'
                    cadenasCohesivo=lineasdefinteraccion(nameinp_k_mas_1+'.inp', '*Surface Interaction',3)  #buscamos la cadena a cambiar
                    #reemplazamos la primera linea en los dos archivos de inicio
                    replacement(nameinp_k_mas_1+'.inp', cadenasCohesivo[0], cadenasUinter[0])
                    #reemplazamos la segunda linea en los dos archivos de inicio
                    replacement(nameinp_k_mas_1+'.inp', cadenasCohesivo[1], cadenasUinter[1])
                    #reemplazamos la tercera linea en los dos archivos de inicio. Esta es distinta
                    replacement(nameinp_k_mas_1+'.inp', cadenasCohesivo[2], cadenasUinterN1[0][:])
                    #reemplazamos la linea que pide las salida de las superficies de contacto    
                    cadenaFOcontact=lineasdefinteraccion(nameinp_k_mas_1+'.inp', 'CDISP, CSTRESS,',1)  
                    replacement(nameinp_k_mas_1+'.inp', cadenaFOcontact[0], cadenaSDVcontact[0])
                    # replacement(nameinp_k_mas_1+'.inp', cadenaFOcontact[0], cadenaSDVcontact[0])
                    
                    # job preparation:
                    start = newload + float(0.025*newload)
                    end = 5.25*newload
                    num = 25
                    power = float(1)  # Bias towards the start
                    # 1. Create a normalized linear space (0 to 1)
                    if num == 1:
                        linear_norm = [start]
                    else:
                        linear_norm = [float(i) / (num - 1) for i in range(num)]
        
                    # Apply the power and scale, all inside the comprehension
                    load_vector = [start + (end - start) * (norm_val ** power) for norm_val in linear_norm]
                    print("DEBUG VECTOR:", load_vector)
                    print("--- Launching Abaqus job with Abaqus/CAE (mdb) ---")
                    worker_script_bat = 'run_job.bat'
                    iternum = 1
                    # *************************************************************************
                    # NOTE: for any crack propagation step inside the loop, the load amplitude
                    # values must be the same across all increments of the step, i.e.,
                    # ((0,value),(1.0,value)), consistent with crack initiation steps,
                    # which will enable a most precise crack advance load prediction.
                    # *************************************************************************
                    def PMTESCcritEneSubrV2(name_files,control,cadenanombreSDV, workdir):
                        chdir(workdir)
                        odb = openOdb(name_files + '.odb')    
                        odbv = session.openOdb(name_files + '.odb')
                        myAssembly = odb.rootAssembly
                        #--------------------------------------------------------------------------
                        #||||||||||||| DEFINICION DE CAMINOS DE VARIABLES EN ABAQUS         |||||||||||||||
                        #--------------------------------------------------------------------------  
                        key_step = odb.steps.keys()
                        lastFrame = odb.steps[key_step[0]].frames[-1]
                        #path_despl=odb.steps[key_step[0]].frames[-1].fieldOutputs['U'].\
                        #    getSubset(region=InterNodeSet).values
                        path_SDV1=odb.steps[key_step[0]].frames[-1].fieldOutputs['SDV1     '+cadenanombreSDV].values
                        nodeglob=odb.steps[key_step[0]].frames[-1].fieldOutputs['SDV2     '+cadenanombreSDV].values
                        # path_SDV7=odb.steps[key_step[0]].frames[-1].fieldOutputs['SDV7     '+cadenanombreSDV].values  
                        # path_SDV8=odb.steps[key_step[0]].frames[-1].fieldOutputs ['SDV8     '+cadenanombreSDV].values 
                        path_SDV11=odb.steps[key_step[0]].frames[-1].fieldOutputs['SDV11    '+cadenanombreSDV].values
                        path_SDV12=odb.steps[key_step[0]].frames[-1].fieldOutputs['SDV12    '+cadenanombreSDV].values
                        t2tc=odb.steps[key_step[0]].frames[-1].fieldOutputs['SDV21    '+cadenanombreSDV].values
                        tcrit=odb.steps[key_step[0]].frames[-1].fieldOutputs['SDV22    '+cadenanombreSDV].values
                        GtotE=odb.steps[key_step[0]].frames[-1].fieldOutputs['SDV7     '+cadenanombreSDV].values
                        GcE=odb.steps[key_step[0]].frames[-1].fieldOutputs['SDV8     '+cadenanombreSDV].values
                        # damT=odb.steps[key_step[0]].frames[-1].fieldOutputs['SDV9     '+cadenanombreSDV].values
                        
                        # nodecnt=0 #node counter of nodes which define stress condition
                        Asigma_lst=[] #empty list where to store the interface vars. that define Asigma
                        # for loop of interface nodes:
                        for node in range(len(nodeglob)):
                            if t2tc[node].data>1.0:
                                # nodect+=1
                                Asigma_lst.append([nodeglob[node].data,t2tc[node].data,tcrit[node].data,GtotE[node].data])
                        if len(Asigma_lst)!=0:
                            # print(np.array(Asigma_lst))
                            # ***RANDOM SHUFFLING OF ROWS IN THE Asigma_lst LIST***:
                            # Asetcrit=random.sample()
                            # convert the list into a numpy data array:
                            Asetcritres=np.array(Asigma_lst)
                            # numpy random normal distribution:
                            p=0.0 # 0, 0.1, 0.2, 0.3,...constant of multiplication
                            sigma=1 #standard deviation: 0.25, 0.5, 1,...
                                
                            chirandomv2 = np.random.normal(0,sigma,len(Asigma_lst))
                            tc_rand = Asetcritres[:,2]*(np.ones(len(Asetcritres[:,2])) + float(p)*chirandomv2)
                            tc_rand[tc_rand < 0]=1E-5
                            t2tcx = (Asetcritres[:,1]*Asetcritres[:,2])/(tc_rand)
                            # print(t2tcx, len(t2tcx))
                            # Make a new Asetcrit array and concatenate the t2tcx 1D vector on the last
                            # column; then convert the array into a list and sort it using the new t/tc data
                            Asetcrit_new = np.concatenate([Asetcritres, t2tcx.reshape(len(t2tcx),1)], axis=1)
                            Asetcrit=list(Asetcrit_new)
                            # SORT THE UPDATED LIST IN DESCENDING ORDER OF T_X/TC_RAND
                            Asetcrit.sort(key=lambda Asetcrit:Asetcrit[:][-1],reverse=True)
                            # Save the array as a text file:
                                # Specify the format for each column
                            fmt = ['%d', '%.15e', '%.15e', '%.6e', '%.15e']
                            # header=''
                            # np.savetxt(name_files+'Asetcrit.txt', Asetcrit, fmt=fmt, delimiter='\t')
                        else: print('A_sigma set is empty')
                        
                        # concentrated force OR disp. extraction:
                        dcbtoparm = odb.rootAssembly.nodeSets['TOPN']
                        # control point dispalcement extraction:
                        if control==0:
                            tf_field=lastFrame.fieldOutputs['U'].getScalarField(componentLabel='U2',)
                            dcbtf2_field=tf_field.getSubset(region=dcbtoparm,position=NODAL).values[0].data
                        elif control!=0:
                            # concentrated force extraction:
                            tf_field=lastFrame.fieldOutputs['TF'].getScalarField(componentLabel='TF2',)
                            tf_fieldlst = [fv.data for fv in tf_field.values]
                            # tf_field=tf_field.getSubset(region=dcbtoparm,position=NODAL).values[0].data
                            dcbtf2_field = sum(tf_fieldlst)
                        #para sacar la energia	
                        key_HisReg = odb.steps[key_step[0]].historyRegions.keys()
                        path_HisRegEne = odb.steps[key_step[0]].historyRegions[key_HisReg[0]]
                        
                        sumaenerinter=0.00
                        NeL2damagetotal=[]
                        NeL2damageK=[]
                        
                        for k in range(len(path_SDV1)):
                            damage=path_SDV1[k].data
                            NeG=path_SDV1[k].nodeLabel
                            damageKmenos1=path_SDV11[k].data
                            GtE=GtotE[k].data
                            GcrE=GcE[k].data
                            areacontacto=path_SDV12[k].data
                            sumaenerinter=sumaenerinter+(GtE*areacontacto)+(GcrE*areacontacto)
                            if damage==1.0:
                                NeL2damagetotal.append(NeG)
                                if damageKmenos1==0.0:          
                                    NeL2damageK.append(NeG)
                                
                        #--------------------------------------------------------------------------
                        #Energia de todo el sistema en el setp=j
                        whole_model_region = odb.steps[key_step[0]].historyRegions['Assembly ASSEMBLY']
                        allwk_history = whole_model_region.historyOutputs['ALLWK']
                        allse_history = whole_model_region.historyOutputs['ALLSE']
                        # Work = 2*(path_HisRegEne.historyOutputs['ALLWK'].data[-1][1])
                        # eneDefSolidos = path_HisRegEne.historyOutputs['ALLSE'].data[-1][1]
                        #tomamos 'ALLSE' en lugar de 'ALLIE' porque estrictamente es la energia de deformacion
                        #pero podemos replantearnos poner la otra por el hourglassing
                        # if control==1:
                        #     eneHtotal = eneDefSolidos+sumaenerinter-Work
                        # else:
                        #     eneHtotal = eneDefSolidos+sumaenerinter
                    #
                        # Loop over all steps in the odb
                        for step_name in odb.steps.keys():
                            # Initialize a list to store the results
                            energy_data = []
                            # print(f"Processing Step: {step_name}...")
                            step = odb.steps[step_name]
                            # Get the history region for the whole model from the STEP object
                            # This is where whole-model energy variables are stored.
                            whole_model_region = step.historyRegions['Assembly ASSEMBLY']
                            allwk_history = whole_model_region.historyOutputs['ALLWK']
                            allse_history = whole_model_region.historyOutputs['ALLSE']
                            
                            Uel_frame=[]
                            PI_frame=[]
                            deltaR_frame=[]
                            Ut_frame=[]
                            
                            # Extract total potential energy at the end of first increment:
                            for frame in step.frames:
                                currentFrame = step.frames[frame.frameId]
                                if frame.frameId==1:
                                    # initialize interface energy to zero at each increment:
                                    interf_strnenergy=float(0)
                                    interf_energdiss=float(0)
                                    # initialize fracture surface defined by stress condition:
                                    area_stressc=float(0)
                                    atotal=float(0)
                                    # initialize fracture surface fulfilling the coupled criterion FFM:
                                    area_ccffm=float(0)
                                    
                                    path_SDV1=currentFrame.fieldOutputs['SDV1     '+cadenanombreSDV].values
                                    GtotE=currentFrame.fieldOutputs['SDV7     '+cadenanombreSDV].values
                                    GcE=currentFrame.fieldOutputs['SDV8     '+cadenanombreSDV].values
                                    damT=currentFrame.fieldOutputs['SDV9     '+cadenanombreSDV].values
                                    damE=currentFrame.fieldOutputs['SDV10    '+cadenanombreSDV].values
                                    allse_value = allse_history.data[frame.frameId][1]
                                    ut_value = float(0)
                                    tforce_field = []
                                    ut_field=currentFrame.fieldOutputs['U'].getScalarField(componentLabel='U2',)
                                    ut_value=abs(ut_field.getSubset(region=dcbtoparm,position=NODAL).values[0].data)
                                    Ut_frame.append(ut_value)
                                    for k in range(len(path_SDV1)):
                                        GtE=GtotE[k].data
                                        GcrE=GcE[k].data
                                        damTval=damT[k].data
                                        damEval=damE[k].data
                                        areacontacto=path_SDV12[k].data
                                        atotal=areacontacto+atotal
                                        # if path_SDV1[k].data==float(0): #active spring
                                        interf_strnenergy=interf_strnenergy + (GtE*areacontacto)
                                        # else: #broken spring
                                        #     interf_strnenergy=interf_strnenergy
                                        # Stress criterion fulfilled:
                                        if damTval==1.:
                                            area_stressc = area_stressc + areacontacto
                                        # CCFFM criterion:
                                        if damTval==1. and damEval==1.:
                                            area_ccffm = area_ccffm + areacontacto
                                            interf_energdiss=interf_energdiss + (GcrE*areacontacto)
                                    # concentrated force extraction:
                                    tforce=currentFrame.fieldOutputs['TF'].getScalarField(componentLabel='TF2',)
                                    tforce_field = [fv.data for fv in tforce.getSubset(region=dcbtoparm,position=NODAL).values]
                                    totalf_field = 2*abs(sum(tforce_field))
                                    # at frame=1, elastic strain energy='ALLWK', which is also valid
                                        # for a geometric nonlinear static analysis
                                    if control==1:
                                        PI0 = (allwk_history.data[frame.frameId][1]) - (totalf_field*ut_value)
                                    else:
                                        PI0 = (allwk_history.data[frame.frameId][1])
                                    print('PI_0: ',PI0)
                                    
                            # Loop over all frames (increments) in the current step:
                            for frame in step.frames:
                                currentFrame = step.frames[frame.frameId]
                                # initialize interface energy to zero at each increment:
                                interf_strnenergy=float(0)
                                interf_energdiss=float(0)
                                # initialize fracture surface defined by stress condition:
                                area_stressc=float(0)
                                atotal=float(0)
                                # initialize fracture surface fulfilling the coupled criterion FFM:
                                area_ccffm=float(0)
                                path_SDV1=currentFrame.fieldOutputs['SDV1     '+cadenanombreSDV].values
                                # nodeglob=frame.fieldOutputs['SDV2     '+cadenanombreSDV].values
                                # path_SDV11=frame.fieldOutputs['SDV11    '+cadenanombreSDV].values
                                path_SDV12=currentFrame.fieldOutputs['SDV12    '+cadenanombreSDV].values
                                # t2tc=frame.fieldOutputs['SDV21    '+cadenanombreSDV].values
                                # tcrit=frame.fieldOutputs['SDV22    '+cadenanombreSDV].values
                                GtotE=currentFrame.fieldOutputs['SDV7     '+cadenanombreSDV].values
                                GcE=currentFrame.fieldOutputs['SDV8     '+cadenanombreSDV].values
                                damT=currentFrame.fieldOutputs['SDV9     '+cadenanombreSDV].values
                                damE=currentFrame.fieldOutputs['SDV10    '+cadenanombreSDV].values
                                for k in range(len(path_SDV1)):
                                    GtE=GtotE[k].data
                                    GcrE=GcE[k].data
                                    damTval=damT[k].data
                                    damEval=damE[k].data
                                    areacontacto=path_SDV12[k].data
                                    atotal=areacontacto+atotal
                                    # if path_SDV1[k].data==float(0): #active spring
                                    interf_strnenergy=interf_strnenergy + (GtE*areacontacto)
                                    # else: #broken spring
                                    #     interf_strnenergy=interf_strnenergy
                                        
                                    # Stress criterion fulfilled:
                                    if damTval==1.:
                                        area_stressc = area_stressc + areacontacto
                                    # CCFFM criterion:
                                    if damTval==1. and damEval==1.:
                                        area_ccffm = area_ccffm + areacontacto
                                        interf_energdiss=interf_energdiss + (GcrE*areacontacto)
                                        # interf_strnenergy=interf_strnenergy - (GcrE*areacontacto)
                                        # print(path_SDV1[k].data,interf_strnenergy,interf_energdiss)
                                #
                                allse_value = allse_history.data[frame.frameId][1]
                                ut_value = float(0)
                                tforce_field = []
                                ut_field=currentFrame.fieldOutputs['U'].getScalarField(componentLabel='U2',)
                                ut_value=abs(ut_field.getSubset(region=dcbtoparm,position=NODAL).values[0].data)
                                Ut_frame.append(ut_value)
                                
                                # concentrated force extraction:
                                tforce=currentFrame.fieldOutputs['TF'].getScalarField(componentLabel='TF2',)
                                tforce_field = [fv.data for fv in tforce.getSubset(region=dcbtoparm,position=NODAL).values]
                                totalf_field = 2*abs(sum(tforce_field))
                                
                                deltaR_frame.append(interf_energdiss)
                                delR=deltaR_frame[-1]
                                if control==1:
                                    whole_model_region = step.historyRegions['Assembly ASSEMBLY']
                                    allwk_history = whole_model_region.historyOutputs['ALLWK']
                                    allse_history = whole_model_region.historyOutputs['ALLSE']
                                    # interf_strnenergy = allwk_history.data[frame.frameId][1] - allse_history.data[frame.frameId][1]
                                    # PIf = (allse_history.data[frame.frameId][1]+interf_strnenergy) - (totalf_field*ut_value)
                                    if frame.frameId<=1:
                                        wpot = totalf_field*ut_value
                                        totalPotenergy = (allse_history.data[frame.frameId][1]+interf_strnenergy)-wpot
                                        # totalPotenergy = -0.5*(wpot)
                                        # totalPot0 = (allwk_history.data[frame.frameId][1])-wpot
                                        PI_frame.append(totalPotenergy)
                                        if frame.frameId==1 and PI_frame[0]!=float(0):
                                            PI0 = PI_frame[0]
                                            delPI=PI_frame[-1]-PI_frame[0]
                                        else: delPI = float(0); delR=float(0)
                                        
                                    else:
                                        wpot = totalf_field*ut_value
                                        Uel_f = (allse_history.data[frame.frameId][1]+interf_strnenergy)
                                        totalPotenergy = Uel_f - wpot
                                        # totalPotenergy = -0.5*(wpot)
                                        PI_frame.append(totalPotenergy)
                                        delPI = totalPotenergy - PI0
                                    Uel_frame.append(allse_history.data[frame.frameId][1]+interf_strnenergy)
                                    
                                    # print(PI_frame[-1],(allse_history.data[frame.frameId][1]+interf_strnenergy),-wpot)
                                    energy_data.append([frame.frameId, totalPotenergy, allse_value, delPI, delR, delPI+delR, area_stressc, area_ccffm, 
                                                        sum(tforce_field), ut_value])
                                        
                                else: #fixed-grips control:
                                    # interf_strnenergy = allwk_history.data[frame.frameId][1] - allse_history.data[frame.frameId][1]
                                    PIf = (allse_history.data[frame.frameId][1]+interf_strnenergy)
                                    if frame.frameId==1: PI_frame.append(PIf)
                                    else:   PI_frame.append(PIf)
                                    
                                    Uel_frame.append(allwk_history.data[frame.frameId][1])
                                    if frame.frameId<=1:
                                        if frame.frameId==1 and PI_frame[0]!=float(0):
                                            delPI=PI_frame[-1]-PI_frame[0]
                                        else: delPI = float(0); delR = float(0)
                                        # else: delPI=PIf-PI0
                                    elif frame.frameId>1 and PI_frame[0]!=float(0):
                                        delPI=PI_frame[-1]-PI_frame[0]
                                    elif frame.frameId>1 and PI_frame[0]==float(0):
                                        delPI=PIf-PI0
                                    energy_data.append([frame.frameId, PI_frame[-1], allse_value, delPI, delR, delPI+delR, area_stressc, area_ccffm,
                                                        sum(tforce_field), ut_value])
                                
                                # The .data attribute for these objects is a tuple, e.g., (time, value)
                                # We are interested in the value, which is at index 1.
                                if frame.frameId<=1: 
                                    # interf_strnenergy = allwk_history.data[frame.frameId][1] - allse_history.data[frame.frameId][1]
                                    allwk0 = 1*allwk_history.data[frame.frameId][1]
                                    allse0 = allse_history.data[frame.frameId][1]
                                    if control==1:
                                        totalPot0 = (allwk0) - (totalf_field*ut_value)
                                        # totalPot0 = -0.5*(wpot)
                                    else:
                                        totalPot0 = allwk0
                                    
                                    # Append the extracted data as a list to our main list
                                    # energy_data.append([frame.frameId, totalPot0, allse0, delPI, deltaR_frame[frame.frameId], 0.0, area_stressc, area_ccffm])
                                    # print(frame.frameId,totalPot0,atotal,area_stressc,area_ccffm,interf_energdiss,interf_strnenergy)
                                # else:
                                #     ut_field=currentFrame.fieldOutputs['U'].getScalarField(componentLabel='U2',)
                                #     ut_value=abs(ut_field.getSubset(region=dcbtoparm,position=NODAL).values[0].data)
                                #     # concentrated force extraction:
                                #     tforce=currentFrame.fieldOutputs['TF'].getScalarField(componentLabel='TF2',)
                                #     tforce_field = [fv.data for fv in tforce.getSubset(region=dcbtoparm,position=NODAL).values]
                                #     totalf_field = 2*abs(sum(tforce_field))
                                    
                                #     allwk_value = 1*allwk_history.data[frame.frameId][1]
                                #     allse_value = allse_history.data[frame.frameId][1]
                                #     
                                
                                    # print(frame.frameId, totalf_field, ut_value, PI_frame[-1], delPI)
                                    # Append the extracted data as a list
                                    # energy_data.append([frame.frameId, PI_frame[-1], allse_value, delPI, delR, delPI+delR, area_stressc, area_ccffm])
                                    # print(frame.frameId,PI_frame[-1],atotal,area_stressc,area_ccffm,interf_energdiss,interf_strnenergy)
                            # save the numpy energy evolution as a .txt file:
                            tbl_format=['%d', '%.10e', '%.10e', '%.12e', '%.12e', '%.12e', '%.6e', '%.6e', '%.8e', '%.8e']
                            np.savetxt(name_files+'_delPIdelRevol.txt', np.array(energy_data), delimiter='\t', fmt=tbl_format)
                            # print('Total energy evolution array:')
                            # print(np.array(energy_data))
                        import displayGroupOdbToolset as dgo
                        # Get the viewport
                        viewport = session.viewports[session.currentViewportName]
                        viewport.setValues(displayedObject=odbv)
                        session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
                        session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
                            CONTOURS_ON_DEF, ))
                        session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(
                            curveRefinementLevel=EXTRA_FINE)
                        session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(
                            deformationScaling=UNIFORM, uniformScaleFactor=1)
                        session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                            variableLabel='SDV1     ASSEMBLY_TOP_SURF/ASSEMBLY_BOTTOM_SURF', 
                            outputPosition=NODAL, )
                        # session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                        #     variableLabel='SDV10    ASSEMBLY_TOP_SURF/ASSEMBLY_BOTTOM_SURF', 
                        #     outputPosition=NODAL, )
                        leaf = dgo.LeafFromPartInstance(partInstanceName=('BOTTOM ARM-1', ))
                        session.viewports['Viewport: 1'].odbDisplay.displayGroup.remove(leaf=leaf)
                        session.viewports['Viewport: 1'].view.setProjection(projection=PARALLEL)
                        session.viewports['Viewport: 1'].view.setValues(session.views['Bottom'])
                        session.viewports['Viewport: 1'].viewportAnnotationOptions.setValues(
                            compass=OFF)
                        session.viewports['Viewport: 1'].odbDisplay.contourOptions.setValues(
                            numIntervals=2, maxValue=0, minValue=0)
                        session.viewports['Viewport: 1'].view.setValues(nearPlane=493.329, 
                            farPlane=508.291, width=71.6097, height=24.2584, cameraPosition=(
                            219.901, -499.005, 5.99311), cameraTarget=(219.901, 1.80515, 5.99311))
                        session.pngOptions.setValues(imageSize=(1600, 712))
                        session.printOptions.setValues(vpDecorations=OFF, vpBackground=ON)
                        session.printToFile(
                            fileName=str(name_files)+'.png', 
                            format=PNG, canvasObjects=(session.viewports['Viewport: 1'], ))
                        odb.close()
                        return energy_data[-1][5],interf_strnenergy,NeL2damageK,NeL2damagetotal,dcbtf2_field,energy_data
                    for index, val in enumerate(load_vector):
                        os.chdir(actual_directory)
                        # if index==0:
                        #     load_vector = [start + (end - start) * (norm_val ** power) for norm_val in linear_norm]
                        print("Index: {0}, Value: {1}".format(index, val))
                        # define the new .inp and directory name
                        nameinp_k_mas_1=odbdef+'_crackprop'+'_m'+str(iternum)
                        # modify the job name variable at each iteration:
                        job_name = nameinp_k_mas_1
                        # 1. Generate a temporary name for this specific iteration's model
                        temp_model_name = "{}_Iter_{}".format(nameinpK, index)
                        
                        # 2. DELETE temp model if it was left over from a crash
                        if temp_model_name in mdb.models:
                            del mdb.models[temp_model_name]
                            
                        # 3. COPY the Master Model to the Temp Model
                        # This creates a perfect clone that we can modify safely
                        new_model = mdb.Model(name=temp_model_name, objectToCopy=mdb.models[nameinpK])
                        if index==0:
                            prev_load=newload
                        elif index>0:
                            prev_load=load_vector[index-1]
                        
                        new_model.amplitudes[amplname].setValues(
                            timeSpan=STEP, 
                            smooth=SOLVER_DEFAULT, 
                            data=((0.0, load_vector[index]), (1.0, load_vector[index]))
                        )
                        # CREATE JOB using the TEMP MODEL
                        # 
                        # (or use the standard name, but point it to the new model)
                        # use the exact same job name:
                        if nameinp_k_mas_1 in mdb.jobs: del mdb.jobs[nameinp_k_mas_1]
                        job_name = nameinp_k_mas_1
                        
                        if index == 0:
                            print ("--- Starting Iteration: {} ---".format(job_name))                       
                            # write the job for the first crack propagation step: 
                            mdb.Job(name=nameinp_k_mas_1, model=temp_model_name, description='crack propagation analysis, m1', 
                                    type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, 
                                    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, memory=memoria,
                                    explicitPrecision=SINGLE, nodalOutputPrecision=FULL, echoPrint=OFF, 
                                    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine=nameSUBR, 
                                    scratch=work_dirkm, resultsFormat=ODB, multiprocessingMode=DEFAULT, numDomains=cpus, numCpus=cpus, 
                                    numGPUs=0)
                            # write the input .inp file:
                            new_desc = 'crack propagation analysis, m{}, load={}'.format(iternum, load_vector[index])
                            mdb.jobs[nameinp_k_mas_1].setValues(description=new_desc)
                            mdb.jobs[nameinp_k_mas_1].writeInput(consistencyChecking=OFF)
                            
                            #Buscamos uinter='*Surface Interaction'
                            cadenasCohesivo=lineasdefinteraccion(nameinp_k_mas_1+'.inp', '*Surface Interaction',3)  #buscamos la cadena a cambiar
                            #reemplazamos la primera linea en los dos archivos de inicio
                            replacement(nameinp_k_mas_1+'.inp', cadenasCohesivo[0], cadenasUinter[0])
                            #reemplazamos la segunda linea en los dos archivos de inicio
                            replacement(nameinp_k_mas_1+'.inp', cadenasCohesivo[1], cadenasUinter[1])
                            #reemplazamos la tercera linea en los dos archivos de inicio. Esta es distinta
                            replacement(nameinp_k_mas_1+'.inp', cadenasCohesivo[2], cadenasUinterN11[0][:])
                            #reemplazamos la linea que pide las salida de las superficies de contacto    
                            cadenaFOcontact=lineasdefinteraccion(nameinp_k_mas_1+'.inp', 'CDISP, CSTRESS',1)  
                            replacement(nameinp_k_mas_1+'.inp', cadenaFOcontact[0], cadenaSDVcontact[0])
                            # Add .fil output requests for the input file:
                            add_urdfill_output(nameinp_k_mas_1+'.inp')
                            # copy .inp file into the propagation directory
                            shutil.copy(nameinp_k_mas_1+'.inp',work_dirkm)
                            if not os.path.exists(work_dirkm):
                                os.makedirs(work_dirkm)
                                # print ("Created directory: {}".format(work_dirkm))
                                work_dirkm = os.path.normpath(work_dirkm)
                                print ("Normalized path: {}".format(work_dirkm))
                            else:
                                print ("Directory already exists: {}".format(work_dirkm))
                                work_dirkm = os.path.normpath(work_dirkm)
                                print ("Normalized path: {}".format(work_dirkm))
                            # --- DELETE THE STALE JOB (CRITICAL STEP) ---
                            # otherwise, mdb.Job(...) returns the OLD snapshot.
                            if nameinp_k_mas_1 in mdb.jobs:
                                del mdb.jobs[nameinp_k_mas_1]
                        else:
                            
                            # copy the original .inp file
                            copy(nameinp+'.inp', nameinp_k_mas_1+'.inp')
                            #cargamos el input original de ABAQUS
                            mdb.ModelFromInputFile(inputFileName=nameinp_k_mas_1+'.inp', name=temp_model_name)
                            
                            # -----------------------------------------------------
                            # 1. 
                            mod_name = mdb.models[temp_model_name]
                            # 2. Get a list of all existing Predefined Fields (Initial State, Temp, etc.)
                            # We use .keys() to get a static list so we don't modify the container while iterating
                            p_field_keys = mod_name.predefinedFields.keys()
                            # 3. Loop through and delete them
                            for key in p_field_keys:
                                print ("Cleaning up old Predefined Field: {}".format(key))
                                del mod_name.predefinedFields[key]
                            # update the amplitude field
                            # mdb.models[nameinpK].amplitudes[amplname].setValues(timeSpan=STEP, 
                            #             smooth=SOLVER_DEFAULT, data=((0.0, load_vector[index-1]), (1.0, val)))
                            
                            # DELETE the old job from the last iteration
                            #    (This is necessary if you re-create it)
                            if nameinp_k_mas_1 in mdb.jobs:
                                print("Deletion of job: {}".format(nameinp_k_mas_1))
                                del mdb.jobs[nameinp_k_mas_1]
                            
                            # -----------------------------------------------------
                            #imponemos el estado inicial del job ganador entre las n para que sea el inicio del paso anterior
                            instancesname=mdb.models[temp_model_name].rootAssembly.instances.keys()
                            instanceslista=[]
                            # Populate the list with the model instances
                            for inst in range(len(instancesname)):
                                instanceslista.append(mdb.models[temp_model_name].rootAssembly.instances[instancesname[inst]])
                            # apply the updated initial state to the model database 
                            mdb.models[temp_model_name].InitialState(updateReferenceConfiguration=OFF, 
                                      fileName=odbdef+'_crackprop'+'_m'+str(iternum-1), 
                                      endStep=LAST_STEP, endIncrement=STEP_END, name='Predefined Field-'+str(iternum), 
                                      createStepName='Initial', instances=tuple(instanceslista))
                            # 
                            mdb.models[temp_model_name].steps.changeKey(fromName='Step-'+str(1), toName='Step-'+str(1)+'_m'+str(iternum))
                            mdb.Job(name=nameinp_k_mas_1, model=temp_model_name, description='crack propagation analysis, m{}, load={}'.format(iternum, load_vector[index]),
                                    type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, 
                                    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, memory=memoria,
                                    explicitPrecision=SINGLE, nodalOutputPrecision=FULL, echoPrint=OFF, 
                                    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine=nameSUBR, 
                                    scratch=work_dirkm, resultsFormat=ODB, multiprocessingMode=DEFAULT, numDomains=cpus, numCpus=cpus, 
                                    numGPUs=0)
                            
                            # write the input .inp file:
                            mdb.jobs[nameinp_k_mas_1].writeInput(consistencyChecking=OFF)
                            
                            cadenasUinter=['*Surface Interaction, name=IntProp-1, user, depvar=22, properties=9\n','\n']
                            kn = 300; kt = kn*1
                            cadenasUinterN=[kn,kt,kt,0.09375,0.0,0.0,8.]
                            # 11 starting damage configurations:
                            cadenasUinterN1=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(1)+", "+str(control)+'\n']))]
                            cadenasUinterN2=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(2)+", "+str(control)+'\n']))]
                            cadenasUinterN3=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(3)+", "+str(control)+'\n']))]
                            cadenasUinterN4=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(4)+", "+str(control)+'\n']))]
                            cadenasUinterN5=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(5)+", "+str(control)+'\n']))]
                            cadenasUinterN6=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(6)+", "+str(control)+'\n']))]
                            cadenasUinterN7=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(7)+", "+str(control)+'\n']))]
                            cadenasUinterN8=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(8)+", "+str(control)+'\n']))]
                            cadenasUinterN9=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(9)+", "+str(control)+'\n']))]
                            cadenasUinterN10=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(10)+", "+str(control)+'\n']))]
                            cadenasUinterN11=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(11)+", "+str(control)+'\n']))]
                            
                            Uinterprops=[cadenasUinterN1,cadenasUinterN11]    
                            cadenaSDVcontact=['CDISP, CSTRESS, SDV\n']
                            cadenanombreSDV='ASSEMBLY_TOP_SURF/ASSEMBLY_BOTTOM_SURF'
                            
                            #Buscamos uinter='*Surface Interaction'
                            cadenasCohesivo=lineasdefinteraccion(nameinp_k_mas_1+'.inp', '*Surface Interaction',3)  #buscamos la cadena a cambiar
                            #reemplazamos la primera linea en los dos archivos de inicio
                            replacement(nameinp_k_mas_1+'.inp', cadenasCohesivo[0], cadenasUinter[0])
                            #reemplazamos la segunda linea en los dos archivos de inicio
                            replacement(nameinp_k_mas_1+'.inp', cadenasCohesivo[1], cadenasUinter[1])
                            #reemplazamos la tercera linea en los dos archivos de inicio. Esta es distinta
                            replacement(nameinp_k_mas_1+'.inp', cadenasCohesivo[2], cadenasUinterN1[0][:])
                            #reemplazamos la linea que pide las salida de las superficies de contacto    
                            cadenaFOcontact=lineasdefinteraccion(nameinp_k_mas_1+'.inp', 'CDISP, CSTRESS',1)  
                            replacement(nameinp_k_mas_1+'.inp', cadenaFOcontact[0], cadenaSDVcontact[0])
                            
                            # Add .fil output requests for the input file:
                            add_urdfill_output(nameinp_k_mas_1+'.inp')
                            # Call the function to modify the .inp file with the new load values
                            modify_amplitude_dataV2(nameinp_k_mas_1+'.inp', load_vector[index], load_vector[index])
                            # copy the .inp file into the propagation dir.
                            shutil.copy(nameinp_k_mas_1+'.inp', work_dirkm)
                            os.chdir(actual_directory)
                        # CLEANUP
                        # Delete the temp model to free up memory
                        del mdb.models[temp_model_name]
                        
                        # Define files to monitor
                        lock_file = os.path.join(work_dirkm, job_name + '.lck')
                        
                        # ... (your pre-check to remove old .lck files) ...
                        if os.path.exists(lock_file):
                            os.remove(lock_file)
                                
                        # --- THIS IS THE MODIFIED PART ---
                        # Call the .bat file and pass the job_name as an argument
                        command = [
                            worker_script_bat,
                            job_name,
                            work_dirkm # work directory argument
                        ]
                        
                        print ("Submitting job via Worker .bat file...")
                        # We use .call() because the .bat file is fast
                        # and will exit as soon as the job is submitted.
                        subprocess.call(command, cwd=work_dirkm)
                        # --- END OF MODIFIED PART ---
                        
                        # --- Rest of your script ---
                        
                        # 1. Wait for the job to start
                        # print ("Waiting for job to start (for .lck file)...")
                        time_to_start = 0
                        # ... (rest of your .lck polling logic) ...
                        # while not os.path.exists(lock_file):
                        #     time.sleep(1)
                        #     # ... (etc) ...
                            
                        # 2. Wait for the job to finish
                        # print ("Job is running. Polling {}.lck...".format(job_name))
                        # while os.path.exists(lock_file):
                        #     time.sleep(10)
                        print ("--- Iteration {} complete. ---".format(job_name))
                        # ADVANCE ITERATION POSTPROCESS
                        # change to the crack advance working dir.
                        chdir(work_dirkm)
                        enerHtotalN,sumaenerinterN,NeL2damageKN,NeL2damagetotalN,tf,energy_evol=PMTESCcritEneSubrV2(job_name,control,cadenanombreSDV,work_dirkm)
                        model_outputfile(actual_directory,sim_file,k,11,val,tf,enerHtotalN,sumaenerinterN,len(NeL2damageKN),len(NeL2damagetotalN),energy_evol)
                        datos_salida.append([float(val),float(tf),float(enerHtotalN),sumaenerinterN,len(NeL2damageKN),job_name,len(NeL2damagetotalN)])
                        stepdata.append([k,val,tf,enerHtotalN,sumaenerinterN,len(NeL2damageKN),len(NeL2damagetotalN),energy_evol[-1][5]])
                        KDn_outputfile(actual_directory,Dirname,k,11,enerHtotalN,sumaenerinterN,len(NeL2damageKN),tf,energy_evol)
                        # Delete unnecessary simulation files:
                        # 
                        inp_file = os.path.join(work_dirkm, job_name + '.inp')
                        mdl_file = os.path.join(work_dirkm, job_name + '.mdl')
                        prt_file = os.path.join(work_dirkm, job_name + '.prt')
                        res_file = os.path.join(work_dirkm, job_name + '.res')
                        odb_file = os.path.join(work_dirkm, job_name + '.odb')
                        stt_file = os.path.join(work_dirkm, job_name + '.stt')
                        delete_list = [inp_file]
                        for item in delete_list:
                            os.remove(item)
                        iternum += 1
                        
                    print("--- Sequential analysis complete. ---")
                    
                    chdir(work_dirkm)
                    # List of file extensions to delete
                    extensions_to_delete = ['*.dat', '*.sim', '*.com', '*.jnl', '*.stt', '*.res','*.mdl'\
                                    '*.mtx','*.pes','*.par','*.pmg','*.ipm','*.pyc','*.fil','*.prt']
                    for ext in extensions_to_delete:
                        # glob.glob finds all files matching the pattern
                        files_to_remove = glob.glob(ext)
                        for f in files_to_remove:
                            try:
                                os.remove(f)
                                "  - Deleted: {}".format(f)
                            except OSError as e:
                                "  - Error deleting file {}: {}".format(f, e)
                                
                    chdir(actual_directory)
                    if godb==1:
                        list_odbs=(['*_K'+str(k-1)+'.odb','*_K'+str(k-1)+'.inp','*_K'+str(k)+'.odb',
                                    '*_K'+str(k)+'.sta','*_K'+str(k-1)+'.sta','*_K'+str(k)+'.inp','*_crackprop'+'_m*'+'.inp'])
                        mover_archivos(list_odbs,Dirname)
                        #borramos el resto de archivos de ese paso                                          
                        list_delete = ('*.sta', '*.dat', '*.sim', '*.prt', '*.msg', '*.com', '*.jnl',\
                                       '*.mtx','*.pes','*.par','*.pmg','*.odb','*.ipm','*.pyc','*.res','*.mdl','*.stt')
                        borrar_archivos(list_delete)
                        
                    end = time.time()
                    print('Simulation has terminated in %.2f' %((end-start)/60)+' min')
                    sys.exit()

            # FORCE-CONTROL ANALYSIS
            else:
                # CRACK PROPAGATION AND ARREST STAGE: FORCE CONTROL
                # -------------------------------------------------------------
                def create_crack_restart_inpFctrl(file_name, n, displacement_mag):
                    """
                    Generates an Abaqus restart .inp file for crack propagation step 'n'.
                    
                    Args:
                        file_name (str): Name of the new .inp file (without extension).
                        n (int): The current loop iterator (step number). Starts at 2.
                        displacement_mag (float): The target displacement value for this step (e.g., 0.0395).
                    """
                    
                    # Calculate the step to read from
                    prev_step = n
                    
                    # Create unique names for this iteration to avoid conflicts
                    amp_name = "AMP-BIS-{}".format(n)
                    step_name = "Step-{}".format(n)
                    
                    full_filename = file_name + '.inp'
                    
                    print ("Writing restart file: {}".format(full_filename))
                    
                    with open(full_filename, 'w') as f:
                        # --- HEADER ---
                        f.write("*HEADING\n")
                        f.write("Crack Propagation Restart - Step {}\n".format(n))
                        
                        # --- RESTART READ ---
                        # The 'END STEP' is crucial here to reset KINC to 1
                        f.write("*RESTART, READ, STEP={}, END STEP\n".format(prev_step))
                        
                        # --- AMPLITUDE DEFINITION ---
                        # We define a NEW amplitude for this specific step.
                        # Format: Time, Value, Time, Value
                        # This creates a constant hold at 'displacement_mag' across the step (0.0 to 1.0)
                        f.write("**\n")
                        f.write("** New Amplitude for Step {}\n".format(n))
                        f.write("*AMPLITUDE, NAME={}\n".format(amp_name))
                        f.write("0.0, {}, 1.0, {}\n".format(displacement_mag, displacement_mag))
                        
                        # --- STEP DEFINITION ---
                        f.write("**\n")
                        f.write("** STEP: {}\n".format(step_name))
                        f.write("**\n")
                        f.write("*Step, name={}, nlgeom=NO, inc=25, unsymm=NO\n".format(step_name))
                        f.write("*Static, direct\n")
                        f.write("0.04, 1., \n")
                        
                        # --- TRACTION BOUNDARY CONDITIONS ---
                        f.write("**\n")
                        f.write("** LOADS\n")
                        f.write("**\n")
                        # Apply Positive Force F to TOPN
                        f.write("** Name: Load-1   Type: Concentrated force\n")
                        f.write("*Cload, amplitude={}\n".format(amp_name))
                        f.write("Topcn, 2, 1.\n")
                        f.write("** Name: Load-2   Type: Concentrated force\n")
                        f.write("*Cload, amplitude={}\n".format(amp_name))
                        f.write("Topen, 2, 0.5\n")
                        
                        # Apply Negative Force F to BOTN
                        f.write("** Name: Load-3   Type: Concentrated force\n")
                        f.write("*Cload, amplitude={}\n".format(amp_name))
                        f.write("Botcn, 2, -1.\n")
                        f.write("** Name: Load-4   Type: Concentrated force\n")
                        f.write("*Cload, amplitude={}\n".format(amp_name))
                        f.write("Boten, 2, -0.5\n")
                        
                        # --- OUTPUT REQUESTS ---
                        f.write("**\n")
                        f.write("** OUTPUT REQUESTS\n")
                        f.write("** This block is for the URDFIL subroutine\n")
                        f.write("*ENERGY FILE, FREQUENCY=1\n")
                        
                        # Important: Write restart for the NEXT loop iteration
                        f.write("*Restart, write, frequency=1\n")
                        
                        f.write("**\n")
                        f.write("** FIELD OUTPUT: F-Output-1\n")
                        f.write("**\n")
                        f.write("*Output, field\n")
                        f.write("*Node Output, nset=INITIALFRONT\n")
                        f.write("COORD, \n")
                        
                        f.write("**\n")
                        f.write("** FIELD OUTPUT: F-Output-2\n")
                        f.write("**\n")
                        f.write("*Contact Output, nset=\"TOP ARM-1\"\n")
                        f.write("CDISP, CSTRESS, SDV\n")
                        
                        f.write("**\n")
                        f.write("** FIELD OUTPUT: F-Output-3\n")
                        f.write("**\n")
                        f.write("*Element Output, elset=\"TOP ARM-1\", directions=NO\n")
                        f.write("SDV, \n")
                        
                        f.write("**\n")
                        f.write("** FIELD OUTPUT: F-Output-4\n")
                        f.write("**\n")
                        f.write("*Node Output, nset=\"TOP ARM-1\"\n")
                        f.write("U, \n")
                        
                        f.write("**\n")
                        f.write("** FIELD OUTPUT: F-Output-5\n")
                        f.write("**\n")
                        f.write("*Node Output, nset=TOPN\n")
                        f.write("CF, TF, UT\n")
                        
                        f.write("**\n")
                        f.write("** HISTORY OUTPUT: H-Output-1\n")
                        f.write("**\n")
                        f.write("*Output, history\n")
                        f.write("*Energy Output\n")
                        f.write("ALLIE, ALLSE, ALLWK\n")
                        
                        f.write("**\n")
                        f.write("** HISTORY OUTPUT: H-Output-2\n")
                        f.write("**\n")
                        f.write("*Node Output, nset=TOPN\n")
                        f.write("CF2, TF2, U2\n")
                        
                        f.write("*End Step\n")
                    return full_filename           
                chdir(actual_directory)
                # name of the Abaqus worker script
                worker_script = 'runrestart_job.bat' 
                # define the new .inp and directory name
                nameinp_k_mas_1=odbdef+'_crackprop'+'_m'+str(1)
                work_dirkm = os.path.join(actual_directory, 
                             os.path.basename(odbdef+'_crackprop'+'_m'+str(0) + '_scratch'))
                # create a new crack propagation working directory
                os.mkdir(work_dirkm)
                # copy the worker script into the new work directory
                shutil.copy(worker_script, work_dirkm)
                
                # files we need to export
                inp_file = os.path.join(odbdef + '.inp')
                mdl_file = os.path.join(work_dir+'\\'+odbdef + '.mdl')
                prt_file = os.path.join(work_dir+'\\'+odbdef + '.prt')
                res_file = os.path.join(work_dir+'\\'+odbdef + '.res')
                odb_file = os.path.join(work_dir+'\\'+odbdef + '.odb')
                stt_file = os.path.join(work_dir+'\\'+odbdef + '.stt')
                export_fil = [mdl_file,prt_file,odb_file,nameSUBR,res_file,
                              stt_file,common_file]
                # export the files into the new directory
                for f in export_fil:
                    shutil.copy(f, work_dirkm)
                """
                # *************************************************************************
                # NOTE: for any crack propagation step inside the loop, the load amplitude
                # values must be the same across all increments of the step, i.e.,
                # ((0,value),(1.0,value)), consistent with crack initiation steps,
                # which will enable a most precise crack advance load prediction;
                # should be invariant to the bisection algorithm scheme.
                # *************************************************************************
                """
                # propagation m-substeps, for loop:
                for m in range(1,11):
                    # change to the script parent dir.
                    chdir(actual_directory)
                    if m==1:
                        oldjob_name=odbdef
                    else: 
                        oldjob_name=odbdef+'_m'+str(m-1)
                    inpfil_name = odbdef+'_m'+str(m)
                    inp_filename = create_crack_restart_inpFctrl(inpfil_name, m, newload)
                    shutil.copy(inp_filename, work_dirkm)
                    os.remove(inp_filename)
                    # Change to the crack propagation dir.
                    chdir(work_dirkm)
                    # Define files to monitor
                    lock_file = os.path.join(inpfil_name + '.lck')
                    
                    # ... (your pre-check to remove old .lck files) ...
                    # if os.path.exists(lock_file):
                    #     os.remove(lock_file)
                            
                    # --- THIS IS THE MODIFIED PART ---
                    # Call the .bat file and pass the job_name as an argument
                    command = [
                        worker_script,
                        inpfil_name,
                        oldjob_name,
                        work_dirkm # work directory argument
                    ]
                    
                    print("Submitting job via Worker .bat file...")
                    # We use .call() because the .bat file is fast
                    # and will exit as soon as the job is submitted.
                    subprocess.call(command, cwd=work_dirkm)
                    # --- Rest of your script ---
                    
                    # 1. Wait for the job to start
                    # print ("Waiting for job to start (for .lck file)...")
                    # time_to_start = 0
                    # ... (rest of your .lck polling logic) ...
                    # while not os.path.exists(lock_file):
                    #     time.sleep(1)
                        # ... (etc) ...
                    
                    # 2. Wait for the job to finish
                    # print ("Job is running. Polling {}.lck...".format(inpfil_name))
                    # while os.path.exists(lock_file):
                    #     time.sleep(10)
                    print ("--- Iteration {} complete. ---".format(inpfil_name))
                    odb_file = inpfil_name
                    enerHtotalN,sumaenerinterN,NeL2damageKN,NeL2damagetotalN,tf,energy_evol=PMTESCabaqusV0.PMTESCcritEneSubr(odb_file,control,cadenanombreSDV,work_dirkm)
                    # delete redundant files                                        
                    list_delete = ('*.sta', '*.dat', '*.sim', '*.msg', '*.com', '*.jnl',\
                                   '*.mtx','*.pes','*.pmg','*.ipm','*.pyc','*.fil')
                        
                print('*CRACK PROPAGATION STOPPED AT LOAD SUBSTEP NUMBER: {} *'.format(m))
                # delete redundant files                                        
                list_delete = ('*.sta', '*.dat', '*.sim', '*.prt', '*.msg', '*.com', '*.jnl',\
                               '*.mtx','*.pes','*.par','*.pmg','*.ipm','*.pyc','*.res','*.mdl','*.stt','*.fil')
                borrar_archivos(list_delete)
                chdir(actual_directory)
                for Dn in range(1,len(Uinterprops)+1):
                    nameinpDNk=nameinp+str(k-1)+'_N'+str(Dn)+'_K'+str(k)
                    job_name = '{}'.format(nameinpDNk)
                    work_dir = job_name + '_scratch'
                    shutil.rmtree(work_dir)
                end = time.time()
                print('Program has terminated in %.2f' %((end-start)/60)+' min')
                sys.exit()
                    # if m==1:
                    #     odbdef=Pmtesc_postp[indx][-1]
                    #     # shutil.copy(odbdef+'.odb',work_dirkm)
                    #     # shutil.copy(odbdef+'.res',work_dirkm)
                    #     # shutil.copy(odbdef+'.prt',work_dirkm)
                    #     # shutil.copy(common_file,work_dirkm)
                    # elif m>1:
                    #     nameinp_k_mf=nameinp+'_crackp'+'_m'+str(m-1)
                    #     job_namef = '{}'.format(nameinp_k_mf)
                    #     work_dirkmf = os.path.join(work_dir, os.path.basename(job_namef + '_scratch'))
                        
                    #     odbdef=work_dirkmf+'\\'+nameinp_k_mf
                    #     # shutil.copy(odbdef+'.odb',work_dirkm)
                    #     # shutil.copy(odbdef+'.res',work_dirkm)
                    #     # shutil.copy(odbdef+'.prt',work_dirkm)
                    #     # shutil.copy(common_file,work_dirkm)
                    # # load original ABAQUS input file
                    # mdb.ModelFromInputFile(inputFileName=nameinp+'.inp', name=nameinp_k_m)
                    # instancesname=mdb.models[nameinp_k_m].rootAssembly.instances.keys()
                    # instanceslista=[]
                    # for inst in range(len(instancesname)):
                    #     instanceslista.append(mdb.models[nameinp_k_m].rootAssembly.instances[instancesname[inst]])
                        
                    # mdb.models[nameinp_k_m].InitialState(updateReferenceConfiguration=OFF, fileName=odbdef, 
                    #           endStep=LAST_STEP, endIncrement=STEP_END, name='Predefined Field-1', 
                    #           createStepName='Initial', instances=tuple(instanceslista))
                    # # change step name:
                    # mdb.models[nameinp_k_m].steps.changeKey(fromName='Step-1', toName='Step-'+str(2)+'m'+str(m))
                    
                    # cadenasUinter=['*Surface Interaction, name=IntProp-1, user, depvar=22, properties=9\n','\n']
                    # kn = 300; kt = kn*1
                    # cadenasUinterN=[kn,kt,kt,0.09375,0.0,0.0,8.]
                    # # 11 starting damage configurations:
                    # cadenasUinterN1=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(1)+", "+str(control)+'\n']))]
                    # cadenasUinterN2=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(2)+", "+str(control)+'\n']))]
                    # cadenasUinterN3=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(3)+", "+str(control)+'\n']))]
                    # cadenasUinterN4=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(4)+", "+str(control)+'\n']))]
                    # cadenasUinterN5=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(5)+", "+str(control)+'\n']))]
                    # cadenasUinterN6=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(6)+", "+str(control)+'\n']))]
                    # cadenasUinterN7=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(7)+", "+str(control)+'\n']))]
                    # cadenasUinterN8=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(8)+", "+str(control)+'\n']))]
                    # cadenasUinterN9=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(9)+", "+str(control)+'\n']))]
                    # cadenasUinterN10=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(10)+", "+str(control)+'\n']))]
                    # cadenasUinterN11=[", ".join(map(str, [str(val) for val in cadenasUinterN]+[str(11)+", "+str(control)+'\n']))]
                    # Uinterprops=[cadenasUinterN1,cadenasUinterN11]    
                    # cadenaSDVcontact=['CDISP, CSTRESS, SDV\n']
                    # cadenanombreSDV='ASSEMBLY_TOP_SURF/ASSEMBLY_BOTTOM_SURF'
                    # nameodbcargainicial=nameinp+'_crackp'+'_m'+str(m)
                    # #Buscamos uinter='*Surface Interaction'
                    # cadenasCohesivo=lineasdefinteraccion(nameinp_k_m+'.inp', '*Surface Interaction',3)  #buscamos la cadena a cambiar
                    # #reemplazamos la primera linea en los dos archivos de inicio
                    # replacement(nameinp_k_m+'.inp', cadenasCohesivo[0], cadenasUinter[0])
                    # #reemplazamos la segunda linea en los dos archivos de inicio
                    # replacement(nameinp_k_m+'.inp', cadenasCohesivo[1], cadenasUinter[1])
                    # #reemplazamos la tercera linea en los dos archivos de inicio. Esta es distinta
                    # replacement(nameinp_k_m+'.inp', cadenasCohesivo[2], cadenasUinterN1[0][:])
                    # #reemplazamos la linea que pide las salida de las superficies de contacto    
                    # cadenaFOcontact=lineasdefinteraccion(nameinp_k_m+'.inp', 'CDISP, CSTRESS',1)  #buscamos la cadena a cambiar
                    # replacement(nameinp_k_m+'.inp', cadenaFOcontact[0], cadenaSDVcontact[0])
                    # # replacement(nameinp_k_m+'.inp', cadenaFOcontact[0], cadenaSDVcontact[0])
                    # mdb.models[nameinp_k_m].amplitudes[amplname].setValues(timeSpan=STEP, smooth=SOLVER_DEFAULT, data=((0.0, Fk), (1.0, Fk)))
                    
                    # # write the job for the crack propagation step: 
                    # mdb.Job(name=nameinp_k_m, model=nameinp_k_m, description='crack propagation analysis', 
                    #         type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, 
                    #         memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, memory=memoria,
                    #         explicitPrecision=SINGLE, nodalOutputPrecision=FULL, echoPrint=OFF, 
                    #         modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine=nameSUBR, 
                    #         scratch=work_dir, resultsFormat=ODB, multiprocessingMode=DEFAULT, numDomains=cpus, numCpus=cpus, 
                    #         numGPUs=0)
                    # # write the input .inp file:
                    # mdb.jobs[nameinp_k_m].writeInput(consistencyChecking=OFF)
                    
                    # # Add .fil output requests for the input file:
                    # add_urdfill_output(nameinp_k_m+'.inp')
                    # #Buscamos uinter='*Surface Interaction'
                    # cadenasCohesivo=lineasdefinteraccion(nameinp_k_m+'.inp', '*Surface Interaction',3)  #buscamos la cadena a cambiar
                    # #reemplazamos la primera linea en los dos archivos de inicio
                    # replacement(nameinp_k_m+'.inp', cadenasCohesivo[0], cadenasUinter[0])
                    # #reemplazamos la segunda linea en los dos archivos de inicio
                    # replacement(nameinp_k_m+'.inp', cadenasCohesivo[1], cadenasUinter[1])
                    # #reemplazamos la tercera linea en los dos archivos de inicio. Esta es distinta
                    # replacement(nameinp_k_m+'.inp', cadenasCohesivo[2], cadenasUinterN1[0][:])
                    # #reemplazamos la linea que pide las salida de las superficies de contacto    
                    # cadenaFOcontact=lineasdefinteraccion(nameinp_k_m+'.inp', 'CDISP, CSTRESS',1)  #buscamos la cadena a cambiar
                    # replacement(nameinp_k_m+'.inp', cadenaFOcontact[0], cadenaSDVcontact[0])
                    # # replacement(nameinp_k_m+'.inp', cadenaFOcontact[0], cadenaSDVcontact[0])
                    # mdb.models[nameinp_k_m].amplitudes[amplname].setValues(timeSpan=STEP, smooth=SOLVER_DEFAULT, data=((0.0, Fk), (1.0, Fk)))
                    
                    # # shutil.copy(nameSUBR, work_dirkm)
                    # myJobm = mdb.JobFromInputFile(name=nameinp_k_m, 
                    #     inputFileName=nameinp_k_m, type=ANALYSIS, 
                    #     atTime=None, waitMinutes=0, waitHours=0, queue=None,
                    #     memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, memory=memoria, 
                    #     explicitPrecision=SINGLE, nodalOutputPrecision=FULL, userSubroutine=nameSUBR, 
                    #     scratch=work_dir, resultsFormat=ODB, multiprocessingMode=DEFAULT, numDomains=cpus, numCpus=cpus, 
                    #     numGPUs=0)
                    # myJobm.submit(consistencyChecking=OFF)
                    # myJobm.waitForCompletion()
                    # # chdir(work_dirkm)
                    
                    # # Define the lock file path
                    # lock_file = os.path.join(work_dir, myJobm.name + '.lck')
            
        elif k>2 and NeL2damagek!=0 and stepdata[-1][-1]<0 and eps0>bisectol and k<=numiter:
            print(k, Fk, incre, eps0, NeL2damagek, energy_evol[-1][5])
            incre = (Fk-Fkm1)
            if NeL2damagekm1!=0 and stepdata[-2][-1]>0:
                newload=Fkm1+(0.5*abs(incre))
            elif NeL2damagekm1!=0 and stepdata[-2][-1]<0:
                matchdam = np.where((steps_data[:,5] != 0) & (steps_data[:,-1] > 0))[0]
                if matchdam.size > 0:
                    print("Last matching index:", km1_stepindx)
                    km1_indx = matchdam[-1]
                    Fkm1dam0 = steps_data[km1_indx,1]
                    incre = Fk-Fkm1dam0
                    newload=Fkm1dam0 + 0.5*abs(incre)
                else:
                    matchdam = np.where((steps_data[:,5] == 0))[0]
                    if matchdam.size>0:
                        km1_indx = matchdam[-1]
                        Fkm1dam0 = steps_data[km1_indx,1]
                        incre = Fk-Fkm1dam0
                        newload=Fkm1dam0 + 0.5*abs(incre)
                    else:
                        newload=0.5*Fk
                    
            elif NeL2damagekm1==0:
                newload=Fkm1+(0.5*abs(incre))
            print("Updated applied load: ",newload)
            
            job_name = '{}'.format(odbdef)
            work_dirodb=job_name + '_scratch'
            odb_file=odbdef
            #creamos los inp nuevos
            nameinpK=nameinp+str(k)
            shutil.copy(nameinp+str(k-1)+'.inp',nameinpK+'.inp')
            mdb.ModelFromInputFile(inputFileName=nameinpK+'.inp', name=nameinpK)
            mdb.models[nameinpK].amplitudes[amplname].setValues(timeSpan=STEP, smooth=SOLVER_DEFAULT, data=((0.0, newload), (1.0, newload)))
            # Create the job from the model
            mdb.Job(name=nameinpK, model=nameinpK)
            # write the .inp with newload
            mdb.jobs[nameinpK].writeInput(consistencyChecking=OFF)
            if godb==1:
                list_odbs=(['*_K'+str(k-1)+'.odb','*_K'+str(k-1)+'.inp','*_K'+str(k)+'.odb',
                            '*_K'+str(k)+'.sta','*_K'+str(k-1)+'.sta','*_K'+str(k)+'.inp'])
                mover_archivos(list_odbs,Dirname)
                #borramos el resto de archivos de ese paso                                          
                list_delete = ('*.sta', '*.dat', '*.sim', '*.prt', '*.msg', '*.com', '*.jnl',\
                               '*.mtx','*.pes','*.par','*.pmg','*.odb','*.ipm','*.pyc','*.res','*.mdl','*.stt')
                borrar_archivos(list_delete)
            for Dn in range(1,len(Uinterprops)+1):
                nameinpDNk=nameinp+str(k-1)+'_N'+str(Dn)+'_K'+str(k)
                job_name = '{}'.format(nameinpDNk)
                work_dir = job_name + '_scratch'
                shutil.rmtree(work_dir)
        elif k>2 and NeL2damagek!=0 and stepdata[-1][-1]>0 and k<=numiter:
            print(k, Fk, incre, eps0, NeL2damagek, stepdata[-1][-1])
            incre = (Fk-Fkm1)
            if NeL2damagekm1!=0 and stepdata[-2][-1]>0:
                # matchdam = np.where(steps_data[:,5] != 0 and steps_data[:,-1] < 0)[0]
                matchdam = np.where((steps_data[:,5] != 0) & (steps_data[:,-1] < 0))[0]
                km1_indx = matchdam[-1]
                Fkm1dam1 = steps_data[km1_indx,1]
                incre = Fkm1dam1 - Fk
                newload=Fkm1dam1 - 0.5*abs(incre)
            elif NeL2damagekm1!=0 and stepdata[-2][-1]<0:
                newload=Fk + 0.5*abs(incre)
            elif NeL2damagekm1==0:
                # matchdam1 = np.where(steps_data[:,5] != 0 and steps_data[:,-1]<0)[0]
                matchdam1 = np.where((steps_data[:,5] != 0) & (steps_data[:,-1] < 0))[0]
                if matchdam1.size > 0:
                    km1_indx = matchdam1[-1]
                    # print("Last matching index:", km1_indx)
                    Fkm1dam1 = steps_data[km1_indx,1]
                    incre = Fkm1dam1 - Fk
                    newload=Fk + abs(0.5*incre)
                else:
                    newload = 2*Fk
                    
            job_name = '{}'.format(odbdef)
            work_dirodb=job_name + '_scratch'
            odb_file=odbdef
            #creamos los inp nuevos
            nameinpK=nameinp+str(k)
            shutil.copy(nameinp+str(k-1)+'.inp',nameinpK+'.inp')
            mdb.ModelFromInputFile(inputFileName=nameinpK+'.inp', name=nameinpK)
            mdb.models[nameinpK].amplitudes[amplname].setValues(timeSpan=STEP, smooth=SOLVER_DEFAULT, data=((0.0, newload), (1.0, newload)))
            # Create the job from the model
            mdb.Job(name=nameinpK, model=nameinpK)
            # write the .inp with newload
            mdb.jobs[nameinpK].writeInput(consistencyChecking=OFF)
            if godb==1:
                list_odbs=(['*_K'+str(k-1)+'.odb','*_K'+str(k-1)+'.inp','*_K'+str(k)+'.odb',
                            '*_K'+str(k)+'.sta','*_K'+str(k-1)+'.sta','*_K'+str(k)+'.inp'])
                mover_archivos(list_odbs,Dirname)
                #borramos el resto de archivos de ese paso                                          
                list_delete = ('*.sta', '*.dat', '*.sim', '*.prt', '*.msg', '*.com', '*.jnl',\
                               '*.mtx','*.pes','*.par','*.pmg','*.odb','*.ipm','*.pyc','*.res','*.mdl','*.stt')
                borrar_archivos(list_delete)
            for Dn in range(1,len(Uinterprops)+1):
                nameinpDNk=nameinp+str(k-1)+'_N'+str(Dn)+'_K'+str(k)
                job_name = '{}'.format(nameinpDNk)
                work_dir = job_name + '_scratch'
                shutil.rmtree(work_dir)
        elif k>2 and NeL2damagek==0 and (k<=numiter):
            print(k, Fk, incre, eps0, NeL2damagek, NeL2damagekm1)
            incre = abs(Fk-Fkm1)
            if incre == 0:
                incre = 0.5*Fk
                print('Increment updated from delta=0')
            if NeL2damagekm1!=0 and stepdata[km1_stepindx][-1]>0:
                newload=2*Fkm1
            if NeL2damagekm1!=0 and stepdata[km1_stepindx][-1]<0:
                newload=Fkm1 - (0.5*incre)
            if NeL2damagekm1==0:
                matchdam1 = np.where((steps_data[:,5] != 0) & (steps_data[:,-1] < 0))[0]
                if matchdam1.size > 0:
                    km1_indx = matchdam1[-1]
                    # print("Last matching index:", km1_indx)
                    Fkm1dam1 = steps_data[km1_indx,1]
                    incre = Fkm1dam1 - Fk
                    newload=Fk + abs(0.5*incre)
                else:
                    newload = 2*Fk
            
            job_name = '{}'.format(odbdef)
            work_dirodb=job_name + '_scratch'
            odb_file=odbdef
            #creamos los inp nuevos
            nameinpK=nameinp+str(k)
            shutil.copy(nameinp+str(k-1)+'.inp',nameinpK+'.inp')
            mdb.ModelFromInputFile(inputFileName=nameinpK+'.inp', name=nameinpK)
            # Create the job from the model
            mdb.Job(name=nameinpK, model=nameinpK)
            mdb.models[nameinpK].amplitudes[amplname].setValues(timeSpan=STEP, smooth=SOLVER_DEFAULT, data=((0.0, newload), (1.0, newload)))
            
            # write the .inp with newload
            mdb.jobs[nameinpK].writeInput(consistencyChecking=OFF)
            if godb==1:
                list_odbs=(['*_K'+str(k-1)+'.odb','*_K'+str(k-1)+'.inp','*_K'+str(k)+'.odb',
                            '*_K'+str(k)+'.sta','*_K'+str(k-1)+'.sta','*_K'+str(k)+'.inp'])
                mover_archivos(list_odbs,Dirname)
                #borramos el resto de archivos de ese paso                                          
                list_delete = ('*.sta', '*.dat', '*.sim', '*.prt', '*.msg', '*.com', '*.jnl',\
                               '*.mtx','*.pes','*.par','*.pmg','*.odb','*.ipm','*.pyc','*.res','*.mdl','*.stt')
                borrar_archivos(list_delete)
            for Dn in range(1,len(Uinterprops)+1):
                nameinpDNk=nameinp+str(k-1)+'_N'+str(Dn)+'_K'+str(k)
                job_name = '{}'.format(nameinpDNk)
                work_dir = job_name + '_scratch'
                shutil.rmtree(work_dir)
        if k>numiter:
            if godb==1:
                list_odbs=(['*_K'+str(k-1)+'.odb','*_K'+str(k-1)+'.inp','*_K'+str(k)+'.odb',
                            '*_K'+str(k)+'.sta','*_K'+str(k-1)+'.sta','*_K'+str(k)+'.inp'])
                mover_archivos(list_odbs,Dirname)
                #borramos el resto de archivos de ese paso                                          
                list_delete = (['*_K'+str(k-1)+'.*','*'+str(k-1)+'.inp'])
                borrar_archivos(list_delete)
            for Dn in range(1,len(Uinterprops)+1):
                nameinpDNk=nameinp+str(k-1)+'_N'+str(Dn)+'_K'+str(k)
                job_name = '{}'.format(nameinpDNk)
                work_dir = job_name + '_scratch'
                shutil.rmtree(work_dir)
            end = time.time()
            print('Program aborted, number of steps exceeded the iteration limit')
            print('Program has terminated in %.2f' %((end-start)/60)+' min')
            sys.exit()
        
    list_delete = ('*.dat', '*.sim', '*.prt', '*.msg', '*.com', '*.jnl',\
                   '*.mtx','*.pes','*.par','*.pmg','*.odb','*.ipm','*.pyc','*.res','*.mdl','*.stt','*.fil')
    borrar_archivos(list_delete)
    k+=1    
# LOAD BISECTION ALGORITHM ENDS HERE

#


