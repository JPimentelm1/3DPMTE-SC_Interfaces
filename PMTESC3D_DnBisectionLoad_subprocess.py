# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 15:12:59 2023

@authors: Mar Munoz Reja Moreno;
          Jose M. Pimentel: load bisection algorithm, starting damage dist.
        implementation (uinter), Python's subprocess module for abaqus job cluster
        execution, jobs parallelization approach.
"""

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
# Redirect stdout to the command prompt
sys.stdout = sys.__stdout__
import PMTESCabaqus

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
    outfile.write(str('%d\t %.8e\t' %(energy_evol[-1][0],energy_evol[-1][-3]))+'\n')
    
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
            if cleaned_line == target_header.replace(" ", ""):
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
control=1
# initial load increment:
incr=1
factorcarga=1.0

# initial load:
cargainicial=1
# initial step:
k=1
# Bisection algorithm stop tolerance
bisectol=1.E-3
# Initial bisection tolerance
eps0=1
# Bisection iteration counter
bisect_itercnt=1
# Bisection limit of iterations
numiter=45
#
Nmin = round(-mth.log10(bisectol)/mth.log10(2))

nameinp='DCBuinter_L237m1'
UINTER_lst=['UINTERLEBIM3D_kincxit.for','UINTERLEBIMAMA3D.for']
# nameUMAT=UINTER_lst[0]
nameSUBR=UINTER_lst[0]
amplname='AMP-1'
godb=1  #si quiero guardar todos los odb pongo 1, si no 0
#definicion de mi uinter para cambiar
cadenasUinter=['*Surface Interaction, name=IntProp-1, user, depvar=22, properties=9\n','\n']
cadenasUinterN=[300.,600.,600.,0.09375,0.,0.,8.]

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
Dirname='AdaptiveF_ODBs-'+nameinp+'_knn'+str(cadenasUinterN[0])+'_mu'+str(cadenasUinterN[-1])+'_control'+str(control)

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
# simfile_obj=open(sim_file,'w')  
# simfile_obj.close()

nameodbcargainicial=nameinp+'-carga-inicial'
shutil.copy(nameinp+'.inp',nameodbcargainicial+'.inp')
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
# myJob.submit(consistencyChecking=OFF)
# myJob.waitForCompletion()
print("--- Now entering monitoring loop... ---")

# --- 3. MONITOR THE JOB (POLLING LOOP) ---
# The 'process.poll()' method checks if the subprocess has terminated.
# It returns 'None' if it's still running.
# while process.poll() is None:
#     "Waiting for '{}' to finish... (checking again in N seconds)".format(jobName)
#     # 'time.sleep()' makes the script pause, preventing it from
#     # using 100% CPU while it waits.
#     time.sleep(10)
# The return code (0 for success) is captured.
return_code = process.wait()    
"\n--- Abaqus job has finished ---"
"--- The Abaqus process returned exit code: {} ---".format(return_code)

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
newload = PMTESCabaqus.PMTESCcriTen_carga(nameodbcargainicial, cargainicial, cadenanombreSDV)
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
        Joblist.append(myJobN)
        
        myJobN.setValues(
        scratch=work_dir,
        userSubroutine=job_sub_path,  # Point to local subroutine,
        numCpus=1,          # MUST BE 1 FOR PARALLEL EXECUTION
        numDomains=1,        # Disable domain decomposition
        multiprocessingMode=DEFAULT)
        
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
    # last_job_in_list = Joblist[-1]
    # last_job_name = last_job_in_list.name
    # # NOTE: This path construction MUST match the one used in your submission loop.
    # last_job_work_dir = os.path.join(os.getcwd(), last_job_name + '_scratch')
    # expected_odb_path = os.path.join(last_job_work_dir, last_job_name + '.odb')
    # print('The script will proceed after this file is created: {}'.format(expected_odb_path))
    
    # 2. WAITING LOOP: This part remains the same
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
        enerHtotalN,sumaenerinterN,NeL2damageKN,NeL2damagetotalN,tf,energy_evol=PMTESCabaqus.PMTESCcritEneSubr(odb_file,control,cadenanombreSDV,work_dirodb)
        model_outputfile(actual_directory,sim_file,k,Dn,newload,tf,enerHtotalN,sumaenerinterN,len(NeL2damageKN),len(NeL2damagetotalN),energy_evol)
        Pmtesc_postp.append([enerHtotalN,sumaenerinterN,NeL2damageKN,NeL2damagetotalN,tf,nameinpDNk])
        datos_salida.append([k,Dn,float(newload),float(tf),float(enerHtotalN),sumaenerinterN,len(NeL2damageKN),nameinpDNk,len(NeL2damagetotalN)])
        print([k,Dn,tf,enerHtotalN,sumaenerinterN,len(NeL2damageKN),len(NeL2damagetotalN),nameinpDNk,energy_evol[-1][-1]]) 
        # move the job directory to the parent script directory:
        # shutil.move(work_dir, Dirname)
        
    #%% Algorithm flow control for the ODB file with minimum total energi PI+DeltaR:
    MinPIpDeltaR=min(row[0] for row in Pmtesc_postp); print(MinPIpDeltaR)
    min_row = min(Pmtesc_postp, key=lambda x:x[0])
    for row in range(len(Pmtesc_postp)):
        if len(Pmtesc_postp[row][2])>0.0 and Pmtesc_postp[row][0]<=min_row[0]:
                odbdef=Pmtesc_postp[row][-1]
                NeL2damage=Pmtesc_postp[row][2]
                NeL2damagetotal=Pmtesc_postp[row][3]
                job_name = '{}'.format(odbdef)
        
        else:
            if len(Pmtesc_postp[row][2])==0.0:
                # 
                odbdef=Pmtesc_postp[row][-1]
                NeL2damage=Pmtesc_postp[row][2]
                NeL2damagetotal=Pmtesc_postp[row][3]							 
                job_name = '{}'.format(odbdef)
    print(odbdef, len(NeL2damage))
    work_dirodb=Dirname+'//'+job_name + '_scratch'
    odb_file=job_name +'_scratch'+'//'+odbdef
    enerHtotalN,sumaenerinterN,NeL2damageKN,NeL2damagetotalN,tf,energy_evol=PMTESCabaqus.PMTESCcritEneSubr(odb_file,control,cadenanombreSDV,work_dirodb)
    model_outputfile(actual_directory,sim_file,k,Dn,newload,tf,enerHtotalN,sumaenerinterN,len(NeL2damageKN),len(NeL2damagetotalN),energy_evol)
    datos_salida.append([float(newload),float(tf),float(enerHtotalN),sumaenerinterN,len(NeL2damageKN),odbdef,len(NeL2damagetotalN)])
    
    stepdata.append([k,newload,tf,enerHtotalN,sumaenerinterN,len(NeL2damageKN),len(NeL2damagetotalN),energy_evol[-1][-3]])
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
            
        elif k>2 and NeL2damagek!=0 and stepdata[-1][-1]<0 and eps0<bisectol and k<=numiter:
            print(k, Fk, tf, eps0, NeL2damagek)
            if godb==1:
                list_odbs=(['*_K'+str(k-1)+'.odb','*_K'+str(k-1)+'.inp','*_K'+str(k)+'.inp',
                            '*_K'+str(k)+'.sta','*_K'+str(k-1)+'.sta'])
                mover_archivos(list_odbs,Dirname)
            #borramos el resto de archivos de ese paso                                          
            list_delete = ('*.sta', '*.dat', '*.sim', '*.prt', '*.msg', '*.com', '*.jnl',\
                           '*.mtx','*.pes','*.par','*.pmg','*.odb','*.ipm','*.pyc','*.res','*.mdl','*.stt')
            borrar_archivos(list_delete)                                    
            chdir(actual_directory)
            end = time.time()
            print('Simulation has terminated in %.2f' %((end-start)/60)+' min')
            sys.exit()
        elif k>2 and NeL2damagek!=0 and stepdata[-1][-1]<0 and eps0>bisectol and k<=numiter:
            print(k, Fk, incre, eps0, NeL2damagek, energy_evol[-1][-3])
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

        if k>numiter:
            if godb==1:
                list_odbs=(['*_K'+str(k-1)+'.odb','*_K'+str(k-1)+'.inp','*_K'+str(k)+'.odb',
                            '*_K'+str(k)+'.sta','*_K'+str(k-1)+'.sta','*_K'+str(k)+'.inp'])
                mover_archivos(list_odbs,Dirname)
                #borramos el resto de archivos de ese paso                                          
                list_delete = (['*_K'+str(k-1)+'.*','*'+str(k-1)+'.inp'])
                borrar_archivos(list_delete)
            end = time.time()
            print('Program aborted, number of steps exceeded the iteration limit')
            print('Program has terminated in %.2f' %((end-start)/60)+' min')
            sys.exit()
        
    list_delete = ('*.dat', '*.sim', '*.prt', '*.msg', '*.com', '*.jnl',\
                   '*.mtx','*.pes','*.par','*.pmg','*.odb','*.ipm','*.pyc','*.res','*.mdl','*.stt','*.fil')
    borrar_archivos(list_delete)
    k+=1    
    # LOAD BISECTION ALGORITHM ENDS HERE
    
#%%preparamos el siguiente inp desde el odb anterior
#importamos el inp inicial sin CI impuestas para crear el nuevo modelo del paso siguiente 
nameinp_k_mas_1=nameinp+str(k+1)
copy(nameinp+'.inp',nameinp_k_mas_1+'.inp')
#cargamos el input original de ABAQUS
mdb.ModelFromInputFile(inputFileName=nameinp+'.inp', name=nameinp_k_mas_1)
#imponemos el estado inicial del job ganador entre las n para que sea el inicio del paso anterior
instancesname=mdb.models[nameinp_k_mas_1].rootAssembly.instances.keys()
instanceslista=[]
for inst in range(len(instancesname)):
    instanceslista.append(mdb.models[nameinp_k_mas_1].rootAssembly.instances[instancesname[inst]])
    
mdb.models[nameinp_k_mas_1].InitialState(updateReferenceConfiguration=OFF, fileName=odbdef, 
          endStep=LAST_STEP, endIncrement=STEP_END, name='Predefined Field-1', 
          createStepName='Initial', instances=tuple(instanceslista))

# cambiamos el nombre del step  
mdb.models[nameinp_k_mas_1].steps.changeKey(fromName='Step-1', toName='Step-'+str(k+1))
# cambiamos la carga
#%%aumentamos la carga si no hay mas rotura
if(collections.Counter(NeL2damagetotal)==collections.Counter(NeL2damageKm1)):
    newload=incr+newload
mdb.models[nameinp_k_mas_1].amplitudes[amplname].setValues(timeSpan=STEP, smooth=SOLVER_DEFAULT, data=((0.0, newload), (1.0, newload)))
#%%guardamos el dagno anterior
NeL2damageKm1=NeL2damagetotal

#escribimos el inp final de k mas 1 para empezar el siguiente paso  
mdb.Job(name=nameinp_k_mas_1, model=nameinp_k_mas_1, description='', 
        type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, memory=memoria,
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine=nameSUBR, 
        scratch=actual_directory, resultsFormat=ODB, multiprocessingMode=DEFAULT, numDomains=cpus, numCpus=cpus, 
        numGPUs=0)

mdb.jobs[nameinp_k_mas_1].writeInput(consistencyChecking=OFF)
#dejamos el inp basico por donde empezara el bucle con los diferentes inicios
nameinpK=nameinp_k_mas_1
print('Updated applied load: ',newload)
#%% agnadimos el nuevo odb al nuevo odb completo para la salida
if k==1:
    nameodbcompleto=nameinp+'_completo.odb'
    copy(odbdef+'.odb',nameodbcompleto)
else:
    call('abaqus restartjoin originalodb='+nameodbcompleto+' restartodb='+odbdef+'.odb'+' history',shell=True)
### borramos archivos menos el odb que sera opcional. No se puede hacer antes 
    #guardamos lo odb si los queremos
    if godb==1:
        list_odbs=(['*_K'+str(k-1)+'.odb'])
        mover_archivos(list_odbs,Dirname)
    #borramos el resto de archivos de ese paso                                          
    list_delete = (['*_K'+str(k-1)+'.*','*'+str(k-1)+'.inp'])
    #borrar_archivos(list_delete)
#%%final  
filedicc=open(actual_directory+'/datospasoscompletoUINTER.pkl','wb')    
pickle.dump(datos_salida,filedicc)
filedicc.close()
if godb==1:
    list_odbs=(['*_K'+str(k)+'.odb'])
    mover_archivos(list_odbs,Dirname)

move('datospasoscompletoUINTER.pkl',Dirname)
move(nameodbcompleto,Dirname)
list_delete = (['*_K'+str(k)+'.*','*'+str(k)+'.inp','*'+str(k+1)+'.inp'])
borrar_archivos(list_delete)    
list_delete = ('*.sta', '*.dat', '*.sim', '*.prt', '*.msg', '*.com', '*.jnl',\
               '*.mtx','*.pes','*.par','*.pmg','*.odb','*.ipm','*.pyc','*.res','*.mdl','*.stt')
borrar_archivos(list_delete)
end = time.time()
print('La simulacion ha terminado satisfactoriamente en %.2f' %((end-start)/60)+' min')

##%%
##PARA ABRIR UN ARCHIVO PKL:
import pickle
with open(Dirname+'/datospasoscompletoUINTER.pkl','rb') as f:
    datospasoscompletoUINTER = pickle.load(f)
# Convert PKL 2 numpy array:
stepdata_np=np.array(datospasoscompletoUINTER)
np.savetxt(Dirname+'/table_stepdata.txt', stepdata_np, delimiter='\t')