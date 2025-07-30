# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 15:12:59 2023

@authors: Mar Munoz Reja Moreno
            Jose M. Pimentel
"""

from os import chdir, path, mkdir, getcwd, remove
import os
from shutil import rmtree, move,copy
from glob import glob
from subprocess import call
import fileinput
import sys
import time
import collections
import pickle
import numpy as np
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
            
def model_outputfile(sim_file,K,Dn,newload,TF,enerHtotalN1,sumaenerinterN1,lNeL2damageKN1,lNeL2damagetotal):
    outfile=open(sim_file,'a')
    outfile.write(str('%d\t %d\t %.8e\t %.8e\t %.8e\t %.8e\t %d\t'%(K,Dn,newload,TF,enerHtotalN1,sumaenerinterN1,lNeL2damageKN1)))
    outfile.write(str('%d\t'%(lNeL2damagetotal))+'\n')
    
    outfile.close()
    
#%% datos de entrada
call('cls',shell=True)
actual_directory = getcwd()
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
cargainicial=1
# initial step:
k=1
# Bisection algorithm stop tolerance
bisectol=250.0E-4
# Initial bisection tolerance
eps0=1
# Bisection iteration counter
bisect_itercnt=1
# Bisection limit of iterations
numiter=39

nameinp='DCBuinter_m2'
nameSUBR='UINTERLEBIM3D_Mod.f'
amplname='AMP-1'
godb=1  #si quiero guardar todos los odb pongo 1, si no 0
#definicion de mi uinter para cambiar
cadenasUinter=['*Surface Interaction, name=IntProp-1, user, depvar=22, properties=8\n','\n']

cadenasUinterN=[150.,  533.3,  533.3, .75, .25, 8]
cadenasUinterN1='150.,  533.3,  533.3,  0.75,   0.25,  0.25, 8., 1.\n'    
cadenasUinterN2='150.,  533.3,  533.3,  0.75,   0.25,  0.25, 8., 2.\n'
cadenasUinterN3='150.,  533.3,  533.3,  0.75,   0.25,  0.25, 8., 3.\n'
cadenasUinterN4='150.,  533.3,  533.3,  0.75,   0.25,  0.25, 8., 4.\n'
cadenasUinterN5='150.,  533.3,  533.3,  0.75,   0.25,  0.25, 8., 5.\n'
cadenasUinterN6='150.,  533.3,  533.3,  0.75,   0.25,  0.25, 8., 6.\n'
cadenasUinterN7='150.,  533.3,  533.3,  0.75,   0.25,  0.25, 8., 7.\n' 
cadenasUinterN8='150.,  533.3,  533.3,  0.75,   0.25,  0.25, 8., 8.\n' 
cadenasUinterN9='150.,  533.3,  533.3,  0.75,   0.25,  0.25, 8., 9.\n'
cadenasUinterN10='150.,  533.3,  533.3,  0.75,   0.25,  0.25, 8., 10.\n'
cadenasUinterN11='150.,  533.3,  533.3,  0.75,   0.25,  0.25, 8., 11.\n'
Uinterprops=[cadenasUinterN1,cadenasUinterN2,cadenasUinterN3,cadenasUinterN4,
             cadenasUinterN5,cadenasUinterN6,cadenasUinterN7,cadenasUinterN8,
             cadenasUinterN9,cadenasUinterN10,cadenasUinterN11]    
cadenaSDVcontact=['CDISP, CSTRESS, SDV\n']
# cadenanombreSDV='ASSEMBLY_VIGA_INF_SUPER_SUP/ASSEMBLY_VIGA_SUP_SUPER_INF'
cadenanombreSDV='ASSEMBLY_TOP_SURF/ASSEMBLY_BOTTOM_SURF'
Dirname='AdaptiveF_ODBs-'+nameinp+'_knn'+str(cadenasUinterN[0])+'_mu'+str(cadenasUinterN[-1])

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
copy(nameinp+'.inp',nameodbcargainicial+'.inp')
#Buscamos uinter='*Surface Interaction'
cadenasCohesivo=lineasdefinteraccion(nameodbcargainicial+'.inp', '*Surface Interaction',3)  #buscamos la cadena a cambiar
#reemplazamos la primera linea en los dos archivos de inicio
replacement(nameodbcargainicial+'.inp', cadenasCohesivo[0], cadenasUinter[0])
#reemplazamos la segunda linea en los dos archivos de inicio
replacement(nameodbcargainicial+'.inp', cadenasCohesivo[1], cadenasUinter[1])
#reemplazamos la tercera linea en los dos archivos de inicio. Esta es distinta
replacement(nameodbcargainicial+'.inp', cadenasCohesivo[2], cadenasUinterN1[:])
#reemplazamos la linea que pide las salida de las superficies de contacto    
cadenaFOcontact=lineasdefinteraccion(nameodbcargainicial+'.inp', 'CDISP, CSTRESS',1)  #buscamos la cadena a cambiar
#cadenaSDVcontact=['CDISP, CSTRESS, SDV\n']
replacement(nameodbcargainicial+'.inp', cadenaFOcontact[0], cadenaSDVcontact[0])
replacement(nameodbcargainicial+'.inp', cadenaFOcontact[0], cadenaSDVcontact[0])
myJob = mdb.JobFromInputFile(name=nameodbcargainicial, 
        inputFileName=nameodbcargainicial, type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None,
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, memory=memoria, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, userSubroutine=nameSUBR, 
        scratch=actual_directory, parallelizationMethodExplicit=DOMAIN, numDomains=cpus, 
        activateLoadBalancing=False, multiprocessingMode=DEFAULT, numCpus=cpus)

myJob.submit(consistencyChecking=OFF)
myJob.waitForCompletion()
###obtener la nueva carga del odb
# print(nameodbcargainicial, cargainicial, cadenanombreSDV)
newload = PMTESCabaqus.PMTESCcriTen_carga(nameodbcargainicial, cargainicial, cadenanombreSDV)
# 
newload=newload*factorcarga 
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
copy(nameinp+'.inp',nameinpK+'.inp')

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

##

###inicializamos variables:
NeL2damageKm1=[]
datos_salida=[]
stepdata=[]
#%% LOAD BISECTION ALGORITHM BEGINS HERE
# bucle k. Ahora la k y la m estan acopladas. Simplemente cuando hay dagno no aumento carga y cuando termina de dagnar aumento
while (k<=numiter):
    # 
    #%% 
    #%
    # Base configuration
    base_dir = 'parallel_jobs_k{}'.format(k)
    subroutine_file = nameSUBR  # Update this path
    input_file = nameinpK+'.inp'  # Update this path
    # Create base directory
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)
    #SE PREPARAN DESDE EL INP ORIGEN, PARA REVISAR POSTERIORMENTE LOS INP SI QUEREMOS
    Joblist=[]
    Pmtesc_postp=[]
    for Dn in range(1,len(Uinterprops)+1):
        #creamos los inp nuevos
        nameinpDNk=nameinpK+'_N'+str(Dn)+'_K'+str(k)
        # nameinpN2k=nameinpK+'_N2_K'+str(k)
        copy(nameinpK+'.inp',nameinpDNk+'.inp')  
        # copy(nameinpK+'.inp',nameinpN2k+'.inp')  
        # print(Uinterprops[:][Dn-1])
        
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
        replacement(nameinpDNk+'.inp', cadenasCohesivo[2], Uinterprops[:][Dn-1])
        # replacement(nameinpN2k+'.inp', cadenasCohesivo[2], cadenasUinterN2[0])
        #reemplazamos la linea que pide las salida de las superficies de contacto    
        cadenaFOcontact=lineasdefinteraccion(nameinpDNk+'.inp', 'CDISP, CSTRESS',1)  #buscamos la cadena a cambiar
        cadenaSDVcontact=['CDISP, CSTRESS, SDV\n']
        replacement(nameinpDNk+'.inp', cadenaFOcontact[0], cadenaSDVcontact[0])
        # replacement(nameinpN2k+'.inp', cadenaFOcontact[0], cadenaSDVcontact[0])

        #%
        myJobN = mdb.JobFromInputFile(name=nameinpDNk, 
            inputFileName=nameinpDNk, type=ANALYSIS, 
            atTime=None, waitMinutes=0, waitHours=0, queue=None,
            memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, memory=memoria, 
            explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, userSubroutine=nameSUBR, 
            scratch=actual_directory, parallelizationMethodExplicit=DOMAIN, numDomains=cpus, 
            activateLoadBalancing=False, multiprocessingMode=DEFAULT, numCpus=cpus)
        Joblist.append(myJobN)
        # print(Joblist)  
        
    # ***ABAQUS JOB PARALLELIZATION APPROACH:***
    for job in range(len(Joblist)):
        nameinpDNk=nameinpK+'_N'+str(job+1)+'_K'+str(k)
        # Create unique work directory for each job
        job_name = '{}'.format(nameinpDNk)
        work_dir = os.path.join(base_dir, job_name)
        
        if not os.path.exists(work_dir):
            os.makedirs(work_dir)
        # Copy subroutine to job directory (CRITICAL STEP)
        job_sub_path = os.path.join(work_dir, os.path.basename(subroutine_file))
        shutil.copyfile(subroutine_file, job_sub_path)
        
        Joblist[job].setValues(
        scratch=work_dir,
        userSubroutine=job_sub_path,  # Point to local subroutine,
        numCpus=1,          # MUST BE 1 FOR PARALLEL EXECUTION
        numDomains=1,        # Disable domain decomposition
        multiprocessingMode=DEFAULT)
        Joblist[job].submit(consistencyChecking=OFF)
    # Create thread for job monitoring
    def monitor_job(job_obj):
        job_obj.waitForCompletion()
        
    # Create and start monitoring threads
    threads = []
    for job in Joblist:
        t = threading.Thread(target=monitor_job, args=(job,))
        t.daemon = True  # Ensures threads exit if main program exits
        t.start()
        threads.append(t)
        
        # print(Pmtesc_postp)
        
    # Wait for all threads to complete
    while any(t.is_alive() for t in threads):
        time.sleep(5)  # Check every 5 seconds
        # print("Jobs still running: %d", sum(t.is_alive() for t in threads))
    # Check job statuses
    failed_jobs = []
    for i, job in enumerate(Joblist):
        status = job.status
        if status in [ABORTED, TERMINATED]:
            failed_jobs.append((i, job.name, status))
        # print("Job %s (%d) status: %s", (job.name, i, status))
    
    if failed_jobs:
        msg = "Failed jobs detected:\n"
        for idx, name, status in failed_jobs:
            msg += "  - Job %d (%s): %s\n" , (idx, name, status)
        getWarningReply(msg, buttons=(YES, NO)) 
    else:
        print("All jobs completed successfully!")
    # 
    os.chdir(os.getcwd())
    for Dn in range(1,len(Uinterprops)+1):
        nameinpDNk=nameinpK+'_N'+str(Dn)+'_K'+str(k)
        # ODB Post-processing; simulation data output:
        enerHtotalN,sumaenerinterN,NeL2damageKN,NeL2damagetotalN,tf=PMTESCabaqus.PMTESCcritEneUMAT(nameinpDNk,control,cadenanombreSDV)
        model_outputfile(sim_file,k,Dn,newload,tf,enerHtotalN,sumaenerinterN,len(NeL2damageKN),len(NeL2damagetotalN))
        Pmtesc_postp.append([enerHtotalN,sumaenerinterN,NeL2damageKN,NeL2damagetotalN,tf,nameinpDNk])
        datos_salida.append([k,Dn,float(newload),float(tf),float(enerHtotalN),sumaenerinterN,len(NeL2damageKN),nameinpDNk,len(NeL2damagetotalN)])
        print([k,Dn,tf,enerHtotalN,sumaenerinterN,len(NeL2damageKN),len(NeL2damagetotalN),nameinpDNk])        
    #%% Algorithm flow control for the ODB file with minimum total energi PI+DeltaR:
    MinPIpDeltaR=min(row[0] for row in Pmtesc_postp); 
    print(MinPIpDeltaR)
    min_row = min(Pmtesc_postp, key=lambda x:x[0])
    for row in range(len(Pmtesc_postp)):
        if len(Pmtesc_postp[row][2])>0.0 and Pmtesc_postp[row][0]<=min_row[0]:
                odbdef=Pmtesc_postp[row][-1]
                NeL2damage=Pmtesc_postp[row][2]
                NeL2damagetotal=Pmtesc_postp[row][3]
        
        else:
            if len(Pmtesc_postp[row][2])==0.0:
                # 
                odbdef=Pmtesc_postp[row][-1]
                NeL2damage=Pmtesc_postp[row][2]
                NeL2damagetotal=Pmtesc_postp[row][3]							 
    
    enerHtotalN,sumaenerinterN,NeL2damageKN,NeL2damagetotalN,tf=PMTESCabaqus.PMTESCcritEneUMAT(odbdef,control,cadenanombreSDV)
    model_outputfile(sim_file,k,Dn,newload,tf,enerHtotalN,sumaenerinterN,len(NeL2damageKN),len(NeL2damagetotalN))
    datos_salida.append([float(newload),float(tf),float(enerHtotalN),sumaenerinterN,len(NeL2damageKN),odbdef,len(NeL2damagetotalN)])
    
    stepdata.append([k,newload,tf,enerHtotalN,sumaenerinterN,len(NeL2damageKN),len(NeL2damagetotalN)])
    sim_filenp=np.array(stepdata)
    np.save('sim_file.npy', sim_filenp)
    
    if k==1:
        steps_data=np.load('sim_file.npy')
        print(steps_data)
        # kstepindx=steps_data[-1,0]==k
        # Fk = (steps_data[kstepindx,1])
        eps0=1.0; 
        newload=newload+incr
        # print(k,newload)
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
        
        k+=1
    else:
        steps_data=np.load('sim_file.npy')
        print(steps_data)
        kstepindx=steps_data[-1,0]==k
        km1stepindx=steps_data[-2,0]==k-1
        Fk = newload
        Fkm1 = steps_data[-2,1]
        incre = Fk-Fkm1
        eps0 = abs(incre)/abs(Fk)
        NeL2damagek=steps_data[-1,5]
        NeL2damagekm1=steps_data[k-2,5]
        print(k, Fk, incre, eps0, NeL2damagek, NeL2damagekm1)
        if incre==0:
            k+=1
            incre = 0.25*Fk
            eps0 = abs(incre)/abs(Fk)
            #creamos los inp nuevos
            nameinpK=nameinp+str(k)
            copy(nameinp+str(k-1)+'.inp',nameinpK+'.inp') 
            mdb.ModelFromInputFile(inputFileName=nameinpK+'.inp', name=nameinpK)
            mdb.models[nameinpK].amplitudes[amplname].setValues(timeSpan=STEP, smooth=SOLVER_DEFAULT, data=((0.0, newload), (1.0, newload)))
            # Create the job from the model
            mdb.Job(name=nameinpK, model=nameinpK)
            # write the .inp with newload
            mdb.jobs[nameinpK].writeInput(consistencyChecking=OFF)
            
        if NeL2damagek!=0 and eps0<bisectol and (k<=numiter):
            print(k, Fk, tf, eps0, NeL2damagek)
            if godb==1:
                list_odbs=(['*_K'+str(k)+'.odb','*_K'+str(k-1)+'.odb','*_K'+str(k-1)+'.inp',
                            '*_K'+str(k)+'.inp','*_K'+str(k)+'.sta','*_K'+str(k-1)+'.sta'])
                mover_archivos(list_odbs,Dirname)
            #borramos el resto de archivos de ese paso                                          
            list_delete = ('*.sta', '*.dat', '*.sim', '*.prt', '*.msg', '*.com', '*.jnl',\
                           '*.mtx','*.pes','*.par','*.pmg','*.odb','*.ipm','*.pyc','*.res','*.mdl','*.stt')
            borrar_archivos(list_delete)
            end = time.time()
            print('Simulation has terminated in %.2f' %((end-start)/60)+' min')
            sys.exit() 
        if NeL2damagek==0 and eps0<bisectol and (k<=numiter):
            incre = (Fk-Fkm1)
            print(k, Fk, incre, eps0, NeL2damagek, NeL2damagekm1)
            if NeL2damagekm1!=0:
                newload=newload+abs(incre*0.5)
            elif NeL2damagekm1==0:
                newload=newload+(1.25*abs(incre))
            
            #creamos los inp nuevos
            nameinpK=nameinp+str(k)
            copy(nameinp+str(k-1)+'.inp',nameinpK+'.inp') 
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
            k+=1
        if NeL2damagek!=0 and eps0>bisectol and (k<=numiter):
            print(k, Fk, incre, eps0, NeL2damagek, NeL2damagekm1)
            incre = (Fk-Fkm1)
            if NeL2damagekm1!=0:
                newload=abs(newload-(1.25*abs(incre)))
            elif NeL2damagekm1==0:
                newload=abs(newload-abs(incre*0.5))
            
            #creamos los inp nuevos
            nameinpK=nameinp+str(k)
            copy(nameinp+str(k-1)+'.inp',nameinpK+'.inp') 
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
            k+=1
        if NeL2damagek==0 and eps0>bisectol and (k<=numiter):
            print(k, Fk, incre, eps0, NeL2damagek, NeL2damagekm1)
            incre = (Fk-Fkm1)
            if NeL2damagekm1!=0:
                newload=abs(newload+abs(incre*0.5))
            elif NeL2damagekm1==0:
                newload=abs(newload+(1.25*abs(incre)))
            
            #creamos los inp nuevos
            nameinpK=nameinp+str(k)
            copy(nameinp+str(k-1)+'.inp',nameinpK+'.inp') 
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
            k+=1
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
        
        
    if godb==1:
        list_odbs=(['*_K'+str(k-1)+'.odb','*_K'+str(k-1)+'.inp','*_K'+str(k)+'.odb',
                    '*_K'+str(k)+'.sta','*_K'+str(k-1)+'.sta','*_K'+str(k)+'.inp'])
        mover_archivos(list_odbs,Dirname)
        #borramos el resto de archivos de ese paso                                                                               
        list_delete = (['*_K'+str(k-1)+'.*','*'+str(k-1)+'.inp'])
        # borrar_archivos(list_delete)
        list_delete = ('*.sta', '*.dat', '*.sim', '*.prt', '*.msg', '*.com', '*.jnl',\
                       '*.mtx','*.pes','*.par','*.pmg','*.odb','*.ipm','*.pyc','*.res','*.mdl','*.stt')
        borrar_archivos(list_delete)
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
print(stepdata_np)