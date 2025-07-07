# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 15:12:59 2023

@author: Mar Munoz Reja Moreno
"""

from os import chdir, path, mkdir, getcwd, remove
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

def model_outputfile(sim_file,K,newload,TF,enerHtotalN1,sumaenerinterN1,lNeL2damageKN1,enerHtotalN2,sumaenerinterN2,lNeL2damageKN2,lNeL2damagetotal):
    outfile=open(sim_file,'a')
    outfile.write(str('%d\t %.8e\t %.8e\t %.8e\t %.8e\t %d\t'%(K,newload,TF,enerHtotalN1,sumaenerinterN1,lNeL2damageKN1)))
    outfile.write(str('%.8e\t %.8e\t %d\t %d\t'%(enerHtotalN2,sumaenerinterN2,lNeL2damageKN2,lNeL2damagetotal))+'\n')
    
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
incr=0
factorcarga=1.0
# initial load:
cargainicial=1
pasosK=1
nameinp='DCBuinter_m2'
UINTER_lst=['UINTERLEBIM3D_Mod.for','UINTERLEBIM3D_BKMod.for']
nameUMAT=UINTER_lst[0]
amplname='AMP-1'
godb=1  #si quiero guardar todos los odb pongo 1, si no 0
#definicion de mi uinter para cambiar
cadenasUinter=['*Surface Interaction, name=IntProp-1, user, depvar=20, properties=8\n','\n']
# cadenasUinterN1=['150.,  600.,  600.,  0.75,   0.25,  0.25, 8., 1.\n']  
cadenasUinterN=[150.,  600.,  600., .75, .25, 8]
cadenasUinterN1=['150.,  600.,  600.,  0.75,   0.25,  0.25, 8., 1.\n']    
cadenasUinterN2=['150.,  600.,  600.,  0.75,   0.25,  0.25, 8., 2.\n']     
cadenaSDVcontact=['CDISP, CSTRESS, SDV\n']
# cadenanombreSDV='ASSEMBLY_VIGA_INF_SUPER_SUP/ASSEMBLY_VIGA_SUP_SUPER_INF'
cadenanombreSDV='ASSEMBLY_TOP_SURF/ASSEMBLY_BOTTOM_SURF'
Dirname='ODBs-'+nameinp+'_knn'+str(cadenasUinterN[0])+'_mu'+str(cadenasUinterN[-1])

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
replacement(nameodbcargainicial+'.inp', cadenasCohesivo[2], cadenasUinterN1[0])
#reemplazamos la linea que pide las salida de las superficies de contacto    
cadenaFOcontact=lineasdefinteraccion(nameodbcargainicial+'.inp', 'CDISP, CSTRESS',1)  #buscamos la cadena a cambiar
#cadenaSDVcontact=['CDISP, CSTRESS, SDV\n']
replacement(nameodbcargainicial+'.inp', cadenaFOcontact[0], cadenaSDVcontact[0])
replacement(nameodbcargainicial+'.inp', cadenaFOcontact[0], cadenaSDVcontact[0])
myJob = mdb.JobFromInputFile(name=nameodbcargainicial, 
        inputFileName=nameodbcargainicial, type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None,
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, memory=memoria, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, userSubroutine=nameUMAT, 
        scratch=actual_directory, parallelizationMethodExplicit=DOMAIN, numDomains=cpus, 
        activateLoadBalancing=False, multiprocessingMode=DEFAULT, numCpus=cpus)

myJob.submit(consistencyChecking=OFF)
myJob.waitForCompletion()
###obtener la nueva carga del odb
# print(nameodbcargainicial, cargainicial, cadenanombreSDV)
newload = PMTESCabaqus.PMTESCcriTen_carga(nameodbcargainicial, cargainicial, cadenanombreSDV)
# How is the newload variable calculated?
newload=newload*factorcarga 
###borrar archivos
if godb==1:
    list_odbs=([nameodbcargainicial+'*.odb'])
    mover_archivos(list_odbs,Dirname)
#borramos el resto de archivos de ese paso                                          
list_delete = ([nameodbcargainicial+'.*'])
#borrar_archivos(list_delete)
###preparamos nuestro inp k=1
nameinpK=nameinp+str(1)
copy(nameinp+'.inp',nameinpK+'.inp')

###importo inp y impongo nueva carga
print(newload)
mdb.ModelFromInputFile(inputFileName=nameinpK+'.inp', name=nameinpK)
mdb.models[nameinpK].amplitudes[amplname].setValues(timeSpan=STEP, smooth=SOLVER_DEFAULT, data=((0.0, newload), (1.0, newload)))
mdb.Job(name=nameinpK, model=nameinpK, description='', 
            type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, 
            memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, memory=memoria,
            explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
            modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine=nameUMAT, 
            scratch=actual_directory, resultsFormat=ODB, multiprocessingMode=DEFAULT, numDomains=cpus, numCpus=cpus, 
            numGPUs=0)

####escribo el inp de k=1 con la nueva carga
mdb.jobs[nameinpK].writeInput(consistencyChecking=OFF)

##

###inicializamos variables:
NeL2damageKm1=[]
datos_salida=[]
#%%bucle k. Ahora la k y la m estan acopladas. Simplemente cuando hay dagno no aumento carga y cuando termina de dagnar aumento
for k in range(1,pasosK+1):
    #%% hacemos dos inp: uno con N=1 y otro con N=2 para el primer paso
    #%%PREPARAMOS LOS DOS INP PARA CADA N. 
    #SE PREPARAN DESDE EL INP ORIGEN, PARA REVISAR POSTERIORMENTE LOS INP SI QUEREMOS
    #creamos los inp nuevos
    nameinpN1k=nameinpK+'_N1_K'+str(k)
    nameinpN2k=nameinpK+'_N2_K'+str(k)
    copy(nameinpK+'.inp',nameinpN1k+'.inp')  
    copy(nameinpK+'.inp',nameinpN2k+'.inp')  
    #abrimos el inp y buscamos la frase con las pripiedades de la UMAT o Uinter
    #Buscamos uinter='*Surface Interaction'
    cadenasCohesivo=lineasdefinteraccion(nameinpK+'.inp', '*Surface Interaction',3)  #buscamos la cadena a cambiar
    #reemplazamos la primera linea en los dos archivos de inicio
    replacement(nameinpN1k+'.inp', cadenasCohesivo[0], cadenasUinter[0])
    replacement(nameinpN2k+'.inp', cadenasCohesivo[0], cadenasUinter[0])
    #reemplazamos la segunda linea en los dos archivos de inicio
    replacement(nameinpN1k+'.inp', cadenasCohesivo[1], cadenasUinter[1])
    replacement(nameinpN2k+'.inp', cadenasCohesivo[1], cadenasUinter[1])
    #reemplazamos la tercera linea en los dos archivos de inicio. Esta es distinta
    replacement(nameinpN1k+'.inp', cadenasCohesivo[2], cadenasUinterN1[0])
    replacement(nameinpN2k+'.inp', cadenasCohesivo[2], cadenasUinterN2[0])
    #reemplazamos la linea que pide las salida de las superficies de contacto    
    cadenaFOcontact=lineasdefinteraccion(nameinpK+'.inp', 'CDISP, CSTRESS',1)  #buscamos la cadena a cambiar
    cadenaSDVcontact=['CDISP, CSTRESS, SDV\n']
    replacement(nameinpN1k+'.inp', cadenaFOcontact[0], cadenaSDVcontact[0])
    replacement(nameinpN2k+'.inp', cadenaFOcontact[0], cadenaSDVcontact[0])

    #%%CORREMO EN PARALELO LOS DOS INP DESDE EL ARCHIVO PARA PODER HACERLO CON UINTER
    myJobN1 = mdb.JobFromInputFile(name=nameinpN1k, 
        inputFileName=nameinpN1k, type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None,
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, memory=memoria, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, userSubroutine=nameUMAT, 
        scratch=actual_directory, parallelizationMethodExplicit=DOMAIN, numDomains=cpus, 
        activateLoadBalancing=False, multiprocessingMode=DEFAULT, numCpus=cpus)
    myJobN2 = mdb.JobFromInputFile(name=nameinpN2k, 
        inputFileName=nameinpN2k, type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None,
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, memory=memoria, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, userSubroutine=nameUMAT, 
        scratch=actual_directory, parallelizationMethodExplicit=DOMAIN, numDomains=cpus, 
        activateLoadBalancing=False, multiprocessingMode=DEFAULT, numCpus=cpus)

    myJobN1.submit(consistencyChecking=OFF)
    myJobN2.submit(consistencyChecking=OFF)
    myJobN1.waitForCompletion()
    myJobN2.waitForCompletion()

    #%%ABRIMOS LOS ODB 
    enerHtotalN1,sumaenerinterN1,NeL2damageKN1,NeL2damagetotalN1,tf=PMTESCabaqus.PMTESCcritEneUMAT(nameinpN1k, control,cadenanombreSDV)
    enerHtotalN2,sumaenerinterN2,NeL2damageKN2,NeL2damagetotalN2,tf=PMTESCabaqus.PMTESCcritEneUMAT(nameinpN2k, control,cadenanombreSDV)
    
    #%%decidimos que hacer en paso siguiente
    if len(NeL2damageKN1)>0.0 or len(NeL2damageKN2)>0.0:
        if enerHtotalN1<enerHtotalN2:
            odbdef=nameinpN1k 
            NeL2damage=NeL2damageKN1
            NeL2damagetotal=NeL2damagetotalN1
        else:
            odbdef=nameinpN2k
            NeL2damage=NeL2damageKN2
            NeL2damagetotal=NeL2damagetotalN2
    else:
        #si no hay nada dagnado podemos tomar cualquiera de los dos. 
        odbdef=nameinpN1k
        NeL2damage=NeL2damageKN1
        NeL2damagetotal=NeL2damagetotalN1								 
    ###salida de datos, solo para revisar o graficar
    enerHtotalN,sumaenerinterN,NeL2damageKN,NeL2damagetotalN,tf=PMTESCabaqus.PMTESCcritEneUMAT(odbdef,control,cadenanombreSDV)
    model_outputfile(sim_file,k,newload,tf,enerHtotalN,sumaenerinterN,len(NeL2damageKN),enerHtotalN2,sumaenerinterN2,len(NeL2damageKN2),len(NeL2damagetotalN))
    datos_salida.append([float(newload),float(tf),float(enerHtotalN),sumaenerinterN,len(NeL2damageKN),enerHtotalN2,sumaenerinterN2,len(NeL2damageKN2),odbdef,len(NeL2damagetotalN)])
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
        
    # print(tuple(instanceslista))
    mdb.models[nameinp_k_mas_1].InitialState(updateReferenceConfiguration=OFF, fileName=odbdef, 
              endStep=LAST_STEP, endIncrement=STEP_END, name='Predefined Field-1', 
              createStepName='Initial', instances=tuple(instanceslista))
    
    # cambiamos el nombre del step  
    mdb.models[nameinp_k_mas_1].steps.changeKey(fromName='Step-1', toName='Step-'+str(k+1))
    # cambiamos la carga
    #%%aumentamos la carga si no hay mas rotura
    if(collections.Counter(NeL2damagetotalN)==collections.Counter(NeL2damageKm1)):
        newload=incr+newload
    mdb.models[nameinp_k_mas_1].amplitudes[amplname].setValues(timeSpan=STEP, smooth=SOLVER_DEFAULT, data=((0.0, newload), (1.0, newload)))
    #%%guardamos el dagno anterior
    NeL2damageKm1=NeL2damagetotalN
    
    #escribimos el inp final de k mas 1 para empezar el siguiente paso  
    mdb.Job(name=nameinp_k_mas_1, model=nameinp_k_mas_1, description='', 
            type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, 
            memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, memory=memoria,
            explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
            modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine=nameUMAT, 
            scratch=actual_directory, resultsFormat=ODB, multiprocessingMode=DEFAULT, numDomains=cpus, numCpus=cpus, 
            numGPUs=0)

    mdb.jobs[nameinp_k_mas_1].writeInput(consistencyChecking=OFF)
    #dejamos el inp basico por donde empezara el bucle con los diferentes inicios
    nameinpK=nameinp_k_mas_1
    print('Updated applied load: ',newload)
    #%% agnadimos el nuevo odb al nuevo odb completo para la salida
    if k==1:
        nameodbcompleto=nameinp+'_completo.odb'
        copy (odbdef+'.odb',nameodbcompleto)
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
# Save the numpy array as a text file with tab delimiters:
    # Specify the format for each column
# fmtbl = ['%.6f', '%.6f', '%.6f', '%d', '%.6f', '%.6f', '%d']
# np.savetxt(Dirname+'.txt', stepdata_np, fmt=fmtbl, delimiter='\t')