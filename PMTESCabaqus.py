# ********************************************************************************
# POSTPROCESADO ODB PARA ORDENAR NODOS

# ********************************************************************************

#-------------------------------------------------------
#||||	 IMPORTACION DE MODULOS Y APERTURA DE ODB	||||		   
#-------------------------------------------------------
from os import chdir,path,mkdir
import os
from shutil import move, rmtree
from abaqus import *
from abaqusConstants import *
from odbAccess import*
import visualization
#from sys import argv
from  math import sqrt #para raiz cuadrada
from string import replace
import sys
import numpy as np
import random
#from pickle import dump
#import pickle
# Redirect stdout to the command prompt
sys.stdout = sys.__stdout__
# np.random.seed() #the sequence of generated random numbers will be different 
                # each time you run the code, as it will depend on the exact moment the seed is set,
                # i.e. it is dependent on the machine time.
# to obtain the same random sequence for each program run:
np.random.seed(1357111317)

def PMTESCcriTen_carga(name_files,carga,cadenanombreSDV):
    odb = openOdb(name_files + '.odb')
    myAssembly = odb.rootAssembly
    old_load=float(carga)
	#-------------------------------------------------------

	#--------------------------------------------------------------------------
	#|||||||||||||    		 DEFINICION DE VARIABLES            |||||||||||||||
	#--------------------------------------------------------------------------

    #InterElementSet = odb.rootAssembly.elementSets[Eset]
    key_step = odb.steps.keys()
    lastFrame = odb.steps[key_step[0]].frames[-1]
    #se toma el frames[1] para poder hacer el calculo sin nada de interfase rota
    path_SDV5=odb.steps[key_step[0]].frames[-1].fieldOutputs['SDV5     '+cadenanombreSDV].values 
    path_SDV6=odb.steps[key_step[0]].frames[-1].fieldOutputs['SDV6     '+cadenanombreSDV].values 
	#-----------------------------------------------------------------------------------------------------------------------
	#||||||||||||||||||CREACION DE TABLA CON LOS PUNTOS DE INTEGRACION A DANAR ORDENADOS|||||||||||||||||||
	#-----------------------------------------------------------------------------------------------------------------------
    nodosdamT=[]
    for k in range(len(path_SDV5)):
        GtotT=path_SDV5[k].data #la G en cada nodo segun el criterio tensional
        GcT=path_SDV6[k].data #la Gc en cada nodo segun el criterio tensional
        if GcT!=0:
            fc=sqrt(GtotT/GcT)
        else:
            fc=0
        nodosdamT.append(fc)
	#    prueba.append([fc,GtotT,GcT])
    nodosdamT.sort(key=lambda nodosdamT:nodosdamT,reverse=True)
	#prueba.sort(key=lambda prueba:prueba[0],reverse=True)
	###tomo la del ultimo nodo que sera el mas tensionado
    new_load=old_load/nodosdamT[0]
    # -------------------------------------------------------------------------
# 	EXTRACTION OF NODE COORDINATES AND SDV ALONG THE INITIAL CRACK FRONT:
    # -------------------------------------------------------------------------
    instanceobj=[]
    for instanceName in myAssembly.instances.keys():
        instanceobj=myAssembly.instances.keys(instanceName)
    crackfrontnodes = odb.rootAssembly.nodeSets['INITIALFRONT']
    crackfrontelems = odb.rootAssembly.elementSets['INITIALFRONT']
    crackfront_coords=lastFrame.fieldOutputs['COORD'].getScalarField(componentLabel='COOR3',)
    # crackfront_sdvi=lastFrame.fieldOutputs['SDV17']
    crackfront_nfield=crackfront_coords.getSubset(region=crackfrontnodes,position=NODAL)
    # crackfront_G1field=crackfront_sdvi.getSubset(region=crackfrontelems, position=ELEMENT_NODAL)
    crackf_coordvals=crackfront_nfield.values
    # crackf_G1vals=crackfront_G1field.values
    initialfront_coords=[]
    initialfront_G1data=[]
    for node in crackf_coordvals:
        initialfront_coords.append(node.data)
    initialfront_coords.sort(key=lambda initialfront_coords:initialfront_coords,reverse=False)
    # for elem in crackf_G1vals:
    #     initialfront_G1data.append(elem.data)
    # initialfront_G1data.sort(key=lambda initialfront_G1data:initialfront_G1data[:],reverse=False)
    
    # print(initialfront_coords)
    # print(initialfront_G1data)
    odb.close()
    return new_load
# %%
def PMTESCcritEneUMAT(name_files, control,cadenanombreSDV):
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
    
    nodecnt=0 #node counter of nodes which define stress condition
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
        p=0.2 # 0, 0.1, 0.2, 0.3,... constant of multiplication
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
    
    # concentrated force extraction:
    dcbtoparm = odb.rootAssembly.nodeSets['TOPRP']
    tf_field=lastFrame.fieldOutputs['TF'].getScalarField(componentLabel='TF2',)
    dcbtf2_field=tf_field.getSubset(region=dcbtoparm,position=NODAL).values[0].data
    
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
    Work = 2*(path_HisRegEne.historyOutputs['ALLWK'].data[-1][1])
    eneDefSolidos = path_HisRegEne.historyOutputs['ALLSE'].data[-1][1]
    #tomamos 'ALLSE' en lugar de 'ALLIE' porque estrictamente es la energia de deformacion
    #pero podemos replantearnos poner la otra por el hourglassing
    if control==1:
        eneHtotal = eneDefSolidos+sumaenerinter-Work
    else:
        eneHtotal = eneDefSolidos+sumaenerinter
#
    # Initialize a list to store the results
    energy_data = []
    # Loop over all steps in the odb
    for step_name in odb.steps.keys():
        # print(f"Processing Step: {step_name}...")
        step = odb.steps[step_name]
        # Get the history region for the whole model from the STEP object
        # This is where whole-model energy variables are stored.
        whole_model_region = step.historyRegions['Assembly ASSEMBLY']
        allwk_history = whole_model_region.historyOutputs['ALLWK']
        allse_history = whole_model_region.historyOutputs['ALLSE']
        
        PI_frame=[]
        deltaR_frame=[]
        # Loop over all frames (increments) in the current step:
        for frame in step.frames:
            # initialize interface energy to zero at each increment:
            interf_strnenergy=float(0)
            interf_energdiss=float(0)
            path_SDV1=frame.fieldOutputs['SDV1     '+cadenanombreSDV].values
            # nodeglob=frame.fieldOutputs['SDV2     '+cadenanombreSDV].values
            # path_SDV11=frame.fieldOutputs['SDV11    '+cadenanombreSDV].values
            path_SDV12=frame.fieldOutputs['SDV12    '+cadenanombreSDV].values
            # t2tc=frame.fieldOutputs['SDV21    '+cadenanombreSDV].values
            # tcrit=frame.fieldOutputs['SDV22    '+cadenanombreSDV].values
            GtotE=frame.fieldOutputs['SDV7     '+cadenanombreSDV].values
            GcE=frame.fieldOutputs['SDV8     '+cadenanombreSDV].values
            for k in range(len(path_SDV12)):
                GtE=GtotE[k].data
                GcrE=GcE[k].data
                areacontacto=path_SDV12[k].data
                interf_strnenergy=interf_strnenergy+(GtE*areacontacto)+(GcrE*areacontacto)
                if path_SDV1[k].data==1.:
                    interf_energdiss=interf_energdiss + (GcrE*areacontacto)
            # 
            if control==1:
                PIf = allse_history.data[frame.frameId][1]+interf_strnenergy-2*allwk_history.data[frame.frameId][1]
                PI_frame.append(PIf)
            else:
                PIf = allse_history.data[frame.frameId][1]+interf_strnenergy
                PI_frame.append(PIf)
            deltaR_frame.append(interf_energdiss)
            # The .data attribute for these objects is a tuple, e.g., (time, value)
            # We are interested in the value, which is at index 1.
            if frame.frameId<1: 
                allwk0 = 2*allwk_history.data[frame.frameId][1]
                allse0 = allse_history.data[frame.frameId][1]
                if control==1:
                    totalPot0 = allse0+interf_strnenergy-allwk0
                else:
                    totalPot0 = allse0+interf_strnenergy
                
                # Append the extracted data as a list to our main list
                energy_data.append([frame.frameId, totalPot0, allse0, 0.0, deltaR_frame[0], 0.0])
            
            else:        
                allwk_value = 2*allwk_history.data[frame.frameId][1]
                allse_value = allse_history.data[frame.frameId][1]
                if control==1:
                    totalPotenergy = allse_value+interf_strnenergy-allwk_value
                else:
                    totalPotenergy = allse_value+interf_strnenergy
                
                delPI=totalPotenergy-PI_frame[1]
                delR=deltaR_frame[-1]
                # Append the extracted data as a list
                energy_data.append([frame.frameId, totalPotenergy, allse_value, delPI, delR, delPI+delR])
        # save the numpy energy evolution as a .txt file:
        tbl_format=['%d', '%.10e', '%.10e', '%.12e', '%.12e', '%.12e']
        np.savetxt(name_files+'_delPIdelRevol.txt', np.array(energy_data), delimiter='\t', fmt=tbl_format)
        print('Total energy evolution array:')
        print(np.array(energy_data))
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
    return eneHtotal,interf_strnenergy,NeL2damageK,NeL2damagetotal,dcbtf2_field,energy_data
#
def PMTESCcritEneSubr(name_files, control,cadenanombreSDV, workdir):
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
    
    nodecnt=0 #node counter of nodes which define stress condition
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
        p=0.0 # 0, 0.1, 0.2, 0.3,... constant of multiplication
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
    # Initialize a list to store the results
    energy_data = []
    # Loop over all steps in the odb
    for step_name in odb.steps.keys():
        # print(f"Processing Step: {step_name}...")
        step = odb.steps[step_name]
        # Get the history region for the whole model from the STEP object
        # This is where whole-model energy variables are stored.
        whole_model_region = step.historyRegions['Assembly ASSEMBLY']
        allwk_history = whole_model_region.historyOutputs['ALLWK']
        allse_history = whole_model_region.historyOutputs['ALLSE']
        
        PI_frame=[]
        deltaR_frame=[]
        # Loop over all frames (increments) in the current step:
        for frame in step.frames:
            # initialize interface energy to zero at each increment:
            interf_strnenergy=float(0)
            interf_energdiss=float(0)
            # initialize fracture surface defined by stress condition:
            area_stressc=float(0)
            atotal=float(0)
            # initialize fracture surface fulfilling the coupled criterion FFM:
            area_ccffm=float(0)
            path_SDV1=frame.fieldOutputs['SDV1     '+cadenanombreSDV].values
            # nodeglob=frame.fieldOutputs['SDV2     '+cadenanombreSDV].values
            # path_SDV11=frame.fieldOutputs['SDV11    '+cadenanombreSDV].values
            path_SDV12=frame.fieldOutputs['SDV12    '+cadenanombreSDV].values
            # t2tc=frame.fieldOutputs['SDV21    '+cadenanombreSDV].values
            # tcrit=frame.fieldOutputs['SDV22    '+cadenanombreSDV].values
            GtotE=frame.fieldOutputs['SDV7     '+cadenanombreSDV].values
            GcE=frame.fieldOutputs['SDV8     '+cadenanombreSDV].values
            damT=frame.fieldOutputs['SDV9     '+cadenanombreSDV].values
            damE=frame.fieldOutputs['SDV10    '+cadenanombreSDV].values
            for k in range(len(path_SDV12)):
                GtE=GtotE[k].data
                GcrE=GcE[k].data
                damTval=damT[k].data
                damEval=damE[k].data
                areacontacto=path_SDV12[k].data
                atotal=areacontacto+atotal
                interf_strnenergy=interf_strnenergy+(GtE*areacontacto)+(GcrE*areacontacto)
                # if damTval==1. and damEval==1.:
                    
                if damTval==1.:
                    area_stressc = area_stressc + areacontacto
                if damTval==1. and damEval==1.:
                    area_ccffm = area_ccffm + areacontacto
                    interf_energdiss=interf_energdiss + (GcrE*areacontacto)
                    print(path_SDV1[k].data,interf_strnenergy,interf_energdiss)
            # 
            if control==1:
                PIf = allse_history.data[frame.frameId][1]+interf_strnenergy-allwk_history.data[frame.frameId][1]
                PI_frame.append(PIf)
            else:
                PIf = allse_history.data[frame.frameId][1]+interf_strnenergy
                PI_frame.append(PIf)
            deltaR_frame.append(interf_energdiss)
            # The .data attribute for these objects is a tuple, e.g., (time, value)
            # We are interested in the value, which is at index 1.
            if frame.frameId<1: 
                allwk0 = 1*allwk_history.data[frame.frameId][1]
                allse0 = allse_history.data[frame.frameId][1]
                if control==1:
                    totalPot0 = allse0+interf_strnenergy-allwk0
                else:
                    totalPot0 = allse0+interf_strnenergy
                
                # Append the extracted data as a list to our main list
                energy_data.append([frame.frameId, totalPot0, allse0, 0.0, deltaR_frame[0], 0.0, area_stressc, area_ccffm])
            
            else:        
                allwk_value = 1*allwk_history.data[frame.frameId][1]
                allse_value = allse_history.data[frame.frameId][1]
                if control==1:
                    totalPotenergy = allse_value+interf_strnenergy-allwk_value
                else:
                    totalPotenergy = allse_value+interf_strnenergy
                
                delPI=totalPotenergy-PI_frame[1]
                delR=deltaR_frame[-1]
                # Append the extracted data as a list
                energy_data.append([frame.frameId, totalPotenergy, allse_value, delPI, delR, delPI+delR, area_stressc, area_ccffm])
            # print(frame.frameId,totalPotenergy,atotal,area_stressc,area_ccffm,interf_energdiss,interf_strnenergy)
        # save the numpy energy evolution as a .txt file:
        tbl_format=['%d', '%.10e', '%.10e', '%.12e', '%.12e', '%.12e', '%.6e', '%.6e']
        np.savetxt(name_files+'_delPIdelRevol.txt', np.array(energy_data), delimiter='\t', fmt=tbl_format)
        print('Total energy evolution array:')
        print(np.array(energy_data))
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
    return energy_data[-1][1],interf_strnenergy,NeL2damageK,NeL2damagetotal,dcbtf2_field,energy_data
   