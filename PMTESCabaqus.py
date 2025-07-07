# ********************************************************************************
# POSTPROCESADO ODB PARA ORDENAR NODOS

# ********************************************************************************

#-------------------------------------------------------
#||||	 IMPORTACION DE MODULOS Y APERTURA DE ODB	||||		   
#-------------------------------------------------------
from os import chdir,path,mkdir
from shutil import move, rmtree
from odbAccess import*
#from sys import argv
from  math import sqrt #para raiz cuadrada
from string import replace
import sys
#from pickle import dump
#import pickle
# Redirect stdout to the command prompt
sys.stdout = sys.__stdout__

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
    print(initialfront_coords)
    # print(initialfront_G1data)
    odb.close()
    return new_load
# %%
def PMTESCcritEneUMAT(name_files, control,cadenanombreSDV):
    odb = openOdb(name_files + '.odb')    
    myAssembly = odb.rootAssembly
    #--------------------------------------------------------------------------
    #||||||||||||| DEFINICION DE CAMINOS DE VARIABLES EN ABAQUS         |||||||||||||||
    #--------------------------------------------------------------------------  
    key_step = odb.steps.keys()
    lastFrame = odb.steps[key_step[0]].frames[-1]
    #InterElementSet = odb.rootAssembly.elementSets[Eset]    
    #path_despl=odb.steps[key_step[0]].frames[-1].fieldOutputs['U'].\
    #    getSubset(region=InterNodeSet).values
    path_SDV1=odb.steps[key_step[0]].frames[-1].fieldOutputs['SDV1     '+cadenanombreSDV].values
    path_SDV7=odb.steps[key_step[0]].frames[-1].fieldOutputs['SDV7     '+cadenanombreSDV].values  
    path_SDV8=odb.steps[key_step[0]].frames[-1].fieldOutputs['SDV8     '+cadenanombreSDV].values 
    path_SDV11=odb.steps[key_step[0]].frames[-1].fieldOutputs['SDV11    '+cadenanombreSDV].values
    path_SDV12=odb.steps[key_step[0]].frames[-1].fieldOutputs['SDV12    '+cadenanombreSDV].values
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
        GtotE=path_SDV7[k].data
        GcE=path_SDV8[k].data
        areacontacto=path_SDV12[k].data
        sumaenerinter=sumaenerinter+GtotE*areacontacto+GcE*areacontacto
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
    odb.close()
    return eneHtotal,sumaenerinter,NeL2damageK,NeL2damagetotal,dcbtf2_field
#
#    