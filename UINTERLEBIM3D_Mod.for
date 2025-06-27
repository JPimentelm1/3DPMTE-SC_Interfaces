CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                           C
C          SUBROUTINE  LEBIM with AMA                       C 
C          Authors:    Mar Munoz Reja           
C                                                           C
C          University Seville, Spain,  Agoust  2024         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 

      SUBROUTINE UINTER(STRESS,DDSDDR,AMKI,AMSKI,FLUX,DDFDDT,DDSDDT,
     1     DDFDDR,STATEV,SED,SFD,SPD,SVD,SCD,PNEWDT,RDISP,DRDISP,
     2     TEMP,DTEMP,PREDEF,DPRED,TIME,DTIME,FREQR,CINAME,SLNAME,
     3     MSNAME,PROPS,COORDS,ALOCALDIR,DROT,AREA,CHRLNGTH,NODE,NDIR,
     4     NSTATV,NPRED,NPROPS,MCRD,KSTEP,KINC,KIT,LINPER,LOPENCLOSE,
     5     LSTATE,LSDI,LPRINT)
C 
      INCLUDE 'ABA_PARAM.INC'
C 
      DIMENSION STRESS(NDIR),DDSDDR(NDIR,NDIR),FLUX(2),DDFDDT(2,2),
     $     DDSDDT(NDIR,2),DDFDDR(2,NDIR),STATEV(NSTATV),RDISP(NDIR),
     $     DRDISP(NDIR),TEMP(2),DTEMP(2),PREDEF(2,NPRED),DPRED(2,NPRED),
     $     TIME(2),PROPS(NPROPS),COORDS(MCRD),ALOCALDIR(MCRD,MCRD),
     $     DROT(2,2),AMKI(NDIR,NDIR),AMSKI(NDIR,NDIR)
      CHARACTER*80 CINAME,SLNAME,MSNAME
      PARAMETER(TOLER=1.D-12, ZERO=0.D0, ONE=1.D0, TWO=2.D0,
     $     HALF=ONE/TWO)
c
      INTEGER ND,PT,ELEM,NODEL,NODECT
      REAL*8 factor
      INTEGER error_ap,error_lec
      
      real*8 DDS(ndir,ndir),psig1,psig2,lambda1,lambda2, Cets, Lets, Le
      real*8 Pi,psiGcrit1,psiGcrit2,mu,Gtot,GcT,GcE
      real*8 GiT,GiiT,GiiiT,GiE,GiiE,GiiiE,t,taun,tc,tauc,sigmac
      real*8 Knn,Ktt1,Ktt2, GIct,GIIct,GIIIct, Nini
      real*8 KnnE,KttE1,KttE2,signoN
	integer k,j,i,nprops,node
	integer kstep,kinc,nstatv,damage
	real*8 sse,spd,scd,rpl,drpldt,dtime,temp,dtemp,celent,pnewdt
	character*200 fullname
      
C      real*8 DDS(ndir,ndir), psig,lambda,GIcb, Cets, Lets, Le
C      real*8 sigmacb,Pi,Gi,Gii,Gtot,ktkn,Gc,psiGcrit, mu,GcE
C      real*8 xnue,ebulk3,eg2,elam,trval,Gs,Gn,Gt,Knn,Ktt,Kss,K33
C      real*8 GiE,GiiE,KnnE,KttE,signoN
C	integer k1,k2,k,j,i,ndi,nshr,ntens,nprops,node
C	integer layer,kspt,kstep,kinc,nstatv,damage
C	real*8 sse,spd,scd,rpl,drpldt,dtime,temp,dtemp,celent,pnewdt
C	character*200 fullname

!      write(*,*) ddsddr, node, kstep,kinc,kit,statev(1)
      dds =0.d0
	psig1=0.d0
      psig2=0.d0
	stress=0.d0
	Cets = 0.d0 
	Lets = 0.d0	
      
cccccccccccccccccccccccccc definicion de NUEVAS variables CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      !READ THE MATERIAL PROPERTIES FROM THE INPUT FILE      
      Knn = PROPS(1)
      Ktt1 = PROPS(2)
      Ktt2 = PROPS(3)
      GIct = PROPS(4)
!      GIIct = PROPS(5) !estos no se utilizan si tomamos el criterio de HS
!      GIIIct = PROPS(6) !estos no se utilizan si tomamos el criterio de HS

      !Lambda de HS. No se utililiza si tomamos el criterio cuadratico
      lambda1 = PROPS(5) 
      !Lambda de HS. No se utililiza si tomamos el criterio cuadratico
      lambda2 = PROPS(6) 
      
      mu = PROPS(7)
      Nini = PROPS(8)
!      Le = PROPS(11)     !solo se utiliza si querremos hacer una Interaccion Nodo-Super. Si la hacemos sup-sup no es necesario
!cccccccccccccccccccccccccc inicializacion de variables para quitar basura:
      Pi=ACOS(-1d0)
      GiT=0.0d0
      GiiT=0.0d0
      GiiiT=0.0d0
      GiE=0.0d0
      GiiE=0.0d0 
      GiiiE=0.0d0       
!cccccccccccccccccccccccccc INICIALIZAMOS EL DAGNO:
c     Si estamos en k=1 y j=1, el damage y el damagemenosK debe ser cero
c     Si estamos en k=cualquiera y j=1, el damage por CT, por CT y resto de variables deben ser cero
c     pero tanto el dagno como el damagemenosK deben mantenerse del k anterior, por eso la tomamos de statev(1) 
c     Para el resto de casos:
c     damagemenosK debe ser NO reversible, por eso lo debemos guardar de un paso k a otro
c     se actualiza al pricipio del AMA (j=1) para un paso k cualquiera, pero en el resto de pasos j no cambia
c     solo pasa de un paso j a otro mediante la statev(13)   
c     Y el damageT debe permanecer cte dentro de cada j, por eso se toma desde statev(11) que no cambia
c     pero damageE si se actualiza en cada j, y se hace a traves de statev(11) que se decide al final de la UMAT
C     Por eso se empieza a calcular con el damageE de la iteracion anterior, j-1    
      IF ((kinc*kstep).eq.1) THEN
          damage=0.00d0 
          damagemenosK=0.00d0
      ELSEIF ((kinc*kstep).eq.kstep) THEN
        damagemenosK=statev(1)
        damage=statev(1) 
        damageT=0.00d0 
        damageE=0.00d0
      ELSE         
        damage=statev(1) ! ==0 desde k=1 y j=1 hasta que cambie
        damageT=statev(9) ! ==0, desde k=1 y j=1 hasta que cambie
        damageE=statev(10) ! ==0, desde k=1 y j=1 hasta que cambie
        damagemenosK=statev(11)
      ENDIF
cccccccccccccccccccccccccc CALCULO DE LA MATRIZ TANGENTE Y DE LAS TENSIONES 
c     La UINTER es llamada, como minimo, 2 veces en cada paso j. Esto es j.1 y j.2
c     algunas veces se llama mas veces si no encuentra convergencia
c     en j.1, j.2,..., j.n, RDISP(j)= RDISP(j-1)+ DRDISP(j-1), siempre, es decir 
c     el desplazamiento relativo entre superficies RDISP(j) que convergio en la ultima j 
c     Sin embargo, en j.1, DRDISP(j)=0, y en j.2,..., j.n, RDISP(j) se actualiza hasta el subpaso n.
c     por eso todos estos calculos se deben hacer con RDISP(j) a diferencia de la UMAT que debe actualizarse STRAN(j)+DSTRAN(j). 

c     funcion dagno para normales dependiendo del signo del desplazamiento relativo RDISP
c     se debe recordar que para la UINTER el signo positivo indica contacto entre superficies
      signoN=SIGN(1.d0,(-RDISP(1)))
      if (signoN.le.0.00d0) then
      !cuando es compresion siempre hay rigidez normal entre superficies
        funN=1.0d0 
      else 
      !cuando es traccion no hay rigidez normal entre superficies
        funN=1.0d0-damage 
      end if
      !no hay nunca rigidez tangencial entre superficies cuando esta dagnado
      funT=1.0d0-damage 
c     Rigideces del resorte dependiendo del dagno. Si esta dagnado es cero.
      KnnE=Knn*funN
      KttE1=Ktt1*funT
      KttE2=Ktt2*funT
      
c	Actualizacion de la matriz rigidez de la interfase
c     Esta matriz es la misma dentro de cada j, y dependera del dagno del j-1.
c     Es como las dos partes del AMA: i)Se calcula FEM con el dagno j-1 en la primera llamada de la UINTER 
c     en la segunda llamada de la UINTER, dentro del incremento j, se miniminiza la energia respecto al dagno.
c     Hay que notar que el dagno, como el resto de variables de estado NO cambia dentro de cada j
c     solo cambia al final del j
      DDSDDR(1,1)=KnnE
      DDSDDR(2,2)=KttE1
      DDSDDR(3,3)=KttE2
     
c     Calculo del vector STRESS en funcion de la matriz rigidez de la interfase
c     sigma_nn:
      DO i=1,3
        STRESS(i)=DDSDDR(i,i)*(RDISP(i))
      ENDDO
      taun=DSQRT(STRESS(2)**2+STRESS(3)**2)
c     sigma_tt:
C      STRESS(2)=DDSDDR(2,2)*(RDISP(2)) 
C      STRESS(3)=DDSDDR(3,3)*(RDISP(3)) 
 
cccccccccccccccccccccccccc CALCULO de energias (tensional y energetica) para todos los nodos de las interfases       
c     calculo de energia del criterio tensional
      GiT=(STRESS(1))**2.d0/(2.d0*Knn)
      GiiT=(STRESS(2))**2.d0/(2.d0*Ktt1)
      GiiiT=(STRESS(3))**2.d0/(2.d0*Ktt2)
c     calculo de energia del criterio energetico 
C     en comparacion con la UMAT aqui el RDISP es desplazamiento relativo
c     al final de a iteración
      GiE=Knn*((RDISP(1))**2.d0)/(2.d0)
      GiiE=Ktt1*((RDISP(2))**2.d0)/(2.d0)
      GiiiE=Ktt2*((RDISP(3))**2.d0)/(2.d0)
      
cccccccccccccccccccc     calculo de la energia critica en la primera iteracion para el CT y el CE     
      if (kinc.eq.1) then
c     en la UINTER el desplazamiento de despegue es negativo y contacto positivo   
        psig1=datan2(taun*dsqrt(Knn/Ktt1),-STRESS(1))
C        psig2=datan2(STRESS(3)*dsqrt(Knn/Ktt2),-STRESS(1))
	  psiGcrit1=pi/(2.d0*(1.d0-lambda1))
C        psiGcrit2=pi/(2.d0*(1.d0-lambda2))
        !No damage condition !!!!!!REVISAR SI ESTA BIEN
        IF (abs(psig1).gt.psiGcrit1) then
        GcT = GIct*1d8 
        else
        GcT = GIct*(1.d0 + (dtan(psig1*(1.d0-lambda1)))**2.d0)
        endif
	  GcE=Gct*mu
c     calculo de energia total en la primera iteracion para el CT y ese debe ser el mismo en tomo el AMA
        IF (signoN.GT.0.D0)THEN
            GtotT=GiT+GiiT+GiiiT
        ELSE
            GtotT=GiiT+GiiiT
        ENDIF
      else 
c     si NO es la primera iteracion se debe tomar las energias criticas y 
c     la energia total de CT (sustituye al vector tension) de las variables de estado
c     porque es igual que tomarlo de la primera ITR.
c     Tened en cuenta que tanto la Gct como GcE deben ser tomadas de la primera ITR
c     sume a la energia disipada cuando abramos el odb en el script de python
        psig1=statev(3)
        psig2=statev(4)
        GtotT=statev(5)
        Gct=statev(6)
	  GcE=Gct*mu
	  
	endif
C     traction vector norm:
      t = DSQRT(STRESS(1)**2.d0 + taun**2.d0)
C     normal critical strength:     
      sigmac=(DSQRT(2*GIct*KnnE))*
     @DSQRT((1.d0+(dtan(psig1*(1.d0-lambda1)))**2.d0))*DCOS(psig1)
C     shear critical strength:
      tauc=dsqrt(Ktt1/Knn)*(DSQRT(2*GIct*KnnE))*
     @DSQRT((1.d0+(dtan(psig1*(1.d0-lambda1)))**2.d0))*DSIN(psig1)
C     critical strength:   
      tc = DSQRT(sigmac**2.d0 + tauc**2.d0)
      
c     calculo de energia total para cualquier iteracion para el CE con el dagno decidido al final de la ITR anterior	
	IF (signoN.GT.0.D0)THEN
            GtotE=GiE+GiiE+GiiiE
      ELSE
            GtotE=GiiE+GiiiE
      ENDIF


cccccccccccccccccccccccccc CALCULO DEL DAGNO DEL CT, CE y DEFINITIVO     
c     aunque hayamos calculados todas las variables en cada PI, para 
c     el calculo del CT y CE unicamente se debe contar con los PI NO dagnados en los pasos anteriores
c     por eso solo nos metemos en el bucle si damagemenosK=0
      IF (damagemenosK.le.1d-18) THEN !cumple irreversibilidad
c     En el primer j, al final de la UMAT, se evalua el CT         
c     Y en el damageE, unicamente se decide los difernetes inicios dependiendo del CT 
      !si estamos al final en el primer incremento            
	  if (kinc.eq.1) then
	      NODECT=0
	      if (t.GE.tc) then !cumple CT
	      WRITE(*,*) t/tc 
C	      CALL GETNORMAL('TOP_SURF',NODE,TIME(1),COORDS,NORMAL,STATUS)
	          damageT=1.d0      !entonces el PI esta dagando por CT
	          !inicio todo dagnado del CT     
                if (abs(Nini-2.d0).le.1d-18) damageE=1.d0
                !inicio nada dagnado del CE     
                if (abs(Nini-1.d0).le.1d-18) damageE=0.d0     
c               write(*,*)Nini, abs(Nini-2.d0)
            end if !ya se termina el CT y se fija para 
                   !el resto de pasos j con la statev(11)
c     al final de cualquier j diferente de 1, se evalua el damageE porque CT ya esta evaluado
 	  else
 	      if (GtotE.ge.GcE) then
 	          damageE=1.d0
 	      else 
 	          damageE=0.d0
 	      end if
        end if !salimos de la evalucaion de damageT y del damageE

c     Despues de hacer el CE se debe combinar con el CT
 	  if ((abs(damageT).le.1d-18).or.(abs(damageE).le.1d-18)) then
 	      damage=0.d0
 	  else
 	      damage=1.0
        end if
      !el final del bucle de damagemenosK. Garantiza la irreversibilidad
      ENDIF 

C     fuera del bucle que garantiza la irreversibilidad preparamos la energia y energia citica por PI
c     para calcular despues el funcional de energia en cada N. Para todos los PI de las interfases
c     esten dagnados o no. Notar que algun PI dagnado puede cambiar de traccion a compresion en algun k
c     por eso debe hacerse fuera del bucle anterior
C     ver apuntes de goodnote       
      IF (damage.le.1d-18) THEN
      !se hace cero para que no sume en el funcional de energia
        GcE=0.d0 
        GtotE=GiE+GiiE+GiiiE
      ELSEIF (damagemenosK.GT.1d-5) THEN 
      !si no rompe en este paso solo suma energia a compresion
        GcE=0.d0 !no sume en el funcional de energia
        GtotE=GiE*((1.d0-signoN)/2.d0) !solo se sumaria energia si esta 
        !a compresion porque ya esta roto del paso anterior
      ELSE !si rompe en este paso suma energia a compresion 
           !y energia disipada
        GtotE=GiE*((1.d0-signoN)/2.d0)
      ENDIF
      
cccccccccccccccccccccccccc END NEW MARCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
c    ***********************************************************************************************************
c    ***********************************************************************************************************  
c	salida de datos. Todos los datos han de ser utilizado con anterioridad
      statev(1)=damage
      statev(2)=node
      statev(3)=psig1
      statev(4)=psig2
      statev(5)=GtotT
      statev(6)=Gct
      statev(7)=GtotE
      statev(8)=GcE
	statev(9)=damageT
	statev(10)=damageE
      statev(11)=damagemenosK
      statev(12)=area
      statev(13)=STRESS(1)
      statev(14)=STRESS(2)
      statev(15)=STRESS(3)
      statev(16)=DSQRT(STRESS(1)**2 + STRESS(2)**2 + STRESS(3)**2)
C	IF (KSTEP*KINC.EQ.KSTEP) then
C	  CALL GETINTERNAL("TOP ARM-1", 7, 0, NODEL, JRCD)
C	  WRITE INSTANCE NAME, LOCAL NODE NUMBER
C	  WRITE(*,*) CMNAME, CPNAME
C	  WRITE(*,*) KSTEP,damage,RDISP,STRESS(1)
C      ENDIF
      
	RETURN
	END 