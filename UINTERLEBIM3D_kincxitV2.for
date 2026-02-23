CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                           C
C          SUBROUTINE  LEBIM with AMA                       C 
C          Authors:    Mar Munoz Reja           
C                      Jose M. Pimentel                     C
C          University of Seville, Spain,   2025             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      MODULE random_r8
        IMPLICIT NONE
        PRIVATE ! Makes everything in the module private
        PUBLIC :: r8_normal_01 ! Makes only this function public
        
        REAL(KIND=8), PARAMETER :: PI = ACOS(-1.0D0)
      
      CONTAINS
      function r8_normal_01 ( )
    !*****************************************************************************80
    !
    !! r8_normal_01() returns a unit pseudonormal R8.
    !
    !  Discussion:
    !
    !    The standard normal probability distribution function (PDF) has
    !    mean 0 and standard deviation 1.
    !
    !  Licensing:
    !
    !    This code is distributed under the MIT license.
    !
    !  Modified:
    !
    !    06 August 2013
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Output, real ( kind = rk ) R8_NORMAL_01, a normally distributed
    !    random value.
    !
          implicit none

          integer, parameter :: rk = kind ( 1.0D+00 )
            
          real ( kind = rk ) r1
          real ( kind = rk ) r2
          real ( kind = rk ) r8_normal_01
          real ( kind = rk ), parameter :: r8_pi = ACOS(-1.0d0)
          real ( kind = rk ) x
          ! Variables for seeding
          LOGICAL, SAVE :: first_call = .TRUE.
          INTEGER :: i, n, clock_val
          INTEGER, ALLOCATABLE :: seed_array(:)
          
            ! --- SEEDING LOGIC (Runs only once) ---
            IF (first_call) THEN
            ! Get size of seed array for this compiler
                CALL RANDOM_SEED(SIZE=n)      
                ALLOCATE(seed_array(n))
                CALL SYSTEM_CLOCK(COUNT=clock_val) ! Get machine time
                
                DO i = 1, n
                ! Create a unique seed array
                    seed_array(i) = clock_val + i * 37 
                END DO
                
                CALL RANDOM_SEED(PUT=seed_array)
                DEALLOCATE(seed_array)
                first_call = .FALSE.
            END IF
        
          call random_number ( harvest = r1 )
          call random_number ( harvest = r2 )
          x = sqrt ( - 2.0D+00 * log ( r1 ) ) *
     &        cos ( 2.0D+00 * r8_pi * r2 )

          r8_normal_01 = x

          return
          end function r8_normal_01
      !*****************************************************************************
      END MODULE random_r8
      
C     |||||||START OF UINTER SUBROUTINE|||||||
      SUBROUTINE UINTER(STRESS,DDSDDR,AMKI,AMSKI,FLUX,DDFDDT,DDSDDT,
     1     DDFDDR,STATEV,SED,SFD,SPD,SVD,SCD,PNEWDT,RDISP,DRDISP,
     2     TEMP,DTEMP,PREDEF,DPRED,TIME,DTIME,FREQR,CINAME,SLNAME,
     3     MSNAME,PROPS,COORDS,ALOCALDIR,DROT,AREA,CHRLNGTH,NODE,NDIR,
     4     NSTATV,NPRED,NPROPS,MCRD,KSTEP,KINC,KIT,LINPER,LOPENCLOSE,
     5     LSTATE,LSDI,LPRINT)
C     
      USE random_r8
      INCLUDE 'ABA_PARAM.INC'
C     Include the definition of the common block
      INCLUDE 'my_common.inc'
C 
      DIMENSION STRESS(NDIR),DDSDDR(NDIR,NDIR),FLUX(2),DDFDDT(2,2),
     $     DDSDDT(NDIR,2),DDFDDR(2,NDIR),STATEV(NSTATV),RDISP(NDIR),
     $     DRDISP(NDIR),TEMP(2),DTEMP(2),PREDEF(2,NPRED),DPRED(2,NPRED),
     $     TIME(2),PROPS(NPROPS),COORDS(MCRD),ALOCALDIR(MCRD,MCRD),
     $     DROT(2,2),AMKI(NDIR,NDIR),AMSKI(NDIR,NDIR)
      CHARACTER*80 CINAME,SLNAME,MSNAME
      PARAMETER (ONE=1.0d0, TWO=2.0d0, HALF=ONE/TWO,
     &          h=0.05d0, ZERO=0.0d0,
     &          three=3.0d0, four=4.0d0, five=5.0d0, six=6.0d0,
     &          seven=7.0d0, eight=8.0d0, nine=9.0d0, ten=10.0d0,
     &          eleven=11.0d0)
      COMMON /CRACK_DATA/ NUM_BROKEN_INC
      INTEGER NUM_BROKEN_INC
      INTEGER ND,PT,ELEM,NODEL,NODECT
      REAL*8 factor
      INTEGER error_ap,error_lec

      ! User-defined penalty stiffness (adjust based on material properties)
      REAL*8 :: penetration, max_penetration, K_penalty, adaptive_k
      COMMON /PENETRATION_STATE/ max_penetration
      
      real*8 DDS(ndir,ndir),psig1,psig2,lambda1,lambda2, Cets, Lets, Le
      real*8 PI,psiGcrit1,psiGcrit2,mu,Gtot,GcT,GcE
      real*8 GiT,GiiT,GiiiT,GiE,GiiE,GiiiE,t,taun,tc,tauc,sigmac,sigmaeq
      real*8 Knn,Ktt1,Ktt2, GIct,GIIct,GIIIct,Nini,Ninip,mu_s,eps0
      real*8 KnnE,KttE1,KttE2,signoN
	integer k,j,i,nprops,node,control
	integer kstep,kinc,nstatv,damage
	real*8 sse,spd,scd,rpl,drpldt,dtime,temp,dtemp,celent,pnewdt
	character*255 fullname, program_path
      COMMON /directory/ program_path
      ! Common block for global node count and initialization flag
      INTEGER INIT_FLAG, cshapes
      REAL*8 NODNUM, NODECC
      COMMON /NODE_COUNT/ INIT_FLAG
      COMMON /uext_var/ NODNUM, NODECC, endflag
      
C --- Use a dynamically allocatable array to store failed node IDs
      INTEGER, ALLOCATABLE, SAVE :: failedNodeIDs(:)
      INTEGER, PARAMETER :: MAX_FAILED = 500000
      INTEGER, SAVE :: numFailedNodes = 0
      INTEGER, SAVE :: capacity = 0
      LOGICAL :: isAlreadyFailed
      ! maximum of 500 step increments:
      REAL*8 KINC_failn(500,2), KINC_failnm1(500,2) 
      REAL*8 :: FAILED_DATA(MAX_FAILED,10), mint2tc, maxt2tc, mingte2gc,
     &          maxgte2gc, deltagt2gc, tcrand, randr8, prand, randnr8
      LOGICAL :: CALLED
      REAL*8  :: TARGET_TIME, deltat2tc, filesignl
      COMMON /FAILED_N/ FAILED_DATA, KINC_failn, mint2tc, maxt2tc, 
     &        mingte2gc, maxgte2gc, deltagt2gc, deltat2tc, tcrand, 
     &        KINC_failnm1
      DATA CALLED /.FALSE./
      
C     END OF DATA DECLARATION BLOCK      
C     ---------------------------------------------------------------
      TOLER = 1.0D-5
      ! Initialize on first call
      IF (INIT_FLAG .EQ. 0) THEN
          NODNUM = 0
          INIT_FLAG = 1
      END IF
      !
C      IF (endflag == .TRUE.) THEN
C        NODECC = ZERO
C      END IF
      
C --- Initialize array on first entry
      IF (.NOT. ALLOCATED(failedNodeIDs)) THEN
        capacity = MAX_FAILED ! Initial guess for number of failed nodes
        ALLOCATE(failedNodeIDs(capacity))
      END IF

      dds =0.d0
	psig1=0.d0
      psig2=0.d0
	stress=0.d0
	Cets = 0.d0 
	Lets = 0.d0	
      
c      !READ THE MATERIAL PROPERTIES FROM THE INPUT FILE      
      Knn = PROPS(1)
      Ktt1 = PROPS(2)
      Ktt2 = PROPS(3)
      GIct = PROPS(4) !fracture energy associated to inf. crack growth

      !Lambda de HS. No se utililiza si tomamos el criterio cuadratico
      lambda1 = PROPS(5) 
      !Lambda de HS. No se utililiza si tomamos el criterio cuadratico
      lambda2 = PROPS(6) 
      mu = PROPS(7)
      Nini = PROPS(8)
C     Force:(1) or displacement:(0) control
      control = PROPS(9)
C      WRITE(*,*) 'PROPS(9) = ', control 
C      WRITE(*,'(A,I0)') 'NINT(PROPS(9)) = ', NINT(control)
!cccccccccccccccccccccccccc inicializacion de variables
      prand = ZERO
      Pi=ACOS(-ONE)
      GiT=0.0d0
      GiiT=0.0d0
      GiiiT=0.0d0
      GiE=0.0d0
      GiiE=0.0d0 
      GiiiE=0.0d0
!      Control of stress or energy-based debonding shape:
!      cshapes=0 : stress-based
!      cshapes=1 : energy-based
      cshapes=0 
!cccccccccccccccccccccccccc INICIALIZAMOS EL DAGNO:
c     Si estamos en k=1 y j=1, el damage y el damagemenosK debe ser cero
c     Si estamos en k=cualquiera y j=1, el damage por CT, por CT y resto de variables deben ser cero
c     pero tanto el dagno como el damagemenosK deben mantenerse del k anterior, por eso la tomamos de statev(1) 
c     Para el resto de casos:
c     damagemenosK debe ser NO reversible, por eso lo debemos guardar de un paso k a otro
c     se actualiza al pricipio del AMA (j=1) para un paso k cualquiera, pero en el resto de pasos j no cambia
c     solo pasa de un paso j a otro mediante la statev(13)   
c     Y el damageT debe permanecer cte dentro de cada j, por eso se toma desde statev(11) que no cambia
c     pero damageE si se actualiza en cada j, y se hace a traves de statev(11) que se decide al final de la uinter
C     Por eso se empieza a calcular con el damageE de la iteracion anterior, j-1    
      IF ((kinc*kstep).eq.1) THEN
          damage=0.00d0 
          damagemenosK=0.00d0
      ELSE IF ((kinc*kstep).eq.kstep) THEN
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
cccccc CALCULO DE LA MATRIZ TANGENTE Y DE LAS TENSIONES 
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
      K_penalty = TWO*KnnE  ! x Knn (in consistent units)
C     *****************************************************************      
      ! Check for penetration (negative gap = penetration)
!      IF (RDISP(1).GT.ZERO) THEN
!          penetration = RDISP(1) 
!          
!          ! stiffness normal contact:
!          adaptive_k = K_penalty
!          STRESS(1) = adaptive_k * RDISP(1) 
!          DDSDDR(1,1) = adapive_k   !Derivative: d(stress(1))/d(rdisp)
!          
!      ELSE
!        ! Separation: apply full stiffness for numerical stability
!          DDSDDR(1,1) = KnnE
!      ENDIF
C     *****************************************************************

c	Actualizacion de la matriz rigidez de la interfase
c     Esta matriz es la misma dentro de cada j, y dependera del dagno del j-1.
c     Es como las dos partes del AMA: i)Se calcula FEM con el dagno j-1 en la primera llamada de la UINTER 
c     en la segunda llamada de la UINTER, dentro del incremento j, se miniminiza la energia respecto al dagno.
c     Hay que notar que el dagno, como el resto de variables de estado NO cambia dentro de cada j
c     solo cambia al final del j
      DDSDDR(1,1)=KnnE
      DDSDDR(2,2)=KttE1
      DDSDDR(3,3)=KttE2 
      
c     Calculo del vector tension:
      DO i=1,3
        STRESS(i)=DDSDDR(i,i)*(RDISP(i))
      ENDDO
      taun=DSQRT(STRESS(2)**2+STRESS(3)**2)
 
cccccccccccccccccccccccccc CALCULO de energias (tensional y energetica) para todos los nodos de las interfases       
c     calculo de energia del criterio tensional
      IF (RDISP(1)<0.d0) THEN
        GiT=(STRESS(1))**2.d0/(2.d0*Knn)
      ELSE
        GiT=0.d0
      ENDIF
      GiiT=(STRESS(2))**2.d0/(2.d0*Ktt1)
      GiiiT=(STRESS(3))**2.d0/(2.d0*Ktt2)
c     calculo de energia del criterio energetico 
C     en comparacion con la UMAT aqui el RDISP es desplazamiento relativo
c     al final de la iteración
      IF (RDISP(1)<0.d0) THEN
        GiE=Knn*((RDISP(1))**2.d0)/(2.d0)
      ELSE
        GiE=ZERO
      ENDIF
      GiiE=Ktt1*((RDISP(2))**2.d0)/(2.d0)
      GiiiE=Ktt2*((RDISP(3))**2.d0)/(2.d0)
      
cccccc calculo de la energia critica en la primera iteracion para el CT y el CE     
      if (KINC.EQ.1) then
c     en la UINTER el desplazamiento de despegue es negativo y contacto positivo   
        psig1=datan2(taun*dsqrt(Knn/Ktt1),ABS(STRESS(1)))
C        psig2=datan2(STRESS(3)*dsqrt(Knn/Ktt2),-STRESS(1))
	  psiGcrit1=pi/(2.d0*(1.d0-lambda1))
C        psiGcrit2=pi/(2.d0*(1.d0-lambda2))
        !No damage condition !!!!!!REVISAR SI ESTA BIEN
        IF (psig1.gt.psiGcrit1) then
        GcT = GIct*1d8 
        else
        GcT = GIct*(1.d0 + (dtan(psig1*(1.d0-lambda1)))**2.d0)
        endif
	  GcE=Gct*mu
c     calculo de energia total en la primera iteracion para el CT y ese debe ser el mismo en tomo el AMA
        IF (signoN.GT.0.D0) THEN
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
C     strength ratio:
      mu_s=tauc*sigmac**-ONE
C     equivalent traction module:
      IF (signoN.LT.ZERO) THEN
        sigmaeq=DSQRT((mu_s**-2.0d0)*taun**2.0d0)
      ELSE
        sigmaeq=DSQRT(STRESS(1)**2.d0 + ((mu_s**-2.0d0)*taun**2.0d0))
      ENDIF
      
c     calculo de energia total para cualquier iteracion para el CE con el dagno decidido al final de la ITR anterior	
	IF (signoN.GT.0.D0) THEN
            GtotE=GiE+GiiE+GiiiE
      ELSE
            GtotE=GiiE+GiiiE
      ENDIF


cccccc CALCULO DEL DAGNO DEL CT, CE y DEFINITIVO     
c     aunque hayamos calculados todas las variables en cada PI, para 
c     el calculo del CT y CE unicamente se debe contar con los PI NO dagnados en los pasos anteriores
c     por eso solo nos metemos en el bucle si damagemenosK=0
      IF (damagemenosK.LE.1d-18) THEN !cumple irreversibilidad
c     En el primer j, al final de la uinter, se evalua el CT         
c     Y en el damageE, unicamente se decide los difernetes inicios dependiendo del CT 
      !si estamos al final en el primer incremento
      
        ! Estimate when the first increment ends            
	  IF ((KINC*KSTEP.EQ.KSTEP)) THEN
	      TARGET_TIME = TIME(1) + DTIME
	      IF ((GtotT.GE.Gct)) THEN !cumple CT
                damageT=ONE !entonces el nodo esta dagando por CT
                NODNUM = NODNUM + damageT
	          CALL RANDOM_SEED(NODE)
!	          CALL RANDOM_NUMBER(randr8)
	          randnr8 = r8_normal_01( )
	          ! Store vars in array
                IF (NODNUM < MAX_FAILED) THEN
                  FAILED_DATA(NODNUM, 1) = NODNUM
                  FAILED_DATA(NODNUM, 2) = NODE
                  FAILED_DATA(NODNUM, 3) = t/tc
                  FAILED_DATA(NODNUM, 4:6) = COORDS(1:3)
                  FAILED_DATA(NODNUM, 7) = GtotT
                  FAILED_DATA(NODNUM, 8) = GtotT/Gct
                  FAILED_DATA(NODNUM, 9) = 1+(randnr8*prand)
                  tcrand = Gct*FAILED_DATA(NODNUM, 9)
                  tc = tcrand
                  if (tcrand.LT.ZERO) tc = 1.0D-4
                  FAILED_DATA(NODNUM, 10) = GtotT/tc
                ENDIF
	          
	          mint2tc=MINVAL(FAILED_DATA(:,10))
	          maxt2tc=MAXVAL(FAILED_DATA(:,10))
	          deltat2tc=(maxt2tc-mint2tc)
	          mingte2gtc=MINVAL(FAILED_DATA(:, 8))
	          maxgte2gtc=MAXVAL(FAILED_DATA(:, 8))
	          deltagt2gc=maxgte2gtc-mingte2gtc
!	          CALL SORT_FAILED_DATA(FAILED_DATA, NODNUM, 6)

!	          IF (NODNUM.GT.0) THEN
                
                IF (GtotT/tcrand.GE.(mint2tc+(ZERO*(deltat2tc)))
     &          .and.abs(Nini-ONE).le.1d-18) THEN
                damageE=ONE
!                    if (NODE.EQ.FAILED_DATA(i,2)) then
!                write(*,*) NODE, t/tc, GtotE/tcrand, NINT(damageE)
!                    WRITE(*,*) int(FAILED_DATA(i,1)), 
!     &          int(FAILED_DATA(i,2)), int(damageE), FAILED_DATA(i,3)

                ELSE IF (GtotT/tcrand.GE.(mint2tc+(0.10d0*(deltat2tc)))
     &          .and.    abs(Nini-TWO).le.1d-18) THEN
                damageE=ONE
!                    if (NODE.EQ.FAILED_DATA(i,2)) then
!                write(*,*) NODE, t/tc, NINT(damageE)
                
                ELSE IF (GtotT/tcrand.GE.(mint2tc+(0.20d0*(deltat2tc)))
     &          .and.    abs(Nini-three).le.1d-18) THEN
                damageE=ONE
!                    if (NODE.EQ.FAILED_DATA(i,2)) then
!                write(*,*) NODE, t/tc, NINT(damageE)
                
                ELSE IF (GtotT/tcrand.GE.(mint2tc+(0.30d0*(deltat2tc)))
     &          .and.    abs(Nini-four).le.1d-18) THEN
                damageE=ONE
!                    if (NODE.EQ.FAILED_DATA(i,2)) then
!                write(*,*) NODE, t/tc, NINT(damageE)
                
                ELSE IF (GtotT/tcrand.GE.(mint2tc+(0.40d0*(deltat2tc)))
     &          .and.    abs(Nini-five).le.1d-18) THEN
                damageE=ONE
!                    if (NODE.EQ.FAILED_DATA(i,2)) then
!                write(*,*) NODE, t/tc, NINT(damageE)
                
                ELSE IF (GtotT/tcrand.GE.(mint2tc+(0.50d0*(deltat2tc)))
     &          .and.    abs(Nini-six).le.1d-18) THEN
                damageE=ONE
!                    if (NODE.EQ.FAILED_DATA(i,2)) then
!                write(*,*) NODE, t/tc, NINT(damageE)
                
                ELSE IF (GtotT/tcrand.GE.(mint2tc+(0.60d0*(deltat2tc)))
     &          .and.    abs(Nini-seven).le.1d-18) THEN
                damageE=ONE
!                    if (NODE.EQ.FAILED_DATA(i,2)) then
!                write(*,*) NODE, t/tc, NINT(damageE)
                
                ELSE IF (GtotT/tcrand.GE.(mint2tc+(0.70d0*(deltat2tc)))
     &          .and.    abs(Nini-eight).le.1d-18) THEN
                damageE=ONE
!                    if (NODE.EQ.FAILED_DATA(i,2)) then
!                write(*,*) NODE, t/tc, NINT(damageE)
                
                ELSE IF (GtotT/tcrand.GE.(mint2tc+(0.80d0*(deltat2tc)))
     &          .and.    abs(Nini-nine).le.1d-18) THEN
                damageE=ONE
!                    if (NODE.EQ.FAILED_DATA(i,2)) then
!                write(*,*) NODE, t/tc, NINT(damageE)
                
                ELSE IF (GtotT/tcrand.GE.(mint2tc+(0.90d0*(deltat2tc)))
     &          .and.    abs(Nini-ten).le.1d-18) THEN
                damageE=ONE
!                    if (NODE.EQ.FAILED_DATA(i,2)) then
!                write(*,*) NODE, t/tc, NINT(damageE)
                
                ELSE IF (abs(Nini-eleven).le.1d-18) THEN
                damageE=ZERO
!                    if (NODE.EQ.FAILED_DATA(i,2)) then
!                write(*,*) NODE, t/tc, NINT(damageE)
!                WRITE(*,*) int(FAILED_DATA(i,1)), 
!     &          int(FAILED_DATA(i,2)), int(damageE), FAILED_DATA(i,3)
!                else
!                    write(*,*) NODE, t/tc, NINT(damageE)

            ENDIF
            IF (control.EQ.1) THEN
                IF ((abs(damageT).GT.1d-18).AND.(abs(damageE).GT.1d-18))
     &          THEN
                damage=ONE
C     	          WRITE(*,*) KINC, NODE, damage
                ELSE
                damage=ZERO
                ENDIF      
                !ya se termina el CT y se fija para 
                !el resto de pasos j con la statev(11)
            ENDIF
            
        ENDIF            
      
c     al final de cualquier j diferente de 1, se evalua el damageE porque CT ya esta evaluado
 	  ELSE IF (KINC.GT.ONE) THEN
 	      if ((GtotE.GE.GcE)) then
 	          damageE=ONE
 	      else
 	          damageE=ZERO
 	      end if
 	      
c           CC-FFM:
 	      IF ((abs(damageT).GT.1d-18).AND.(abs(damageE).LT.1d-18)) THEN
 	      damage=ZERO
 	      NUM_BROKEN_INC = NUM_BROKEN_INC
 	      NODECC=NODECC+ONE
            KINC_failn(KINC,1)=KINC !
            IF (control.EQ.0) THEN
                KINC_failn(KINC,2)=NUM_BROKEN_INC
            ELSE
                KINC_failn(KINC,2)=NUM_BROKEN_INC
            END IF
            KINC_failnm1(KINC,1)=KINC
            KINC_failnm1(KINC,2)=NODECC
            IF (KINC.GT.3) THEN
                eps0=ABS(KINC_failn(kinc,2)-KINC_failn(kinc-1,2))
C                WRITE(*, *) 'ALLIE at increment ', KINC, ':', ALLSEv
C               KINC_failn(KINC,2).EQ.KINC_failn(KINC-1,2) termination check:
                IF ((eps0.LT.toler).AND.(endflag.EQ.TRUE)) THEN
!                OPEN(UNIT=30, ACTION='READ',
!     +          file=TRIM(program_path)//'\'//'eps0_value.txt',
!     +          STATUS='OLD')
!                READ(30,*) filesignl
!                CLOSE(30)
                
!                WRITE(*,*) 'ALLSE at increment',KINC,':',
!     &          ALLWKv, ALLSEv, control
                WRITE(*,*) damage,KINC_failn(kinc,1),KINC_failn(kinc,2)
                WRITE(*,*)'UINTER: dam0, CONVERGENCE MET; EXIT.'
                CALL XIT()
                
                ELSE IF ((KINC.GT.10).AND.(endflag.EQ.TRUE)) THEN
                WRITE(*,*) 'UINTER: dam0, MAX. ITER. REACHED:',KINC
                WRITE(7,*) 'UINTER: dam0, MAX. ITER. REACHED:',KINC,':',
     &          ALLWKv, ALLSEv
                CALL XIT()
                ENDIF
            ENDIF
C       simultaneous fulfillment of nodal stress- and energy- criterion:            
 	  ELSE IF ((abs(damageT).GT.1d-18).AND.(abs(damageE).GT.1d-18)) THEN
 	      damage=ONE
! 	      GcE=Gct*mu
! 	      NODECC=NODECC+damage
            !$OMP ATOMIC
            NUM_BROKEN_INC = NUM_BROKEN_INC + 1
            KINC_failn(KINC,1)=KINC !
            IF (control.EQ.0) THEN
                KINC_failn(KINC,2)=NUM_BROKEN_INC
            ELSE
                KINC_failn(KINC,2)=NUM_BROKEN_INC
            END IF
            IF (KINC.GT.3) THEN
                eps0=ABS(KINC_failn(kinc,2)-KINC_failn(kinc-1,2))
C               KINC_failn(KINC,2).EQ.KINC_failn(KINC-1,2) termination check:
                IF ((eps0.LT.toler).AND.(endflag.EQ.TRUE)) THEN
!                WRITE(*,*) 'ALLSE at increment',KINC,':',
!     &          ALLWKv, ALLSEv, control
                WRITE(*,*) damage,KINC_failn(kinc,1),KINC_failn(kinc,2)
                WRITE(*,*) 'UINTER: dam1, CONVERGENCE MET; CRACK ONSET.'
                WRITE(*,*) 
     &          'CONVERGENCE ACHIEVED: No binary delta damage detected.'
                CALL XIT()
!                OPEN(UNIT=30, ACTION='READ',
!     +          file=TRIM(program_path)//'\'//'eps0_value.txt',
!     +          STATUS='OLD')
!                READ(30,*) filesignl
!                CLOSE(30)
                
!                IF (filesignl.EQ.ONE) THEN
!                Ninip = ONE
                
!                IF ((GtotT.GE.Gct)) THEN !cumple CT
!                damageT=ONE !entonces el nodo esta dagando por CT
!                NODNUM = NODNUM + damageT
!	          ! Store vars in array
!	          CALL RANDOM_SEED(NODE)
!!	          CALL RANDOM_NUMBER(randr8)
!	          randnr8 = r8_normal_01( )
!                IF (NODNUM < MAX_FAILED) THEN
!                  FAILED_DATA(NODNUM, 1) = NODNUM
!                  FAILED_DATA(NODNUM, 2) = NODE
!                  FAILED_DATA(NODNUM, 3) = t/tc
!                  FAILED_DATA(NODNUM, 4:6) = COORDS(1:3)
!                  FAILED_DATA(NODNUM, 7) = GtotT
!                  FAILED_DATA(NODNUM, 8) = GtotT/Gct
!                  FAILED_DATA(NODNUM, 9) = 1+(randnr8*prand)
!                  tcrand = Gct*FAILED_DATA(NODNUM, 9)
!                  if (tcrand.LT.ZERO) tcrand = 1.0D-4
!                  FAILED_DATA(NODNUM, 10) = GtotT/tcrand
!                ENDIF
!	          mint2tc=MINVAL(FAILED_DATA(:,10))
!	          maxt2tc=MAXVAL(FAILED_DATA(:,10))
!	          deltat2tc=(maxt2tc-mint2tc)
!	          mingte2gtc=MINVAL(FAILED_DATA(:, 8))
!	          maxgte2gtc=MAXVAL(FAILED_DATA(:, 8))
!	          deltagt2gc=maxgte2gtc-mingte2gtc                
!                    IF (GtotT/tcrand.GE.(mint2tc+(ZERO*(deltat2tc)))
!     &              .and.abs(Ninip-ONE).le.1d-18) THEN
!                damageE=ONE
!!                    if (NODE.EQ.FAILED_DATA(i,2)) then
!!                write(*,*) NODE, t/tc, GtotE/tcrand, NINT(damageE)
!!                    WRITE(*,*) int(FAILED_DATA(i,1)), 
!!     &          int(FAILED_DATA(i,2)), int(damageE), FAILED_DATA(i,3)
!
!                ELSE IF (GtotT/tcrand.GE.(mint2tc+(0.10d0*(deltat2tc)))
!     &          .and.    abs(Ninip-TWO).le.1d-18) THEN
!                damageE=ONE
!!                    if (NODE.EQ.FAILED_DATA(i,2)) then
!!                write(*,*) NODE, t/tc, NINT(damageE)
!                
!                ELSE IF (GtotT/tcrand.GE.(mint2tc+(0.20d0*(deltat2tc)))
!     &          .and.    abs(Ninip-three).le.1d-18) THEN
!                damageE=ONE
!!                    if (NODE.EQ.FAILED_DATA(i,2)) then
!!                write(*,*) NODE, t/tc, NINT(damageE)
!                
!                ELSE IF (GtotT/tcrand.GE.(mint2tc+(0.30d0*(deltat2tc)))
!     &          .and.    abs(Ninip-four).le.1d-18) THEN
!                damageE=ONE
!!                    if (NODE.EQ.FAILED_DATA(i,2)) then
!!                write(*,*) NODE, t/tc, NINT(damageE)
!                
!                ELSE IF (GtotT/tcrand.GE.(mint2tc+(0.40d0*(deltat2tc)))
!     &          .and.    abs(Ninip-five).le.1d-18) THEN
!                damageE=ONE
!!                    if (NODE.EQ.FAILED_DATA(i,2)) then
!!                write(*,*) NODE, t/tc, NINT(damageE)
!                
!                ELSE IF (GtotT/tcrand.GE.(mint2tc+(0.50d0*(deltat2tc)))
!     &          .and.    abs(Ninip-six).le.1d-18) THEN
!                damageE=ONE
!!                    if (NODE.EQ.FAILED_DATA(i,2)) then
!!                write(*,*) NODE, t/tc, NINT(damageE)
!                
!                ELSE IF (GtotT/tcrand.GE.(mint2tc+(0.60d0*(deltat2tc)))
!     &          .and.    abs(Ninip-seven).le.1d-18) THEN
!                damageE=ONE
!!                    if (NODE.EQ.FAILED_DATA(i,2)) then
!!                write(*,*) NODE, t/tc, NINT(damageE)
!                
!                ELSE IF (GtotT/tcrand.GE.(mint2tc+(0.70d0*(deltat2tc)))
!     &          .and.    abs(Ninip-eight).le.1d-18) THEN
!                damageE=ONE
!!                    if (NODE.EQ.FAILED_DATA(i,2)) then
!!                write(*,*) NODE, t/tc, NINT(damageE)
!                
!                ELSE IF (GtotT/tcrand.GE.(mint2tc+(0.80d0*(deltat2tc)))
!     &          .and.    abs(Ninip-nine).le.1d-18) THEN
!                damageE=ONE
!!                    if (NODE.EQ.FAILED_DATA(i,2)) then
!!                write(*,*) NODE, t/tc, NINT(damageE)
!                
!                ELSE IF (GtotT/tcrand.GE.(mint2tc+(0.90d0*(deltat2tc)))
!     &          .and.    abs(Ninip-ten).le.1d-18) THEN
!                damageE=ONE
!!                    if (NODE.EQ.FAILED_DATA(i,2)) then
!!                write(*,*) NODE, t/tc, NINT(damageE)
!                
!                ELSE IF (abs(Ninip-eleven).le.1d-18) THEN
!                damageE=ZERO
!!                    if (NODE.EQ.FAILED_DATA(i,2)) then
!!                write(*,*) NODE, t/tc, NINT(damageE)
!!                WRITE(*,*) int(FAILED_DATA(i,1)), 
!!     &          int(FAILED_DATA(i,2)), int(damageE), FAILED_DATA(i,3)
!!                else
!!                    write(*,*) NODE, t/tc, NINT(damageE)
!
!                ENDIF
                
!            ENDIF  
 	      
! 	      ELSE
!                WRITE(*,*)'UINTER: dam1, CONVERGENCE MET; EXIT.'
!                WRITE(7,*)'UINTER: dam1, CONVERGENCE MET; EXIT.',KINC
!                CALL XIT()
!           ENDIF

                KINC_failn(KINC,1)=KINC !
                IF (control.EQ.0) THEN
                    KINC_failn(KINC,2)=NUM_BROKEN_INC
                ELSE
                    KINC_failn(KINC,2)=NUM_BROKEN_INC
                END IF
                eps0=ABS(KINC_failn(kinc,2)-KINC_failn(kinc-1,2))
                
            ENDIF
            
            ENDIF
            
            IF ((KINC.GT.4).AND.(control.EQ.1)) THEN
                WRITE(*,*) 'UINTER: force control max. iteration, exit.'
                CALL XIT()
            ENDIF
        ENDIF
        ENDIF !salimos de la evaluacion de damageT y del damageE

        IF ((KINC.GT.3).AND.(all(FAILED_DATA(:,10).EQ.ZERO))) THEN
            ! FAILED_DATA array has not been updated
            WRITE(*,*) 'UINTER: stress cond. array is empty, exit.'
            CALL XIT()
        ENDIF
             
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
      ELSE IF (damagemenosK.GT.1d-5) THEN 
      !si no rompe en este paso solo suma energia a compresion
        GcE=0.d0 !no sume en el funcional de energia
        GtotE=GiE*((1.d0-signoN)/2.d0) !solo se sumaria energia si esta 
        !a compresion porque ya esta roto del paso anterior
      ELSE !si rompe en este paso suma energia a compresion 
           !y energia disipada
        GtotE=GiE*((1.d0-signoN)/2.d0)
        
      ENDIF
      
cccccccccccccccccccccccccc END NEW MAR CCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
c    ***********************************************************************************************************
c    ***********************************************************************************************************  
c	
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
      IF (RDISP(1).GT.0.d0) THEN
        statev(13)=-1.d0*STRESS(1)
      ELSE
        statev(13)=ABS(STRESS(1))
      ENDIF
      statev(14)=STRESS(2)
      statev(15)=STRESS(3)
      statev(16)=taun
      statev(17)=GiT
      statev(18)=GiiT
      statev(19)=GiiiT
      statev(20)=RDISP(1)
      statev(21)=GtotT/tc
      statev(22)=tc
	
	RETURN
	END
C     |||||||||END OF UINTER SUBROUTINE|||||||||
	
C     |||||||START OF UEXTERNALDB SUBROUTINE|||||||
C     _____________________________________________
      SUBROUTINE UEXTERNALDB(LOP, LRESTART, TIME, DTIME, KSTEP, KINC)
      INCLUDE 'ABA_PARAM.INC'
      INCLUDE 'my_common.inc'
      DIMENSION TIME(2)

      LOGICAL endflag
      COMMON /uext_var/ NODNUM, NODECC, endflag
      CHARACTER*256 OUTDIR, program_path
      COMMON /directory/ program_path
      ! 
      COMMON /CRACK_DATA/ NUM_BROKEN_INC
      INTEGER NUM_BROKEN_INC
      
      IF (LOP.EQ.1) THEN
        numRecorded = 0
      END IF
      IF (LOP.EQ.0) THEN
        ! Start of the analysis job
        ! Get the program file directory
        CALL GETOUTDIR(OUTDIR, LENOUTDIR)
        program_path=OUTDIR
        
      ENDIF
      IF (LOP .EQ. 1) THEN
         NUM_BROKEN_INC = 0  ! Clear the tally for the new increment
      END IF
      
      IF ((LOP .EQ. 1).OR.(LOP .EQ. 2)) THEN
        ! Start/End of the increment, reset the node count:
        NODECC = 0.0d0
      ENDIF
      IF ((LOP == 5)) THEN
        ! Start of the analysis step
        ! reset stress crit. nodal count to zero
        NODNUM = 0.d0
      ENDIF
      ! End of the analysis increment, endflag=TRUE 
      IF (LOP == 2) THEN 
        endflag = .TRUE.
        CALL GETOUTDIR(OUTDIR, LENOUTDIR)
        program_path=OUTDIR
        
      ELSE
        endflag = .FALSE.  
      ENDIF
      
      IF ((KINC.GE.1) .AND. (LOP .EQ. 2)) THEN
         WRITE(*,*) 'Increment ', KINC, ' concluded.'
         WRITE(*,*) 'Total springs broken (D=1):', 
     &              NUM_BROKEN_INC

!         ! TERMINATION LOGIC
!         IF (NUM_BROKEN_INC .EQ. 0) THEN
!            WRITE(*,*) 
!     &      'CONVERGENCE ACHIEVED: No binary damage detected.'
!            ! To stop Abaqus gracefully from here, we call the error utility
!            ! with a specific "Exit" code or use a terminal signal.
!            CALL XIT()
!         END IF
      END IF
      
      RETURN
      END
C     ___________________________________________
C     |||||||END OF UEXTERNALDB SUBROUTINE|||||||