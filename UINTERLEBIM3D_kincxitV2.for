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
      !Lambda HS. 
      lambda1 = PROPS(5) 
      !Lambda HS.
      lambda2 = PROPS(6) 
      mu = PROPS(7)
      Nini = PROPS(8)
C     Force:(1) or displacement:(0) control
      control = PROPS(9)
C      WRITE(*,*) 'PROPS(9) = ', control 
C      WRITE(*,'(A,I0)') 'NINT(PROPS(9)) = ', NINT(control)
!cccccccccccccccccccccccccc inicializacion de variables
      prand = 0.25d0
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
!cccccccccccccccccccccccccc DAMAGE INITIALIZATION:
c     If we are at k=1 and j=1, damage and damagemenosK must be zero.               
c     If we are at any k and j=1, the CT damage, the CT damage and the rest        
c     of the variables must be zero.                                               
c     But both damage and damagemenosK must be preserved from the previous         
c     k, which is why we take them from statev(1).                                 
c     For the remaining cases:                                                     
c     damagemenosK must be NON‑reversible, so we must store it from one            
c     k‑step to the next.                                                          
c     It is updated at the beginning of the AMA (j=1) for any k‑step, but          
c     in the remaining j‑steps it does not change.                                 
c     It only transfers from one j‑step to the next through statev(13).            
c     And damageT must remain constant within each j, which is why it is           
c     taken from statev(11), which does not change.                                
c     But damageE is updated in each j, and this is done through statev(11),       
c     which is determined at the end of the iteration.                             
c     That is why the calculation starts using the damageE from the previous       
c     iteration, j-1.                                                               

      IF ((kinc*kstep).eq.1) THEN
          damage=0.00d0 
          damagemenosK=0.00d0
      ELSE IF ((kinc*kstep).eq.kstep) THEN
        damagemenosK=statev(1)
        damage=statev(1) 
        damageT=0.00d0 
        damageE=0.00d0
      ELSE         
        damage=statev(1) ! ==0 from k=1 and j=1 until change
        damageT=statev(9) ! ==0, from k=1 and j=1 until change
        damageE=statev(10) ! ==0, from k=1 and j=1 until change
        damagemenosK=statev(11)
      ENDIF
cccccc TANGENT STIFFNESS MATRIX CALCULATION 
c     The UINTER is called at least 2 times in each j‑step. These are j.1 and j.2. 
c     Sometimes it is called more times if convergence is not achieved.             
c     In j.1, j.2, ..., j.n, RDISP(j) = RDISP(j‑1) + DRDISP(j‑1), always; that is,  
c     the relative displacement between surfaces RDISP(j) that converged in the     
c     previous j.                                                                  
c     However, in j.1, DRDISP(j) = 0, and in j.2, ..., j.n, RDISP(j) is updated    
c     up to substep n.                                                             
c     Therefore, all these calculations must be done using RDISP(j), unlike in     
c     UMAT where STRAN(j) + DSTRAN(j) must be updated.                             
                                                                                   
c     Damage function for normal components depending on the sign of the           
c     relative displacement RDISP.                                                 
c     It must be remembered that for UINTER a positive sign indicates contact      
c     between surfaces.                                                            

      signoN=SIGN(1.d0,(-RDISP(1)))
      if (signoN.le.0.00d0) then
      !when in compression there is always normal stiffness
        funN=1.0d0 
      else 
      !when in traction there is no normal stiffness between surfaces
        funN=1.0d0-damage 
      end if
      !no tangential stiffness when the spring is damaged
      funT=1.0d0-damage 
c     Spring stiffness in terms of the damage variable
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
c     Update of the interface stiffness matrix.                                     
c     This matrix is the same within each j, and it depends on the damage           
c      from j-1                                                                  
c     It works like the two parts of the AMA: (i) FEM is computed with the          
c     damage from j‑1 in the first call to UINTER;                                   
c     in the second call to UINTER, within increment j, the energy is               
c     minimized with respect to the damage.                                         
c     It must be noted that the damage, like the rest of the state variables,       
c     does NOT change within each j;                                                
c     it only changes at the end of j.                                              


      DDSDDR(1,1)=KnnE
      DDSDDR(2,2)=KttE1
      DDSDDR(3,3)=KttE2 
      
c     stress vector calculation:
      DO i=1,3
        STRESS(i)=DDSDDR(i,i)*(RDISP(i))
      ENDDO
      taun=DSQRT(STRESS(2)**2+STRESS(3)**2)
 
c     COMPUTATION of energies (tensile and energetic) for all nodes of the interfaces.                                                                   
c     Computation of the energy for the tensile criterion.                          
      IF (RDISP(1)<0.d0) THEN
        GiT=(STRESS(1))**2.d0/(2.d0*Knn)
      ELSE
        GiT=0.d0
      ENDIF
      GiiT=(STRESS(2))**2.d0/(2.d0*Ktt1)
      GiiiT=(STRESS(3))**2.d0/(2.d0*Ktt2)
c     Computation of the energy for the energetic criterion.                         
c     RDISP represents the relative displacement at the end of the iteration.                                                  
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
C        psig1=datan2(STRESS(2)*dsqrt(Knn/Ktt1),-STRESS(1))
        psig1=datan2(ABS(taun),-STRESS(1))
	  psiGcrit1=pi/(2.d0*(1.d0-lambda1))
C        psiGcrit2=pi/(2.d0*(1.d0-lambda2))
        !No damage condition !!!!!!REVISAR SI ESTA BIEN
        IF (psig1.gt.psiGcrit1) then
        GcT = GIct*1d8 
        else
        GcT = GIct*(1.d0 + (dtan(psig1*(1.d0-lambda1)))**2.d0)
        endif
	  GcE=GcT*mu
c     calculo de energia total en la primera iteracion para el CT y ese debe ser el mismo en tomo el AMA
        IF (signoN.GT.0.D0) THEN
            GtotT=GiT+GiiT+GiiiT
        ELSE
            GtotT=GiiT+GiiiT
        ENDIF
      else 
c     If it is NOT the first iteration, the critical energies and the total SC        
c     energy (which replaces the stress vector) must be taken from the state         
c     variables, because this is equivalent to taking them from the first ITR.       
c     Keep in mind that both Gct and GcE must be taken from the first ITR            
c     and added to the dissipated energy when we open the ODB in the Python script.  

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
      
c     Computation of the total energy for any iteration for the EC, using the        
c     damage value determined at the end of the previous iteration.                  
	IF (signoN.GT.0.D0) THEN
            GtotE=GiE+GiiE+GiiiE
      ELSE
            GtotE=GiiE+GiiiE
      ENDIF


c     ***** CALCULATION of SC, EC, and final damage *****                                      
c     Although all variables are computed in each node, for the calculation            
c     of SC and EC only the springs that were NOT damaged in previous steps              
c     must be considered.                                                            
c     Therefore, we only enter the loop if damagemenosK = 0.                         

      IF (damagemenosK.LE.1d-18) THEN 
c     Irreversibility condition is satisfied.                                        
c     In the first j, at the end of the UINTER, the SC is evaluated.                 
c     And in damageE, only the different initial values are decided                  
c     depending on the SC, if we are at the end of the first increment.              
      
        ! Estimate when the first increment ends            
	  if ((KINC*KSTEP.EQ.KSTEP)) then
	      TARGET_TIME = TIME(1) + DTIME
	      IF ((GtotT.GE.Gct)) THEN !cumple CT
                damageT=ONE 
                NODNUM = NODNUM + damageT
	          ! Store vars in array
	          randnr8 = r8_normal_01( )
                IF (NODNUM < MAX_FAILED) THEN
                  FAILED_DATA(NODNUM, 1) = NODNUM
                  FAILED_DATA(NODNUM, 2) = NODE
                  FAILED_DATA(NODNUM, 3) = t/tc
                  FAILED_DATA(NODNUM, 4:6) = COORDS(1:3)
                  FAILED_DATA(NODNUM, 7) = GtotT
                  FAILED_DATA(NODNUM, 8) = GtotT/Gct
                  FAILED_DATA(NODNUM, 9) = 1.0d0+(randnr8*prand)
                  tcrand = Gct*FAILED_DATA(NODNUM, 9)

                  if (tcrand.LT.ZERO) tcrand = 1.0D-4
                  FAILED_DATA(NODNUM, 10) = GtotT/tcrand
                ENDIF
	          
	          mint2tc=MINVAL(FAILED_DATA(:,10))
	          maxt2tc=MAXVAL(FAILED_DATA(:,10))
	          deltat2tc=(maxt2tc-mint2tc)
	          mingte2gtc=MINVAL(FAILED_DATA(:, 8))
	          maxgte2gtc=MAXVAL(FAILED_DATA(:, 8))
	          deltagt2gc=maxgte2gtc-mingte2gtc
                
                IF (GtotT/tcrand.GE.(mint2tc+(ZERO*(deltat2tc)))
     &          .and.(abs(Nini-ONE).le.1d-18))
     &              THEN
                damageE=ONE

                ELSE IF (GtotT/tcrand.GE.(mint2tc+(0.10d0*(deltat2tc)))
     &          .and.(abs(Nini-two).le.1d-18))
     &          THEN
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

            ENDIF
C           End of stress condition for jAMA=0 or KINC=1
            IF (control.EQ.1) THEN
                IF ((abs(damageT).GT.1d-18).AND.(abs(damageE).GT.1d-18))
     &          THEN
                damage=ONE
C     	          WRITE(*,*) KINC, NODE, damage
                ELSE
                damage=ZERO
                ENDIF      

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
            KINC_failn(KINC,2)=NUM_BROKEN_INC
            KINC_failnm1(KINC,1)=KINC
            KINC_failnm1(KINC,2)=NODECC
            IF (KINC.GE.3) THEN
                eps0=ABS(KINC_failn(kinc,2)-KINC_failn(kinc-1,2))
C               KINC_failn(KINC,2).EQ.KINC_failn(KINC-1,2) termination check:
                IF ((eps0.LT.toler).AND.(endflag.EQ.TRUE)) THEN
                WRITE(*,*) KINC_failn(kinc,1),INT(KINC_failn(kinc,2))
                WRITE(*,*)'UINTER: dam0, CONVERGENCE MET; EXIT.'
                CALL XIT()
                
                ENDIF
            ENDIF
C       simultaneous fulfillment of nodal stress- and energy- criterion:            
 	  ELSE IF ((abs(damageT).GT.1d-18).AND.(abs(damageE).GT.1d-18)) THEN
 	      damage=ONE
            NUM_BROKEN_INC = NUM_BROKEN_INC + 1
            KINC_failn(KINC,1)=KINC !
            KINC_failn(KINC,2)=NUM_BROKEN_INC
            IF (KINC.GE.3) THEN
                eps0=ABS(KINC_failn(kinc,2)-KINC_failn(kinc-1,2))
C               KINC_failn(KINC,2).EQ.KINC_failn(KINC-1,2) termination check:
                IF ((eps0.LT.toler).AND.(endflag.EQ.TRUE)) THEN
                WRITE(*,*) KINC_failn(kinc,1),INT(KINC_failn(kinc,2))
                WRITE(*,*) 'UINTER: dam1, CONVERGENCE MET; CRACK ONSET.'
                WRITE(*,*) 
     &          'CONVERGENCE ACHIEVED: No binary delta damage detected.'
                CALL XIT()
                
                    if (filesignl.EQ.ONE) then
                    Ninip = ONE
                    IF ((GtotT.GE.Gct)) THEN !cumple CT
                damageT=ONE !entonces el nodo esta dagando por CT
                NODNUM = NODNUM + damageT
	          ! Store vars in array
!	          CALL RANDOM_NUMBER(randr8)
	          randnr8 = r8_normal_01( )
                IF (NODNUM < MAX_FAILED) THEN
                  FAILED_DATA(NODNUM, 1) = NODNUM
                  FAILED_DATA(NODNUM, 2) = NODE
                  FAILED_DATA(NODNUM, 3) = t/tc
                  FAILED_DATA(NODNUM, 4:6) = COORDS(1:3)
                  FAILED_DATA(NODNUM, 7) = GtotT
                  FAILED_DATA(NODNUM, 8) = GtotT/Gct
                  FAILED_DATA(NODNUM, 9) = 1+(randnr8*prand)
                  tcrand = Gct*FAILED_DATA(NODNUM, 9)
                  if (tcrand.LT.ZERO) tcrand = 1.0D-4
                  FAILED_DATA(NODNUM, 10) = GtotT/tcrand
                ENDIF
	          
	          mint2tc=MINVAL(FAILED_DATA(:,10))
	          maxt2tc=MAXVAL(FAILED_DATA(:,10))
	          deltat2tc=(maxt2tc-mint2tc)
	          mingte2gtc=MINVAL(FAILED_DATA(:, 8))
	          maxgte2gtc=MAXVAL(FAILED_DATA(:, 8))
	          deltagt2gc=maxgte2gtc-mingte2gtc               
            ENDIF  
             
            IF ((abs(damageT).GT.1d-18).AND.(abs(damageE).GT.1d-18)) 
     &      THEN
 	          damage=ONE
 	      ELSE
 	          damage=ZERO
 	      ENDIF
 	      
 	      else
 	              WRITE(*,*)'UINTER: dam1, CONVERGENCE MET; EXIT.'
 	              WRITE(7,*)'UINTER: dam1, CONVERGENCE MET; EXIT.',KINC
                    CALL XIT()
                    endif
                KINC_failn(KINC,1)=KINC !
                KINC_failn(KINC,2)=NUM_BROKEN_INC
                eps0=ABS(KINC_failn(kinc,2)-KINC_failn(kinc-1,2))
                
                ENDIF
            ENDIF

        ENDIF
        ENDIF !salimos de la evaluacion de damageT y del damageE

        IF ((KINC.GE.3).AND.(all(FAILED_DATA(:,10).EQ.ZERO))) THEN
            ! FAILED_DATA array has not been updated
            WRITE(*,*) 'UINTER: stress cond. array is empty, exit.'
            CALL XIT()
        ENDIF
    
      !el final del bucle de damagemenosK. Garantiza la irreversibilidad
      ENDIF

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
      statev(21)=GtotT/GcT
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
         IF ((NUM_BROKEN_INC.EQ.0) .AND. (KINC.GT.2)) THEN
         WRITE(*,*) 'No springs were broken (D=0): solver exit'
         CALL XIT()
         END IF
      END IF
      
      RETURN
      END
C     ___________________________________________
C     |||||||END OF UEXTERNALDB SUBROUTINE|||||||
