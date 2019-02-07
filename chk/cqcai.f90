!*==CQCAI.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CQCAI
SUBROUTINE CQCAI(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--CQCAI5
  !***BEGIN PROLOGUE  CQCAI
  !***SUBSIDIARY
  !***PURPOSE  Quick check for SLATEC subroutines
  !            CAIRY, CBIRY
  !***LIBRARY   SLATEC
  !***CATEGORY  C10D
  !***TYPE      COMPLEX (CQCAI-C, ZQCAI-Z)
  !***KEYWORDS  QUICK CHECK, CAIRY, CBIRY
  !***AUTHOR  Amos, Don, (SNL)
  !           Goudy, Sue, (SNL)
  !           Walton, Lee, (SNL)
  !***DESCRIPTION
  !
  ! *Usage:
  !
  !        INTEGER  LUN, KPRINT, IPASS
  !
  !        CALL CQCAI (LUN, KPRINT, IPASS)
  !
  ! *Arguments:
  !
  !     LUN    :IN  is the unit number to which output is to be written.
  !
  !     KPRINT :IN  controls the amount of output, as specified in the
  !                 SLATEC Guidelines.
  !
  !     IPASS  :OUT indicates whether the test passed or failed.
  !                 A value of one is good, indicating no failures.
  !
  ! *Description:
  !
  !   CQCAI is a quick check routine for the complex Airy functions
  !    generated by subroutines CAIRY and CBIRY.
  !
  !   CQCAI generates Airy functions and their derivatives from CAIRY
  !    and CBIRY and checks them against the Wronskian evaluation
  !    in the Z plane.
  !
  !***REFERENCES  Abramowitz, M. and Stegun, I. A., Handbook
  !                 of Mathematical Functions, Dover Publications,
  !                 New York, 1964.
  !               Amos, D. E., A Subroutine Package for Bessel
  !                 Functions of a Complex Argument and Nonnegative
  !                 Order, SAND85-1018, May, 1985.
  !***ROUTINES CALLED  CAIRY, CBIRY, I1MACH, R1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   890831  Revised to meet new SLATEC standards
  !
  !***END PROLOGUE  CQCAI
  !
  !*Internal Notes:
  !   Machine constants are defined by functions I1MACH and R1MACH.
  !
  !   The parameter MQC can have values 1 (the default) for a faster,
  !   less definitive test or 2 for a slower, more definitive test.
  !
  !**End
  !
  !  Set test complexity parameter.
  !
  INTEGER MQC
  PARAMETER (MQC=1)
  !
  !  Declare arguments.
  !
  INTEGER Lun , Kprint , Ipass
  !
  !  Declare external functions.
  !
  INTEGER I1MACH
  REAL R1MACH
  EXTERNAL I1MACH , R1MACH
  !
  !  Declare local variables.
  !
  COMPLEX con1 , con2 , con3 , cv , cw , cy , w , y , z , zr
  REAL aa , ab , acw , acy , alim , arg , arzr , atol , av , a1 , a2 , ct , &
    c23 , dig , elim , eps , er , ertol , film , fnul , fpi , hpi , pi , &
    pi3 , r , rl , rm , rpi , rtpi , rzr , r1m4 , r1m5 , slak , spi , &
    st , t , tol , tpi , tpi3 , ts
  INTEGER i , icase , icl , ierr , il , ir , irb , irset , it , itl , k , &
    kdo , keps , kode , k1 , k2 , lflg , nz1 , nz2 , nz3 , nz4
  DIMENSION kdo(20) , keps(20) , t(20) , w(20) , y(20)
  !
  !***FIRST EXECUTABLE STATEMENT  CQCAI
  IF ( Kprint>=2 ) THEN
    WRITE (Lun,99001)
    99001   FORMAT (' QUICK CHECK ROUTINE FOR THE AIRY FUNCTIONS FROM ',&
      'CAIRY AND CBIRY'/)
  ENDIF
  !-----------------------------------------------------------------------
  !     Set parameters related to machine constants.
  !     TOL is the approximate unit roundoff limited to 1.0E-18.
  !     ELIM is the approximate exponential over- and underflow limit.
  !     exp(-ELIM).lt.exp(-ALIM)=exp(-ELIM)/TOL    and
  !     exp(ELIM).gt.exp(ALIM)=exp(ELIM)*TOL       are intervals near
  !     underflow and overflow limits where scaled arithmetic is done.
  !     RL is the lower boundary of the asymptotic expansion for large Z.
  !     DIG = number of base 10 digits in TOL = 10**(-DIG).
  !     FNUL is the lower boundary of the asymptotic series for large FNU.
  !-----------------------------------------------------------------------
  r1m4 = R1MACH(4)
  tol = MAX(r1m4,1.0E-18)
  atol = 100.0E0*tol
  aa = -LOG10(r1m4)
  k1 = I1MACH(12)
  k2 = I1MACH(13)
  r1m5 = R1MACH(5)
  k = MIN(ABS(k1),ABS(k2))
  elim = 2.303E0*(k*r1m5-3.0E0)
  ab = aa*2.303E0
  alim = elim + MAX(-ab,-41.45E0)
  dig = MIN(aa,18.0E0)
  slak = 3.0E0 + 4.0E0*(-LOG10(tol)-7.0E0)/11.0E0
  slak = MAX(slak,3.0E0)
  ertol = tol*10.0E0**slak
  rl = 1.2E0*dig + 3.0E0
  rm = 0.5E0*(alim+elim)
  rm = MIN(rm,200.0E0)
  rm = MAX(rm,rl+10.0E0)
  fnul = 10.0E0 + 6.0E0*(dig-3.0E0)
  IF ( Kprint>=2 ) THEN
    WRITE (Lun,99002)
    99002   FORMAT (' PARAMETERS'/5X,'TOL ',8X,'ELIM',8X,'ALIM',8X,'RL  ',8X,'FNUL',&
      8X,'DIG')
    WRITE (Lun,99003) tol , elim , alim , rl , fnul , dig
    99003   FORMAT (6E12.4/)
  ENDIF
  !-----------------------------------------------------------------------
  !     Generate angles for construction of complex Z to be used in tests.
  !-----------------------------------------------------------------------
  fpi = ATAN(1.0E0)
  hpi = fpi + fpi
  pi = hpi + hpi
  tpi = pi + pi
  rpi = 1.0E0/pi
  tpi3 = tpi/3.0E0
  spi = pi/6.0E0
  pi3 = spi + spi
  rtpi = 1.0E0/tpi
  a1 = rtpi*COS(spi)
  a2 = rtpi*SIN(spi)
  con1 = CMPLX(COS(tpi3),SIN(tpi3))
  con2 = CMPLX(a1,-a2)
  con3 = CMPLX(rpi,0.0E0)
  c23 = 2.0E0/3.0E0
  !-----------------------------------------------------------------------
  !     KDO(K), K = 1,IL  determines which of the IL angles in -PI to PI
  !     are used to compute values of Z.
  !       KDO(K) = 0  means that the index K will be used for one or two
  !                   values of Z, depending on the choice of KEPS(K)
  !              = 1  means that the index K and the corresponding angle
  !                   will be skipped
  !     KEPS(K), K = 1,IL determines which of the angles get incremented
  !     up and down to put values of Z in regions where different
  !     formulae are used.
  !       KEPS(K)  = 0  means that the angle will be used without change
  !                = 1  means that the angle will be incremented up and
  !                   down by EPS
  !     The angles to be used are stored in the T(I) array, I = 1,ITL.
  !-----------------------------------------------------------------------
  IF ( MQC/=2 ) THEN
    icl = 1
    il = 5
    DO i = 1 , il
      kdo(i) = 0
      keps(i) = 0
    ENDDO
  ELSE
    icl = 2
    il = 7
    DO i = 1 , il
      kdo(i) = 0
      keps(i) = 0
    ENDDO
    keps(2) = 1
    keps(3) = 1
    keps(5) = 1
    keps(6) = 1
  ENDIF
  i = 2
  eps = 0.01E0
  film = il - 1
  t(1) = -pi + eps
  DO k = 2 , il
    IF ( kdo(k)==0 ) THEN
      t(i) = pi*(-il+2*k-1)/film
      IF ( keps(k)/=0 ) THEN
        ts = t(i)
        t(i) = ts - eps
        i = i + 1
        t(i) = ts + eps
      ENDIF
      i = i + 1
    ENDIF
  ENDDO
  itl = i - 1
  !-----------------------------------------------------------------------
  !     Test values of Z in -PI.lt.arg(Z).le.PI near formula boundaries.
  !-----------------------------------------------------------------------
  IF ( Kprint>=2 ) THEN
    WRITE (Lun,99004)
    99004   FORMAT (' CHECKS IN THE Z PLANE'/)
  ENDIF
  lflg = 0
  DO icase = 1 , icl
    !-----------------------------------------------------------------------
    !     ICASE = 1 computes wron(AI(Z),BI(Z))     =CON3
    !     ICASE = 2 computes wron(AI(Z),AI(Z*CON1))=CON2
    !-----------------------------------------------------------------------
    DO kode = 1 , 2
      DO irset = 1 , 3
        irb = MIN(irset,2)
        DO ir = irb , 4
          !------------ switch (irset)
          SELECT CASE (irset)
            CASE (2)
              r = (2.0E0*(4-ir)+rl*(ir-1))/3.0E0
            CASE (3)
              r = (rl*(4-ir)+rm*(ir-1))/3.0E0
            CASE DEFAULT
              r = 2.0E0*(ir-1)/3.0E0
          END SELECT
          !------------ end switch
          DO it = 1 , itl
            !----------------------------------------------------------------------
            !     The following values are set before the DO 30 loop:
            !            C23 = 2/3
            !           CON1 = cmplx(cos(2PI/3),sin(2PI/3))
            !           CON2 = cmplx(cos(PI/6),-sin(PI/6)/2PI
            !           CON3 = cmplx(1/PI,0)
            !----------------------------------------------------------------------
            ct = COS(t(it))
            st = SIN(t(it))
            IF ( ABS(ct)<atol ) ct = 0.0E0
            IF ( ABS(st)<atol ) st = 0.0E0
            z = CMPLX(r*ct,r*st)
            zr = CMPLX(c23,0.0E0)*z*SQRT(z)
            rzr = REAL(zr)
            arzr = ABS(rzr)
            !-------------- Check for possible underflow or overflow
            IF ( arzr/=0.0E0 ) THEN
              arg = -arzr - 0.5E0*LOG(arzr) + 0.226E0
              arg = arg + arg
              !---------------- Skip test for this case?
              IF ( arg<(-elim) ) CYCLE
            ENDIF
            CALL CAIRY(z,0,kode,y(1),nz1,ierr)
            CALL CAIRY(z,1,kode,y(2),nz2,ierr)
            IF ( icase==1 ) THEN
              !---------------- Compare 1/PI with Wronskian of CAIRY(Z) and CBIRY(Z).
              CALL CBIRY(z,0,kode,w(1),ierr)
              CALL CBIRY(z,1,kode,w(2),ierr)
              IF ( kode==2 ) THEN
                !-----------------------------------------------------------------------
                !     When KODE = 2, the scaling factor exp(-zeta1-zeta2) is 1.0 for
                !     -PI.lt.arg(Z).le.PI/3 and exp(-2.0*zeta1) for PI/3.lt.arg(Z)
                !     .le.PI where zeta1 = zeta2 in this range. This is due to the fact
                !     that arg(Z*CON1) is taken to be in (-PI,PI) by the principal
                !     square root.
                !-----------------------------------------------------------------------
                !------------------ Adjust scaling factor.
                cv = CMPLX(arzr,0.0E0) - zr
                cv = EXP(cv)
                w(1) = w(1)*cv
                w(2) = w(2)*cv
              ENDIF
              cv = con3
            ELSE
              !---------------- Compare exp(-i*PI/6)/2PI with Wronskian of CAIRY(Z)
              !                 and CAIRY(Z*exp(2i*PI/3)).
              cv = z*con1
              CALL CAIRY(cv,0,kode,w(1),nz3,ierr)
              CALL CAIRY(cv,1,kode,w(2),nz4,ierr)
              IF ( kode==2 ) THEN
                IF ( t(it)>=pi3 ) THEN
                  !-------------------- Adjust scaling factor.
                  cv = zr + zr
                  cv = EXP(-cv)
                  w(1) = w(1)*cv
                  w(2) = w(2)*cv
                ENDIF
              ENDIF
              w(2) = w(2)*con1
              cv = con2
            ENDIF
            !-----------------------------------------------------------------------
            !     Error relative to maximum term
            !-----------------------------------------------------------------------
            av = ABS(cv)
            cw = y(1)*w(2)
            cy = y(2)*w(1)
            cy = cw - cy - cv
            acy = ABS(y(1))*ABS(w(2))
            acw = ABS(w(1))*ABS(y(2))
            av = MAX(acw,acy,av)
            er = ABS(cy)/av
            IF ( er>=ertol ) THEN
              IF ( lflg==0 ) THEN
                IF ( Kprint>=2 ) THEN
                  WRITE (Lun,99005) ertol
                  99005                 FORMAT (' CASES WHICH VIOLATE THE RELATIVE ERROR',&
                    ' TEST WITH ERTOL = ',E12.4/)
                  WRITE (Lun,99006)
                  99006                 FORMAT (' INPUT TO CAIRY AND ERROR')
                ENDIF
                IF ( Kprint>=3 ) THEN
                  WRITE (Lun,99007)
                  99007                 FORMAT (' COMPARISON VALUE AND WRONSKIAN')
                  WRITE (Lun,99008)
                  99008                 FORMAT (' RESULTS FROM CAIRY AND/OR CBIRY')
                  WRITE (Lun,99009)
                  99009                 FORMAT (' TEST CASE INDICES'/)
                ENDIF
              ENDIF
              lflg = 1
              IF ( Kprint>=2 ) THEN
                WRITE (Lun,99010) z , er
                99010               FORMAT (12X,'INPUT:    Z=',2E12.4,5X,'ERROR:   ER=',E12.4)
              ENDIF
              IF ( Kprint>=3 ) THEN
                WRITE (Lun,99011) cv , cy
                99011               FORMAT (' COMPARISON VALUE:   CV=',2E12.4/8X,&
                  'WRONSKIAN:   CY=',2E12.4)
                WRITE (Lun,99012) nz1 , y(1) , nz2 , y(2)
                99012               FORMAT (10X,'RESULTS:  NZ1=',I3,4X,'Y(1)=',2E12.4/20X,&
                  'NZ2=',I3,4X,'Y(2)=',2E12.4)
                IF ( icase==1 ) THEN
                  WRITE (Lun,99013) w(1) , w(2)
                  99013                 FORMAT (31X,'W(1)=',2E12.4/31X,'W(2)=',2E12.4)
                ELSE
                  WRITE (Lun,99014) nz3 , w(1) , nz4 , w(2)
                  99014                 FORMAT (20X,'NZ3=',I3,4X,'W(1)=',2E12.4/20X,'NZ4=',I3,4X,&
                    'W(2)=',2E12.4)
                ENDIF
                WRITE (Lun,99015) it , ir , irset , icase
                99015               FORMAT (13X,'CASE:   IT=',I3,4X,'IR=',I3,4X,'IRSET=',I3,4X,&
                  'ICASE=',I3,4X/)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  IF ( Kprint>=2 ) THEN
    IF ( lflg==0 ) THEN
      WRITE (Lun,99016)
      99016     FORMAT (' QUICK CHECKS OK')
    ELSE
      WRITE (Lun,99017)
      99017     FORMAT (' ***',5X,'FAILURE(S) FOR CAIRY IN THE Z PLANE')
    ENDIF
  ENDIF
  Ipass = 0
  IF ( lflg==0 ) Ipass = 1
  IF ( Ipass==1.AND.Kprint>=2 ) THEN
    WRITE (Lun,99018)
    99018   FORMAT (/' ****** CAIRY  PASSED ALL TESTS  ******'/)
  ENDIF
  IF ( Ipass==0.AND.Kprint>=1 ) THEN
    WRITE (Lun,99019)
    99019   FORMAT (/' ****** CAIRY  FAILED SOME TESTS ******'/)
  ENDIF
END SUBROUTINE CQCAI
