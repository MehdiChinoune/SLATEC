!*==SDASQC.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK SDASQC
SUBROUTINE SDASQC(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--SDASQC5
  !*** Start of declarations inserted by SPAG
  REAL SDJAC1 , SDJAC2 , SDRES1
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  SDASQC
  !***PURPOSE  Quick check for SLATEC routine SDASSL.
  !***LIBRARY   SLATEC (DASSL)
  !***CATEGORY  I1A2
  !***TYPE      SINGLE PRECISION (SDASQC-S, DDASQC-D)
  !***KEYWORDS  QUICK CHECK, SDASSL
  !***AUTHOR  PETZOLD, LINDA R., (LLNL)
  !             COMPUTING AND MATHEMATICS RESEARCH DIVISION
  !             LAWRENCE LIVERMORE NATIONAL LABORATORY
  !             L - 316, P.O. BOX 808,
  !             LIVERMORE, CA.    94550
  !***DESCRIPTION
  !       Demonstration program for SDASSL.
  !
  !       SDASSL is used to solve two simple problems,
  !       one with a full Jacobian, the other with a banded Jacobian,
  !       with the 2 appropriate values of info(5) in each case.
  !       If the errors are too large, or other difficulty occurs,
  !       a warning message is printed.  All output is on unit LUN.
  !
  !       To run the demonstration problems with full print, call
  !       SDASQC with KPRINT = 3.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  EDIT2, SDASSL, SDJAC1, SDJAC2, SDRES1, SDRES2
  !***REVISION HISTORY  (YYMMDD)
  !   851113  DATE WRITTEN
  !   880608  Revised to meet new LDOC standards.
  !   880615  Revised to meet new prologue standards.
  !   891016  Converted to 4.0 format (FNF).
  !   901001  Minor improvements to prologue.  (FNF)
  !   901003  Added some error prints and improved output formats.  (FNF)
  !   901009  Corrected GAMS classification code.  (FNF)
  !   901009  Changed AMAX1 to MAX.  (FNF)
  !   901030  Made all declarations explicit; added 1P's to formats. (FNF)
  !***END PROLOGUE  SDASQC
  !
  INTEGER Lun , Kprint , Ipass
  !
  EXTERNAL EDIT2 , SDASSL , SDJAC1 , SDJAC2 , SDRES1 , SDRES2
  !
  INTEGER i , idid , info(15) , iout , ipar(1) , ires , iwork(45) , j190 , &
    j290 , liw , lrw , ml , mu , neq , nerr , nfe , nje , nout , nqu , &
    nst
  REAL atol , delta(25) , dtout , er , er1 , er2 , erm , ero , hu , rpar(1)&
    , rtol , rwork(550) , t , tout , tout1 , y(25) , yprime(25) , yt1 , &
    yt2
  !
  DATA tout1/1.0E0/ , dtout/1.0E0/
  !
  !***FIRST EXECUTABLE STATEMENT  SDASQC
  Ipass = 1
  nerr = 0
  rtol = 0.0E0
  atol = 1.0E-3
  lrw = 550
  liw = 45
  !
  ! FIRST PROBLEM
  !
  neq = 2
  nout = 10
  IF ( Kprint>=2 ) WRITE (Lun,99001) neq , rtol , atol
  99001 FORMAT ('1'/1X,' DEMONSTRATION PROGRAM FOR SDASSL',///1X,&
    ' PROBLEM 1..   LINEAR DIFFERENTIAL/ALGEBRAIC SYSTEM..',/1X,&
    ' X1DOT + 10.0*X1 = 0,  X1 + X2 = 1',/1X,&
    ' X1(0) = 1.0, X2(0) = 0.0',/1X,' NEQ =',I2/1X,' RTOL =',1P,E10.1,&
    '   ATOL =',E10.1)
  !
  DO j190 = 1 , 2
    DO i = 1 , 15
      info(i) = 0
    ENDDO
    IF ( j190==2 ) info(5) = 1
    !
    IF ( Kprint>2 ) WRITE (Lun,99002) info(5)
    99002   FORMAT (////1X,' INFO(5) =',I3//6X,'T',15X,'X1',14X,'X2',12X,'NQ',6X,&
      'H',12X/)
    !
    t = 0.0E0
    y(1) = 1.0E0
    y(2) = 0.0E0
    yprime(1) = -10.0E0
    yprime(2) = 10.0E0
    tout = tout1
    ero = 0.0E0
    DO iout = 1 , nout
      CALL SDASSL(SDRES1,neq,t,y,yprime,tout,info,rtol,atol,idid,rwork,lrw,&
        iwork,liw,rpar,ipar,SDJAC1)
      hu = rwork(7)
      nqu = iwork(8)
      IF ( Kprint>2 ) THEN
        WRITE (Lun,99003) t , y(1) , y(2) , nqu , hu
        99003       FORMAT (1X,1P,E15.5,E16.5,E16.5,I6,E14.3)
      ENDIF
      !
      IF ( idid<0 ) EXIT
      yt1 = EXP(-10.0E0*t)
      yt2 = 1.0E0 - yt1
      er1 = ABS(yt1-y(1))
      er2 = ABS(yt2-y(2))
      er = MAX(er1,er2)/atol
      ero = MAX(ero,er)
      IF ( er>1000.0E0 ) THEN
        IF ( Kprint>=2 ) WRITE (Lun,99010) t
        !
        nerr = nerr + 1
      ENDIF
      tout = tout + dtout
    ENDDO
    IF ( idid<0 ) THEN
      IF ( Kprint>=2 ) WRITE (Lun,99011) idid , t
      !
      nerr = nerr + 1
    ENDIF
    nst = iwork(11)
    nfe = iwork(12)
    nje = iwork(13)
    IF ( Kprint>2 ) WRITE (Lun,99012) nst , nfe , nje , ero
    !
  ENDDO
  !
  ! SECOND PROBLEM
  !
  neq = 25
  ml = 5
  mu = 0
  iwork(1) = ml
  iwork(2) = mu
  nout = 5
  IF ( Kprint>=2 ) WRITE (Lun,99004) neq , ml , mu , rtol , atol
  99004 FORMAT ('1'/1X,' DEMONSTRATION PROGRAM FOR SDASSL',///1X,&
    ' PROBLEM 2.. YDOT = A * Y , WHERE ',&
    ' A IS A BANDED LOWER TRIANGULAR MATRIX',/1X,&
    '  DERIVED FROM 2-D ADVECTION PDE',/1X,' NEQ =',I3,'   ML =',I2,&
    '   MU =',I2/1X,' RTOL =',1P,E10.1,'   ATOL =',E10.1)
  !
  DO j290 = 1 , 2
    DO i = 1 , 15
      info(i) = 0
    ENDDO
    info(6) = 1
    IF ( j290==2 ) info(5) = 1
    !
    IF ( Kprint>2 ) WRITE (Lun,99005) info(5)
    99005   FORMAT (////1X,' INFO(5) =',I3//6X,'T',14X,'MAX.ERR.',5X,'NQ',6X,'H'/)
    !
    t = 0.0E0
    DO i = 2 , neq
      y(i) = 0.0E0
    ENDDO
    y(1) = 1.0E0
    DO i = 1 , neq
      delta(i) = 0.0E0
    ENDDO
    !        Following is to initialize YPRIME.
    CALL SDRES2(t,y,delta,yprime,ires,rpar,ipar)
    tout = 0.01E0
    ero = 0.0E0
    DO iout = 1 , nout
      CALL SDASSL(SDRES2,neq,t,y,yprime,tout,info,rtol,atol,idid,rwork,lrw,&
        iwork,liw,rpar,ipar,SDJAC2)
      CALL EDIT2(y,t,erm)
      hu = rwork(7)
      nqu = iwork(8)
      IF ( Kprint>2 ) WRITE (Lun,99006) t , erm , nqu , hu
      99006     FORMAT (1X,1P,E15.5,E14.3,I6,E14.3)
      !
      IF ( idid<0 ) EXIT
      er = erm/atol
      ero = MAX(ero,er)
      IF ( er>1000.0E0 ) THEN
        IF ( Kprint>=2 ) WRITE (Lun,99010) t
        nerr = nerr + 1
      ENDIF
      tout = tout*10.0E0
    ENDDO
    IF ( idid<0 ) THEN
      IF ( Kprint>=2 ) WRITE (Lun,99011) idid , t
      nerr = nerr + 1
    ENDIF
    nst = iwork(11)
    nfe = iwork(12)
    nje = iwork(13)
    IF ( Kprint>2 ) WRITE (Lun,99012) nst , nfe , nje , ero
    !
  ENDDO
  IF ( Kprint>=2 ) WRITE (Lun,99007) nerr
  !
  99007 FORMAT (////1X,' NUMBER OF ERRORS ENCOUNTERED =',I3)
  !
  IF ( nerr>0 ) Ipass = 0
  IF ( (Ipass==1).AND.(Kprint>1) ) WRITE (Lun,99008)
  99008 FORMAT (/,' ----------SDASSL PASSED ALL TESTS----------')
  IF ( (Ipass==0).AND.(Kprint/=0) ) WRITE (Lun,99009)
  99009 FORMAT (/,' **********SDASSL FAILED SOME TESTS*********')
  99010 FORMAT (//' WARNING.. ERROR EXCEEDS 1000 * TOLERANCE','  WHEN  T =',1P,&
    E13.5//)
  99011 FORMAT (//'TROUBLE..  SDASSL RETURNED  IDID =',I4,'  WHEN  T =',1P,E13.5)
  99012 FORMAT (//1X,' FINAL STATISTICS FOR THIS RUN..',/1X,' NUMBER OF STEPS =',&
    I5/1X,' NUMBER OF F-S   =',I5/1X,' NUMBER OF J-S   =',I5/1X,&
    ' ERROR OVERRUN =',1P,E10.2)
END SUBROUTINE SDASQC
