MODULE TEST48_MOD
  IMPLICIT NONE

CONTAINS
  !** EDIT2
  SUBROUTINE EDIT2(Y,T,Erm)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to SDASQC.
    !***
    ! **Library:**   SLATEC (DASSL)
    !***
    ! **Type:**      SINGLE PRECISION (EDIT2-S, DEDIT2-D)
    !***
    ! **Author:**  PETZOLD, LINDA R., (LLNL)
    !***
    ! **See also:**  SDASQC
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   891013  DATE WRITTEN
    !   901001  Converted prologue to 4.0 format and made all argument
    !           declarations explicit.  (FNF)
    !   901009  Changed AMAX1 to MAX.  (FNF)
    !   901030  Removed FLOAT's; made all local declarations explicit. (FNF)

    REAL Y(*), T, Erm
    INTEGER i, j, k
    REAL a1, a2, er, ex, yt
    REAL, PARAMETER :: alph1 = 1.0E0, alph2 = 1.0E0
    INTEGER, PARAMETER :: ng = 5
    !* FIRST EXECUTABLE STATEMENT  EDIT2
    Erm = 0.0E0
    IF ( T==0.0E0 ) RETURN
    ex = 0.0E0
    IF ( T<=30.0E0 ) ex = EXP(-2.0E0*T)
    a2 = 1.0E0
    DO j = 1, ng
      a1 = 1.0E0
      DO i = 1, ng
        k = i + (j-1)*ng
        yt = T**(i+j-2)*ex*a1*a2
        er = ABS(Y(k)-yt)
        Erm = MAX(Erm,er)
        a1 = a1*alph1/i
      END DO
      a2 = a2*alph2/j
    END DO
  END SUBROUTINE EDIT2
  !** SDASQC
  SUBROUTINE SDASQC(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for SLATEC routine SDASSL.
    !***
    ! **Library:**   SLATEC (DASSL)
    !***
    ! **Category:**  I1A2
    !***
    ! **Type:**      SINGLE PRECISION (SDASQC-S, DDASQC-D)
    !***
    ! **Keywords:**  QUICK CHECK, SDASSL
    !***
    ! **Author:**  PETZOLD, LINDA R., (LLNL)
    !             COMPUTING AND MATHEMATICS RESEARCH DIVISION
    !             LAWRENCE LIVERMORE NATIONAL LABORATORY
    !             L - 316, P.O. BOX 808,
    !             LIVERMORE, CA.    94550
    !***
    ! **Description:**
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
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  EDIT2, SDASSL, SDJAC1, SDJAC2, SDRES1, SDRES2

    !* REVISION HISTORY  (YYMMDD)
    !   851113  DATE WRITTEN
    !   880608  Revised to meet new LDOC standards.
    !   880615  Revised to meet new prologue standards.
    !   891016  Converted to 4.0 format (FNF).
    !   901001  Minor improvements to prologue.  (FNF)
    !   901003  Added some error prints and improved output formats.  (FNF)
    !   901009  Corrected GAMS classification code.  (FNF)
    !   901009  Changed AMAX1 to MAX.  (FNF)
    !   901030  Made all declarations explicit; added 1P's to formats. (FNF)

    !
    INTEGER Lun, Kprint, Ipass
    !
    INTEGER i, idid, info(15), iout, ipar(1), ires, iwork(45), j190, &
      j290, liw, lrw, ml, mu, neq, nerr, nfe, nje, nout, nqu, nst
    REAL atol, delta(25), er, er1, er2, erm, ero, hu, rpar(1), &
      rtol, rwork(550), t, tout, y(25), yprime(25), yt1, yt2
    !
    REAL, PARAMETER :: tout1 = 1.0E0, dtout = 1.0E0
    !
    !* FIRST EXECUTABLE STATEMENT  SDASQC
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
    IF ( Kprint>=2 ) WRITE (Lun,99001) neq, rtol, atol
    99001 FORMAT ('1'/1X,' DEMONSTRATION PROGRAM FOR SDASSL',///1X,&
      ' PROBLEM 1..   LINEAR DIFFERENTIAL/ALGEBRAIC SYSTEM..',/1X,&
      ' X1DOT + 10.0*X1 = 0,  X1 + X2 = 1',/1X,&
      ' X1(0) = 1.0, X2(0) = 0.0',/1X,' NEQ =',I2/1X,' RTOL =',1P,E10.1,&
      '   ATOL =',E10.1)
    !
    DO j190 = 1, 2
      DO i = 1, 15
        info(i) = 0
      END DO
      IF ( j190==2 ) info(5) = 1
      !
      IF ( Kprint>2 ) WRITE (Lun,99002) info(5)
      99002 FORMAT (////1X,' INFO(5) =',I3//6X,'T',15X,'X1',14X,'X2',12X,'NQ',6X,&
        'H',12X/)
      !
      t = 0.0E0
      y(1) = 1.0E0
      y(2) = 0.0E0
      yprime(1) = -10.0E0
      yprime(2) = 10.0E0
      tout = tout1
      ero = 0.0E0
      DO iout = 1, nout
        CALL SDASSL(SDRES1,neq,t,y,yprime,tout,info,rtol,atol,idid,rwork,lrw,&
          iwork,liw,rpar,ipar,SDJAC1)
        hu = rwork(7)
        nqu = iwork(8)
        IF ( Kprint>2 ) THEN
          WRITE (Lun,99003) t, y(1), y(2), nqu, hu
          99003 FORMAT (1X,1P,E15.5,E16.5,E16.5,I6,E14.3)
        END IF
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
        END IF
        tout = tout + dtout
      END DO
      IF ( idid<0 ) THEN
        IF ( Kprint>=2 ) WRITE (Lun,99011) idid, t
        !
        nerr = nerr + 1
      END IF
      nst = iwork(11)
      nfe = iwork(12)
      nje = iwork(13)
      IF ( Kprint>2 ) WRITE (Lun,99012) nst, nfe, nje, ero
      !
    END DO
    !
    ! SECOND PROBLEM
    !
    neq = 25
    ml = 5
    mu = 0
    iwork(1) = ml
    iwork(2) = mu
    nout = 5
    IF ( Kprint>=2 ) WRITE (Lun,99004) neq, ml, mu, rtol, atol
    99004 FORMAT ('1'/1X,' DEMONSTRATION PROGRAM FOR SDASSL',///1X,&
      ' PROBLEM 2.. YDOT = A * Y, WHERE ',&
      ' A IS A BANDED LOWER TRIANGULAR MATRIX',/1X,&
      '  DERIVED FROM 2-D ADVECTION PDE',/1X,' NEQ =',I3,'   ML =',I2,&
      '   MU =',I2/1X,' RTOL =',1P,E10.1,'   ATOL =',E10.1)
    !
    DO j290 = 1, 2
      DO i = 1, 15
        info(i) = 0
      END DO
      info(6) = 1
      IF ( j290==2 ) info(5) = 1
      !
      IF ( Kprint>2 ) WRITE (Lun,99005) info(5)
      99005 FORMAT (////1X,' INFO(5) =',I3//6X,'T',14X,'MAX.ERR.',5X,'NQ',6X,'H'/)
      !
      t = 0.0E0
      DO i = 2, neq
        y(i) = 0.0E0
      END DO
      y(1) = 1.0E0
      DO i = 1, neq
        delta(i) = 0.0E0
      END DO
      !        Following is to initialize YPRIME.
      CALL SDRES2(t,y,delta,yprime,ires,rpar,ipar)
      tout = 0.01E0
      ero = 0.0E0
      DO iout = 1, nout
        CALL SDASSL(SDRES2,neq,t,y,yprime,tout,info,rtol,atol,idid,rwork,lrw,&
          iwork,liw,rpar,ipar,SDJAC2)
        CALL EDIT2(y,t,erm)
        hu = rwork(7)
        nqu = iwork(8)
        IF ( Kprint>2 ) WRITE (Lun,99006) t, erm, nqu, hu
        99006 FORMAT (1X,1P,E15.5,E14.3,I6,E14.3)
        !
        IF ( idid<0 ) EXIT
        er = erm/atol
        ero = MAX(ero,er)
        IF ( er>1000.0E0 ) THEN
          IF ( Kprint>=2 ) WRITE (Lun,99010) t
          nerr = nerr + 1
        END IF
        tout = tout*10.0E0
      END DO
      IF ( idid<0 ) THEN
        IF ( Kprint>=2 ) WRITE (Lun,99011) idid, t
        nerr = nerr + 1
      END IF
      nst = iwork(11)
      nfe = iwork(12)
      nje = iwork(13)
      IF ( Kprint>2 ) WRITE (Lun,99012) nst, nfe, nje, ero
      !
    END DO
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
  !** SDJAC1
  SUBROUTINE SDJAC1(T,Y,Yprime,Pd,Cj,Rpar,Ipar)
    IMPLICIT NONE
    !>
    !***
    !  First Jacobian evaluator for SDASQC.
    !***
    ! **Library:**   SLATEC (DASSL)
    !***
    ! **Type:**      SINGLE PRECISION (SDJAC1-S, DDJAC1-D)
    !***
    ! **Author:**  PETZOLD, LINDA R., (LLNL)
    !***
    ! **See also:**  SDASQC
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   891013  DATE WRITTEN
    !   901001  Converted prologue to 4.0 format and made all argument
    !           declarations explicit.  (FNF)

    INTEGER Ipar(*)
    REAL T, Y(*), Yprime(*), Pd(2,2), Cj, Rpar(*)
    !* FIRST EXECUTABLE STATEMENT  SDJAC1
    Pd(1,1) = Cj + 10.0E0
    Pd(1,2) = 0.0E0
    Pd(2,1) = 1.0E0
    Pd(2,2) = 1.0E0
  END SUBROUTINE SDJAC1
  !** SDJAC2
  SUBROUTINE SDJAC2(T,Y,Yprime,Pd,Cj,Rpar,Ipar)
    IMPLICIT NONE
    !>
    !***
    !  Second Jacobian evaluator for SDASQC.
    !***
    ! **Library:**   SLATEC (DASSL)
    !***
    ! **Type:**      SINGLE PRECISION (SDJAC2-S, DDJAC2-D)
    !***
    ! **Author:**  PETZOLD, LINDA R., (LLNL)
    !***
    ! **See also:**  SDASQC
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   891013  DATE WRITTEN
    !   901001  Converted prologue to 4.0 format and made all argument
    !           declarations explicit.  (FNF)
    !   901001  Eliminated 7-character variable names MBANDPn by explicitly
    !           including MBAND+n in expressions.  (FNF)
    !   901030  Made all local declarations explicit.  (FNF)

    INTEGER Ipar(*)
    REAL T, Y(*), Yprime(*), Pd(11,25), Cj, Rpar(*)
    INTEGER j, mband
    REAL, PARAMETER :: alph1 = 1.0E0, alph2 = 1.0E0
    INTEGER, PARAMETER :: ng = 5
    INTEGER, PARAMETER :: ml = 5, mu = 0, neq = 25
    !* FIRST EXECUTABLE STATEMENT  SDJAC2
    mband = ml + mu + 1
    DO j = 1, neq
      Pd(mband,j) = -2.0E0 - Cj
      Pd(mband+1,j) = alph1
      Pd(mband+2,j) = 0.0E0
      Pd(mband+3,j) = 0.0E0
      Pd(mband+4,j) = 0.0E0
      Pd(mband+5,j) = alph2
    END DO
    DO j = 1, neq, ng
      Pd(mband+1,j) = 0.0E0
    END DO
  END SUBROUTINE SDJAC2
  !** SDRES1
  SUBROUTINE SDRES1(T,Y,Yprime,Delta,Ires,Rpar,Ipar)
    IMPLICIT NONE
    !>
    !***
    !  First residual evaluator for SDASQC.
    !***
    ! **Library:**   SLATEC (DASSL)
    !***
    ! **Type:**      SINGLE PRECISION (SDRES1-S, DDRES1-D)
    !***
    ! **Author:**  PETZOLD, LINDA R., (LLNL)
    !***
    ! **See also:**  SDASQC
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   891013  DATE WRITTEN
    !   901001  Converted prologue to 4.0 format and made all argument
    !           declarations explicit.  (FNF)

    INTEGER Ires, Ipar(*)
    REAL T, Y(*), Yprime(*), Delta(*), Rpar(*)
    !* FIRST EXECUTABLE STATEMENT  SDRES1
    Delta(1) = Yprime(1) + 10.0E0*Y(1)
    Delta(2) = Y(2) + Y(1) - 1.0E0
  END SUBROUTINE SDRES1
  !** SDRES2
  SUBROUTINE SDRES2(T,Y,Yprime,Delta,Ires,Rpar,Ipar)
    IMPLICIT NONE
    !>
    !***
    !  Second residual evaluator for SDASQC.
    !***
    ! **Library:**   SLATEC (DASSL)
    !***
    ! **Type:**      SINGLE PRECISION (SDRES2-S, DDRES2-D)
    !***
    ! **Author:**  PETZOLD, LINDA R., (LLNL)
    !***
    ! **See also:**  SDASQC
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   891013  DATE WRITTEN
    !   901001  Converted prologue to 4.0 format and made all argument
    !           declarations explicit.  (FNF)
    !   901030  Made all local declarations explicit.  (FNF)

    INTEGER Ires, Ipar(*)
    REAL T, Y(*), Yprime(*), Delta(*), Rpar(*)
    INTEGER i, j, k
    REAL d
    REAL, PARAMETER :: alph1 = 1.0E0, alph2 = 1.0E0
    INTEGER, PARAMETER :: ng = 5
    !* FIRST EXECUTABLE STATEMENT  SDRES2
    DO j = 1, ng
      DO i = 1, ng
        k = i + (j-1)*ng
        d = -2.0E0*Y(k)
        IF ( i/=1 ) d = d + Y(k-1)*alph1
        IF ( j/=1 ) d = d + Y(k-ng)*alph2
        Delta(k) = d - Yprime(k)
      END DO
    END DO
  END SUBROUTINE SDRES2
END MODULE TEST48_MOD
!** TEST48
PROGRAM TEST48
  USE TEST48_MOD
  IMPLICIT NONE
  !>
  !***
  !  Driver for testing SLATEC subprograms
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  I1A2
  !***
  ! **Type:**      SINGLE PRECISION (TEST48-S, TEST49-D)
  !***
  ! **Keywords:**  QUICK CHECK DRIVER
  !***
  ! **Author:**  SLATEC Common Mathematical Library Committee
  !***
  ! **Description:**
  !
  !- Usage:
  !     One input data record is required
  !         READ (LIN, '(I1)') KPRINT
  !
  !- Arguments:
  !     KPRINT = 0  Quick checks - No printing.
  !                 Driver       - Short pass or fail message printed.
  !              1  Quick checks - No message printed for passed tests,
  !                                short message printed for failed tests.
  !                 Driver       - Short pass or fail message printed.
  !              2  Quick checks - Print short message for passed tests,
  !                                fuller information for failed tests.
  !                 Driver       - Pass or fail message printed.
  !              3  Quick checks - Print complete quick check results.
  !                 Driver       - Pass or fail message printed.
  !
  !- Description:
  !     Driver for testing SLATEC subprograms
  !       SDASSL
  !
  !***
  ! **References:**  Fong, Kirby W., Jefferson, Thomas H., Suyehiro,
  !                 Tokihiko, Walton, Lee, Guidelines to the SLATEC Common
  !                 Mathematical Library, March 21, 1989.
  !***
  ! **Routines called:**  I1MACH, SDASQC, XERMAX, XSETF, XSETUN

  !* REVISION HISTORY  (YYMMDD)
  !   891013  DATE WRITTEN
  !   901001  Converted prologue to 4.0 format.  (FNF)
  !   901009  Corrected GAMS classification code.  (FNF)

  INTEGER I1MACH
  INTEGER ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST48
  lun = I1MACH(2)
  lin = I1MACH(1)
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  CALL GET_ARGUMENT(kprint)
  CALL XSETUN(lun)
  IF ( kprint<=1 ) THEN
    CALL XSETF(0)
  ELSE
    CALL XSETF(1)
  END IF
  CALL XERMAX(1000)
  !
  !     Test SDASSL
  !
  CALL SDASQC(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST48 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST48 *************')
  END IF
  STOP
END PROGRAM TEST48
