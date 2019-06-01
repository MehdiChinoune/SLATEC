MODULE TEST49_MOD
  IMPLICIT NONE

CONTAINS
  !** DEDIT2
  SUBROUTINE DEDIT2(Y,T,Erm)
    !>
    !  Subsidiary to DDASQC.
    !***
    ! **Library:**   SLATEC (DASSL)
    !***
    ! **Type:**      DOUBLE PRECISION (EDIT2-S, DEDIT2-D)
    !***
    ! **Author:**  PETZOLD, LINDA R., (LLNL)
    !***
    ! **See also:**  DDASQC
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   891013  DATE WRITTEN
    !   901001  Converted prologue to 4.0 format and made all argument
    !           declarations explicit.  (FNF)
    !   901009  Changed AMAX1 to MAX.  (FNF)
    !   901030  Removed FLOAT's; made all local declarations explicit. (FNF)

    REAL(8) :: Y(*), T, Erm
    INTEGER i, j, k
    REAL(8) :: a1, a2, er, ex, yt
    REAL(8), PARAMETER :: alph1 = 1.0D0, alph2 = 1.0D0
    INTEGER, PARAMETER :: ng = 5
    !* FIRST EXECUTABLE STATEMENT  DEDIT2
    Erm = 0.0D0
    IF ( T==0.0D0 ) RETURN
    ex = 0.0D0
    IF ( T<=30.0D0 ) ex = EXP(-2.0D0*T)
    a2 = 1.0D0
    DO j = 1, ng
      a1 = 1.0D0
      DO i = 1, ng
        k = i + (j-1)*ng
        yt = T**(i+j-2)*ex*a1*a2
        er = ABS(Y(k)-yt)
        Erm = MAX(Erm,er)
        a1 = a1*alph1/i
      END DO
      a2 = a2*alph2/j
    END DO
  END SUBROUTINE DEDIT2
  !** DDASQC
  SUBROUTINE DDASQC(Lun,Kprint,Ipass)
    !>
    !  Quick check for SLATEC routine DDASSL.
    !***
    ! **Library:**   SLATEC (DASSL)
    !***
    ! **Category:**  I1A2
    !***
    ! **Type:**      DOUBLE PRECISION (SDASQC-S, DDASQC-D)
    !***
    ! **Keywords:**  DDASSL, QUICK CHECK
    !***
    ! **Author:**  PETZOLD, LINDA R., (LLNL)
    !             COMPUTING AND MATHEMATICS RESEARCH DIVISION
    !             LAWRENCE LIVERMORE NATIONAL LABORATORY
    !             L - 316, P.O. BOX 808,
    !             LIVERMORE, CA.    94550
    !***
    ! **Description:**
    !       Demonstration program for DDASSL.
    !
    !       DDASSL is used to solve two simple problems,
    !       one with a full Jacobian, the other with a banded Jacobian,
    !       with the 2 appropriate values of info(5) in each case.
    !       If the errors are too large, or other difficulty occurs,
    !       a warning message is printed.  All output is on unit LUN.
    !
    !       To run the demonstration problems with full print, call
    !       DDASQC with KPRINT = 3.
    !
    !***
    ! **References:**  (NONE)
    !***
    ! **Routines called:**  DDASSL, DDJAC1, DDJAC2, DDRES1, DDRES2, DEDIT2

    !* REVISION HISTORY  (YYMMDD)
    !   851113  DATE WRITTEN
    !   880608  Revised to meet new LDOC standards.
    !   880615  Revised to meet new prologue standards.
    !   891016  Converted to 4.0 format (FNF).
    !   901001  Minor improvements to prologue.  (FNF)
    !   901003  Added some error prints and improved output formats.  (FNF)
    !   901009  Corrected GAMS classification code.  (FNF)
    !   901009  Changed AMAX1 to MAX.  (FNF)
    !   901009  Constructed double precision version.  (FNF)
    !   901030  Made all declarations explicit; added 1P's to formats. (FNF)
    USE slatec, ONLY : DDASSL
    !
    INTEGER Lun, Kprint, Ipass
    !
    INTEGER i, idid, info(15), iout, ires, iwork(45), j190, &
      j290, liw, lrw, ml, mu, neq, nerr, nfe, nje, nout, nqu, nst
    REAL(8) :: atol(1), delta(25), er, er1, er2, erm, ero, hu, rtol(1), &
      rwork(550), t, tout, y(25), yprime(25), yt1, yt2
    !
    REAL(8), PARAMETER :: tout1 = 1.0D0, dtout = 1.0D0
    !
    !* FIRST EXECUTABLE STATEMENT  DDASQC
    Ipass = 1
    nerr = 0
    rtol = 0.0D0
    atol = 1.0D-3
    lrw = 550
    liw = 45
    !
    ! FIRST PROBLEM
    !
    neq = 2
    nout = 10
    IF ( Kprint>=2 ) WRITE (Lun,99001) neq, rtol, atol
    99001 FORMAT ('1'/1X,' DEMONSTRATION PROGRAM FOR DDASSL',///1X,&
      ' PROBLEM 1..   LINEAR DIFFERENTIAL/ALGEBRAIC SYSTEM..',/1X,&
      ' X1DOT + 10.0*X1 = 0,  X1 + X2 = 1',/1X,&
      ' X1(0) = 1.0, X2(0) = 0.0',/1X,' NEQ =',I2/1X,' RTOL =',1P,D10.1,&
      '   ATOL =',D10.1)
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
      t = 0.0D0
      y(1) = 1.0D0
      y(2) = 0.0D0
      yprime(1) = -10.0D0
      yprime(2) = 10.0D0
      tout = tout1
      ero = 0.0D0
      DO iout = 1, nout
        CALL DDASSL(DDRES1,neq,t,y,yprime,tout,info,rtol,atol,idid,rwork,lrw,&
          iwork,liw,DDJAC1)
        hu = rwork(7)
        nqu = iwork(8)
        IF ( Kprint>2 ) THEN
          WRITE (Lun,99003) t, y(1), y(2), nqu, hu
          99003 FORMAT (1X,1P,D15.5,D16.5,D16.5,I6,D14.3)
        END IF
        !
        IF ( idid<0 ) EXIT
        yt1 = EXP(-10.0D0*t)
        yt2 = 1.0D0 - yt1
        er1 = ABS(yt1-y(1))
        er2 = ABS(yt2-y(2))
        er = MAX(er1,er2)/atol(1)
        ero = MAX(ero,er)
        IF ( er>1000.0D0 ) THEN
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
    99004 FORMAT ('1'/1X,' DEMONSTRATION PROGRAM FOR DDASSL',///1X,&
      ' PROBLEM 2.. YDOT = A * Y, WHERE ',&
      ' A IS A BANDED LOWER TRIANGULAR MATRIX',/1X,&
      '  DERIVED FROM 2-D ADVECTION PDE',/1X,' NEQ =',I3,'   ML =',I2,&
      '   MU =',I2/1X,' RTOL =',1P,D10.1,'   ATOL =',D10.1)
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
      t = 0.0D0
      DO i = 2, neq
        y(i) = 0.0D0
      END DO
      y(1) = 1.0D0
      DO i = 1, neq
        delta(i) = 0.0D0
      END DO
      !        Following is to initialize YPRIME.
      CALL DDRES2(t,y,delta,yprime,ires)
      tout = 0.01D0
      ero = 0.0D0
      DO iout = 1, nout
        CALL DDASSL(DDRES2,neq,t,y,yprime,tout,info,rtol,atol,idid,rwork,lrw,&
          iwork,liw,DDJAC2)
        CALL DEDIT2(y,t,erm)
        hu = rwork(7)
        nqu = iwork(8)
        IF ( Kprint>2 ) WRITE (Lun,99006) t, erm, nqu, hu
        99006 FORMAT (1X,1P,D15.5,D14.3,I6,D14.3)
        !
        IF ( idid<0 ) EXIT
        er = erm/atol(1)
        ero = MAX(ero,er)
        IF ( er>1000.0D0 ) THEN
          IF ( Kprint>=2 ) WRITE (Lun,99010) t
          nerr = nerr + 1
        END IF
        tout = tout*10.0D0
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
    99008 FORMAT (/,' ----------DDASSL PASSED ALL TESTS----------')
    IF ( (Ipass==0).AND.(Kprint/=0) ) WRITE (Lun,99009)
    99009 FORMAT (/,' **********DDASSL FAILED SOME TESTS*********')
    99010 FORMAT (//' WARNING.. ERROR EXCEEDS 1000 * TOLERANCE','  WHEN  T =',1P,&
      D13.5//)
    99011 FORMAT (//'TROUBLE..  DDASSL RETURNED  IDID =',I4,'  WHEN  T =',1P,D13.5)
    99012 FORMAT (//1X,' FINAL STATISTICS FOR THIS RUN..',/1X,' NUMBER OF STEPS =',&
      I5/1X,' NUMBER OF F-S   =',I5/1X,' NUMBER OF J-S   =',I5/1X,&
      ' ERROR OVERRUN =',1P,D10.2)
  END SUBROUTINE DDASQC
  !** DDJAC1
  SUBROUTINE DDJAC1(T,Y,Yprime,Pd,Cj)
    !>
    !  First Jacobian evaluator for DDASQC.
    !***
    ! **Library:**   SLATEC (DASSL)
    !***
    ! **Type:**      DOUBLE PRECISION (SDJAC1-S, DDJAC1-D)
    !***
    ! **Author:**  PETZOLD, LINDA R., (LLNL)
    !***
    ! **See also:**  DDASQC
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   891013  DATE WRITTEN
    !   901001  Converted prologue to 4.0 format and made all argument
    !           declarations explicit.  (FNF)

    REAL(8) :: T, Y(:), Yprime(:), Pd(:,:), Cj
    !* FIRST EXECUTABLE STATEMENT  DDJAC1
    Pd(1,1) = Cj + 10.0D0
    Pd(1,2) = 0.0D0
    Pd(2,1) = 1.0D0
    Pd(2,2) = 1.0D0
  END SUBROUTINE DDJAC1
  !** DDJAC2
  SUBROUTINE DDJAC2(T,Y,Yprime,Pd,Cj)
    !>
    !  Second Jacobian evaluator for DDASQC.
    !***
    ! **Library:**   SLATEC (DASSL)
    !***
    ! **Type:**      DOUBLE PRECISION (SDJAC2-S, DDJAC2-D)
    !***
    ! **Author:**  PETZOLD, LINDA R., (LLNL)
    !***
    ! **See also:**  DDASQC
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   891013  DATE WRITTEN
    !   901001  Converted prologue to 4.0 format and made all argument
    !           declarations explicit.  (FNF)
    !   901001  Eliminated 7-character variable names MBANDPn by explicitly
    !           including MBAND+n in expressions.  (FNF)
    !   901030  Made all local declarations explicit.  (FNF)

    REAL(8) :: T, Y(:), Yprime(:), Pd(:,:), Cj
    INTEGER :: j, mband
    REAL(8), PARAMETER :: alph1 = 1.0D0, alph2 = 1.0D0
    INTEGER, PARAMETER :: ng = 5
    INTEGER, PARAMETER :: ml = 5, mu = 0, neq = 25
    !* FIRST EXECUTABLE STATEMENT  DDJAC2
    mband = ml + mu + 1
    DO j = 1, neq
      Pd(mband,j) = -2.0D0 - Cj
      Pd(mband+1,j) = alph1
      Pd(mband+2,j) = 0.0D0
      Pd(mband+3,j) = 0.0D0
      Pd(mband+4,j) = 0.0D0
      Pd(mband+5,j) = alph2
    END DO
    DO j = 1, neq, ng
      Pd(mband+1,j) = 0.0D0
    END DO
  END SUBROUTINE DDJAC2
  !** DDRES1
  SUBROUTINE DDRES1(T,Y,Yprime,Delta,Ires)
    !>
    !  First residual evaluator for DDASQC.
    !***
    ! **Library:**   SLATEC (DASSL)
    !***
    ! **Type:**      DOUBLE PRECISION (SDRES1-S, DDRES1-D)
    !***
    ! **Author:**  PETZOLD, LINDA R., (LLNL)
    !***
    ! **See also:**  DDASQC
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   891013  DATE WRITTEN
    !   901001  Converted prologue to 4.0 format and made all argument
    !           declarations explicit.  (FNF)

    INTEGER :: Ires
    REAL(8) :: T, Y(:), Yprime(:), Delta(:)
    !* FIRST EXECUTABLE STATEMENT  DDRES1
    Delta(1) = Yprime(1) + 10.0D0*Y(1)
    Delta(2) = Y(2) + Y(1) - 1.0D0
  END SUBROUTINE DDRES1
  !** DDRES2
  SUBROUTINE DDRES2(T,Y,Yprime,Delta,Ires)
    !>
    !  Second residual evaluator for DDASQC.
    !***
    ! **Library:**   SLATEC (DASSL)
    !***
    ! **Type:**      DOUBLE PRECISION (SDRES2-S, DDRES2-D)
    !***
    ! **Author:**  PETZOLD, LINDA R., (LLNL)
    !***
    ! **See also:**  DDASQC
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   891013  DATE WRITTEN
    !   901001  Converted prologue to 4.0 format and made all argument
    !           declarations explicit.  (FNF)
    !   901030  Made all local declarations explicit.  (FNF)

    INTEGER :: Ires
    REAL(8) :: T, Y(:), Yprime(:), Delta(:)
    INTEGER :: i, j, k
    REAL(8) :: d
    REAL(8), PARAMETER :: alph1 = 1.0D0, alph2 = 1.0D0
    INTEGER, PARAMETER :: ng = 5
    !* FIRST EXECUTABLE STATEMENT  DDRES2
    DO j = 1, ng
      DO i = 1, ng
        k = i + (j-1)*ng
        d = -2.0D0*Y(k)
        IF ( i/=1 ) d = d + Y(k-1)*alph1
        IF ( j/=1 ) d = d + Y(k-ng)*alph2
        Delta(k) = d - Yprime(k)
      END DO
    END DO
  END SUBROUTINE DDRES2
END MODULE TEST49_MOD
!** TEST49
PROGRAM TEST49
  USE TEST49_MOD, ONLY : DDASQC
  USE slatec, ONLY : I1MACH, XSETF, XSETUN, XERMAX
  USE common_mod, ONLY : GET_ARGUMENT
  IMPLICIT NONE
  !>
  !  Driver for testing SLATEC subprograms
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  I1A2
  !***
  ! **Type:**      DOUBLE PRECISION (TEST48-S, TEST49-D)
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
  !       DDASSL
  !
  !***
  ! **References:**  Fong, Kirby W., Jefferson, Thomas H., Suyehiro,
  !                 Tokihiko, Walton, Lee, Guidelines to the SLATEC Common
  !                 Mathematical Library, March 21, 1989.
  !***
  ! **Routines called:**  DDASQC, I1MACH, XERMAX, XSETF, XSETUN

  !* REVISION HISTORY  (YYMMDD)
  !   891013  DATE WRITTEN
  !   901001  Converted prologue to 4.0 format.  (FNF)
  !   901009  Corrected GAMS classification code.  (FNF)
  !   901009  Constructed double precision version.  (FNF)
  INTEGER ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST49
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
  !     Test DDASSL
  !
  CALL DDASQC(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST49 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST49 *************')
  END IF
  STOP
END PROGRAM TEST49
