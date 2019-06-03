MODULE TEST13_MOD
  USE service, ONLY : SP, DP
  IMPLICIT NONE

CONTAINS
  !** QCRC
  SUBROUTINE QCRC(Lun,Kprint,Ipass)
    !>
    !  Quick check for RC.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Pexton, R. L., (LLNL)
    !***
    ! **Description:**
    !
    !            QUICK TEST FOR CARLSON INTEGRAL RC
    !
    !***
    ! **Routines called:**  NUMXER, R1MACH, RC, XERCLR, XGETF, XSETF

    !* REVISION HISTORY  (YYMMDD)
    !   790801  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   910708  Minor modifications in use of KPRINT.  (WRB)
    USE slatec, ONLY : NUMXER, R1MACH, RC, XERCLR, XGETF, XSETF
    INTEGER Kprint, Ipass, contrl, kontrl, Lun, ier
    INTEGER ipass1, ipass2, ipass3, ipass4
    REAL(SP) pi, trc, dif
    !* FIRST EXECUTABLE STATEMENT  QCRC
    CALL XERCLR
    CALL XGETF(contrl)
    IF ( Kprint>=3 ) THEN
      kontrl = +1
    ELSE
      kontrl = 0
    END IF
    CALL XSETF(kontrl)
    !
    !  FORCE ERROR 1
    !
    IF ( Kprint>=3 ) WRITE (Lun,99001)
    99001 FORMAT (' RC - FORCE ERROR 1 TO OCCUR')
    trc = RC(-1.0E0,-1.0E0,ier)
    ier = NUMXER(ier)
    IF ( ier==1 ) THEN
      ipass1 = 1
    ELSE
      ipass1 = 0
    END IF
    CALL XERCLR
    !
    !  FORCE ERROR 2
    !
    IF ( Kprint>=3 ) WRITE (Lun,99002)
    99002 FORMAT (' RC - FORCE ERROR 2 TO OCCUR')
    trc = RC(R1MACH(1),R1MACH(1),ier)
    ier = NUMXER(ier)
    IF ( ier==2 ) THEN
      ipass2 = 1
    ELSE
      ipass2 = 0
    END IF
    CALL XERCLR
    !
    !  FORCE ERROR 3
    !
    IF ( Kprint>=3 ) WRITE (Lun,99003)
    99003 FORMAT (' RC - FORCE ERROR 3 TO OCCUR')
    trc = RC(R1MACH(2),R1MACH(2),ier)
    ier = NUMXER(ier)
    IF ( ier==3 ) THEN
      ipass3 = 1
    ELSE
      ipass3 = 0
    END IF
    CALL XERCLR
    !
    !  ARGUMENTS IN RANGE
    !
    pi = 3.1415926535897932E0
    trc = RC(0.0E0,0.25E0,ier)
    CALL XERCLR
    dif = trc - pi
    IF ( (ABS(dif/pi)<1000.0E0*R1MACH(4)).AND.(ier==0) ) THEN
      ipass4 = 1
    ELSE
      ipass4 = 0
    END IF
    Ipass = MIN(ipass1,ipass2,ipass3,ipass4)
    IF ( Kprint<=0 ) THEN
    ELSEIF ( Kprint==1 ) THEN
      IF ( Ipass/=1 ) WRITE (Lun,99006)
    ELSEIF ( Ipass==1 ) THEN
      WRITE (Lun,99004)
      99004 FORMAT (' RC - PASSED')
    ELSE
      WRITE (Lun,99006)
      IF ( ipass4==0 ) WRITE (Lun,99005) pi, trc, dif
      99005 FORMAT (' CORRECT ANSWER =',1PE14.6/'COMPUTED ANSWER =',&
        E14.6/'     DIFFERENCE =',E14.6)
    END IF
    CALL XSETF(contrl)
    99006 FORMAT (' RC - FAILED')
  END SUBROUTINE QCRC
  !** QCRD
  SUBROUTINE QCRD(Lun,Kprint,Ipass)
    !>
    !  Quick check for RD.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Pexton, R. L., (LLNL)
    !***
    ! **Description:**
    !
    !            QUICK TEST FOR CARLSON INTEGRAL RD
    !
    !***
    ! **Routines called:**  NUMXER, R1MACH, RD, XERCLR, XGETF, XSETF

    !* REVISION HISTORY  (YYMMDD)
    !   790801  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   910708  Minor modifications in use of KPRINT.  (WRB)
    USE slatec, ONLY : NUMXER, R1MACH, RD, XERCLR, XGETF, XSETF
    INTEGER Kprint, Ipass, contrl, kontrl, Lun, ier
    INTEGER ipass1, ipass2, ipass3, ipass4
    REAL(SP) blem, trd, dif
    !* FIRST EXECUTABLE STATEMENT  QCRD
    CALL XERCLR
    CALL XGETF(contrl)
    IF ( Kprint>=3 ) THEN
      kontrl = +1
    ELSE
      kontrl = 0
    END IF
    CALL XSETF(kontrl)
    !
    !  FORCE ERROR 1
    !
    IF ( Kprint>=3 ) WRITE (Lun,99001)
    99001 FORMAT (' RD - FORCE ERROR 1 TO OCCUR')
    trd = RD(-1.0E0,-1.0E0,-1.0E0,ier)
    ier = NUMXER(ier)
    IF ( ier==1 ) THEN
      ipass1 = 1
    ELSE
      ipass1 = 0
    END IF
    CALL XERCLR
    !
    !  FORCE ERROR 2
    !
    IF ( Kprint>=3 ) WRITE (Lun,99002)
    99002 FORMAT (' RD - FORCE ERROR 2 TO OCCUR')
    trd = RD(1.0E0,1.0E0,-1.0E0,ier)
    ier = NUMXER(ier)
    IF ( ier==2 ) THEN
      ipass2 = 1
    ELSE
      ipass2 = 0
    END IF
    CALL XERCLR
    !
    !  FORCE ERROR 3
    !
    IF ( Kprint>=3 ) WRITE (Lun,99003)
    99003 FORMAT (' RD - FORCE ERROR 3 TO OCCUR')
    trd = RD(R1MACH(2),R1MACH(2),R1MACH(2),ier)
    ier = NUMXER(ier)
    IF ( ier==3 ) THEN
      ipass3 = 1
    ELSE
      ipass3 = 0
    END IF
    CALL XERCLR
    !
    !  ARGUMENTS IN RANGE
    !  BLEM=3 * LEMNISCATE CONSTANT B
    !
    blem = 1.79721035210338831E0
    trd = RD(0.0E0,2.0E0,1.0E0,ier)
    CALL XERCLR
    dif = trd - blem
    IF ( (ABS(dif/blem)<1000.0E0*R1MACH(4)).AND.(ier==0) ) THEN
      ipass4 = 1
    ELSE
      Ipass = 0
    END IF
    Ipass = MIN(ipass1,ipass2,ipass3,ipass4)
    IF ( Kprint<=0 ) THEN
    ELSEIF ( Kprint==1 ) THEN
      IF ( Ipass/=1 ) WRITE (Lun,99006)
    ELSEIF ( Ipass==1 ) THEN
      WRITE (Lun,99004)
      99004 FORMAT (' RD - PASSED')
    ELSE
      WRITE (Lun,99006)
      IF ( ipass4==0 ) WRITE (Lun,99005) blem, trd, dif
      99005 FORMAT (' CORRECT ANSWER =',1PE14.6/'COMPUTED ANSWER =',&
        E14.6/'     DIFFERENCE =',E14.6)
    END IF
    CALL XSETF(contrl)
    99006 FORMAT (' RD - FAILED')
  END SUBROUTINE QCRD
  !** QCRF
  SUBROUTINE QCRF(Lun,Kprint,Ipass)
    !>
    !  Quick check for RF.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Pexton, R. L., (LLNL)
    !***
    ! **Description:**
    !
    !            QUICK TEST FOR CARLSON INTEGRAL RF
    !
    !***
    ! **Routines called:**  NUMXER, R1MACH, RF, XERCLR, XGETF, XSETF

    !* REVISION HISTORY  (YYMMDD)
    !   790801  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   910708  Minor modifications in use of KPRINT.  (WRB)
    USE slatec, ONLY : NUMXER, R1MACH, RF, XERCLR, XGETF, XSETF
    INTEGER Kprint, Ipass, contrl, kontrl, Lun, ier
    INTEGER ipass1, ipass2, ipass3, ipass4
    REAL(SP) alem, trf, dif
    !* FIRST EXECUTABLE STATEMENT  QCRF
    CALL XERCLR
    CALL XGETF(contrl)
    IF ( Kprint>=3 ) THEN
      kontrl = +1
    ELSE
      kontrl = 0
    END IF
    CALL XSETF(kontrl)
    !
    !  FORCE ERROR 1
    !
    IF ( Kprint>=3 ) WRITE (Lun,99001)
    99001 FORMAT (' RF - FORCE ERROR 1 TO OCCUR')
    trf = RF(-1.0E0,-1.0E0,-1.0E0,ier)
    ier = NUMXER(ier)
    IF ( ier==1 ) THEN
      ipass1 = 1
    ELSE
      ipass1 = 0
    END IF
    CALL XERCLR
    !
    !  FORCE ERROR 2
    !
    IF ( Kprint>=3 ) WRITE (Lun,99002)
    99002 FORMAT (' RF - FORCE ERROR 2 TO OCCUR')
    trf = RF(R1MACH(1),R1MACH(1),R1MACH(1),ier)
    ier = NUMXER(ier)
    IF ( ier==2 ) THEN
      ipass2 = 1
    ELSE
      ipass2 = 0
    END IF
    CALL XERCLR
    !
    !  FORCE ERROR 3
    !
    IF ( Kprint>=3 ) WRITE (Lun,99003)
    99003 FORMAT (' RF - FORCE ERROR 3 TO OCCUR')
    trf = RF(R1MACH(2),R1MACH(2),R1MACH(2),ier)
    ier = NUMXER(ier)
    IF ( ier==3 ) THEN
      ipass3 = 1
    ELSE
      ipass3 = 0
    END IF
    CALL XERCLR
    !
    !  ARGUMENTS IN RANGE
    !  ALEM=LEMNISCATE CONSTANT A
    !
    alem = 1.311028777146059905E0
    trf = RF(0.0E0,1.0E0,2.0E0,ier)
    CALL XERCLR
    dif = trf - alem
    IF ( (ABS(dif/alem)<1000.0E0*R1MACH(4)).AND.(ier==0) ) THEN
      ipass4 = 1
    ELSE
      ipass4 = 0
    END IF
    Ipass = MIN(ipass1,ipass2,ipass3,ipass4)
    IF ( Kprint==0 ) THEN
    ELSEIF ( Kprint==1 ) THEN
      IF ( Ipass/=1 ) WRITE (Lun,99006)
    ELSEIF ( Ipass==1 ) THEN
      WRITE (Lun,99004)
      99004 FORMAT (' RF - PASSED')
    ELSE
      WRITE (Lun,99006)
      IF ( ipass4==0 ) WRITE (Lun,99005) alem, trf, dif
      99005 FORMAT (' CORRECT ANSWER =',1PE14.6/'COMPUTED ANSWER =',&
        E14.6/'     DIFFERENCE =',E14.6)
    END IF
    CALL XSETF(contrl)
    99006 FORMAT (' RF - FAILED')
  END SUBROUTINE QCRF
  !** QCRJ
  SUBROUTINE QCRJ(Lun,Kprint,Ipass)
    !>
    !  Quick check for RJ.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Pexton, R. L., (LLNL)
    !***
    ! **Description:**
    !
    !            QUICK TEST FOR CARLSON INTEGRAL RJ
    !
    !***
    ! **Routines called:**  NUMXER, R1MACH, RJ, XERCLR, XGETF, XSETF

    !* REVISION HISTORY  (YYMMDD)
    !   790801  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   910708  Minor modifications in use of KPRINT.  (WRB)
    USE slatec, ONLY : NUMXER, R1MACH, RJ, XERCLR, XGETF, XSETF
    INTEGER Kprint, Ipass, contrl, kontrl, Lun, ier
    INTEGER ipass1, ipass2, ipass3, ipass4
    REAL(SP) consj, trj, dif
    !* FIRST EXECUTABLE STATEMENT  QCRJ
    CALL XERCLR
    CALL XGETF(contrl)
    IF ( Kprint>=3 ) THEN
      kontrl = +1
    ELSE
      kontrl = 0
    END IF
    CALL XSETF(kontrl)
    !
    !  FORCE ERROR 1
    !
    IF ( Kprint>=3 ) WRITE (Lun,99001)
    99001 FORMAT (' RJ - FORCE ERROR 1 TO OCCUR')
    trj = RJ(-1.0E0,-1.0E0,-1.0E0,-1.0E0,ier)
    ier = NUMXER(ier)
    IF ( ier==1 ) THEN
      ipass1 = 1
    ELSE
      ipass1 = 0
    END IF
    CALL XERCLR
    !
    !  FORCE ERROR 2
    !
    IF ( Kprint>=3 ) WRITE (Lun,99002)
    99002 FORMAT (' RJ - FORCE ERROR 2 TO OCCUR')
    trj = RJ(R1MACH(1),R1MACH(1),R1MACH(1),R1MACH(1),ier)
    ier = NUMXER(ier)
    IF ( ier==2 ) THEN
      ipass2 = 1
    ELSE
      ipass2 = 0
    END IF
    CALL XERCLR
    !
    !  FORCE ERROR 3
    !
    IF ( Kprint>=3 ) WRITE (Lun,99003)
    99003 FORMAT (' RJ - FORCE ERROR 3 TO OCCUR')
    trj = RJ(R1MACH(2),R1MACH(2),R1MACH(2),R1MACH(2),ier)
    ier = NUMXER(ier)
    IF ( ier==3 ) THEN
      ipass3 = 1
    ELSE
      ipass3 = 0
    END IF
    CALL XERCLR
    !
    !  ARGUMENTS IN RANGE
    !
    consj = 0.142975796671567538E0
    trj = RJ(2.0E0,3.0E0,4.0E0,5.0E0,ier)
    CALL XERCLR
    dif = trj - consj
    IF ( (ABS(dif/consj)<1000.0E0*R1MACH(4)).AND.(ier==0) ) THEN
      ipass4 = 1
    ELSE
      ipass4 = 0
    END IF
    Ipass = MIN(ipass1,ipass2,ipass3,ipass4)
    IF ( Kprint<=0 ) THEN
    ELSEIF ( Kprint==1 ) THEN
      IF ( Ipass/=1 ) WRITE (Lun,99006)
    ELSEIF ( Ipass==1 ) THEN
      WRITE (Lun,99004)
      99004 FORMAT (' RJ - PASSED')
    ELSE
      WRITE (Lun,99006)
      IF ( ipass4==0 ) WRITE (Lun,99005) consj, trj, dif
      99005 FORMAT (' CORRECT ANSWER =',1PE14.6/'COMPUTED ANSWER =',&
        E14.6/'     DIFFERENCE =',E14.6)
    END IF
    CALL XSETF(contrl)
    99006 FORMAT (' RJ - FAILED')
  END SUBROUTINE QCRJ
END MODULE TEST13_MOD
!** TEST13
PROGRAM TEST13
  USE TEST13_MOD, ONLY : QCRC, QCRD, QCRF, QCRJ
  USE slatec, ONLY : I1MACH, XSETF, XSETUN, XERMAX
  USE common_mod, ONLY : GET_ARGUMENT
  IMPLICIT NONE
  !>
  !  Driver for testing SLATEC subprograms
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  C14
  !***
  ! **Type:**      SINGLE PRECISION (TEST13-S, TEST14-D)
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
  !        RC       RD       RF       RJ
  !
  !***
  ! **References:**  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  I1MACH, QCRC, QCRD, QCRF, QCRJ, XERMAX, XSETF,
  !                    XSETUN

  !* REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)
  INTEGER ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST13
  lun = I1MACH(2)
  lin = I1MACH(1)
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  CALL GET_ARGUMENT(kprint)
  CALL XERMAX(1000)
  CALL XSETUN(lun)
  IF ( kprint<=1 ) THEN
    CALL XSETF(0)
  ELSE
    CALL XSETF(1)
  END IF
  !
  !     Test single precision Carlson elliptic routines
  !
  CALL QCRC(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  CALL QCRD(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  CALL QCRF(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  CALL QCRJ(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST13 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST13  *************')
  END IF
  STOP
END PROGRAM TEST13
