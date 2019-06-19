MODULE TEST14_MOD
  USE service, ONLY : SP, DP
  IMPLICIT NONE

CONTAINS
  !** QCDRC
  SUBROUTINE QCDRC(Lun,Kprint,Ipass)
    !> Quick check for DRC.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Pexton, R. L., (LLNL)
    !***
    ! **Description:**
    !
    !            QUICK TEST FOR CARLSON INTEGRAL DRC
    !
    !***
    ! **Routines called:**  D1MACH, DRC, NUMXER, XERCLR, XGETF, XSETF

    !* REVISION HISTORY  (YYMMDD)
    !   790801  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   910708  Minor modifications in use of KPRINT.  (WRB)
    USE slatec, ONLY : D1MACH, DRC, NUMXER, XERCLR, XGETF, XSETF
    INTEGER :: Kprint, Ipass, contrl, kontrl, Lun, ier
    INTEGER :: ipass1, ipass2, ipass3, ipass4
    REAL(DP) :: pi, trc, dif
    !* FIRST EXECUTABLE STATEMENT  QCDRC
    CALL XERCLR
    CALL XGETF(contrl)
    IF( Kprint>=3 ) THEN
      kontrl = +1
    ELSE
      kontrl = 0
    END IF
    CALL XSETF(kontrl)
    !
    !  FORCE ERROR 1
    !
    IF( Kprint>=3 ) WRITE (Lun,99001)
    99001 FORMAT (' DRC - FORCE ERROR 1 TO OCCUR')
    trc = DRC(-1.0D0,-1.0D0,ier)
    ier = NUMXER(ier)
    IF( ier==1 ) THEN
      ipass1 = 1
    ELSE
      ipass1 = 0
    END IF
    CALL XERCLR
    !
    !  FORCE ERROR 2
    !
    IF( Kprint>=3 ) WRITE (Lun,99002)
    99002 FORMAT (' DRC - FORCE ERROR 2 TO OCCUR')
    trc = DRC(D1MACH(1),D1MACH(1),ier)
    ier = NUMXER(ier)
    IF( ier==2 ) THEN
      ipass2 = 1
    ELSE
      ipass2 = 0
    END IF
    CALL XERCLR
    !
    !  FORCE ERROR 3
    !
    IF( Kprint>=3 ) WRITE (Lun,99003)
    99003 FORMAT (' DRC - FORCE ERROR 3 TO OCCUR')
    trc = DRC(D1MACH(2),D1MACH(2),ier)
    ier = NUMXER(ier)
    IF( ier==3 ) THEN
      ipass3 = 1
    ELSE
      ipass3 = 0
    END IF
    CALL XERCLR
    !
    !  ARGUMENTS IN RANGE
    !
    pi = 3.141592653589793238462643383279D0
    trc = DRC(0.0D0,0.25D0,ier)
    CALL XERCLR
    dif = trc - pi
    IF( (ABS(dif/pi)<1000.0D0*D1MACH(4)) .AND. (ier==0) ) THEN
      ipass4 = 1
    ELSE
      ipass4 = 0
    END IF
    Ipass = MIN(ipass1,ipass2,ipass3,ipass4)
    IF( Kprint<=0 ) THEN
    ELSEIF( Kprint==1 ) THEN
      IF( Ipass/=1 ) WRITE (Lun,99006)
    ELSEIF( Ipass==1 ) THEN
      WRITE (Lun,99004)
      99004 FORMAT (' DRC - PASSED')
    ELSE
      WRITE (Lun,99006)
      IF( ipass4==0 ) WRITE (Lun,99005) pi, trc, dif
      99005 FORMAT (' CORRECT ANSWER =',1PD21.14/'COMPUTED ANSWER =',&
        D21.14/'     DIFFERENCE =',D21.14)
    END IF
    CALL XSETF(contrl)
    99006 FORMAT (' DRC - FAILED')
  END SUBROUTINE QCDRC
  !** QCDRD
  SUBROUTINE QCDRD(Lun,Kprint,Ipass)
    !> Quick check for DRD.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Pexton, R. L., (LLNL)
    !***
    ! **Description:**
    !
    !            QUICK TEST FOR CARLSON INTEGRAL DRD
    !
    !***
    ! **Routines called:**  D1MACH, DRD, NUMXER, XERCLR, XGETF, XSETF

    !* REVISION HISTORY  (YYMMDD)
    !   790801  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   910708  Minor modifications in use of KPRINT.  (WRB)
    !   930214  Added more digits to BLEM.  (WRB)
    USE slatec, ONLY : D1MACH, DRD, NUMXER, XERCLR, XGETF, XSETF
    INTEGER :: Kprint, Ipass, contrl, kontrl, Lun, ier
    INTEGER :: ipass1, ipass2, ipass3, ipass4
    REAL(DP) :: blem, trd, dif
    !* FIRST EXECUTABLE STATEMENT  QCDRD
    CALL XERCLR
    CALL XGETF(contrl)
    IF( Kprint>=3 ) THEN
      kontrl = +1
    ELSE
      kontrl = 0
    END IF
    CALL XSETF(kontrl)
    !
    !  FORCE ERROR 1
    !
    IF( Kprint>=3 ) WRITE (Lun,99001)
    99001 FORMAT (' DRD - FORCE ERROR 1 TO OCCUR')
    trd = DRD(-1.0D0,-1.0D0,-1.0D0,ier)
    ier = NUMXER(ier)
    IF( ier==1 ) THEN
      ipass1 = 1
    ELSE
      ipass1 = 0
    END IF
    CALL XERCLR
    !
    !  FORCE ERROR 2
    !
    IF( Kprint>=3 ) WRITE (Lun,99002)
    99002 FORMAT (' DRD - FORCE ERROR 2 TO OCCUR')
    trd = DRD(1.0D0,1.0D0,-1.0D0,ier)
    ier = NUMXER(ier)
    IF( ier==2 ) THEN
      ipass2 = 1
    ELSE
      ipass2 = 0
    END IF
    CALL XERCLR
    !
    !  FORCE ERROR 3
    !
    IF( Kprint>=3 ) WRITE (Lun,99003)
    99003 FORMAT (' DRD - FORCE ERROR 3 TO OCCUR')
    trd = DRD(D1MACH(2),D1MACH(2),D1MACH(2),ier)
    ier = NUMXER(ier)
    IF( ier==3 ) THEN
      ipass3 = 1
    ELSE
      ipass3 = 0
    END IF
    CALL XERCLR
    !
    !  ARGUMENTS IN RANGE
    !  BLEM=3 * LEMNISCATE CONSTANT B
    !
    blem = 1.797210352103388311159883738420485817341D0
    trd = DRD(0.0D0,2.0D0,1.0D0,ier)
    CALL XERCLR
    dif = trd - blem
    IF( (ABS(dif/blem)<1000.0D0*D1MACH(4)) .AND. (ier==0) ) THEN
      ipass4 = 1
    ELSE
      Ipass = 0
    END IF
    Ipass = MIN(ipass1,ipass2,ipass3,ipass4)
    IF( Kprint<=0 ) THEN
    ELSEIF( Kprint==1 ) THEN
      IF( Ipass/=1 ) WRITE (Lun,99006)
    ELSEIF( Ipass==1 ) THEN
      WRITE (Lun,99004)
      99004 FORMAT (' DRD - PASSED')
    ELSE
      WRITE (Lun,99006)
      IF( ipass4==0 ) WRITE (Lun,99005) blem, trd, dif
      99005 FORMAT (' CORRECT ANSWER =',1PD21.14/'COMPUTED ANSWER =',&
        D21.14/'     DIFFERENCE =',D21.14)
    END IF
    CALL XSETF(contrl)
    99006 FORMAT (' DRD - FAILED')
  END SUBROUTINE QCDRD
  !** QCDRF
  SUBROUTINE QCDRF(Lun,Kprint,Ipass)
    !> Quick check for DRF.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Pexton, R. L., (LLNL)
    !***
    ! **Description:**
    !
    !            QUICK TEST FOR CARLSON INTEGRAL DRF
    !
    !***
    ! **Routines called:**  D1MACH, DRF, NUMXER, XERCLR, XGETF, XSETF

    !* REVISION HISTORY  (YYMMDD)
    !   790801  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   910708  Minor modifications in use of KPRINT.  (WRB)
    !   930214  Added more digits to ALEM.  (WRB)
    USE slatec, ONLY : D1MACH, DRF, NUMXER, XERCLR, XGETF, XSETF
    INTEGER :: Kprint, Ipass, contrl, kontrl, Lun, ier
    INTEGER :: ipass1, ipass2, ipass3, ipass4
    REAL(DP) :: alem, trf, dif
    !* FIRST EXECUTABLE STATEMENT  QCDRF
    CALL XERCLR
    CALL XGETF(contrl)
    IF( Kprint>=3 ) THEN
      kontrl = +1
    ELSE
      kontrl = 0
    END IF
    CALL XSETF(kontrl)
    !
    !  FORCE ERROR 1
    !
    IF( Kprint>=3 ) WRITE (Lun,99001)
    99001 FORMAT (' DRF - FORCE ERROR 1 TO OCCUR')
    trf = DRF(-1.0D0,-1.0D0,-1.0D0,ier)
    ier = NUMXER(ier)
    IF( ier==1 ) THEN
      ipass1 = 1
    ELSE
      ipass1 = 0
    END IF
    CALL XERCLR
    !
    !  FORCE ERROR 2
    !
    IF( Kprint>=3 ) WRITE (Lun,99002)
    99002 FORMAT (' DRF - FORCE ERROR 2 TO OCCUR')
    trf = DRF(D1MACH(1),D1MACH(1),D1MACH(1),ier)
    ier = NUMXER(ier)
    IF( ier==2 ) THEN
      ipass2 = 1
    ELSE
      ipass2 = 0
    END IF
    CALL XERCLR
    !
    !  FORCE ERROR 3
    !
    IF( Kprint>=3 ) WRITE (Lun,99003)
    99003 FORMAT (' DRF - FORCE ERROR 3 TO OCCUR')
    trf = DRF(D1MACH(2),D1MACH(2),D1MACH(2),ier)
    ier = NUMXER(ier)
    IF( ier==3 ) THEN
      ipass3 = 1
    ELSE
      ipass3 = 0
    END IF
    CALL XERCLR
    !
    !  ARGUMENTS IN RANGE
    !  ALEM=LEMNISCATE CONSTANT A
    !
    alem = 1.3110287771460599052324197949455597068D0
    trf = DRF(0.0D0,1.0D0,2.0D0,ier)
    CALL XERCLR
    dif = trf - alem
    IF( (ABS(dif/alem)<1000.0D0*D1MACH(4)) .AND. (ier==0) ) THEN
      ipass4 = 1
    ELSE
      ipass4 = 0
    END IF
    Ipass = MIN(ipass1,ipass2,ipass3,ipass4)
    IF( Kprint==0 ) THEN
    ELSEIF( Kprint==1 ) THEN
      IF( Ipass/=1 ) WRITE (Lun,99006)
    ELSEIF( Ipass==1 ) THEN
      WRITE (Lun,99004)
      99004 FORMAT (' DRF - PASSED')
    ELSE
      WRITE (Lun,99006)
      IF( ipass4==0 ) WRITE (Lun,99005) alem, trf, dif
      99005 FORMAT (' CORRECT ANSWER =',1PD21.14/'COMPUTED ANSWER =',&
        D21.14/'     DIFFERENCE =',D21.14)
    END IF
    CALL XSETF(contrl)
    99006 FORMAT (' DRF - FAILED')
  END SUBROUTINE QCDRF
  !** QCDRJ
  SUBROUTINE QCDRJ(Lun,Kprint,Ipass)
    !> Quick check for DRJ.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Pexton, R. L., (LLNL)
    !***
    ! **Description:**
    !
    !            QUICK TEST FOR CARLSON INTEGRAL DRJ
    !
    !***
    ! **Routines called:**  D1MACH, DRJ, NUMXER, XERCLR, XGETF, XSETF

    !* REVISION HISTORY  (YYMMDD)
    !   790801  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   910708  Minor modifications in use of KPRINT.  (WRB)
    !   930214  Added more digits to CONSJ.  (WRB)
    USE slatec, ONLY : D1MACH, DRJ, NUMXER, XERCLR, XGETF, XSETF
    INTEGER :: Kprint, Ipass, contrl, kontrl, Lun, ier
    INTEGER :: ipass1, ipass2, ipass3, ipass4
    REAL(DP) :: consj, trj, dif
    !* FIRST EXECUTABLE STATEMENT  QCDRJ
    CALL XERCLR
    CALL XGETF(contrl)
    IF( Kprint>=3 ) THEN
      kontrl = +1
    ELSE
      kontrl = 0
    END IF
    CALL XSETF(kontrl)
    !
    !  FORCE ERROR 1
    !
    IF( Kprint>=3 ) WRITE (Lun,99001)
    99001 FORMAT (' DRJ - FORCE ERROR 1 TO OCCUR')
    trj = DRJ(-1.0D0,-1.0D0,-1.0D0,-1.0D0,ier)
    ier = NUMXER(ier)
    IF( ier==1 ) THEN
      ipass1 = 1
    ELSE
      ipass1 = 0
    END IF
    CALL XERCLR
    !
    !  FORCE ERROR 2
    !
    IF( Kprint>=3 ) WRITE (Lun,99002)
    99002 FORMAT (' DRJ - FORCE ERROR 2 TO OCCUR')
    trj = DRJ(D1MACH(1),D1MACH(1),D1MACH(1),D1MACH(1),ier)
    ier = NUMXER(ier)
    IF( ier==2 ) THEN
      ipass2 = 1
    ELSE
      ipass2 = 0
    END IF
    CALL XERCLR
    !
    !  FORCE ERROR 3
    !
    IF( Kprint>=3 ) WRITE (Lun,99003)
    99003 FORMAT (' DRJ - FORCE ERROR 3 TO OCCUR')
    trj = DRJ(D1MACH(2),D1MACH(2),D1MACH(2),D1MACH(2),ier)
    ier = NUMXER(ier)
    IF( ier==3 ) THEN
      ipass3 = 1
    ELSE
      ipass3 = 0
    END IF
    CALL XERCLR
    !
    !  ARGUMENTS IN RANGE
    !
    consj = 0.14297579667156753833233879421985774801D0
    trj = DRJ(2.0D0,3.0D0,4.0D0,5.0D0,ier)
    CALL XERCLR
    dif = trj - consj
    IF( (ABS(dif/consj)<1000.0D0*D1MACH(4)) .AND. (ier==0) ) THEN
      ipass4 = 1
    ELSE
      ipass4 = 0
    END IF
    Ipass = MIN(ipass1,ipass2,ipass3,ipass4)
    IF( Kprint<=0 ) THEN
    ELSEIF( Kprint==1 ) THEN
      IF( Ipass/=1 ) WRITE (Lun,99006)
    ELSEIF( Ipass==1 ) THEN
      WRITE (Lun,99004)
      99004 FORMAT (' DRJ - PASSED')
    ELSE
      WRITE (Lun,99006)
      IF( ipass4==0 ) WRITE (Lun,99005) consj, trj, dif
      99005 FORMAT (' CORRECT ANSWER =',1PD21.14/'COMPUTED ANSWER =',&
        D21.14/'     DIFFERENCE =',D21.14)
    END IF
    CALL XSETF(contrl)
    99006 FORMAT (' DRJ - FAILED')
  END SUBROUTINE QCDRJ
END MODULE TEST14_MOD
!** TEST14
PROGRAM TEST14
  USE TEST14_MOD, ONLY : QCDRC, QCDRD, QCDRF, QCDRJ
  USE slatec, ONLY : I1MACH, XSETF, XSETUN, XERMAX
  USE common_mod, ONLY : GET_ARGUMENT
  IMPLICIT NONE
  !> Driver for testing SLATEC subprograms
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  C14
  !***
  ! **Type:**      DOUBLE PRECISION (TEST13-S, TEST14-D)
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
  !        DRC      DRD      DRF      DRJ
  !
  !***
  ! **References:**  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  I1MACH, QCDRC, QCDRD, QCDRF, QCDRJ, XERMAX, XSETF,
  !                    XSETUN

  !* REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)
  INTEGER :: ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST14
  lun = I1MACH(2)
  lin = I1MACH(1)
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  CALL GET_ARGUMENT(kprint)
  CALL XERMAX(1000)
  CALL XSETUN(lun)
  IF( kprint<=1 ) THEN
    CALL XSETF(0)
  ELSE
    CALL XSETF(1)
  END IF
  !
  !     Test double precision Carlson elliptic routines
  !
  CALL QCDRC(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  CALL QCDRD(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  CALL QCDRF(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  CALL QCDRJ(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST14 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST14  *************')
  END IF
  STOP
END PROGRAM TEST14
