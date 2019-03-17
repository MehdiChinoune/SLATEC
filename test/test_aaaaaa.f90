MODULE TEST01_MOD
  IMPLICIT NONE

CONTAINS
  !DECK QC6A
  SUBROUTINE QC6A(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !***BEGIN PROLOGUE  QC6A
    !***PURPOSE  Test subroutine AAAAAA.
    !***LIBRARY   SLATEC
    !***TYPE      ALL (QC6A-A)
    !***AUTHOR  Boland, W. Robert, (LANL)
    !***DESCRIPTION
    !
    ! *Usage:
    !
    !        INTEGER  LUN, KPRINT, IPASS
    !
    !        CALL  QC6A (LUN, KPRINT, IPASS)
    !
    ! *Arguments:
    !
    !     LUN   :IN  is the unit number to which output is to be written.
    !
    !     KPRINT:IN  controls the amount of output, as specified in the
    !                SLATEC Guidelines.
    !
    !     IPASS:OUT  indicates whether the test passed or failed.
    !                A value of one is good, indicating no failures.
    !
    ! *Description:
    !
    !   This routine tests the SLATEC routine AAAAAA to see if the version
    !   number in the SLATEC library source is the same as the quick check
    !   version number.
    !
    !***ROUTINES CALLED  AAAAAA
    !***REVISION HISTORY  (YYMMDD)
    !   890713  DATE WRITTEN
    !   921215  Updated for Version 4.0.  (WRB)
    !   930701  Updated for Version 4.1.  (WRB)
    !***END PROLOGUE  QC6A
    !
    !*Internal Notes:
    !
    !     Data set-up is done via a PARAMETER statement.
    !
    !**End
    !
    !  Declare arguments.
    !
    INTEGER Lun, Kprint, Ipass
    !
    !  DECLARE VARIABLES.
    !
    CHARACTER(16) :: ver, VERSN
    PARAMETER (VERSN=' 4.1')
    !
    !***FIRST EXECUTABLE STATEMENT  QC6A
    IF ( Kprint>=3 ) WRITE (Lun,99001)
    !
    ! FORMATs.
    !
    99001 FORMAT ('1 CODE TO TEST SLATEC ROUTINE AAAAAA')
    CALL AAAAAA(ver)
    IF ( ver==VERSN ) THEN
      Ipass = 1
      IF ( Kprint>=3 ) THEN
        WRITE (Lun,99006)
        WRITE (Lun,99002) ver
        99002 FORMAT (' *** Passed -- version number = ',A16)
      ENDIF
    ELSE
      Ipass = 0
      IF ( Kprint>=3 ) WRITE (Lun,99006)
      IF ( Kprint>=2 ) WRITE (Lun,99003) ver, VERSN
      99003 FORMAT (' *** Failed -- version number from AAAAAA = ',A16,&
        ' but expected version number = ',A16)
    ENDIF
    !
    !     Terminate.
    !
    IF ( Kprint>=2.AND.Ipass==1 ) WRITE (Lun,99004)
    99004 FORMAT (/' ************QC6A   PASSED  ALL TESTS ****************')
    IF ( Kprint>=1.AND.Ipass==0 ) WRITE (Lun,99005)
    99005 FORMAT (/' ************QC6A   FAILED SOME TESTS ****************')
    RETURN
    99006 FORMAT (/' QC6A RESULTS')
    !------------- LAST LINE OF QC6A FOLLOWS -----------------------------
  END SUBROUTINE QC6A
END MODULE TEST01_MOD
!DECK TEST01
PROGRAM TEST01
  USE TEST01_MOD
  IMPLICIT NONE
  !***BEGIN PROLOGUE  TEST01
  !***PURPOSE  Driver for testing SLATEC subprogram
  !            AAAAAA
  !***LIBRARY   SLATEC
  !***CATEGORY  Z
  !***TYPE      ALL (TEST01-A)
  !***KEYWORDS  AAAAAA, QUICK CHECK DRIVER
  !***AUTHOR  SLATEC Common Mathematical Library Committee
  !***DESCRIPTION
  !
  ! *Usage:
  !     One input data record is required
  !         READ (LIN, '(I1)') KPRINT
  !
  ! *Arguments:
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
  ! *Description:
  !     Driver for testing SLATEC subprogram
  !        AAAAAA
  !
  !***REFERENCES  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***ROUTINES CALLED  I1MACH, QC6A, XERMAX, XSETF, XSETUN
  !***REVISION HISTORY  (YYMMDD)
  !   890713  DATE WRITTEN
  !   900524  Cosmetic changes to code.  (WRB)
  !***END PROLOGUE  TEST01
  INTEGER I1MACH
  INTEGER ipass, kprint, lin, lun, nfail
  !***FIRST EXECUTABLE STATEMENT  TEST01
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
  ENDIF
  !
  !     Test AAAAAA
  !
  CALL QC6A(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST01 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST01  *************')
  ENDIF
  STOP
END PROGRAM TEST01
