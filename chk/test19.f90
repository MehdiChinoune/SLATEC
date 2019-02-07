!*==TEST19.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK TEST19
PROGRAM TEST19
  IMPLICIT NONE
  !*--TEST195
  !***BEGIN PROLOGUE  TEST19
  !***PURPOSE  Driver for testing SLATEC subprograms
  !***LIBRARY   SLATEC
  !***CATEGORY  D1B
  !***KEYWORDS  QUICK CHECK DRIVER
  !***TYPE      DOUBLE PRECISION (TEST18-S, TEST19-D, TEST20-C)
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
  !     Driver for testing SLATEC subprograms
  !        double precision Levels 2 and 3 BLAS routines
  !
  !***REFERENCES  Kirby W. Fong,  Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***ROUTINES CALLED  I1MACH, DBLAT2, DBLAT3, XERMAX, XSETF, XSETUN
  !***REVISION HISTORY  (YYMMDD)
  !   920601  DATE WRITTEN
  !***END PROLOGUE  TEST19
  INTEGER ipass , kprint , lin , lun , nfail
  !     .. External Functions ..
  INTEGER I1MACH
  EXTERNAL I1MACH
  !***FIRST EXECUTABLE STATEMENT  TEST19
  lun = I1MACH(2)
  lin = I1MACH(1)
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  READ (lin,'(I1)') kprint
  CALL XERMAX(1000)
  CALL XSETUN(lun)
  IF ( kprint<=1 ) THEN
    CALL XSETF(0)
  ELSE
    CALL XSETF(1)
  ENDIF
  !
  !     Test double precision Level 2 BLAS routines
  !
  CALL DBLAT2(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test double precision Level 3 BLAS routines
  !
  CALL DBLAT3(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001   FORMAT (/' --------------TEST19 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002   FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST19 *************')
  ENDIF
  STOP
END PROGRAM TEST19
