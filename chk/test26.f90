!*==TEST26.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK TEST26
PROGRAM TEST26
  IMPLICIT NONE
  !*--TEST265
  !***BEGIN PROLOGUE  TEST26
  !***PURPOSE  Driver for testing SLATEC subprograms.
  !***LIBRARY   SLATEC
  !***CATEGORY  D2A4, D2B4
  !***TYPE      DOUBLE PRECISION (TEST25-S, TEST26-D)
  !***KEYWORDS  QUICK CHECK DRIVER
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
  !       Double precision SLAP subprograms
  !
  !***REFERENCES  Fong, Kirby W., Jefferson, Thomas H., Suyehiro,
  !                 Tokihiko, Walton, Lee, Guidelines to the SLATEC Common
  !                 Mathematical Library, March 21, 1989.
  !***ROUTINES CALLED  DLAPQC, I1MACH, XERMAX, XSETF, XSETUN
  !***REVISION HISTORY  (YYMMDD)
  !   920401  DATE WRITTEN
  !   920511  Added complete declaration section.  (WRB)
  !***END PROLOGUE  TEST26
  !     .. Local Scalars ..
  INTEGER ipass , kprint , lin , lun , nfail
  !     .. External Functions ..
  INTEGER I1MACH
  EXTERNAL I1MACH
  !     .. External Subroutines ..
  EXTERNAL DLAPQC , XERMAX , XSETF , XSETUN
  !***FIRST EXECUTABLE STATEMENT  TEST26
  lun = I1MACH(2)
  lin = I1MACH(1)
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  READ (lin,'(I1)') kprint
  CALL XSETUN(lun)
  IF ( kprint<=1 ) THEN
    CALL XSETF(0)
  ELSE
    CALL XSETF(1)
  ENDIF
  CALL XERMAX(1000)
  !
  !     Test SLAP (double precision)
  !
  CALL DLAPQC(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001   FORMAT (/' --------------TEST26 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002   FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST26 *************')
  ENDIF
  STOP
END PROGRAM TEST26
