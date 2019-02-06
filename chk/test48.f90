!*==TEST48.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK TEST48
      PROGRAM TEST48
      IMPLICIT NONE
!*--TEST485
!*** Start of declarations inserted by SPAG
      INTEGER I1MACH
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  TEST48
!***PURPOSE  Driver for testing SLATEC subprograms
!***LIBRARY   SLATEC
!***CATEGORY  I1A2
!***TYPE      SINGLE PRECISION (TEST48-S, TEST49-D)
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
!       SDASSL
!
!***REFERENCES  Fong, Kirby W., Jefferson, Thomas H., Suyehiro,
!                 Tokihiko, Walton, Lee, Guidelines to the SLATEC Common
!                 Mathematical Library, March 21, 1989.
!***ROUTINES CALLED  I1MACH, SDASQC, XERMAX, XSETF, XSETUN
!***REVISION HISTORY  (YYMMDD)
!   891013  DATE WRITTEN
!   901001  Converted prologue to 4.0 format.  (FNF)
!   901009  Corrected GAMS classification code.  (FNF)
!***END PROLOGUE  TEST48
      INTEGER ipass , kprint , lin , lun , nfail
!***FIRST EXECUTABLE STATEMENT  TEST48
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
!     Test SDASSL
!
      CALL SDASQC(lun,kprint,ipass)
      IF ( ipass==0 ) nfail = nfail + 1
!
!     Write PASS or FAIL message
!
      IF ( nfail==0 ) THEN
        WRITE (lun,99001)
99001   FORMAT (/' --------------TEST48 PASSED ALL TESTS----------------')
      ELSE
        WRITE (lun,99002) nfail
99002   FORMAT (/' ************* WARNING -- ',I5,
     &          ' TEST(S) FAILED IN PROGRAM TEST48 *************')
      ENDIF
      STOP
      END PROGRAM TEST48
