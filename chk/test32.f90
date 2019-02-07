!*==TEST32.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK TEST32
PROGRAM TEST32
  IMPLICIT NONE
  !*--TEST325
  !*** Start of declarations inserted by SPAG
  INTEGER I1MACH
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  TEST32
  !***PURPOSE  Driver for testing SLATEC subprograms
  !***LIBRARY   SLATEC
  !***CATEGORY  E1A
  !***TYPE      SINGLE PRECISION (TEST32-S, TEST33-D)
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
  !        PCHIP
  !
  !***REFERENCES  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***ROUTINES CALLED  I1MACH, PCHQK1, PCHQK2, PCHQK3, PCHQK4, PCHQK5,
  !                    XERMAX, XSETF, XSETUN
  !***REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900314  Added new quick checks PCHQK3, PCHQK4.  (FNF)
  !   900315  Corrected category record.  (FNF)
  !   900321  Moved IPASS to call sequences for SLATEC standards.  (FNF)
  !   900524  Cosmetic changes to code.  (WRB)
  !   930318  Added new quick check PCHQK5.  (WRB,FNF)
  !***END PROLOGUE  TEST32
  INTEGER ipass , kprint , lin , lun , nfail
  !***FIRST EXECUTABLE STATEMENT  TEST32
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
  !     Test PCHIP evaluators
  !
  CALL PCHQK1(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test PCHIP integrators
  !
  CALL PCHQK2(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test PCHIP interpolators
  !
  CALL PCHQK3(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test PCHIP monotonicity checker
  !
  CALL PCHQK4(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test PCH to B-spline conversion.
  !
  CALL PCHQK5(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001   FORMAT (/' --------------TEST32 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002   FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST32 *************')
  ENDIF
  STOP
END PROGRAM TEST32
