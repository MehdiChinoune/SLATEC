*DECK TEST32
      PROGRAM TEST32
C***BEGIN PROLOGUE  TEST32
C***PURPOSE  Driver for testing SLATEC subprograms
C***LIBRARY   SLATEC
C***CATEGORY  E1A
C***TYPE      SINGLE PRECISION (TEST32-S, TEST33-D)
C***KEYWORDS  QUICK CHECK DRIVER
C***AUTHOR  SLATEC Common Mathematical Library Committee
C***DESCRIPTION
C
C *Usage:
C     One input data record is required
C         READ (LIN, '(I1)') KPRINT
C
C *Arguments:
C     KPRINT = 0  Quick checks - No printing.
C                 Driver       - Short pass or fail message printed.
C              1  Quick checks - No message printed for passed tests,
C                                short message printed for failed tests.
C                 Driver       - Short pass or fail message printed.
C              2  Quick checks - Print short message for passed tests,
C                                fuller information for failed tests.
C                 Driver       - Pass or fail message printed.
C              3  Quick checks - Print complete quick check results.
C                 Driver       - Pass or fail message printed.
C
C *Description:
C     Driver for testing SLATEC subprograms
C        PCHIP
C
C***REFERENCES  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
C                 and Lee Walton, Guide to the SLATEC Common Mathema-
C                 tical Library, April 10, 1990.
C***ROUTINES CALLED  I1MACH, PCHQK1, PCHQK2, PCHQK3, PCHQK4, PCHQK5,
C                    XERMAX, XSETF, XSETUN
C***REVISION HISTORY  (YYMMDD)
C   890618  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900314  Added new quick checks PCHQK3, PCHQK4.  (FNF)
C   900315  Corrected category record.  (FNF)
C   900321  Moved IPASS to call sequences for SLATEC standards.  (FNF)
C   900524  Cosmetic changes to code.  (WRB)
C   930318  Added new quick check PCHQK5.  (WRB,FNF)
C***END PROLOGUE  TEST32
      INTEGER IPASS, KPRINT, LIN, LUN, NFAIL
C***FIRST EXECUTABLE STATEMENT  TEST32
      LUN   = I1MACH(2)
      LIN   = I1MACH(1)
      NFAIL = 0
C
C     Read KPRINT parameter
C
      READ (LIN, '(I1)') KPRINT
      CALL XERMAX(1000)
      CALL XSETUN(LUN)
      IF (KPRINT .LE. 1) THEN
         CALL XSETF(0)
      ELSE
         CALL XSETF(1)
      ENDIF
C
C     Test PCHIP evaluators
C
      CALL PCHQK1(LUN,KPRINT,IPASS)
      IF (IPASS .EQ. 0) NFAIL = NFAIL + 1
C
C     Test PCHIP integrators
C
      CALL PCHQK2(LUN,KPRINT,IPASS)
      IF (IPASS .EQ. 0) NFAIL = NFAIL + 1
C
C     Test PCHIP interpolators
C
      CALL PCHQK3(LUN,KPRINT,IPASS)
      IF (IPASS .EQ. 0) NFAIL = NFAIL + 1
C
C     Test PCHIP monotonicity checker
C
      CALL PCHQK4(LUN,KPRINT,IPASS)
      IF (IPASS .EQ. 0) NFAIL = NFAIL + 1
C
C     Test PCH to B-spline conversion.
C
      CALL PCHQK5(LUN,KPRINT,IPASS)
      IF (IPASS .EQ. 0) NFAIL = NFAIL + 1
C
C     Write PASS or FAIL message
C
      IF (NFAIL .EQ. 0) THEN
         WRITE (LUN, 9000)
      ELSE
         WRITE (LUN, 9010) NFAIL
      ENDIF
      STOP
 9000 FORMAT (/' --------------TEST32 PASSED ALL TESTS----------------')
 9010 FORMAT (/' ************* WARNING -- ', I5,
     1        ' TEST(S) FAILED IN PROGRAM TEST32 *************')
      END
