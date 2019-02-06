*DECK TEST33
      PROGRAM TEST33
C***BEGIN PROLOGUE  TEST33
C***PURPOSE  Driver for testing SLATEC subprograms
C***LIBRARY   SLATEC
C***CATEGORY  E1A
C***TYPE      DOUBLE PRECISION (TEST32-S, TEST33-D)
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
C        DPCHIP
C
C***REFERENCES  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
C                 and Lee Walton, Guide to the SLATEC Common Mathema-
C                 tical Library, April 10, 1990.
C***ROUTINES CALLED  DPCHQ1, DPCHQ2, DPCHQ3, DPCHQ4, DPCHQ5, I1MACH,
C                    XERMAX, XSETF, XSETUN
C***REVISION HISTORY  (YYMMDD)
C   890618  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900319  Corrected category record.  (FNF)
C   900320  Added new quick checks DPCHQ3, DPCHQ4.  (FNF)
C   900321  Moved IPASS to call sequences for SLATEC standards.  (FNF)
C   900322  Corrected list of routines called.  (FNF)
C   900524  Cosmetic changes to code.  (WRB)
C   930318  Added new quick check DPCHQ5.  (WRB,FNF)
C***END PROLOGUE  TEST33
      INTEGER IPASS, KPRINT, LIN, LUN, NFAIL
C***FIRST EXECUTABLE STATEMENT  TEST33
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
C     Test DPCHIP evaluators
C
      CALL DPCHQ1(LUN,KPRINT,IPASS)
      IF (IPASS .EQ. 0) NFAIL = NFAIL + 1
C
C     Test DPCHIP integrators
C
      CALL DPCHQ2(LUN,KPRINT,IPASS)
      IF (IPASS .EQ. 0) NFAIL = NFAIL + 1
C
C     Test DPCHIP interpolators
C
      CALL DPCHQ3(LUN,KPRINT,IPASS)
      IF (IPASS .EQ. 0) NFAIL = NFAIL + 1
C
C     Test DPCHIP monotonicity checker
C
      CALL DPCHQ4(LUN,KPRINT,IPASS)
      IF (IPASS .EQ. 0) NFAIL = NFAIL + 1
C
C     Test PCH to B-spline conversion.
C
      CALL DPCHQ5(LUN,KPRINT,IPASS)
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
 9000 FORMAT (/' --------------TEST33 PASSED ALL TESTS----------------')
 9010 FORMAT (/' ************* WARNING -- ', I5,
     1        ' TEST(S) FAILED IN PROGRAM TEST33 *************')
      END
