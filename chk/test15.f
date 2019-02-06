*DECK TEST15
      PROGRAM TEST15
C***BEGIN PROLOGUE  TEST15
C***PURPOSE  Driver for testing SLATEC subprograms
C            RC3JJ    RC3JM    RC6J
C***LIBRARY   SLATEC
C***CATEGORY  C19
C***TYPE      SINGLE PRECISION (TEST15-S, TEST16-D)
C***KEYWORDS  3J COEFFICIENTS, 3J SYMBOLS, 6J COEFFICIENTS, 6J SYMBOLS,
C             CLEBSCH-GORDAN COEFFICIENTS, QUICK CHECK DRIVER,
C             RACAH COEFFICIENTS, VECTOR ADDITION COEFFICIENTS,
C             WIGNER COEFFICIENTS
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
C        RC3JJ    RC3JM    RC6J
C
C***REFERENCES  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
C                 and Lee Walton, Guide to the SLATEC Common Mathema-
C                 tical Library, April 10, 1990.
C***ROUTINES CALLED  I1MACH, QC36J, XERMAX, XSETF, XSETUN
C***REVISION HISTORY  (YYMMDD)
C   891130  DATE WRITTEN
C***END PROLOGUE  TEST15
      INTEGER IPASS, KPRINT, LIN, LUN, NFAIL
C***FIRST EXECUTABLE STATEMENT  TEST15
      LUN = I1MACH(2)
      LIN = I1MACH(1)
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
C     Test single precision 3J6J routines
C
      CALL QC36J(LUN, KPRINT, IPASS)
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
 9000 FORMAT (/' --------------TEST15 PASSED ALL TESTS----------------')
 9010 FORMAT (/' ************* WARNING -- ', I5,
     1        ' TEST(S) FAILED IN PROGRAM TEST15  *************')
      END
