*DECK TEST23
      PROGRAM TEST23
C***BEGIN PROLOGUE  TEST23
C***PURPOSE  Driver for testing SLATEC subprograms
C***LIBRARY   SLATEC
C***CATEGORY  D2
C***TYPE      COMPLEX (TEST23-S)
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
C        CGECO    CGEDI    CGEFA    CGESL
C        CGBCO    CGBDI    CGBFA    CGBSL
C        CPOCO    CPODI    CPOFA    CPOSL
C        CPPCO    CPPDI    CPPFA    CPPSL
C        CPBCO    CPBDI    CPBFA    CPBSL
C        CSICO    CSIDI    CSIFA    CSISL
C        CSPCO    CSPDI    CSPFA    CSPSL
C        CHICO    CHIDI    CHIFA    CHISL
C        CHPCO    CHPDI    CHPFA    CHPSL
C        CTRCO    CTRDI      -      CTRSL
C        CGTSL
C        CPTSL
C        CCHDC
C        CQRDC    CQRSL
C        CSVDC
C
C***REFERENCES  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
C                 and Lee Walton, Guide to the SLATEC Common Mathema-
C                 tical Library, April 10, 1990.
C***ROUTINES CALLED  CCHQC, CGBQC, CGECK, CGTQC, CHIQC, CHPQC, CPBQC,
C                    CPOQC, CPPQC, CPTQC, CQRQC, CSIQC, CSPQC, CSVQC,
C                    CTRQC, I1MACH, XERMAX, XSETF, XSETUN
C***REVISION HISTORY  (YYMMDD)
C   890618  DATE WRITTEN
C   890618  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900524  Cosmetic changes to code.  (WRB)
C***END PROLOGUE  TEST23
      INTEGER KPRINT, LIN, LUN, NERR, NFAIL
C***FIRST EXECUTABLE STATEMENT  TEST23
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
C     Test LINPACK routines
C
      CALL CGECK(LUN,KPRINT,NERR)
      NFAIL = NFAIL + NERR
      CALL CGBQC(LUN,KPRINT,NERR)
      NFAIL = NFAIL + NERR
      CALL CPOQC(LUN,KPRINT,NERR)
      NFAIL = NFAIL + NERR
      CALL CPPQC(LUN,KPRINT,NERR)
      NFAIL = NFAIL + NERR
      CALL CPBQC(LUN,KPRINT,NERR)
      NFAIL = NFAIL + NERR
      CALL CSIQC(LUN,KPRINT,NERR)
      NFAIL = NFAIL + NERR
      CALL CSPQC(LUN,KPRINT,NERR)
      NFAIL = NFAIL + NERR
      CALL CHIQC(LUN,KPRINT,NERR)
      NFAIL = NFAIL + NERR
      CALL CHPQC(LUN,KPRINT,NERR)
      NFAIL = NFAIL + NERR
      CALL CTRQC(LUN,KPRINT,NERR)
      NFAIL = NFAIL + NERR
      CALL CGTQC(LUN,KPRINT,NERR)
      NFAIL = NFAIL + NERR
      CALL CPTQC(LUN,KPRINT,NERR)
      NFAIL = NFAIL + NERR
      CALL CCHQC(LUN,KPRINT,NERR)
      NFAIL = NFAIL + NERR
      CALL CQRQC(LUN,KPRINT,NERR)
      NFAIL = NFAIL + NERR
      CALL CSVQC(LUN,KPRINT,NERR)
      NFAIL = NFAIL + NERR
C
C
C     Write PASS or FAIL message
C
      IF (NFAIL .EQ. 0) THEN
         WRITE (LUN, 9000)
      ELSE
         WRITE (LUN, 9010) NFAIL
      ENDIF
      STOP
 9000 FORMAT (/' --------------TEST23 PASSED ALL TESTS----------------')
 9010 FORMAT (/' ************* WARNING -- ', I5,
     1        ' TEST(S) FAILED IN PROGRAM TEST23 *************')
      END
