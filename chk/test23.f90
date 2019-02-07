!*==TEST23.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK TEST23
PROGRAM TEST23
  IMPLICIT NONE
  !*--TEST235
  !*** Start of declarations inserted by SPAG
  INTEGER I1MACH
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  TEST23
  !***PURPOSE  Driver for testing SLATEC subprograms
  !***LIBRARY   SLATEC
  !***CATEGORY  D2
  !***TYPE      COMPLEX (TEST23-S)
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
  !        CGECO    CGEDI    CGEFA    CGESL
  !        CGBCO    CGBDI    CGBFA    CGBSL
  !        CPOCO    CPODI    CPOFA    CPOSL
  !        CPPCO    CPPDI    CPPFA    CPPSL
  !        CPBCO    CPBDI    CPBFA    CPBSL
  !        CSICO    CSIDI    CSIFA    CSISL
  !        CSPCO    CSPDI    CSPFA    CSPSL
  !        CHICO    CHIDI    CHIFA    CHISL
  !        CHPCO    CHPDI    CHPFA    CHPSL
  !        CTRCO    CTRDI      -      CTRSL
  !        CGTSL
  !        CPTSL
  !        CCHDC
  !        CQRDC    CQRSL
  !        CSVDC
  !
  !***REFERENCES  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***ROUTINES CALLED  CCHQC, CGBQC, CGECK, CGTQC, CHIQC, CHPQC, CPBQC,
  !                    CPOQC, CPPQC, CPTQC, CQRQC, CSIQC, CSPQC, CSVQC,
  !                    CTRQC, I1MACH, XERMAX, XSETF, XSETUN
  !***REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)
  !***END PROLOGUE  TEST23
  INTEGER kprint, lin, lun, nerr, nfail
  !***FIRST EXECUTABLE STATEMENT  TEST23
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
  !     Test LINPACK routines
  !
  CALL CGECK(lun,kprint,nerr)
  nfail = nfail + nerr
  CALL CGBQC(lun,kprint,nerr)
  nfail = nfail + nerr
  CALL CPOQC(lun,kprint,nerr)
  nfail = nfail + nerr
  CALL CPPQC(lun,kprint,nerr)
  nfail = nfail + nerr
  CALL CPBQC(lun,kprint,nerr)
  nfail = nfail + nerr
  CALL CSIQC(lun,kprint,nerr)
  nfail = nfail + nerr
  CALL CSPQC(lun,kprint,nerr)
  nfail = nfail + nerr
  CALL CHIQC(lun,kprint,nerr)
  nfail = nfail + nerr
  CALL CHPQC(lun,kprint,nerr)
  nfail = nfail + nerr
  CALL CTRQC(lun,kprint,nerr)
  nfail = nfail + nerr
  CALL CGTQC(lun,kprint,nerr)
  nfail = nfail + nerr
  CALL CPTQC(lun,kprint,nerr)
  nfail = nfail + nerr
  CALL CCHQC(lun,kprint,nerr)
  nfail = nfail + nerr
  CALL CQRQC(lun,kprint,nerr)
  nfail = nfail + nerr
  CALL CSVQC(lun,kprint,nerr)
  nfail = nfail + nerr
  !
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001   FORMAT (/' --------------TEST23 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002   FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST23 *************')
  ENDIF
  STOP
END PROGRAM TEST23
