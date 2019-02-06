*DECK HEADER
      SUBROUTINE HEADER (KPRINT)
C***BEGIN PROLOGUE  HEADER
C***PURPOSE  Print header for BLAS quick checks.
C***LIBRARY   SLATEC
C***AUTHOR  Lawson, C. L., (JPL)
C***ROUTINES CALLED  (NONE)
C***COMMON BLOCKS    COMBLA
C***REVISION HISTORY  (YYMMDD)
C   741212  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920210  Minor modifications to prologue and code.  (WRB)
C***END PROLOGUE  HEADER
      COMMON /COMBLA/ NPRINT, ICASE, N, INCX, INCY, MODE, PASS
      LOGICAL PASS
      CHARACTER*6 L(38)
C
      DATA L(1)  /'  SDOT'/
      DATA L(2)  /' DSDOT'/
      DATA L(3)  /'SDSDOT'/
      DATA L(4)  /'  DDOT'/
      DATA L(5)  /'DQDOTI'/
      DATA L(6)  /'DQDOTA'/
      DATA L(7)  /' CDOTC'/
      DATA L(8)  /' CDOTU'/
      DATA L(9)  /' SAXPY'/
      DATA L(10) /' DAXPY'/
      DATA L(11) /' CAXPY'/
      DATA L(12) /' SROTG'/
      DATA L(13) /' DROTG'/
      DATA L(14) /'  SROT'/
      DATA L(15) /'  DROT'/
      DATA L(16) /'SROTMG'/
      DATA L(17) /'DROTMG'/
      DATA L(18) /' SROTM'/
      DATA L(19) /' DROTM'/
      DATA L(20) /' SCOPY'/
      DATA L(21) /' DCOPY'/
      DATA L(22) /' CCOPY'/
      DATA L(23) /' SSWAP'/
      DATA L(24) /' DSWAP'/
      DATA L(25) /' CSWAP'/
      DATA L(26) /' SNRM2'/
      DATA L(27) /' DNRM2'/
      DATA L(28) /'SCNRM2'/
      DATA L(29) /' SASUM'/
      DATA L(30) /' DASUM'/
      DATA L(31) /'SCASUM'/
      DATA L(32) /' SSCAL'/
      DATA L(33) /' DSCAL'/
      DATA L(34) /' CSCAL'/
      DATA L(35) /'CSSCAL'/
      DATA L(36) /'ISAMAX'/
      DATA L(37) /'IDAMAX'/
      DATA L(38) /'ICAMAX'/
C***FIRST EXECUTABLE STATEMENT  HEADER
      IF (KPRINT .GE. 2) WRITE (NPRINT,9000) ICASE,L(ICASE)
      RETURN
C
 9000 FORMAT (' Test of subprogram number', I3, 2X, A)
      END
