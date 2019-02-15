!DECK HEADER
SUBROUTINE HEADER(Kprint)
  IMPLICIT NONE
  INTEGER ICAse, INCx, INCy, Kprint, MODe, N, NPRint
  !***BEGIN PROLOGUE  HEADER
  !***PURPOSE  Print header for BLAS quick checks.
  !***LIBRARY   SLATEC
  !***AUTHOR  Lawson, C. L., (JPL)
  !***ROUTINES CALLED  (NONE)
  !***COMMON BLOCKS    COMBLA
  !***REVISION HISTORY  (YYMMDD)
  !   741212  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920210  Minor modifications to prologue and code.  (WRB)
  !***END PROLOGUE  HEADER
  COMMON /COMBLA/ NPRint, ICAse, N, INCx, INCy, MODe, PASs
  LOGICAL PASs
  CHARACTER(6) :: l(38)
  !
  DATA l(1)/'  SDOT'/
  DATA l(2)/' DSDOT'/
  DATA l(3)/'SDSDOT'/
  DATA l(4)/'  DDOT'/
  DATA l(5)/'DQDOTI'/
  DATA l(6)/'DQDOTA'/
  DATA l(7)/' CDOTC'/
  DATA l(8)/' CDOTU'/
  DATA l(9)/' SAXPY'/
  DATA l(10)/' DAXPY'/
  DATA l(11)/' CAXPY'/
  DATA l(12)/' SROTG'/
  DATA l(13)/' DROTG'/
  DATA l(14)/'  SROT'/
  DATA l(15)/'  DROT'/
  DATA l(16)/'SROTMG'/
  DATA l(17)/'DROTMG'/
  DATA l(18)/' SROTM'/
  DATA l(19)/' DROTM'/
  DATA l(20)/' SCOPY'/
  DATA l(21)/' DCOPY'/
  DATA l(22)/' CCOPY'/
  DATA l(23)/' SSWAP'/
  DATA l(24)/' DSWAP'/
  DATA l(25)/' CSWAP'/
  DATA l(26)/' SNRM2'/
  DATA l(27)/' DNRM2'/
  DATA l(28)/'SCNRM2'/
  DATA l(29)/' SASUM'/
  DATA l(30)/' DASUM'/
  DATA l(31)/'SCASUM'/
  DATA l(32)/' SSCAL'/
  DATA l(33)/' DSCAL'/
  DATA l(34)/' CSCAL'/
  DATA l(35)/'CSSCAL'/
  DATA l(36)/'ISAMAX'/
  DATA l(37)/'IDAMAX'/
  DATA l(38)/'ICAMAX'/
  !***FIRST EXECUTABLE STATEMENT  HEADER
  IF ( Kprint>=2 ) WRITE (NPRint,99001) ICAse, l(ICAse)
  !
  99001 FORMAT (' Test of subprogram number',I3,2X,A)
  RETURN
END SUBROUTINE HEADER
