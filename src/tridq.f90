!DECK TRIDQ
SUBROUTINE TRIDQ(Mr,A,B,C,Y,D)
  IMPLICIT NONE
  REAL A, B, C, D, Y, z
  INTEGER i, ip, m, mm1, Mr
  !***BEGIN PROLOGUE  TRIDQ
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to POIS3D
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (TRIDQ-S)
  !***AUTHOR  (UNKNOWN)
  !***SEE ALSO  POIS3D
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900308  Renamed routine from TRID to TRIDQ.  (WRB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  TRIDQ
  DIMENSION A(*), B(*), C(*), Y(*), D(*)
  !***FIRST EXECUTABLE STATEMENT  TRIDQ
  m = Mr
  mm1 = m - 1
  z = 1./B(1)
  D(1) = C(1)*z
  Y(1) = Y(1)*z
  DO i = 2, mm1
    z = 1./(B(i)-A(i)*D(i-1))
    D(i) = C(i)*z
    Y(i) = (Y(i)-A(i)*Y(i-1))*z
  ENDDO
  z = B(m) - A(m)*D(mm1)
  IF ( z/=0. ) THEN
    Y(m) = (Y(m)-A(m)*Y(mm1))/z
  ELSE
    Y(m) = 0.
  ENDIF
  DO ip = 1, mm1
    i = m - ip
    Y(i) = Y(i) - D(i)*Y(i+1)
  ENDDO
END SUBROUTINE TRIDQ
