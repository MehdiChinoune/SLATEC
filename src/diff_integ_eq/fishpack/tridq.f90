!** TRIDQ
SUBROUTINE TRIDQ(Mr,A,B,C,Y,D)
  !> Subsidiary to POIS3D
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (TRIDQ-S)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **See also:**  POIS3D
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900308  Renamed routine from TRID to TRIDQ.  (WRB)
  !   900402  Added TYPE section.  (WRB)

  INTEGER :: Mr
  REAL(SP) :: A(Mr), B(Mr), C(Mr), D(Mr), Y(Mr)
  INTEGER :: i, ip, m, mm1
  REAL(SP) :: z
  !* FIRST EXECUTABLE STATEMENT  TRIDQ
  m = Mr
  mm1 = m - 1
  z = 1./B(1)
  D(1) = C(1)*z
  Y(1) = Y(1)*z
  DO i = 2, mm1
    z = 1./(B(i)-A(i)*D(i-1))
    D(i) = C(i)*z
    Y(i) = (Y(i)-A(i)*Y(i-1))*z
  END DO
  z = B(m) - A(m)*D(mm1)
  IF( z/=0. ) THEN
    Y(m) = (Y(m)-A(m)*Y(mm1))/z
  ELSE
    Y(m) = 0.
  END IF
  DO ip = 1, mm1
    i = m - ip
    Y(i) = Y(i) - D(i)*Y(i+1)
  END DO
END SUBROUTINE TRIDQ
