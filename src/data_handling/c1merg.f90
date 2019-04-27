!** C1MERG
SUBROUTINE C1MERG(Tcos,I1,M1,I2,M2,I3)
  !>
  !  Merge two strings of complex numbers.  Each string is
  !            ascending by the real part.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      COMPLEX (S1MERG-S, D1MERG-D, C1MERG-C, I1MERG-I)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !   This subroutine merges two ascending strings of numbers in the
  !   array TCOS.  The first string is of length M1 and starts at
  !   TCOS(I1+1).  The second string is of length M2 and starts at
  !   TCOS(I2+1).  The merged string goes into TCOS(I3+1).  The ordering
  !   is on the real part.
  !
  !***
  ! **See also:**  CMGNBN
  !***
  ! **Routines called:**  CCOPY

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   910408  Modified to use IF-THEN-ELSE.  Make it look like MERGE
  !           which was modified earlier due to compiler problems on
  !           the IBM RS6000.  (RWC)
  !   920130  Code name changed from CMPMRG to C1MERG.  (WRB)
  INTEGER I1, I2, I3, M1, M2
  COMPLEX Tcos(*)
  !
  INTEGER j1, j2, j3
  !
  !* FIRST EXECUTABLE STATEMENT  C1MERG
  IF ( M1==0.AND.M2==0 ) RETURN
  !
  IF ( M1==0.AND.M2/=0 ) THEN
    Tcos(I3+1:I3+M2) = Tcos(I2+1:I2+M2)
    RETURN
  END IF
  !
  IF ( M1/=0.AND.M2==0 ) THEN
    Tcos(I3+1:I3+M1) = Tcos(I1+1:I1+M1)
    RETURN
  END IF
  !
  j1 = 1
  j2 = 1
  j3 = 1
  DO
    !
    IF ( REAL(Tcos(j1+I1))<=REAL(Tcos(I2+j2)) ) THEN
      Tcos(I3+j3) = Tcos(I1+j1)
      j1 = j1 + 1
      IF ( j1>M1 ) THEN
        Tcos(I3+j3+1:I3+j3-j2+M2+1) = Tcos(I2+j2:I2+M2)
        RETURN
      END IF
    ELSE
      Tcos(I3+j3) = Tcos(I2+j2)
      j2 = j2 + 1
      IF ( j2>M2 ) THEN
        Tcos(I3+j3+1:I3+j3-j1+M1+1) = Tcos(I1+j1:I1+M1)
        RETURN
      END IF
    END IF
    j3 = j3 + 1
  END DO
END SUBROUTINE C1MERG
