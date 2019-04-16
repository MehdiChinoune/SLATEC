!** S1MERG
SUBROUTINE S1MERG(Tcos,I1,M1,I2,M2,I3)
  !>
  !***
  !  Merge two strings of ascending real numbers.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (S1MERG-S, D1MERG-D, C1MERG-C, I1MERG-I)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !   This subroutine merges two ascending strings of numbers in the
  !   array TCOS.  The first string is of length M1 and starts at
  !   TCOS(I1+1).  The second string is of length M2 and starts at
  !   TCOS(I2+1).  The merged string goes into TCOS(I3+1).
  !
  !***
  ! **See also:**  GENBUN
  !***
  ! **Routines called:**  SCOPY

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   901120  Modified to use IF-THEN-ELSE.  Previous spaghetti code did
  !           not compile correctly with optimization on the IBM RS6000.
  !           (RWC)
  !   920130  Code name changed from MERGE to S1MERG.  (WRB)

  INTEGER I1, I2, I3, M1, M2
  REAL Tcos(*)
  !
  INTEGER j1, j2, j3
  !
  !* FIRST EXECUTABLE STATEMENT  S1MERG
  IF ( M1==0.AND.M2==0 ) RETURN
  !
  IF ( M1==0.AND.M2/=0 ) THEN
    CALL SCOPY(M2,Tcos(I2+1),1,Tcos(I3+1),1)
    RETURN
  END IF
  !
  IF ( M1/=0.AND.M2==0 ) THEN
    CALL SCOPY(M1,Tcos(I1+1),1,Tcos(I3+1),1)
    RETURN
  END IF
  !
  j1 = 1
  j2 = 1
  j3 = 1
  DO
    !
    IF ( Tcos(I1+j1)<=Tcos(I2+j2) ) THEN
      Tcos(I3+j3) = Tcos(I1+j1)
      j1 = j1 + 1
      IF ( j1>M1 ) THEN
        CALL SCOPY(M2-j2+1,Tcos(I2+j2),1,Tcos(I3+j3+1),1)
        RETURN
      END IF
    ELSE
      Tcos(I3+j3) = Tcos(I2+j2)
      j2 = j2 + 1
      IF ( j2>M2 ) THEN
        CALL SCOPY(M1-j1+1,Tcos(I1+j1),1,Tcos(I3+j3+1),1)
        RETURN
      END IF
    END IF
    j3 = j3 + 1
  END DO
END SUBROUTINE S1MERG
