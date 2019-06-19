!** I1MERG
SUBROUTINE I1MERG(Icos,I1,M1,I2,M2,I3)
  !> Merge two strings of ascending integers.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      INTEGER (S1MERG-S, D1MERG-D, C1MERG-C, I1MERG-I)
  !***
  ! **Author:**  Boland, W. Robert, (LANL)
  !           Clemens, Reginald, (PLK)
  !***
  ! **Description:**
  !
  !   This subroutine merges two ascending strings of integers in the
  !   array ICOS.  The first string is of length M1 and starts at
  !   ICOS(I1+1).  The second string is of length M2 and starts at
  !   ICOS(I2+1).  The merged string goes into ICOS(I3+1).
  !
  !***
  ! **Routines called:**  ICOPY

  !* REVISION HISTORY  (YYMMDD)
  !   920202  DATE WRITTEN
  INTEGER :: I1, I2, I3, M1, M2
  INTEGER :: Icos( MAX(I1+M1,I2+M2,I3+M1+M2) )
  !
  INTEGER :: j1, j2, j3
  !
  !* FIRST EXECUTABLE STATEMENT  I1MERG
  IF( M1==0 .AND. M2==0 ) RETURN
  !
  IF( M1==0 .AND. M2/=0 ) THEN
    Icos(I3+1:I3+M2) = Icos(I2+1:I2+M2)
    RETURN
  END IF
  !
  IF( M1/=0 .AND. M2==0 ) THEN
    Icos(I3+1:I3+M1) = Icos(I1+1:I1+M1)
    RETURN
  END IF
  !
  j1 = 1
  j2 = 1
  j3 = 1
  DO
    !
    IF( Icos(I1+j1)<=Icos(I2+j2) ) THEN
      Icos(I3+j3) = Icos(I1+j1)
      j1 = j1 + 1
      IF( j1>M1 ) THEN
        Icos(I3+j3+1:I3+j3-j2+M2+1) = Icos(I2+j2:I2+M2)
        RETURN
      END IF
    ELSE
      Icos(I3+j3) = Icos(I2+j2)
      j2 = j2 + 1
      IF( j2>M2 ) THEN
        Icos(I3+j3+1:I3+j3-j1+M1+1) = Icos(I1+j1:I1+M1)
        RETURN
      END IF
    END IF
    j3 = j3 + 1
  END DO
END SUBROUTINE I1MERG
