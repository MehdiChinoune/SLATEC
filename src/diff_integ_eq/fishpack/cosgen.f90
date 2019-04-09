!** COSGEN
SUBROUTINE COSGEN(N,Ijump,Fnum,Fden,A)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to GENBUN
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (COSGEN-S, CMPCSG-C)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     This subroutine computes required cosine values in ascending
  !     order.  When IJUMP .GT. 1 the routine computes values
  !
  !        2*COS(J*PI/L), J=1,2,...,L and J .NE. 0(MOD N/IJUMP+1)
  !
  !     where L = IJUMP*(N/IJUMP+1).
  !
  !
  !     when IJUMP = 1 it computes
  !
  !            2*COS((J-FNUM)*PI/(N+FDEN)),  J=1, 2, ... ,N
  !
  !     where
  !        FNUM = 0.5, FDEN = 0.0, for regular reduction values.
  !        FNUM = 0.0, FDEN = 1.0, for B-R and C-R when ISTAG = 1
  !        FNUM = 0.0, FDEN = 0.5, for B-R and C-R when ISTAG = 2
  !        FNUM = 0.5, FDEN = 0.5, for B-R and C-R when ISTAG = 2
  !                                in POISN2 only.
  !
  !***
  ! **See also:**  GENBUN
  !***
  ! **Routines called:**  PIMACH

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)

  REAL A(*), dum, Fden, Fnum, pi, pibyn, PIMACH, x, y
  INTEGER i, Ijump, k, k1, k2, k3, k4, k5, N, np1
  !* FIRST EXECUTABLE STATEMENT  COSGEN
  pi = PIMACH(dum)
  IF ( N/=0 ) THEN
    IF ( Ijump==1 ) THEN
      np1 = N + 1
      y = pi/(N+Fden)
      DO i = 1, N
        x = np1 - i - Fnum
        A(i) = 2.*COS(x*y)
      END DO
    ELSE
      k3 = N/Ijump + 1
      k4 = k3 - 1
      pibyn = pi/(N+Ijump)
      DO k = 1, Ijump
        k1 = (k-1)*k3
        k5 = (k-1)*k4
        DO i = 1, k4
          x = k1 + i
          k2 = k5 + i
          A(k2) = -2.*COS(x*pibyn)
        END DO
      END DO
    END IF
  END IF
END SUBROUTINE COSGEN
