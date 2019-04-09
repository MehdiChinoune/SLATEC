!** CDPSC
SUBROUTINE CDPSC(Ksgn,N,Nq,Yh)
  IMPLICIT NONE
  !>
  !***
  !  Subroutine CDPSC computes the predicted YH values by
  !            effectively multiplying the YH array by the Pascal triangle
  !            matrix when KSGN is +1, and performs the inverse function
  !            when KSGN is -1.
  !***
  ! **Library:**   SLATEC (SDRIVE)
  !***
  ! **Type:**      COMPLEX (SDPSC-S, DDPSC-D, CDPSC-C)
  !***
  ! **Author:**  Kahaner, D. K., (NIST)
  !             National Institute of Standards and Technology
  !             Gaithersburg, MD  20899
  !           Sutherland, C. D., (LANL)
  !             Mail Stop D466
  !             Los Alamos National Laboratory
  !             Los Alamos, NM  87545
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   900329  Initial submission to SLATEC.

  INTEGER i, j, j1, j2, Ksgn, N, Nq
  COMPLEX Yh(N,*)
  !* FIRST EXECUTABLE STATEMENT  CDPSC
  IF ( Ksgn>0 ) THEN
    DO j1 = 1, Nq
      DO j2 = j1, Nq
        j = Nq - j2 + j1
        DO i = 1, N
          Yh(i,j) = Yh(i,j) + Yh(i,j+1)
        END DO
      END DO
    END DO
  ELSE
    DO j1 = 1, Nq
      DO j2 = j1, Nq
        j = Nq - j2 + j1
        DO i = 1, N
          Yh(i,j) = Yh(i,j) - Yh(i,j+1)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE CDPSC
