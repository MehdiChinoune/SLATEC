!*==CDPSC.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK CDPSC
SUBROUTINE CDPSC(Ksgn,N,Nq,Yh)
  IMPLICIT NONE
  !*--CDPSC5
  !***BEGIN PROLOGUE  CDPSC
  !***SUBSIDIARY
  !***PURPOSE  Subroutine CDPSC computes the predicted YH values by
  !            effectively multiplying the YH array by the Pascal triangle
  !            matrix when KSGN is +1, and performs the inverse function
  !            when KSGN is -1.
  !***LIBRARY   SLATEC (SDRIVE)
  !***TYPE      COMPLEX (SDPSC-S, DDPSC-D, CDPSC-C)
  !***AUTHOR  Kahaner, D. K., (NIST)
  !             National Institute of Standards and Technology
  !             Gaithersburg, MD  20899
  !           Sutherland, C. D., (LANL)
  !             Mail Stop D466
  !             Los Alamos National Laboratory
  !             Los Alamos, NM  87545
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   900329  Initial submission to SLATEC.
  !***END PROLOGUE  CDPSC
  INTEGER i , j , j1 , j2 , Ksgn , N , Nq
  COMPLEX Yh(N,*)
  !***FIRST EXECUTABLE STATEMENT  CDPSC
  IF ( Ksgn>0 ) THEN
    DO j1 = 1 , Nq
      DO j2 = j1 , Nq
        j = Nq - j2 + j1
        DO i = 1 , N
          Yh(i,j) = Yh(i,j) + Yh(i,j+1)
        ENDDO
      ENDDO
    ENDDO
  ELSE
    DO j1 = 1 , Nq
      DO j2 = j1 , Nq
        j = Nq - j2 + j1
        DO i = 1 , N
          Yh(i,j) = Yh(i,j) - Yh(i,j+1)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
END SUBROUTINE CDPSC
