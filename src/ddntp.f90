!*==DDNTP.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DDNTP
SUBROUTINE DDNTP(H,K,N,Nq,T,Tout,Yh,Y)
  IMPLICIT NONE
  !*--DDNTP5
  !***BEGIN PROLOGUE  DDNTP
  !***SUBSIDIARY
  !***PURPOSE  Subroutine DDNTP interpolates the K-th derivative of Y at
  !            TOUT, using the data in the YH array.  If K has a value
  !            greater than NQ, the NQ-th derivative is calculated.
  !***LIBRARY   SLATEC (SDRIVE)
  !***TYPE      DOUBLE PRECISION (SDNTP-S, DDNTP-D, CDNTP-C)
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
  !***END PROLOGUE  DDNTP
  INTEGER i , j , jj , K , kk , kused , N , Nq
  DOUBLE PRECISION factor , H , r , T , Tout , Y(*) , Yh(N,*)
  !***FIRST EXECUTABLE STATEMENT  DDNTP
  IF ( K==0 ) THEN
    DO i = 1 , N
      Y(i) = Yh(i,Nq+1)
    ENDDO
    r = ((Tout-T)/H)
    DO jj = 1 , Nq
      j = Nq + 1 - jj
      DO i = 1 , N
        Y(i) = Yh(i,j) + r*Y(i)
      ENDDO
    ENDDO
  ELSE
    kused = MIN(K,Nq)
    factor = 1.D0
    DO kk = 1 , kused
      factor = factor*(Nq+1-kk)
    ENDDO
    DO i = 1 , N
      Y(i) = factor*Yh(i,Nq+1)
    ENDDO
    r = ((Tout-T)/H)
    DO jj = kused + 1 , Nq
      j = kused + 1 + Nq - jj
      factor = 1.D0
      DO kk = 1 , kused
        factor = factor*(j-kk)
      ENDDO
      DO i = 1 , N
        Y(i) = factor*Yh(i,j) + r*Y(i)
      ENDDO
    ENDDO
    DO i = 1 , N
      Y(i) = Y(i)*H**(-kused)
    ENDDO
  ENDIF
END SUBROUTINE DDNTP
