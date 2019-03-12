!DECK CMPOSP
SUBROUTINE CMPOSP(M,N,A,Bb,C,Q,Idimq,B,B2,B3,W,W2,W3,D,Tcos,P)
  IMPLICIT NONE
  INTEGER i, Idimq, ipstor, j, lh, M, mr, N, nr, nrm1, nrmj, nrpj
  !***BEGIN PROLOGUE  CMPOSP
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to CMGNBN
  !***LIBRARY   SLATEC
  !***TYPE      COMPLEX (POISP2-S, CMPOSP-C)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !     Subroutine to solve Poisson's equation with periodic boundary
  !     conditions.
  !
  !***SEE ALSO  CMGNBN
  !***ROUTINES CALLED  CMPOSD, CMPOSN
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  CMPOSP
  !
  COMPLEX A, Bb, C, Q, B, B2, B3, W, W2, W3, D, Tcos, P, s, t
  DIMENSION A(*), Bb(*), C(*), Q(Idimq,*), B(*), B2(*), B3(*), W(*), &
    W2(*), W3(*), D(*), Tcos(*), P(*)
  !***FIRST EXECUTABLE STATEMENT  CMPOSP
  mr = M
  nr = (N+1)/2
  nrm1 = nr - 1
  IF ( 2*nr/=N ) THEN
    !
    !     ODD  NUMBER OF UNKNOWNS
    !
    DO j = 1, nrm1
      nrpj = N + 1 - j
      DO i = 1, mr
        s = Q(i,j) - Q(i,nrpj)
        t = Q(i,j) + Q(i,nrpj)
        Q(i,j) = s
        Q(i,nrpj) = t
      ENDDO
    ENDDO
    DO i = 1, mr
      Q(i,nr) = 2.*Q(i,nr)
    ENDDO
    lh = nrm1/2
    DO j = 1, lh
      nrmj = nr - j
      DO i = 1, mr
        s = Q(i,j)
        Q(i,j) = Q(i,nrmj)
        Q(i,nrmj) = s
      ENDDO
    ENDDO
    CALL CMPOSD(mr,nrm1,2,A,Bb,C,Q,Idimq,B,W,D,Tcos,P)
    ipstor = INT(W(1))
    CALL CMPOSN(mr,nr,2,1,A,Bb,C,Q(1,nr),Idimq,B,B2,B3,W,W2,W3,D,Tcos,P)
    ipstor = MAX(ipstor,INT(REAL(W(1))))
    DO j = 1, nrm1
      nrpj = nr + j
      DO i = 1, mr
        s = .5*(Q(i,nrpj)+Q(i,j))
        t = .5*(Q(i,nrpj)-Q(i,j))
        Q(i,nrpj) = t
        Q(i,j) = s
      ENDDO
    ENDDO
    DO i = 1, mr
      Q(i,nr) = .5*Q(i,nr)
    ENDDO
    DO j = 1, lh
      nrmj = nr - j
      DO i = 1, mr
        s = Q(i,j)
        Q(i,j) = Q(i,nrmj)
        Q(i,nrmj) = s
      ENDDO
    ENDDO
  ELSE
    !
    !     EVEN NUMBER OF UNKNOWNS
    !
    DO j = 1, nrm1
      nrmj = nr - j
      nrpj = nr + j
      DO i = 1, mr
        s = Q(i,nrmj) - Q(i,nrpj)
        t = Q(i,nrmj) + Q(i,nrpj)
        Q(i,nrmj) = s
        Q(i,nrpj) = t
      ENDDO
    ENDDO
    DO i = 1, mr
      Q(i,nr) = 2.*Q(i,nr)
      Q(i,N) = 2.*Q(i,N)
    ENDDO
    CALL CMPOSD(mr,nrm1,1,A,Bb,C,Q,Idimq,B,W,D,Tcos,P)
    ipstor = INT(W(1))
    CALL CMPOSN(mr,nr+1,1,1,A,Bb,C,Q(1,nr),Idimq,B,B2,B3,W,W2,W3,D,Tcos,P)
    ipstor = MAX(ipstor,INT(REAL(W(1))))
    DO j = 1, nrm1
      nrmj = nr - j
      nrpj = nr + j
      DO i = 1, mr
        s = .5*(Q(i,nrpj)+Q(i,nrmj))
        t = .5*(Q(i,nrpj)-Q(i,nrmj))
        Q(i,nrmj) = s
        Q(i,nrpj) = t
      ENDDO
    ENDDO
    DO i = 1, mr
      Q(i,nr) = .5*Q(i,nr)
      Q(i,N) = .5*Q(i,N)
    ENDDO
  ENDIF
  !
  !     RETURN STORAGE REQUIREMENTS FOR P VECTORS.
  !
  W(1) = CMPLX(REAL(ipstor),0.)
END SUBROUTINE CMPOSP
