!DECK DPPFA
SUBROUTINE DPPFA(Ap,N,Info)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  DPPFA
  !***PURPOSE  Factor a real symmetric positive definite matrix stored in
  !            packed form.
  !***LIBRARY   SLATEC (LINPACK)
  !***CATEGORY  D2B1B
  !***TYPE      DOUBLE PRECISION (SPPFA-S, DPPFA-D, CPPFA-C)
  !***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION, PACKED,
  !             POSITIVE DEFINITE
  !***AUTHOR  Moler, C. B., (U. of New Mexico)
  !***DESCRIPTION
  !
  !     DPPFA factors a double precision symmetric positive definite
  !     matrix stored in packed form.
  !
  !     DPPFA is usually called by DPPCO, but it can be called
  !     directly with a saving in time if  RCOND  is not needed.
  !     (time for DPPCO) = (1 + 18/N)*(time for DPPFA) .
  !
  !     On Entry
  !
  !        AP      DOUBLE PRECISION (N*(N+1)/2)
  !                the packed form of a symmetric matrix  A .  The
  !                columns of the upper triangle are stored sequentially
  !                in a one-dimensional array of length  N*(N+1)/2 .
  !                See comments below for details.
  !
  !        N       INTEGER
  !                the order of the matrix  A .
  !
  !     On Return
  !
  !        AP      an upper triangular matrix  R, stored in packed
  !                form, so that  A = TRANS(R)*R .
  !
  !        INFO    INTEGER
  !                = 0  for normal return.
  !                = K  if the leading minor of order  K  is not
  !                     positive definite.
  !
  !
  !     Packed Storage
  !
  !          The following program segment will pack the upper
  !          triangle of a symmetric matrix.
  !
  !                K = 0
  !                DO 20 J = 1, N
  !                   DO 10 I = 1, J
  !                      K = K + 1
  !                      AP(K) = A(I,J)
  !             10    CONTINUE
  !             20 CONTINUE
  !
  !***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***ROUTINES CALLED  DDOT
  !***REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  DPPFA
  INTEGER N, Info
  REAL(8) :: Ap(*)
  !
  REAL(8) :: DDOT, t
  REAL(8) :: s
  INTEGER j, jj, jm1, k, kj, kk
  !***FIRST EXECUTABLE STATEMENT  DPPFA
  jj = 0
  DO j = 1, N
    Info = j
    s = 0.0D0
    jm1 = j - 1
    kj = jj
    kk = 0
    IF ( jm1>=1 ) THEN
      DO k = 1, jm1
        kj = kj + 1
        t = Ap(kj) - DDOT(k-1,Ap(kk+1),1,Ap(jj+1),1)
        kk = kk + k
        t = t/Ap(kk)
        Ap(kj) = t
        s = s + t*t
      ENDDO
    ENDIF
    jj = jj + j
    s = Ap(jj) - s
    IF ( s<=0.0D0 ) RETURN
    Ap(jj) = SQRT(s)
  ENDDO
  Info = 0
  RETURN
END SUBROUTINE DPPFA
