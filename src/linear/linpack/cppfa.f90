!** CPPFA
SUBROUTINE CPPFA(Ap,N,Info)
  !>
  !  Factor a complex Hermitian positive definite matrix stored
  !            in packed form.
  !***
  ! **Library:**   SLATEC (LINPACK)
  !***
  ! **Category:**  D2D1B
  !***
  ! **Type:**      COMPLEX (SPPFA-S, DPPFA-D, CPPFA-C)
  !***
  ! **Keywords:**  LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION, PACKED,
  !             POSITIVE DEFINITE
  !***
  ! **Author:**  Moler, C. B., (U. of New Mexico)
  !***
  ! **Description:**
  !
  !     CPPFA factors a complex Hermitian positive definite matrix
  !     stored in packed form.
  !
  !     CPPFA is usually called by CPPCO, but it can be called
  !     directly with a saving in time if  RCOND  is not needed.
  !     (Time for CPPCO) = (1 + 18/N)*(Time for CPPFA) .
  !
  !     On Entry
  !
  !        AP      COMPLEX (N*(N+1)/2)
  !                the packed form of a Hermitian matrix  A .  The
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
  !                form, so that  A = CTRANS(R)*R .
  !
  !        INFO    INTEGER
  !                = 0  for normal return.
  !                = K  If the leading minor of order  K  is not
  !                     positive definite.
  !
  !
  !     Packed Storage
  !
  !          The following program segment will pack the upper
  !          triangle of a Hermitian matrix.
  !
  !                K = 0
  !                DO 20 J = 1, N
  !                   DO 10 I = 1, J
  !                      K = K + 1
  !                      AP(K) = A(I,J)
  !             10    CONTINUE
  !             20 CONTINUE
  !
  !***
  ! **References:**  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***
  ! **Routines called:**  CDOTC

  !* REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER N, Info
  COMPLEX(SP) Ap(*)
  !
  COMPLEX(SP) t
  REAL(SP) s
  INTEGER j, jj, jm1, k, kj, kk
  !* FIRST EXECUTABLE STATEMENT  CPPFA
  jj = 0
  DO j = 1, N
    Info = j
    s = 0.0E0
    jm1 = j - 1
    kj = jj
    kk = 0
    IF ( jm1>=1 ) THEN
      DO k = 1, jm1
        kj = kj + 1
        t = Ap(kj) - CDOTC(k-1,Ap(kk+1),1,Ap(jj+1),1)
        kk = kk + k
        t = t/Ap(kk)
        Ap(kj) = t
        s = s + REAL(t*CONJG(t))
      END DO
    END IF
    jj = jj + j
    s = REAL(Ap(jj)) - s
    IF ( s<=0.0E0.OR.AIMAG(Ap(jj))/=0.0E0 ) RETURN
    Ap(jj) = CMPLX(SQRT(s),0.0E0)
  END DO
  Info = 0
  RETURN
END SUBROUTINE CPPFA
