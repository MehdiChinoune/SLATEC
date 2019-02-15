!DECK DCHDD
SUBROUTINE DCHDD(R,Ldr,P,X,Z,Ldz,Nz,Y,Rho,C,S,Info)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  DCHDD
  !***PURPOSE  Downdate an augmented Cholesky decomposition or the
  !            triangular factor of an augmented QR decomposition.
  !***LIBRARY   SLATEC (LINPACK)
  !***CATEGORY  D7B
  !***TYPE      DOUBLE PRECISION (SCHDD-S, DCHDD-D, CCHDD-C)
  !***KEYWORDS  CHOLESKY DECOMPOSITION, DOWNDATE, LINEAR ALGEBRA, LINPACK,
  !             MATRIX
  !***AUTHOR  Stewart, G. W., (U. of Maryland)
  !***DESCRIPTION
  !
  !     DCHDD downdates an augmented Cholesky decomposition or the
  !     triangular factor of an augmented QR decomposition.
  !     Specifically, given an upper triangular matrix R of order P,  a
  !     row vector X, a column vector Z, and a scalar Y, DCHDD
  !     determines an orthogonal matrix U and a scalar ZETA such that
  !
  !                        (R   Z )     (RR  ZZ)
  !                    U * (      )  =  (      ) ,
  !                        (0 ZETA)     ( X   Y)
  !
  !     where RR is upper triangular.  If R and Z have been obtained
  !     from the factorization of a least squares problem, then
  !     RR and ZZ are the factors corresponding to the problem
  !     with the observation (X,Y) removed.  In this case, if RHO
  !     is the norm of the residual vector, then the norm of
  !     the residual vector of the downdated problem is
  !     SQRT(RHO**2 - ZETA**2).  DCHDD will simultaneously downdate
  !     several triplets (Z,Y,RHO) along with R.
  !     For a less terse description of what DCHDD does and how
  !     it may be applied, see the LINPACK guide.
  !
  !     The matrix U is determined as the product U(1)*...*U(P)
  !     where U(I) is a rotation in the (P+1,I)-plane of the
  !     form
  !
  !                       ( C(I)     -S(I)     )
  !                       (                    ) .
  !                       ( S(I)       C(I)    )
  !
  !     The rotations are chosen so that C(I) is double precision.
  !
  !     The user is warned that a given downdating problem may
  !     be impossible to accomplish or may produce
  !     inaccurate results.  For example, this can happen
  !     if X is near a vector whose removal will reduce the
  !     rank of R.  Beware.
  !
  !     On Entry
  !
  !         R      DOUBLE PRECISION(LDR,P), where LDR .GE. P.
  !                R contains the upper triangular matrix
  !                that is to be downdated.  The part of  R
  !                below the diagonal is not referenced.
  !
  !         LDR    INTEGER.
  !                LDR is the leading dimension of the array R.
  !
  !         P      INTEGER.
  !                P is the order of the matrix R.
  !
  !         X      DOUBLE PRECISION(P).
  !                X contains the row vector that is to
  !                be removed from R.  X is not altered by DCHDD.
  !
  !         Z      DOUBLE PRECISION(LDZ,N)Z), where LDZ .GE. P.
  !                Z is an array of NZ P-vectors which
  !                are to be downdated along with R.
  !
  !         LDZ    INTEGER.
  !                LDZ is the leading dimension of the array Z.
  !
  !         NZ     INTEGER.
  !                NZ is the number of vectors to be downdated
  !                NZ may be zero, in which case Z, Y, and RHO
  !                are not referenced.
  !
  !         Y      DOUBLE PRECISION(NZ).
  !                Y contains the scalars for the downdating
  !                of the vectors Z.  Y is not altered by DCHDD.
  !
  !         RHO    DOUBLE PRECISION(NZ).
  !                RHO contains the norms of the residual
  !                vectors that are to be downdated.
  !
  !     On Return
  !
  !         R
  !         Z      contain the downdated quantities.
  !         RHO
  !
  !         C      DOUBLE PRECISION(P).
  !                C contains the cosines of the transforming
  !                rotations.
  !
  !         S      DOUBLE PRECISION(P).
  !                S contains the sines of the transforming
  !                rotations.
  !
  !         INFO   INTEGER.
  !                INFO is set as follows.
  !
  !                   INFO = 0  if the entire downdating
  !                             was successful.
  !
  !                   INFO =-1  if R could not be downdated.
  !                             in this case, all quantities
  !                             are left unaltered.
  !
  !                   INFO = 1  if some RHO could not be
  !                             downdated.  The offending RHO's are
  !                             set to -1.
  !
  !***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***ROUTINES CALLED  DDOT, DNRM2
  !***REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  DCHDD
  INTEGER Ldr, P, Ldz, Nz, Info
  REAL(8) :: R(Ldr,*), X(*), Z(Ldz,*), Y(*), S(*)
  REAL(8) :: Rho(*), C(*)
  !
  INTEGER i, ii, j
  REAL(8) :: a, alpha, azeta, norm, DNRM2
  REAL(8) :: DDOT, t, zeta, b, xx, scale
  !
  !     SOLVE THE SYSTEM TRANS(R)*A = X, PLACING THE RESULT
  !     IN THE ARRAY S.
  !
  !***FIRST EXECUTABLE STATEMENT  DCHDD
  Info = 0
  S(1) = X(1)/R(1,1)
  IF ( P>=2 ) THEN
    DO j = 2, P
      S(j) = X(j) - DDOT(j-1,R(1,j),1,S,1)
      S(j) = S(j)/R(j,j)
    ENDDO
  ENDIF
  norm = DNRM2(P,S,1)
  IF ( norm<1.0D0 ) THEN
    alpha = SQRT(1.0D0-norm**2)
    !
    !        DETERMINE THE TRANSFORMATIONS.
    !
    DO ii = 1, P
      i = P - ii + 1
      scale = alpha + ABS(S(i))
      a = alpha/scale
      b = S(i)/scale
      norm = SQRT(a**2+b**2)
      C(i) = a/norm
      S(i) = b/norm
      alpha = scale*norm
    ENDDO
    !
    !        APPLY THE TRANSFORMATIONS TO R.
    !
    DO j = 1, P
      xx = 0.0D0
      DO ii = 1, j
        i = j - ii + 1
        t = C(i)*xx + S(i)*R(i,j)
        R(i,j) = C(i)*R(i,j) - S(i)*xx
        xx = t
      ENDDO
    ENDDO
    !
    !        IF REQUIRED, DOWNDATE Z AND RHO.
    !
    IF ( Nz>=1 ) THEN
      DO j = 1, Nz
        zeta = Y(j)
        DO i = 1, P
          Z(i,j) = (Z(i,j)-S(i)*zeta)/C(i)
          zeta = C(i)*zeta - S(i)*Z(i,j)
        ENDDO
        azeta = ABS(zeta)
        IF ( azeta<=Rho(j) ) THEN
          Rho(j) = Rho(j)*SQRT(1.0D0-(azeta/Rho(j))**2)
        ELSE
          Info = 1
          Rho(j) = -1.0D0
        ENDIF
      ENDDO
    ENDIF
  ELSE
    Info = -1
  ENDIF
END SUBROUTINE DCHDD
