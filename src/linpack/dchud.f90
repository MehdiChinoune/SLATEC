!DECK DCHUD
SUBROUTINE DCHUD(R,Ldr,P,X,Z,Ldz,Nz,Y,Rho,C,S)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  DCHUD
  !***PURPOSE  Update an augmented Cholesky decomposition of the
  !            triangular part of an augmented QR decomposition.
  !***LIBRARY   SLATEC (LINPACK)
  !***CATEGORY  D7B
  !***TYPE      DOUBLE PRECISION (SCHUD-S, DCHUD-D, CCHUD-C)
  !***KEYWORDS  CHOLESKY DECOMPOSITION, LINEAR ALGEBRA, LINPACK, MATRIX,
  !             UPDATE
  !***AUTHOR  Stewart, G. W., (U. of Maryland)
  !***DESCRIPTION
  !
  !     DCHUD updates an augmented Cholesky decomposition of the
  !     triangular part of an augmented QR decomposition.  Specifically,
  !     given an upper triangular matrix R of order P, a row vector
  !     X, a column vector Z, and a scalar Y, DCHUD determines a
  !     unitary matrix U and a scalar ZETA such that
  !
  !
  !                              (R  Z)     (RR   ZZ )
  !                         U  * (    )  =  (        ) ,
  !                              (X  Y)     ( 0  ZETA)
  !
  !     where RR is upper triangular.  If R and Z have been
  !     obtained from the factorization of a least squares
  !     problem, then RR and ZZ are the factors corresponding to
  !     the problem with the observation (X,Y) appended.  In this
  !     case, if RHO is the norm of the residual vector, then the
  !     norm of the residual vector of the updated problem is
  !     SQRT(RHO**2 + ZETA**2).  DCHUD will simultaneously update
  !     several triplets (Z,Y,RHO).
  !     For a less terse description of what DCHUD does and how
  !     it may be applied, see the LINPACK guide.
  !
  !     The matrix U is determined as the product U(P)*...*U(1),
  !     where U(I) is a rotation in the (I,P+1) plane of the
  !     form
  !
  !                       (     C(I)      S(I) )
  !                       (                    ) .
  !                       (    -S(I)      C(I) )
  !
  !     The rotations are chosen so that C(I) is double precision.
  !
  !     On Entry
  !
  !         R      DOUBLE PRECISION(LDR,P), where LDR .GE. P.
  !                R contains the upper triangular matrix
  !                that is to be updated.  The part of R
  !                below the diagonal is not referenced.
  !
  !         LDR    INTEGER.
  !                LDR is the leading dimension of the array R.
  !
  !         P      INTEGER.
  !                P is the order of the matrix R.
  !
  !         X      DOUBLE PRECISION(P).
  !                X contains the row to be added to R.  X is
  !                not altered by DCHUD.
  !
  !         Z      DOUBLE PRECISION(LDZ,N)Z), where LDZ .GE. P.
  !                Z is an array containing NZ P-vectors to
  !                be updated with R.
  !
  !         LDZ    INTEGER.
  !                LDZ is the leading dimension of the array Z.
  !
  !         NZ     INTEGER.
  !                NZ is the number of vectors to be updated
  !                NZ may be zero, in which case Z, Y, and RHO
  !                are not referenced.
  !
  !         Y      DOUBLE PRECISION(NZ).
  !                Y contains the scalars for updating the vectors
  !                Z.  Y is not altered by DCHUD.
  !
  !         RHO    DOUBLE PRECISION(NZ).
  !                RHO contains the norms of the residual
  !                vectors that are to be updated.  If RHO(J)
  !                is negative, it is left unaltered.
  !
  !     On Return
  !
  !         RC
  !         RHO    contain the updated quantities.
  !         Z
  !
  !         C      DOUBLE PRECISION(P).
  !                C contains the cosines of the transforming
  !                rotations.
  !
  !         S      DOUBLE PRECISION(P).
  !                S contains the sines of the transforming
  !                rotations.
  !
  !***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***ROUTINES CALLED  DROTG
  !***REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  DCHUD
  INTEGER Ldr, P, Ldz, Nz
  REAL(8) :: Rho(*), C(*)
  REAL(8) :: R(Ldr,*), X(*), Z(Ldz,*), Y(*), S(*)
  !
  INTEGER i, j, jm1
  REAL(8) :: azeta, scale
  REAL(8) :: t, xj, zeta
  !
  !     UPDATE R.
  !
  !***FIRST EXECUTABLE STATEMENT  DCHUD
  DO j = 1, P
    xj = X(j)
    !
    !        APPLY THE PREVIOUS ROTATIONS.
    !
    jm1 = j - 1
    IF ( jm1>=1 ) THEN
      DO i = 1, jm1
        t = C(i)*R(i,j) + S(i)*xj
        xj = C(i)*xj - S(i)*R(i,j)
        R(i,j) = t
      ENDDO
    ENDIF
    !
    !        COMPUTE THE NEXT ROTATION.
    !
    CALL DROTG(R(j,j),xj,C(j),S(j))
  ENDDO
  !
  !     IF REQUIRED, UPDATE Z AND RHO.
  !
  IF ( Nz>=1 ) THEN
    DO j = 1, Nz
      zeta = Y(j)
      DO i = 1, P
        t = C(i)*Z(i,j) + S(i)*zeta
        zeta = C(i)*zeta - S(i)*Z(i,j)
        Z(i,j) = t
      ENDDO
      azeta = ABS(zeta)
      IF ( azeta/=0.0D0.AND.Rho(j)>=0.0D0 ) THEN
        scale = azeta + Rho(j)
        Rho(j) = scale*SQRT((azeta/scale)**2+(Rho(j)/scale)**2)
      ENDIF
    ENDDO
  ENDIF
END SUBROUTINE DCHUD
