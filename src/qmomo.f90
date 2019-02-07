!*==QMOMO.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK QMOMO
SUBROUTINE QMOMO(Alfa,Beta,Ri,Rj,Rg,Rh,Integr)
  IMPLICIT NONE
  !*--QMOMO5
  !***BEGIN PROLOGUE  QMOMO
  !***PURPOSE  This routine computes modified Chebyshev moments.  The K-th
  !            modified Chebyshev moment is defined as the integral over
  !            (-1,1) of W(X)*T(K,X), where T(K,X) is the Chebyshev
  !            polynomial of degree K.
  !***LIBRARY   SLATEC (QUADPACK)
  !***CATEGORY  H2A2A1, C3A2
  !***TYPE      SINGLE PRECISION (QMOMO-S, DQMOMO-D)
  !***KEYWORDS  MODIFIED CHEBYSHEV MOMENTS, QUADPACK, QUADRATURE
  !***AUTHOR  Piessens, Robert
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !           de Doncker, Elise
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !***DESCRIPTION
  !
  !        MODIFIED CHEBYSHEV MOMENTS
  !        STANDARD FORTRAN SUBROUTINE
  !        REAL VERSION
  !
  !        PARAMETERS
  !           ALFA   - Real
  !                    Parameter in the weight function W(X), ALFA.GT.(-1)
  !
  !           BETA   - Real
  !                    Parameter in the weight function W(X), BETA.GT.(-1)
  !
  !           RI     - Real
  !                    Vector of dimension 25
  !                    RI(K) is the integral over (-1,1) of
  !                    (1+X)**ALFA*T(K-1,X), K = 1, ..., 25.
  !
  !           RJ     - Real
  !                    Vector of dimension 25
  !                    RJ(K) is the integral over (-1,1) of
  !                    (1-X)**BETA*T(K-1,X), K = 1, ..., 25.
  !
  !           RG     - Real
  !                    Vector of dimension 25
  !                    RG(K) is the integral over (-1,1) of
  !                    (1+X)**ALFA*LOG((1+X)/2)*T(K-1,X), K = 1, ..., 25.
  !
  !           RH     - Real
  !                    Vector of dimension 25
  !                    RH(K) is the integral over (-1,1) of
  !                    (1-X)**BETA*LOG((1-X)/2)*T(K-1,X), K = 1, ..., 25.
  !
  !           INTEGR - Integer
  !                    Input parameter indicating the modified
  !                    Moments to be computed
  !                    INTEGR = 1 compute RI, RJ
  !                           = 2 compute RI, RJ, RG
  !                           = 3 compute RI, RJ, RH
  !                           = 4 compute RI, RJ, RG, RH
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   810101  DATE WRITTEN
  !   891009  Removed unreferenced statement label.  (WRB)
  !   891009  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  QMOMO
  !
  REAL Alfa, alfp1, alfp2, an, anm1, Beta, betp1, betp2, ralf, &
    rbet, Rg, Rh, Ri, Rj
  INTEGER i, im1, Integr
  !
  DIMENSION Rg(25), Rh(25), Ri(25), Rj(25)
  !
  !
  !***FIRST EXECUTABLE STATEMENT  QMOMO
  alfp1 = Alfa + 0.1E+01
  betp1 = Beta + 0.1E+01
  alfp2 = Alfa + 0.2E+01
  betp2 = Beta + 0.2E+01
  ralf = 0.2E+01**alfp1
  rbet = 0.2E+01**betp1
  !
  !           COMPUTE RI, RJ USING A FORWARD RECURRENCE RELATION.
  !
  Ri(1) = ralf/alfp1
  Rj(1) = rbet/betp1
  Ri(2) = Ri(1)*Alfa/alfp2
  Rj(2) = Rj(1)*Beta/betp2
  an = 0.2E+01
  anm1 = 0.1E+01
  DO i = 3, 25
    Ri(i) = -(ralf+an*(an-alfp2)*Ri(i-1))/(anm1*(an+alfp1))
    Rj(i) = -(rbet+an*(an-betp2)*Rj(i-1))/(anm1*(an+betp1))
    anm1 = an
    an = an + 0.1E+01
  ENDDO
  IF ( Integr/=1 ) THEN
    IF ( Integr/=3 ) THEN
      !
      !           COMPUTE RG USING A FORWARD RECURRENCE RELATION.
      !
      Rg(1) = -Ri(1)/alfp1
      Rg(2) = -(ralf+ralf)/(alfp2*alfp2) - Rg(1)
      an = 0.2E+01
      anm1 = 0.1E+01
      im1 = 2
      DO i = 3, 25
        Rg(i) = -(an*(an-alfp2)*Rg(im1)-an*Ri(im1)+anm1*Ri(i))&
          /(anm1*(an+alfp1))
        anm1 = an
        an = an + 0.1E+01
        im1 = i
      ENDDO
      IF ( Integr==2 ) GOTO 100
    ENDIF
    !
    !           COMPUTE RH USING A FORWARD RECURRENCE RELATION.
    !
    Rh(1) = -Rj(1)/betp1
    Rh(2) = -(rbet+rbet)/(betp2*betp2) - Rh(1)
    an = 0.2E+01
    anm1 = 0.1E+01
    im1 = 2
    DO i = 3, 25
      Rh(i) = -(an*(an-betp2)*Rh(im1)-an*Rj(im1)+anm1*Rj(i))&
        /(anm1*(an+betp1))
      anm1 = an
      an = an + 0.1E+01
      im1 = i
    ENDDO
    DO i = 2, 25, 2
      Rh(i) = -Rh(i)
    ENDDO
  ENDIF
  100 CONTINUE
  DO i = 2, 25, 2
    Rj(i) = -Rj(i)
  ENDDO
END SUBROUTINE QMOMO
