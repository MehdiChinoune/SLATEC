!*==BLKTRI.f90  processed by SPAG 6.72Dc at 10:55 on  6 Feb 2019
!DECK BLKTRI
SUBROUTINE BLKTRI(Iflg,Np,N,An,Bn,Cn,Mp,M,Am,Bm,Cm,Idimy,Y,Ierror,W)
  IMPLICIT NONE
  !*--BLKTRI5
  !*** Start of declarations inserted by SPAG
  REAL Am, An, Bm, Bn, Cm, Cn, CNV, CPROD, CPRODP, EPS, PROD, &
    PRODP, W, Y
  INTEGER Idimy, Ierror, Iflg, IK, iw1, iw2, iw3, iwah, iwbh, iwd, &
    iwu, iww, K, M, Mp, N, NCMplx, nh, nl, NM
  INTEGER Np, NPP
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  BLKTRI
  !***PURPOSE  Solve a block tridiagonal system of linear equations
  !            (usually resulting from the discretization of separable
  !            two-dimensional elliptic equations).
  !***LIBRARY   SLATEC (FISHPACK)
  !***CATEGORY  I2B4B
  !***TYPE      SINGLE PRECISION (BLKTRI-S, CBLKTR-C)
  !***KEYWORDS  ELLIPTIC PDE, FISHPACK, TRIDIAGONAL LINEAR SYSTEM
  !***AUTHOR  Adams, J., (NCAR)
  !           Swarztrauber, P. N., (NCAR)
  !           Sweet, R., (NCAR)
  !***DESCRIPTION
  !
  !     Subroutine BLKTRI Solves a System of Linear Equations of the Form
  !
  !          AN(J)*X(I,J-1) + AM(I)*X(I-1,J) + (BN(J)+BM(I))*X(I,J)
  !
  !          + CN(J)*X(I,J+1) + CM(I)*X(I+1,J) = Y(I,J)
  !
  !               for I = 1,2,...,M  and  J = 1,2,...,N.
  !
  !     I+1 and I-1 are evaluated modulo M and J+1 and J-1 modulo N, i.e.,
  !
  !          X(I,0) = X(I,N),  X(I,N+1) = X(I,1),
  !          X(0,J) = X(M,J),  X(M+1,J) = X(1,J).
  !
  !     These equations usually result from the discretization of
  !     separable elliptic equations.  Boundary conditions may be
  !     Dirichlet, Neumann, or Periodic.
  !
  !
  !     * * * * * * * * * *     ON INPUT     * * * * * * * * * *
  !
  !     IFLG
  !       = 0  Initialization only.  Certain quantities that depend on NP,
  !            N, AN, BN, and CN are computed and stored in the work
  !            array  W.
  !       = 1  The quantities that were computed in the initialization are
  !            used to obtain the solution X(I,J).
  !
  !       NOTE   A call with IFLG=0 takes approximately one half the time
  !              as a call with IFLG = 1  .  However, the
  !              initialization does not have to be repeated unless NP, N,
  !              AN, BN, or CN change.
  !
  !     NP
  !       = 0  If AN(1) and CN(N) are not zero, which corresponds to
  !            periodic boundary conditions.
  !       = 1  If AN(1) and CN(N) are zero.
  !
  !     N
  !       The number of unknowns in the J-direction. N must be greater
  !       than 4. The operation count is proportional to MNlog2(N), hence
  !       N should be selected less than or equal to M.
  !
  !     AN,BN,CN
  !       One-dimensional arrays of length N that specify the coefficients
  !       in the linear equations given above.
  !
  !     MP
  !       = 0  If AM(1) and CM(M) are not zero, which corresponds to
  !            periodic boundary conditions.
  !       = 1  If AM(1) = CM(M) = 0  .
  !
  !     M
  !       The number of unknowns in the I-direction. M must be greater
  !       than 4.
  !
  !     AM,BM,CM
  !       One-dimensional arrays of length M that specify the coefficients
  !       in the linear equations given above.
  !
  !     IDIMY
  !       The row (or first) dimension of the two-dimensional array Y as
  !       it appears in the program calling BLKTRI.  This parameter is
  !       used to specify the variable dimension of Y.  IDIMY must be at
  !       least M.
  !
  !     Y
  !       A two-dimensional array that specifies the values of the right
  !       side of the linear system of equations given above.  Y must be
  !       dimensioned at least M*N.
  !
  !     W
  !       A one-dimensional array that must be provided by the user for
  !       work space.
  !             If NP=1 define K=INT(log2(N))+1 and set L=2**(K+1) then
  !                     W must have dimension (K-2)*L+K+5+MAX(2N,6M)
  !
  !             If NP=0 define K=INT(log2(N-1))+1 and set L=2**(K+1) then
  !                     W must have dimension (K-2)*L+K+5+2N+MAX(2N,6M)
  !
  !       **IMPORTANT** For purposes of checking, the required dimension
  !                     of W is computed by BLKTRI and stored in W(1)
  !                     in floating point format.
  !
  !     * * * * * * * * * *     On Output     * * * * * * * * * *
  !
  !     Y
  !       Contains the solution X.
  !
  !     IERROR
  !       An error flag that indicates invalid input parameters.  Except
  !       for number zero, a solution is not attempted.
  !
  !       = 0  No error.
  !       = 1  M is less than 5.
  !       = 2  N is less than 5.
  !       = 3  IDIMY is less than M.
  !       = 4  BLKTRI failed while computing results that depend on the
  !            coefficient arrays AN, BN, CN.  Check these arrays.
  !       = 5  AN(J)*CN(J-1) is less than 0 for some J. Possible reasons
  !            for this condition are
  !            1. The arrays AN and CN are not correct.
  !            2. Too large a grid spacing was used in the discretization
  !               of the elliptic equation.
  !            3. The linear equations resulted from a partial
  !               differential equation which was not elliptic.
  !
  !     W
  !       Contains intermediate values that must not be destroyed if
  !       BLKTRI will be called again with IFLG=1.  W(1) contains the
  !       number of locations required by W in floating point format.
  !
  ! *Long Description:
  !
  !     * * * * * * *   Program Specifications    * * * * * * * * * * * *
  !
  !     Dimension of   AN(N),BN(N),CN(N),AM(M),BM(M),CM(M),Y(IDIMY,N)
  !     Arguments      W(See argument list)
  !
  !     Latest         June 1979
  !     Revision
  !
  !     Required       BLKTRI,BLKTRI,PROD,PRODP,CPROD,CPRODP,COMPB,INDXA,
  !     Subprograms    INDXB,INDXC,PPADD,PSGF,PPSGF,PPSPF,BSRH,TEVLS,
  !                    R1MACH
  !
  !     Special        The Algorithm may fail if ABS(BM(I)+BN(J)) is less
  !     Conditions     than ABS(AM(I))+ABS(AN(J))+ABS(CM(I))+ABS(CN(J))
  !                    for some I and J. The Algorithm will also fail if
  !                    AN(J)*CN(J-1) is less than zero for some J.
  !                    See the description of the output parameter IERROR.
  !
  !     Common         CBLKT
  !     Blocks
  !
  !     I/O            None
  !
  !     Precision      Single
  !
  !     Specialist     Paul Swarztrauber
  !
  !     Language       FORTRAN
  !
  !     History        Version 1 September 1973
  !                    Version 2 April     1976
  !                    Version 3 June      1979
  !
  !     Algorithm      Generalized Cyclic Reduction (See Reference below)
  !
  !     Space
  !     Required       Control Data 7600
  !
  !     Portability    American National Standards Institute Fortran.
  !                    The machine accuracy is set using function R1MACH.
  !
  !     Required       None
  !     Resident
  !     Routines
  !
  !     References     Swarztrauber,P. and R. Sweet, 'Efficient FORTRAN
  !                    Subprograms For The Solution Of Elliptic Equations'
  !                    NCAR TN/IA-109, July, 1975, 138 PP.
  !
  !                    Swarztrauber P. ,'A Direct Method For The Discrete
  !                    Solution Of Separable Elliptic Equations', S.I.A.M.
  !                    J. Numer. Anal.,11(1974) PP. 1136-1150.
  !
  !***REFERENCES  P. N. Swarztrauber and R. Sweet, Efficient Fortran
  !                 subprograms for the solution of elliptic equations,
  !                 NCAR TN/IA-109, July 1975, 138 pp.
  !               P. N. Swarztrauber, A direct method for the discrete
  !                 solution of separable elliptic equations, SIAM Journal
  !                 on Numerical Analysis 11, (1974), pp. 1136-1150.
  !***ROUTINES CALLED  BLKTR1, COMPB, CPROD, CPRODP, PROD, PRODP
  !***COMMON BLOCKS    CBLKT
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  BLKTRI
  !
  DIMENSION An(*), Bn(*), Cn(*), Am(*), Bm(*), Cm(*), Y(Idimy,*), &
    W(*)
  EXTERNAL PROD, PRODP, CPROD, CPRODP
  COMMON /CBLKT / NPP, K, EPS, CNV, NM, NCMplx, IK
  !***FIRST EXECUTABLE STATEMENT  BLKTRI
  NM = N
  Ierror = 0
  IF ( M<5 ) THEN
    Ierror = 1
  ELSEIF ( NM<3 ) THEN
    Ierror = 2
  ELSEIF ( Idimy<M ) THEN
    Ierror = 3
  ELSE
    nh = N
    NPP = Np
    IF ( NPP/=0 ) nh = nh + 1
    IK = 2
    K = 1
    DO
      IK = IK + IK
      K = K + 1
      IF ( nh<=IK ) THEN
        nl = IK
        IK = IK + IK
        nl = nl - 1
        iwah = (K-2)*IK + K + 6
        IF ( NPP/=0 ) THEN
          !
          !     DIVIDE W INTO WORKING SUB ARRAYS
          !
          iw1 = iwah
          iwbh = iw1 + NM
          W(1) = iw1 - 1 + MAX(2*NM,6*M)
        ELSE
          iwbh = iwah + NM + NM
          iw1 = iwbh
          W(1) = iw1 - 1 + MAX(2*NM,6*M)
          NM = NM - 1
        ENDIF
        !
        ! SUBROUTINE COMP B COMPUTES THE ROOTS OF THE B POLYNOMIALS
        !
        IF ( Ierror==0 ) THEN
          iw2 = iw1 + M
          iw3 = iw2 + M
          iwd = iw3 + M
          iww = iwd + M
          iwu = iww + M
          IF ( Iflg==0 ) THEN
            CALL COMPB(nl,Ierror,An,Bn,Cn,W(2),W(iwah),W(iwbh))
          ELSEIF ( Mp/=0 ) THEN
            !
            ! SUBROUTINE BLKTR1 SOLVES THE LINEAR SYSTEM
            !
            CALL BLKTR1(nl,An,Bn,Cn,M,Am,Bm,Cm,Idimy,Y,W(2),W(iw1),W(iw2),&
              W(iw3),W(iwd),W(iww),W(iwu),PROD,CPROD)
          ELSE
            CALL BLKTR1(nl,An,Bn,Cn,M,Am,Bm,Cm,Idimy,Y,W(2),W(iw1),W(iw2),&
              W(iw3),W(iwd),W(iww),W(iwu),PRODP,CPRODP)
          ENDIF
        ENDIF
        EXIT
      ENDIF
    ENDDO
  ENDIF
END SUBROUTINE BLKTRI
