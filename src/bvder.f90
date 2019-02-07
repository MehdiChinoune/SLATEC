!*==BVDER.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK BVDER
SUBROUTINE BVDER(X,Y,Yp,G,Ipar)
  IMPLICIT NONE
  !*--BVDER5
  !*** Start of declarations inserted by SPAG
  REAL C, G, X, XSAv, Y, Yp
  INTEGER IGOfx, INHomo, Ipar, IVP, j, k, l, na, NCOmp, NFC, NOFst
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  BVDER
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to BVSUP
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (BVDER-S, DBVDER-D)
  !***AUTHOR  Watts, H. A., (SNLA)
  !***DESCRIPTION
  !
  ! **********************************************************************
  !     NFC = Number of base solution vectors
  !
  !     NCOMP = Number of components per solution vector
  !
  !              1 -- Nonzero particular solution
  !     INHOMO =
  !              2 or 3 -- Zero particular solution
  !
  !             0 -- Inhomogeneous vector term G(X) identically zero
  !     IGOFX =
  !             1 -- Inhomogeneous vector term G(X) not identically zero
  !
  !     G = Inhomogeneous vector term G(X)
  !
  !     XSAV = Previous value of X
  !
  !     C = Normalization factor for the particular solution
  !
  !           0   ( if  NEQIVP = 0 )
  !     IVP =
  !           Number of differential equations integrated due to
  !           the original boundary value problem   ( if  NEQIVP .GT. 0 )
  !
  !     NOFST - For problems with auxiliary initial value equations,
  !             NOFST communicates to the routine FMAT how to access
  !             the dependent variables corresponding to this initial
  !             value problem.  For example, during any call to FMAT,
  !             the first dependent variable for the initial value
  !             problem is in position  Y(NOFST + 1).
  !             See example in SAND77-1328.
  ! **********************************************************************
  !
  !***SEE ALSO  BVSUP
  !***ROUTINES CALLED  (NONE)
  !***COMMON BLOCKS    ML8SZ, MLIVP
  !***REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890921  Realigned order of variables in certain COMMON blocks.
  !           (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910701  Corrected ROUTINES CALLED section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !   920618  Minor restructuring of code.  (RWC, WRB)
  !***END PROLOGUE  BVDER
  DIMENSION Y(*), Yp(*), G(*)
  !
  ! **********************************************************************
  !
  COMMON /ML8SZ / C, XSAv, IGOfx, INHomo, IVP, NCOmp, NFC
  !
  ! **********************************************************************
  !     The COMMON block below is used to communicate with the user
  !     supplied subroutine FMAT.  The user should not alter this
  !     COMMON block.
  !
  COMMON /MLIVP / NOFst
  ! **********************************************************************
  !
  !***FIRST EXECUTABLE STATEMENT  BVDER
  IF ( IVP>0 ) CALL UIVP(X,Y(IVP+1),Yp(IVP+1))
  NOFst = IVP
  na = 1
  DO k = 1, NFC
    CALL FMAT(X,Y(na),Yp(na))
    NOFst = NOFst - NCOmp
    na = na + NCOmp
  ENDDO
  !
  IF ( INHomo/=1 ) RETURN
  CALL FMAT(X,Y(na),Yp(na))
  !
  IF ( IGOfx==0 ) RETURN
  IF ( X/=XSAv ) THEN
    IF ( IVP==0 ) CALL GVEC(X,G)
    IF ( IVP>0 ) CALL UVEC(X,Y(IVP+1),G)
    XSAv = X
  ENDIF
  !
  !     If the user has chosen not to normalize the particular
  !     solution, then C is defined in BVPOR to be 1.0
  !
  !     The following loop is just
  !     CALL SAXPY (NCOMP, 1.0E0/C, G, 1, YP(NA), 1)
  !
  DO j = 1, NCOmp
    l = na + j - 1
    Yp(l) = Yp(l) + G(j)/C
  ENDDO
END SUBROUTINE BVDER
