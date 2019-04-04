!** BVDER
SUBROUTINE BVDER(X,Y,Yp,G,Ipar)
  USE ML, ONLY : NOFst, C, XSAv, IGOfx, INHomo, IVP, NCOmp, NFC
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to BVSUP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (BVDER-S, DBVDER-D)
  !***
  ! **Author:**  Watts, H. A., (SNLA)
  !***
  ! **Description:**
  !
  !- *********************************************************************
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
  !- *********************************************************************
  !
  !***
  ! **See also:**  BVSUP
  !***
  ! **Routines called:**  (NONE)
  !***
  ! COMMON BLOCKS    ML8SZ, MLIVP

  !* REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890921  Realigned order of variables in certain COMMON blocks.
  !           (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910701  Corrected ROUTINES CALLED section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !   920618  Minor restructuring of code.  (RWC, WRB)

  REAL G(*), Y(*), Yp(*), X
  INTEGER Ipar, j, k, l, na
  !* FIRST EXECUTABLE STATEMENT  BVDER
  IF ( IVP>0 ) CALL UIVP(X,Y(IVP+1),Yp(IVP+1))
  NOFst = IVP
  na = 1
  DO k = 1, NFC
    CALL FMAT(X,Y(na),Yp(na))
    NOFst = NOFst - NCOmp
    na = na + NCOmp
  END DO
  !
  IF ( INHomo/=1 ) RETURN
  CALL FMAT(X,Y(na),Yp(na))
  !
  IF ( IGOfx==0 ) RETURN
  IF ( X/=XSAv ) THEN
    IF ( IVP==0 ) CALL GVEC(X,G)
    IF ( IVP>0 ) CALL UVEC(X,Y(IVP+1),G)
    XSAv = X
  END IF
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
  END DO
END SUBROUTINE BVDER
