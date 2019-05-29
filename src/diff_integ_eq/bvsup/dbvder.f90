!** DBVDER
SUBROUTINE DBVDER(X,Y,Yp,G)
  !>
  !  Subsidiary to DBVSUP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (BVDER-S, DBVDER-D)
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
  !             NOFST communicates to the routine DFMAT how to access
  !             the dependent variables corresponding to this initial
  !             value problem.  For example, during any call to DFMAT,
  !             the first dependent variable for the initial value
  !             problem is in position  Y(NOFST + 1).
  !             See example in SAND77-1328.
  !- *********************************************************************
  !
  !***
  ! **See also:**  DBVSUP
  !***
  ! **Routines called:**  (NONE)
  !***
  ! COMMON BLOCKS    DML8SZ, DMLIVP

  !* REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890921  Realigned order of variables in certain COMMON blocks.
  !           (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910701  Corrected ROUTINES CALLED section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !   920618  Minor restructuring of code.  (RWC, WRB)
  USE DML, ONLY : NOFst, C, XSAv, IGOfx, INHomo, IVP, NCOmp, NFC
  REAL(8) :: X, G(*), Y(*), Yp(*)
  INTEGER :: j, k, l, na
  !- *********************************************************************
  !
  !* FIRST EXECUTABLE STATEMENT  DBVDER
  IF ( IVP>0 ) STOP
  NOFst = IVP
  na = 1
  DO k = 1, NFC
    CALL DFMAT(X,Y(na),Yp(na))
    NOFst = NOFst - NCOmp
    na = na + NCOmp
  END DO
  !
  IF ( INHomo/=1 ) RETURN
  CALL DFMAT(X,Y(na),Yp(na))
  !
  IF ( IGOfx==0 ) RETURN
  IF ( X/=XSAv ) THEN
    IF ( IVP==0 ) CALL DGVEC(X,G)
    IF ( IVP>0 ) STOP
    XSAv = X
  END IF
  !
  !     If the user has chosen not to normalize the particular
  !     solution, then C is defined in DBVPOR to be 1.0
  !
  !     The following loop is just
  !     CALL DAXPY (NCOMP, 1.0D0/C, G, 1, YP(NA), 1)
  !
  DO j = 1, NCOmp
    l = na + j - 1
    Yp(l) = Yp(l) + G(j)/C
  END DO
END SUBROUTINE DBVDER
