!** DBVDER
SUBROUTINE DBVDER(X,Y,Yp,G)
  !> Subsidiary to DBVSUP
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
  !           the original boundary value problem   ( if  NEQIVP > 0 )
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
  USE DML, ONLY : nofst_com, c_com, xsav_com, igofx_com, inhomo_com, ivp_com, &
    ncomp_com, nfc_com
  REAL(DP) :: X, Y(:), Yp(:), G(:)
  INTEGER :: j, k, l, na
  !- *********************************************************************
  !
  !* FIRST EXECUTABLE STATEMENT  DBVDER
  IF( ivp_com>0 ) STOP
  nofst_com = ivp_com
  na = 1
  DO k = 1, nfc_com
    CALL DFMAT(X,Y(na:),Yp(na:))
    nofst_com = nofst_com - ncomp_com
    na = na + ncomp_com
  END DO
  !
  IF( inhomo_com/=1 ) RETURN
  CALL DFMAT(X,Y(na:),Yp(na:))
  !
  IF( igofx_com==0 ) RETURN
  IF( X/=xsav_com ) THEN
    IF( ivp_com==0 ) CALL DGVEC(X,G)
    IF( ivp_com>0 ) STOP
    xsav_com = X
  END IF
  !
  !     If the user has chosen not to normalize the particular
  !     solution, then C is defined in DBVPOR to be 1.0
  !
  !     The following loop is just
  !     CALL DAXPY (NCOMP, 1.0D0/C, G, 1, YP(NA), 1)
  !
  DO j = 1, ncomp_com
    l = na + j - 1
    Yp(l) = Yp(l) + G(j)/c_com
  END DO
END SUBROUTINE DBVDER
