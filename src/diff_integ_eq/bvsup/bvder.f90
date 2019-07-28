!** BVDER
SUBROUTINE BVDER(X,Y,Yp,G)
  !> Subsidiary to BVSUP
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
  !           the original boundary value problem   ( if  NEQIVP > 0 )
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
  !   890921  Realigned order of variables in certain COMMON blocks.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910701  Corrected ROUTINES CALLED section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !   920618  Minor restructuring of code.  (RWC, WRB)
  USE ML, ONLY : nofst_com, c_com, xsav_com, igofx_com, inhomo_com, ivp_com, &
    ncomp_com, nfc_com
  !
  REAL(SP), INTENT(IN) :: Y(:), X
  REAL(SP), INTENT(OUT) :: G(:), Yp(:)
  !
  INTEGER :: j, k, l, na
  !* FIRST EXECUTABLE STATEMENT  BVDER
  IF( ivp_com>0 ) STOP
  nofst_com = ivp_com
  na = 1
  DO k = 1, nfc_com
    CALL FMAT(X,Y(na:),Yp(na:))
    nofst_com = nofst_com - ncomp_com
    na = na + ncomp_com
  END DO
  !
  IF( inhomo_com/=1 ) RETURN
  CALL FMAT(X,Y(na:),Yp(na:))
  !
  IF( igofx_com==0 ) RETURN
  IF( X/=xsav_com ) THEN
    IF( ivp_com==0 ) CALL GVEC(X,G)
    IF( ivp_com>0 ) STOP
    xsav_com = X
  END IF
  !
  !     If the user has chosen not to normalize the particular
  !     solution, then C is defined in BVPOR to be 1.0
  !
  !     The following loop is just
  !     CALL SAXPY (NCOMP, 1.0E0/C, G, 1, YP(NA), 1)
  !
  DO j = 1, ncomp_com
    l = na + j - 1
    Yp(l) = Yp(l) + G(j)/c_com
  END DO
  !
END SUBROUTINE BVDER