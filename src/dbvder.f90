!*==DBVDER.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DBVDER
SUBROUTINE DBVDER(X,Y,Yp,G,Ipar)
  IMPLICIT NONE
  !*--DBVDER5
  !***BEGIN PROLOGUE  DBVDER
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DBVSUP
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (BVDER-S, DBVDER-D)
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
  !             NOFST communicates to the routine DFMAT how to access
  !             the dependent variables corresponding to this initial
  !             value problem.  For example, during any call to DFMAT,
  !             the first dependent variable for the initial value
  !             problem is in position  Y(NOFST + 1).
  !             See example in SAND77-1328.
  ! **********************************************************************
  !
  !***SEE ALSO  DBVSUP
  !***ROUTINES CALLED  (NONE)
  !***COMMON BLOCKS    DML8SZ, DMLIVP
  !***REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890921  Realigned order of variables in certain COMMON blocks.
  !           (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910701  Corrected ROUTINES CALLED section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !   920618  Minor restructuring of code.  (RWC, WRB)
  !***END PROLOGUE  DBVDER
  INTEGER IGOfx , INHomo , Ipar , IVP , j , k , l , na , NCOmp , NFC , NOFst
  REAL(8) :: C , G(*) , X , XSAv , Y(*) , Yp(*)
  !
  ! **********************************************************************
  !
  COMMON /DML8SZ/ C , XSAv , IGOfx , INHomo , IVP , NCOmp , NFC
  !
  ! **********************************************************************
  !     The COMMON block below is used to communicate with the user
  !     supplied subroutine DFMAT.  The user should not alter this
  !     COMMON block.
  !
  COMMON /DMLIVP/ NOFst
  ! **********************************************************************
  !
  !***FIRST EXECUTABLE STATEMENT  DBVDER
  IF ( IVP>0 ) CALL DUIVP(X,Y(IVP+1),Yp(IVP+1))
  NOFst = IVP
  na = 1
  DO k = 1 , NFC
    CALL DFMAT(X,Y(na),Yp(na))
    NOFst = NOFst - NCOmp
    na = na + NCOmp
  ENDDO
  !
  IF ( INHomo/=1 ) RETURN
  CALL DFMAT(X,Y(na),Yp(na))
  !
  IF ( IGOfx==0 ) RETURN
  IF ( X/=XSAv ) THEN
    IF ( IVP==0 ) CALL DGVEC(X,G)
    IF ( IVP>0 ) CALL DUVEC(X,Y(IVP+1),G)
    XSAv = X
  ENDIF
  !
  !     If the user has chosen not to normalize the particular
  !     solution, then C is defined in DBVPOR to be 1.0
  !
  !     The following loop is just
  !     CALL DAXPY (NCOMP, 1.0D0/C, G, 1, YP(NA), 1)
  !
  DO j = 1 , NCOmp
    l = na + j - 1
    Yp(l) = Yp(l) + G(j)/C
  ENDDO
END SUBROUTINE DBVDER
