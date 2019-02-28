!DECK DSLVS
SUBROUTINE DSLVS(Wm,Iwm,X,Tem)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  DSLVS
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DDEBDF
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (SLVS-S, DSLVS-D)
  !***AUTHOR  Watts, H. A., (SNLA)
  !***DESCRIPTION
  !
  !   DSLVS solves the linear system in the iteration scheme for the
  !   integrator package DDEBDF.
  !
  !***SEE ALSO  DDEBDF
  !***ROUTINES CALLED  DGBSL, DGESL
  !***COMMON BLOCKS    DDEBD1
  !***REVISION HISTORY  (YYMMDD)
  !   820301  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !   920422  Changed DIMENSION statement.  (WRB)
  !***END PROLOGUE  DSLVS
  !
  INTEGER i, IER, IOWnd, IOWns, Iwm, JSTart, KFLag, L, MAXord, &
    meband, METh, MITer, ml, mu, N, NFE, NJE, NQ, NQU, NST
  REAL(8) :: di, EL0, H, hl0, HMIn, HMXi, HU, phl0, r, &
    ROWnd, ROWns, Tem, TN, UROund, Wm, X
  DIMENSION Wm(*), Iwm(*), X(*), Tem(*)
  COMMON /DDEBD1/ ROWnd, ROWns(210), EL0, H, HMIn, HMXi, HU, TN, &
    UROund, IOWnd(14), IOWns(6), IER, JSTart, KFLag, L, &
    METh, MITer, MAXord, N, NQ, NST, NFE, NJE, NQU
  !     ------------------------------------------------------------------
  !      THIS ROUTINE MANAGES THE SOLUTION OF THE LINEAR SYSTEM ARISING
  !      FROM A CHORD ITERATION.  IT IS CALLED BY DSTOD  IF MITER .NE. 0.
  !      IF MITER IS 1 OR 2, IT CALLS DGESL TO ACCOMPLISH THIS.
  !      IF MITER = 3 IT UPDATES THE COEFFICIENT H*EL0 IN THE DIAGONAL
  !      MATRIX, AND THEN COMPUTES THE SOLUTION.
  !      IF MITER IS 4 OR 5, IT CALLS DGBSL.
  !      COMMUNICATION WITH DSLVS USES THE FOLLOWING VARIABLES..
  !      WM  = DOUBLE PRECISION WORK SPACE CONTAINING THE INVERSE DIAGONAL
  !      MATRIX IF MITER
  !            IS 3 AND THE LU DECOMPOSITION OF THE MATRIX OTHERWISE.
  !            STORAGE OF MATRIX ELEMENTS STARTS AT WM(3).
  !            WM ALSO CONTAINS THE FOLLOWING MATRIX-RELATED DATA..
  !            WM(1) = SQRT(UROUND) (NOT USED HERE),
  !            WM(2) = HL0, THE PREVIOUS VALUE OF H*EL0, USED IF MITER =
  !            3.
  !      IWM = INTEGER WORK SPACE CONTAINING PIVOT INFORMATION, STARTING
  !            AT IWM(21), IF MITER IS 1, 2, 4, OR 5.  IWM ALSO CONTAINS
  !            THE BAND PARAMETERS ML = IWM(1) AND MU = IWM(2) IF MITER IS
  !            4 OR 5.
  !      X   = THE RIGHT-HAND SIDE VECTOR ON INPUT, AND THE SOLUTION
  !            VECTOR ON OUTPUT, OF LENGTH N.
  !      TEM = VECTOR OF WORK SPACE OF LENGTH N, NOT USED IN THIS VERSION.
  !      IER = OUTPUT FLAG (IN COMMON).  IER = 0 IF NO TROUBLE OCCURRED.
  !            IER = -1 IF A SINGULAR MATRIX AROSE WITH MITER = 3.
  !      THIS ROUTINE ALSO USES THE COMMON VARIABLES EL0, H, MITER, AND N.
  !-----------------------------------------------------------------------
  !     BEGIN BLOCK PERMITTING ...EXITS TO 80
  !        BEGIN BLOCK PERMITTING ...EXITS TO 60
  !***FIRST EXECUTABLE STATEMENT  DSLVS
  IER = 0
  SELECT CASE (MITer)
    CASE (3)
      !
      phl0 = Wm(2)
      hl0 = H*EL0
      Wm(2) = hl0
      IF ( hl0/=phl0 ) THEN
        r = hl0/phl0
        DO i = 1, N
          di = 1.0D0 - r*(1.0D0-1.0D0/Wm(i+2))
          !        .........EXIT
          IF ( ABS(di)==0.0D0 ) GOTO 100
          Wm(i+2) = 1.0D0/di
        ENDDO
      ENDIF
      DO i = 1, N
        X(i) = Wm(i+2)*X(i)
        !     ......EXIT
      ENDDO
    CASE (4,5)
      !
      ml = Iwm(1)
      mu = Iwm(2)
      meband = 2*ml + mu + 1
      CALL DGBSL(Wm(3),meband,N,ml,mu,Iwm(21),X,0)
    CASE DEFAULT
      !     ......EXIT
      CALL DGESL(Wm(3),N,N,Iwm(21),X,0)
  END SELECT
  RETURN
  !     ...EXIT
  100  IER = -1
  !     ----------------------- END OF SUBROUTINE DSLVS
  !     -----------------------
  RETURN
END SUBROUTINE DSLVS
