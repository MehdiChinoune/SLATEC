!** SCOEF
PURE SUBROUTINE SCOEF(Yh,Yp,Ncomp,Nrowb,Nfc,B,Beta,Coef,Inhomo,Re,Ae,By,&
    Cvec,Work,Iwork,Iflag,Nfcc)
  !> Subsidiary to BVSUP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (SCOEF-S, DCOEF-D)
  !***
  ! **Author:**  Watts, H. A., (SNLA)
  !***
  ! **Description:**
  !
  !- *********************************************************************
  ! INPUT TO SCOEF
  !- *********************************************************************
  !
  !     YH = Matrix of homogeneous solutions.
  !     YP = Vector containing particular solution.
  !     NCOMP = Number of components per solution vector.
  !     NROWB = First dimension of B in calling program.
  !     NFC = Number of base solution vectors.
  !     NFCC = 2*NFC for the special treatment of complex valued
  !            equations. Otherwise, NFCC=NFC.
  !     NIC = Number of specified initial conditions.
  !     B = Boundary condition matrix at X = Xfinal.
  !     BETA = Vector of nonhomogeneous boundary conditions at X = Xfinal.
  !              1 - Nonzero particular solution
  !     INHOMO = 2 - Zero particular solution
  !              3 - Eigenvalue problem
  !     RE = Relative error tolerance
  !     AE = Absolute error tolerance
  !     BY = Storage space for the matrix  B*YH
  !     CVEC = Storage space for the vector  BETA-B*YP
  !     WORK = Real array of internal storage. Dimension must be >=
  !            NFCC*(NFCC+4)
  !     IWORK = Integer array of internal storage. Dimension must be >=
  !             3+NFCC
  !
  !- *********************************************************************
  ! OUTPUT FROM SCOEF
  !- *********************************************************************
  !
  !     COEF = Array containing superposition constants.
  !     IFLAG = Indicator of success from SUDS in solving the
  !             boundary equations
  !           = 0 Boundary equations are solved
  !           = 1 Boundary equations appear to have many solutions
  !           = 2 Boundary equations appear to be inconsistent
  !           = 3 For this value of an eigenparameter, the boundary
  !               equations have only the zero solution.
  !
  !- *********************************************************************
  !
  !     Subroutine SCOEF solves for the superposition constants from the
  !     linear equations defined by the boundary conditions at X = Xfinal.
  !
  !                          B*YP + B*YH*COEF = BETA
  !
  !- *********************************************************************
  !
  !***
  ! **See also:**  BVSUP
  !***
  ! **Routines called:**  SDOT, SUDS, XGETF, XSETF
  !***
  ! COMMON BLOCKS    ML5MCO

  !* REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890921  Realigned order of variables in certain COMMON blocks.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  USE ML, ONLY : eps_com
  !
  INTEGER, INTENT(IN) :: Inhomo, Ncomp, Nfc, Nfcc, Nrowb
  INTEGER, INTENT(OUT) :: Iflag
  INTEGER, INTENT(INOUT) :: Iwork(*)
  REAL(SP), INTENT(IN) :: Ae, Re
  REAL(SP), INTENT(IN) :: B(Nrowb,Ncomp), Beta(Nrowb), Yh(Ncomp,Nfcc), Yp(Ncomp)
  REAL(SP), INTENT(INOUT) :: Work(*)
  REAL(SP), INTENT(OUT) :: By(Nfcc,Ncomp), Cvec(Nrowb), Coef(Nfcc)
  !
  INTEGER :: i, j, k, kflag, ki, l, mlso, ncomp2, nfccm1
  REAL(SP) :: bbn, bn, brn, bykl, bys, cons, gam, un, ypn
  !
  !     SET UP MATRIX  B*YH  AND VECTOR  BETA - B*YP
  !
  !* FIRST EXECUTABLE STATEMENT  SCOEF
  ncomp2 = Ncomp/2
  DO k = 1, Nfcc
    DO j = 1, Nfc
      l = j
      IF( Nfc/=Nfcc ) l = 2*j - 1
      By(k,l) = DOT_PRODUCT(B(k,1:Ncomp),Yh(1:Ncomp,j))
    END DO
    IF( Nfc/=Nfcc ) THEN
      DO j = 1, Nfc
        l = 2*j
        bykl = DOT_PRODUCT(B(k,1:ncomp2),Yh(ncomp2+1:2*ncomp2,j))
        By(k,l) = DOT_PRODUCT(B(k,ncomp2+1:2*ncomp2),Yh(1,j:j+ncomp2-1)) - bykl
      END DO
    END IF
    SELECT CASE (Inhomo)
      CASE (2)
        !     CASE 2
        Cvec(k) = Beta(k)
      CASE (3)
        !     CASE 3
        Cvec(k) = 0._SP
      CASE DEFAULT
        !     CASE 1
        Cvec(k) = Beta(k) - DOT_PRODUCT(B(k,1:Ncomp),Yp(1:Ncomp))
    END SELECT
  END DO
  cons = ABS(Cvec(1))
  bys = ABS(By(1,1))
  !
  !- *********************************************************************
  !     SOLVE LINEAR SYSTEM
  !
  Iflag = 0
  mlso = 0
  IF( Inhomo==3 ) mlso = 1
  kflag = INT( 0.5_SP*LOG10(eps_com) )
  DO
    CALL SUDS(By,Coef,Cvec,Nfcc,Nfcc,Nfcc,kflag,mlso,Work,Iwork)
    IF( kflag/=3 ) THEN
      IF( kflag==4 ) Iflag = 2
      IF( Nfcc==1 ) THEN
        !
        !- *********************************************************************
        !     TESTING FOR EXISTENCE AND UNIQUENESS OF BOUNDARY-VALUE PROBLEM
        !     SOLUTION IN A SCALAR CASE
        !
        bn = 0._SP
        un = 0._SP
        ypn = 0._SP
        DO k = 1, Ncomp
          un = MAX(un,ABS(Yh(k,1)))
          ypn = MAX(ypn,ABS(Yp(k)))
          bn = MAX(bn,ABS(B(1,k)))
        END DO
        bbn = MAX(bn,ABS(Beta(1)))
        IF( bys<=10._SP*(Re*un+Ae)*bn ) EXIT
        IF( Inhomo/=3 ) RETURN
        Iflag = 3
        Coef(1) = 1._SP
        RETURN
      ELSE
        IF( Inhomo/=3 ) RETURN
        IF( Iwork(1)<Nfcc ) THEN
          DO k = 1, Nfcc
            ki = 4*Nfcc + k
            Coef(k) = Work(ki)
          END DO
          RETURN
        ELSE
          Iflag = 3
          DO k = 1, Nfcc
            Coef(k) = 0._SP
          END DO
          Coef(Nfcc) = 1._SP
          nfccm1 = Nfcc - 1
          DO k = 1, nfccm1
            j = Nfcc - k
            gam = DOT_PRODUCT(By(j,j:Nfcc),Coef(j:Nfcc))/(Work(j)*By(j,j))
            DO i = j, Nfcc
              Coef(i) = Coef(i) + gam*By(j,i)
            END DO
          END DO
          RETURN
        END IF
      END IF
    ELSE
      kflag = 1
      Iflag = 1
    END IF
  END DO
  brn = bbn/bn*bys
  IF( cons>=0.1*brn .AND. cons<=10._SP*brn ) Iflag = 1
  IF( cons>10._SP*brn ) Iflag = 2
  IF( cons<=Re*ABS(Beta(1))+Ae+(Re*ypn+Ae)*bn ) Iflag = 1
  IF( Inhomo==3 ) Coef(1) = 1._SP
  !
  RETURN
END SUBROUTINE SCOEF