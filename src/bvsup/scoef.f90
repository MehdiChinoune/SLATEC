!** SCOEF
SUBROUTINE SCOEF(Yh,Yp,Ncomp,Nrowb,Nfc,Nic,B,Beta,Coef,Inhomo,Re,Ae,By,&
    Cvec,Work,Iwork,Iflag,Nfcc)
  USE ML, ONLY : EPS
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to BVSUP
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
  !     WORK = Real array of internal storage. Dimension must be .GE.
  !            NFCC*(NFCC+4)
  !     IWORK = Integer array of internal storage. Dimension must be .GE.
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
  !   890921  Realigned order of variables in certain COMMON blocks.
  !           (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)

  INTEGER i, Iflag, Inhomo, Iwork(*), j, k, kflag, ki, l, mlso, Ncomp, ncomp2, &
    nf, Nfc, Nfcc, nfccm1, Nic, Nrowb
  REAL Ae, B(Nrowb,*), bbn, Beta(*), bn, brn, By(Nfcc,*), bykl, bys, Coef(*), cons, &
    Cvec(*), gam, un, Work(*), Yh(Ncomp,*), Yp(*), ypn, SDOT, Re
  !
  !     SET UP MATRIX  B*YH  AND VECTOR  BETA - B*YP
  !
  !* FIRST EXECUTABLE STATEMENT  SCOEF
  ncomp2 = Ncomp/2
  DO k = 1, Nfcc
    DO j = 1, Nfc
      l = j
      IF ( Nfc/=Nfcc ) l = 2*j - 1
      By(k,l) = SDOT(Ncomp,B(k,1),Nrowb,Yh(1,j),1)
    ENDDO
    IF ( Nfc/=Nfcc ) THEN
      DO j = 1, Nfc
        l = 2*j
        bykl = SDOT(ncomp2,B(k,1),Nrowb,Yh(ncomp2+1,j),1)
        By(k,l) = SDOT(ncomp2,B(k,ncomp2+1),Nrowb,Yh(1,j),1) - bykl
      ENDDO
    ENDIF
    SELECT CASE (Inhomo)
      CASE (2)
        !     CASE 2
        Cvec(k) = Beta(k)
      CASE (3)
        !     CASE 3
        Cvec(k) = 0.
      CASE DEFAULT
        !     CASE 1
        Cvec(k) = Beta(k) - SDOT(Ncomp,B(k,1),Nrowb,Yp,1)
    END SELECT
  ENDDO
  cons = ABS(Cvec(1))
  bys = ABS(By(1,1))
  !
  !- *********************************************************************
  !     SOLVE LINEAR SYSTEM
  !
  Iflag = 0
  mlso = 0
  IF ( Inhomo==3 ) mlso = 1
  kflag = INT( 0.5*LOG10(EPS) )
  CALL XGETF(nf)
  CALL XSETF(0)
  DO
    CALL SUDS(By,Coef,Cvec,Nfcc,Nfcc,Nfcc,kflag,mlso,Work,Iwork)
    IF ( kflag/=3 ) THEN
      IF ( kflag==4 ) Iflag = 2
      CALL XSETF(nf)
      IF ( Nfcc==1 ) THEN
        !
        !- *********************************************************************
        !     TESTING FOR EXISTENCE AND UNIQUENESS OF BOUNDARY-VALUE PROBLEM
        !     SOLUTION IN A SCALAR CASE
        !
        bn = 0.
        un = 0.
        ypn = 0.
        DO k = 1, Ncomp
          un = MAX(un,ABS(Yh(k,1)))
          ypn = MAX(ypn,ABS(Yp(k)))
          bn = MAX(bn,ABS(B(1,k)))
        ENDDO
        bbn = MAX(bn,ABS(Beta(1)))
        IF ( bys<=10.*(Re*un+Ae)*bn ) EXIT
        IF ( Inhomo/=3 ) RETURN
        Iflag = 3
        Coef(1) = 1.
        RETURN
      ELSE
        IF ( Inhomo/=3 ) RETURN
        IF ( Iwork(1)<Nfcc ) THEN
          DO k = 1, Nfcc
            ki = 4*Nfcc + k
            Coef(k) = Work(ki)
          ENDDO
          RETURN
        ELSE
          Iflag = 3
          DO k = 1, Nfcc
            Coef(k) = 0.
          ENDDO
          Coef(Nfcc) = 1.
          nfccm1 = Nfcc - 1
          DO k = 1, nfccm1
            j = Nfcc - k
            l = Nfcc - j + 1
            gam = SDOT(l,By(j,j),Nfcc,Coef(j),1)/(Work(j)*By(j,j))
            DO i = j, Nfcc
              Coef(i) = Coef(i) + gam*By(j,i)
            ENDDO
          ENDDO
          RETURN
        ENDIF
      ENDIF
    ELSE
      kflag = 1
      Iflag = 1
    ENDIF
  ENDDO
  brn = bbn/bn*bys
  IF ( cons>=0.1*brn.AND.cons<=10.*brn ) Iflag = 1
  IF ( cons>10.*brn ) Iflag = 2
  IF ( cons<=Re*ABS(Beta(1))+Ae+(Re*ypn+Ae)*bn ) Iflag = 1
  IF ( Inhomo==3 ) Coef(1) = 1.
  RETURN
END SUBROUTINE SCOEF
