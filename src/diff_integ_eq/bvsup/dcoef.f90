!** DCOEF
SUBROUTINE DCOEF(Yh,Yp,Ncomp,Nrowb,Nfc,Nic,B,Beta,Coef,Inhomo,Re,Ae,By,&
    Cvec,Work,Iwork,Iflag,Nfcc)
  USE DML, ONLY : EPS
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to DBVSUP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (SCOEF-S, DCOEF-D)
  !***
  ! **Author:**  Watts, H. A., (SNLA)
  !***
  ! **Description:**
  !
  !- *********************************************************************
  ! INPUT to DCOEF
  !- *********************************************************************
  !
  !     YH = matrix of homogeneous solutions.
  !     YP = vector containing particular solution.
  !     NCOMP = number of components per solution vector.
  !     NROWB = first dimension of B in calling program.
  !     NFC = number of base solution vectors.
  !     NFCC = 2*NFC for the special treatment of COMPLEX*16 valued
  !            equations. Otherwise, NFCC=NFC.
  !     NIC = number of specified initial conditions.
  !     B = boundary condition matrix at X = XFINAL.
  !     BETA = vector of nonhomogeneous boundary conditions at X = XFINAL.
  !              1 - nonzero particular solution
  !     INHOMO = 2 - zero particular solution
  !              3 - eigenvalue problem
  !     RE = relative error tolerance.
  !     AE = absolute error tolerance.
  !     BY = storage space for the matrix  B*YH
  !     CVEC = storage space for the vector  BETA-B*YP
  !     WORK = double precision array of internal storage. Dimension must
  !     be GE
  !            NFCC*(NFCC+4)
  !     IWORK = integer array of internal storage. Dimension must be GE
  !             3+NFCC
  !
  !- *********************************************************************
  ! OUTPUT from DCOEF
  !- *********************************************************************
  !
  !     COEF = array containing superposition constants.
  !     IFLAG = indicator of success from DSUDS in solving the
  !             boundary equations.
  !           = 0 boundary equations are solved.
  !           = 1 boundary equations appear to have many solutions.
  !           = 2 boundary equations appear to be inconsistent.
  !           = 3 for this value of an eigenparameter, the boundary
  !               equations have only the zero solution.
  !
  !- *********************************************************************
  !
  !     Subroutine DCOEF solves for the superposition constants from the
  !     linear equations defined by the boundary conditions at X = XFINAL.
  !
  !                          B*YP + B*YH*COEF = BETA
  !
  !- *********************************************************************
  !
  !***
  ! **See also:**  DBVSUP
  !***
  ! **Routines called:**  DDOT, DSUDS, XGETF, XSETF
  !***
  ! COMMON BLOCKS    DML5MC

  !* REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890921  Realigned order of variables in certain COMMON blocks.
  !           (WRB)
  !   890921  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)

  !
  REAL(8) :: DDOT
  INTEGER i, Iflag, Inhomo, Iwork(*), j, k, kflag, ki, l, mlso, Ncomp, ncomp2, &
    nf, Nfc, Nfcc, nfccm1, Nic, Nrowb
  REAL(8) :: Ae, B(Nrowb,*), bbn, Beta(*), bn, brn, By(Nfcc,*), bykl, bys, &
    Coef(*), cons, Cvec(*), gam, Re, un, Work(*), Yh(Ncomp,*), Yp(*), ypn
  !* FIRST EXECUTABLE STATEMENT  DCOEF
  !
  !     SET UP MATRIX  B*YH  AND VECTOR  BETA - B*YP
  !
  ncomp2 = Ncomp/2
  DO k = 1, Nfcc
    DO j = 1, Nfc
      l = j
      IF ( Nfc/=Nfcc ) l = 2*j - 1
      By(k,l) = DDOT(Ncomp,B(k,1),Nrowb,Yh(1,j),1)
    END DO
    IF ( Nfc/=Nfcc ) THEN
      DO j = 1, Nfc
        l = 2*j
        bykl = DDOT(ncomp2,B(k,1),Nrowb,Yh(ncomp2+1,j),1)
        By(k,l) = DDOT(ncomp2,B(k,ncomp2+1),Nrowb,Yh(1,j),1) - bykl
      END DO
    END IF
    SELECT CASE (Inhomo)
      CASE (2)
        !        CASE 2
        Cvec(k) = Beta(k)
      CASE (3)
        !        CASE 3
        Cvec(k) = 0.0D0
      CASE DEFAULT
        !        CASE 1
        Cvec(k) = Beta(k) - DDOT(Ncomp,B(k,1),Nrowb,Yp,1)
    END SELECT
  END DO
  cons = ABS(Cvec(1))
  bys = ABS(By(1,1))
  !
  !     ******************************************************************
  !         SOLVE LINEAR SYSTEM
  !
  Iflag = 0
  mlso = 0
  IF ( Inhomo==3 ) mlso = 1
  kflag = INT( 0.5D0*LOG10(EPS) )
  CALL XGETF(nf)
  CALL XSETF(0)
  DO
    CALL DSUDS(By,Coef,Cvec,Nfcc,Nfcc,Nfcc,kflag,mlso,Work,Iwork)
    IF ( kflag/=3 ) THEN
      IF ( kflag==4 ) Iflag = 2
      CALL XSETF(nf)
      IF ( Nfcc==1 ) THEN
        !
        !        ***************************************************************
        !            TESTING FOR EXISTENCE AND UNIQUENESS OF BOUNDARY-VALUE
        !            PROBLEM SOLUTION IN A SCALAR CASE
        !
        bn = 0.0D0
        un = 0.0D0
        ypn = 0.0D0
        DO k = 1, Ncomp
          un = MAX(un,ABS(Yh(k,1)))
          ypn = MAX(ypn,ABS(Yp(k)))
          bn = MAX(bn,ABS(B(1,k)))
        END DO
        bbn = MAX(bn,ABS(Beta(1)))
        IF ( bys<=10.0D0*(Re*un+Ae)*bn ) THEN
          brn = bbn/bn*bys
          IF ( cons>=0.1D0*brn.AND.cons<=10.0D0*brn ) Iflag = 1
          IF ( cons>10.0D0*brn ) Iflag = 2
          IF ( cons<=Re*ABS(Beta(1))+Ae+(Re*ypn+Ae)*bn ) Iflag = 1
          IF ( Inhomo==3 ) Coef(1) = 1.0D0
        ELSEIF ( Inhomo==3 ) THEN
          Iflag = 3
          Coef(1) = 1.0D0
        END IF
      ELSEIF ( Inhomo==3 ) THEN
        IF ( Iwork(1)<Nfcc ) THEN
          DO k = 1, Nfcc
            ki = 4*Nfcc + k
            Coef(k) = Work(ki)
          END DO
        ELSE
          Iflag = 3
          DO k = 1, Nfcc
            Coef(k) = 0.0D0
          END DO
          Coef(Nfcc) = 1.0D0
          nfccm1 = Nfcc - 1
          DO k = 1, nfccm1
            j = Nfcc - k
            l = Nfcc - j + 1
            gam = DDOT(l,By(j,j),Nfcc,Coef(j),1)/(Work(j)*By(j,j))
            DO i = j, Nfcc
              Coef(i) = Coef(i) + gam*By(j,i)
            END DO
          END DO
        END IF
      END IF
      EXIT
    ELSE
      kflag = 1
      Iflag = 1
    END IF
  END DO
END SUBROUTINE DCOEF