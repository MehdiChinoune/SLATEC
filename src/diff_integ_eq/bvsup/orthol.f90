!** ORTHOL
SUBROUTINE ORTHOL(A,M,N,Nrda,Iflag,Irank,Iscale,Diag,Kpivot,Scales,Cols,Cs)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to BVSUP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (ORTHOL-S)
  !***
  ! **Author:**  Watts, H. A., (SNLA)
  !***
  ! **Description:**
  !
  !   Reduction of the matrix A to upper triangular form by a sequence of
  !   orthogonal HOUSEHOLDER transformations pre-multiplying A
  !
  !   Modeled after the ALGOL codes in the articles in the REFERENCES
  !   section.
  !
  !- *********************************************************************
  !   INPUT
  !- *********************************************************************
  !
  !     A -- Contains the matrix to be decomposed, must be dimensioned
  !           NRDA by N
  !     M -- Number of rows in the matrix, M greater or equal to N
  !     N -- Number of columns in the matrix, N greater or equal to 1
  !     IFLAG -- Indicates the uncertainty in the matrix data
  !             = 0 when the data is to be treated as exact
  !             =-K when the data is assumed to be accurate to about
  !                 K digits
  !     ISCALE -- Scaling indicator
  !               =-1 if the matrix A is to be pre-scaled by
  !               columns when appropriate.
  !               Otherwise no scaling will be attempted
  !     NRDA -- Row dimension of A, NRDA greater or equal to M
  !     DIAG,KPIVOT,COLS -- Arrays of length at least n used internally
  !         ,CS,SCALES
  !
  !- *********************************************************************
  !   OUTPUT
  !- *********************************************************************
  !
  !     IFLAG - Status indicator
  !            =1 for successful decomposition
  !            =2 if improper input is detected
  !            =3 if rank of the matrix is less than N
  !     A -- Contains the reduced matrix in the strictly upper triangular
  !          part and transformation information in the lower part
  !     IRANK -- Contains the numerically determined matrix rank
  !     DIAG -- Contains the diagonal elements of the reduced
  !             triangular matrix
  !     KPIVOT -- Contains the pivotal information, the column
  !               interchanges performed on the original matrix are
  !               recorded here.
  !     SCALES -- Contains the column scaling parameters
  !
  !- *********************************************************************
  !
  !***
  ! **See also:**  BVSUP
  !***
  ! **References:**  G. Golub, Numerical methods for solving linear least
  !                 squares problems, Numerische Mathematik 7, (1965),
  !                 pp. 206-216.
  !               P. Businger and G. Golub, Linear least squares
  !                 solutions by Householder transformations, Numerische
  !                 Mathematik  7, (1965), pp. 269-276.
  !***
  ! **Routines called:**  CSCALE, R1MACH, SDOT, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900402  Added TYPE section.  (WRB)
  !   910408  Updated the AUTHOR and REFERENCES sections.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER Iflag, Irank, Iscale, j, jcol, k, kp, Kpivot(*), l, M, mk, N, Nrda
  REAL A(Nrda,*), acc, akk, anorm, as, asave, Cols(*), Cs(*), css, Diag(*), diagk, &
    dum(1), R1MACH, sad, sc, Scales(*), SDOT, sig, sigma, sruro, uro
  !
  !- *********************************************************************
  !
  !     MACHINE PRECISION (COMPUTER UNIT ROUNDOFF VALUE) IS DEFINED
  !     BY THE FUNCTION R1MACH.
  !
  !* FIRST EXECUTABLE STATEMENT  ORTHOL
  uro = R1MACH(3)
  dum = 0.
  !
  !- *********************************************************************
  !
  IF ( M>=N.AND.N>=1.AND.Nrda>=M ) THEN
    !
    acc = 10.*uro
    IF ( Iflag<0 ) acc = MAX(acc,10.**Iflag)
    sruro = SQRT(uro)
    Iflag = 1
    Irank = N
    !
    !     COMPUTE NORM**2 OF JTH COLUMN AND A MATRIX NORM
    !
    anorm = 0.
    DO j = 1, N
      Kpivot(j) = j
      Cols(j) = SDOT(M,A(1,j),1,A(1,j),1)
      Cs(j) = Cols(j)
      anorm = anorm + Cols(j)
    END DO
    !
    !     PERFORM COLUMN SCALING ON A WHEN SPECIFIED
    !
    CALL CSCALE(A,Nrda,M,N,Cols,Cs,dum,dum,anorm,Scales,Iscale,0)
    !
    anorm = SQRT(anorm)
    !
    !
    !     CONSTRUCTION OF UPPER TRIANGULAR MATRIX AND RECORDING OF
    !     ORTHOGONAL TRANSFORMATIONS
    !
    !
    DO k = 1, N
      mk = M - k + 1
      IF ( k/=N ) THEN
        kp = k + 1
        !
        !        SEARCHING FOR PIVOTAL COLUMN
        !
        DO j = k, N
          IF ( Cols(j)<sruro*Cs(j) ) THEN
            Cols(j) = SDOT(mk,A(k,j),1,A(k,j),1)
            Cs(j) = Cols(j)
          END IF
          IF ( j/=k ) THEN
            IF ( sigma>=0.99*Cols(j) ) CYCLE
          END IF
          sigma = Cols(j)
          jcol = j
        END DO
        IF ( jcol/=k ) THEN
          !
          !        PERFORM COLUMN INTERCHANGE
          !
          l = Kpivot(k)
          Kpivot(k) = Kpivot(jcol)
          Kpivot(jcol) = l
          Cols(jcol) = Cols(k)
          Cols(k) = sigma
          css = Cs(k)
          Cs(k) = Cs(jcol)
          Cs(jcol) = css
          sc = Scales(k)
          Scales(k) = Scales(jcol)
          Scales(jcol) = sc
          DO l = 1, M
            asave = A(l,k)
            A(l,k) = A(l,jcol)
            A(l,jcol) = asave
          END DO
        END IF
      END IF
      !
      !        CHECK RANK OF THE MATRIX
      !
      sig = SDOT(mk,A(k,k),1,A(k,k),1)
      diagk = SQRT(sig)
      IF ( diagk>acc*anorm ) THEN
        !
        !        CONSTRUCT AND APPLY TRANSFORMATION TO MATRIX A
        !
        akk = A(k,k)
        IF ( akk>0. ) diagk = -diagk
        Diag(k) = diagk
        A(k,k) = akk - diagk
        IF ( k/=N ) THEN
          sad = diagk*akk - sig
          DO j = kp, N
            as = SDOT(mk,A(k,k),1,A(k,j),1)/sad
            DO l = k, M
              A(l,j) = A(l,j) + as*A(l,k)
            END DO
            Cols(j) = Cols(j) - A(k,j)**2
          END DO
        END IF
      ELSE
        !
        !        RANK DEFICIENT PROBLEM
        Iflag = 3
        Irank = k - 1
        CALL XERMSG('SLATEC','ORTHOL',&
          'RANK OF MATRIX IS LESS THAN THE NUMBER OF COLUMNS.',1,1)
        RETURN
      END IF
    END DO
    RETURN
  END IF
  Iflag = 2
  CALL XERMSG('SLATEC','ORTHOL','INVALID INPUT PARAMETERS.',2,1)
  RETURN
  !
  !
  RETURN
END SUBROUTINE ORTHOL
