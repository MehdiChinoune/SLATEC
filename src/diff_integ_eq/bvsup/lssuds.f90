!** LSSUDS
PURE SUBROUTINE LSSUDS(A,X,B,N,M,Nrda,U,Nrdu,Iflag,Mlso,Irank,Iscale,Q,Diag,&
    Kpivot,S,Div,Td,Isflg,Scales)
  !> Subsidiary to BVSUP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (LSSUDS-S, DLSSUD-D)
  !***
  ! **Author:**  Watts, H. A., (SNLA)
  !***
  ! **Description:**
  !
  !    LSSUDS solves the underdetermined system of equations  A Z = B,
  !    where A is N by M and N <= M.  In particular, if rank A equals
  !    IRA, a vector X and a matrix U are determined such that X is the
  !    UNIQUE solution of smallest length, satisfying A X = B, and the
  !    columns of U form an orthonormal basis for the null space of A,
  !    satisfying A U = 0 .  Then all solutions Z are given by
  !              Z = X + C(1)*U(1) + ..... + C(M-IRA)*U(M-IRA)
  !    where U(J) represents the J-th column of U and the C(J) are
  !    arbitrary constants.
  !    If the system of equations are not compatible, only the least
  !    squares solution of minimal length is computed.
  !
  !- ********************************************************************
  !   INPUT
  !- ********************************************************************
  !
  !     A -- Contains the matrix of N equations in M unknowns, A remains
  !          unchanged, must be dimensioned NRDA by M.
  !     X -- Solution array of length at least M.
  !     B -- Given constant vector of length N, B remains unchanged.
  !     N -- Number of equations, N greater or equal to 1.
  !     M -- Number of unknowns, M greater or equal to N.
  !     NRDA -- Row dimension of A, NRDA greater or equal to N.
  !     U -- Matrix used for solution, must be dimensioned NRDU by
  !          (M - rank of A).
  !          (storage for U may be ignored when only the minimal length
  !           solution X is desired)
  !     NRDU -- Row dimension of U, NRDU greater or equal to M.
  !             (if only the minimal length solution is wanted,
  !              NRDU=0 is acceptable)
  !     IFLAG -- Status indicator
  !           =0  for the first call (and for each new problem defined by
  !               a new matrix A) when the matrix data is treated as exact
  !           =-K for the first call (and for each new problem defined by
  !               a new matrix A) when the matrix data is assumed to be
  !               accurate to about K digits.
  !           =1  for subsequent calls whenever the matrix A has already
  !               been decomposed (problems with new vectors B but
  !               same matrix A can be handled efficiently).
  !     MLSO -- =0 if only the minimal length solution is wanted.
  !             =1 if the complete solution is wanted, includes the
  !                linear space defined by the matrix U.
  !     IRANK -- Variable used for the rank of A, set by the code.
  !     ISCALE -- Scaling indicator
  !               =-1 if the matrix A is to be pre-scaled by
  !               columns when appropriate.
  !               If the scaling indicator is not equal to -1
  !               no scaling will be attempted.
  !            For most problems scaling will probably not be necessary.
  !     Q -- Matrix used for the transformation, must be dimensioned
  !            NRDA by M.
  !     DIAG,KPIVOT,S, -- Arrays of length at least N used for internal
  !      DIV,TD,SCALES    storage (except for SCALES which is M).
  !     ISFLG -- Storage for an internal variable.
  !
  !- ********************************************************************
  !   OUTPUT
  !- ********************************************************************
  !
  !     IFLAG -- Status indicator
  !            =1 if solution was obtained.
  !            =2 if improper input is detected.
  !            =3 if rank of matrix is less than N.
  !               To continue, simply reset IFLAG=1 and call LSSUDS again.
  !            =4 if the system of equations appears to be inconsistent.
  !               However, the least squares solution of minimal length
  !               was obtained.
  !     X -- Minimal length least squares solution of A Z = B
  !     IRANK -- Numerically determined rank of A, must not be altered
  !              on succeeding calls with input values of IFLAG=1.
  !     U -- Matrix whose M-IRANK columns are mutually orthogonal unit
  !          vectors which span the null space of A. This is to be ignored
  !          when MLSO was set to zero or IFLAG=4 on output.
  !     Q -- Contains the strictly upper triangular part of the reduced
  !           matrix and transformation information.
  !     DIAG -- Contains the diagonal elements of the triangular reduced
  !             matrix.
  !     KPIVOT -- Contains the pivotal information.  The row interchanges
  !               performed on the original matrix are recorded here.
  !     S -- Contains the solution of the lower triangular system.
  !     DIV,TD -- Contains transformation information for rank
  !               deficient problems.
  !     SCALES -- Contains the column scaling parameters.
  !
  !- ********************************************************************
  !
  !***
  ! **See also:**  BVSUP
  !***
  ! **References:**  H. A. Watts, Solving linear least squares problems
  !                 using SODS/SUDS/CODS, Sandia Report SAND77-0683,
  !                 Sandia Laboratories, 1977.
  !***
  ! **Routines called:**  J4SAVE, OHTROL, ORTHOR, R1MACH, SDOT, XERMAX,
  !                    XERMSG, XGETF, XSETF

  !* REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900328  Added TYPE section.  (WRB)
  !   900510  Fixed an error message.  (RWC)
  !   910408  Updated the AUTHOR and REFERENCES sections.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  USE service, ONLY : eps_sp
  !
  INTEGER, INTENT(IN) :: Iscale, M, Mlso, N, Nrda, Nrdu
  INTEGER, INTENT(INOUT) :: Iflag, Irank, Isflg
  INTEGER, INTENT(OUT) :: Kpivot(N)
  REAL(SP), INTENT(IN) :: A(Nrda,M), B(N)
  REAL(SP), INTENT(INOUT) :: Q(Nrda,M), U(Nrdu,M), X(M)
  REAL(SP), INTENT(OUT) :: Diag(N), Div(N), S(N), Scales(M), Td(N)
  !
  INTEGER :: i, irp, j, jr, k, kp, l, nmir, nu, nfat
  REAL(SP) :: gam, gama, res, ss, uro
  !
  !- *********************************************************************
  !
  !     MACHINE PRECISION (COMPUTER UNIT ROUNDOFF VALUE) IS DEFINED
  !     BY THE FUNCTION R1MACH.
  !
  !* FIRST EXECUTABLE STATEMENT  LSSUDS
  uro = eps_sp
  !
  !- *********************************************************************
  !
  IF( N>=1 .AND. M>=N .AND. Nrda>=N ) THEN
    IF( Nrdu==0 .OR. Nrdu>=M ) THEN
      IF( Iflag<=0 ) THEN
        !
        Isflg = -15
        IF( Iflag/=0 ) THEN
          Isflg = Iflag
          nfat = 0
        END IF
        !
        !     COPY MATRIX A INTO MATRIX Q
        !
        DO k = 1, M
          DO j = 1, N
            Q(j,k) = A(j,k)
          END DO
        END DO
        !
        !     USE ORTHOGONAL TRANSFORMATIONS TO REDUCE Q TO LOWER
        !     TRIANGULAR FORM
        !
        CALL ORTHOR(Q,N,M,Nrda,Iflag,Irank,Iscale,Diag,Kpivot,Scales,Div,Td)
        !
        IF( Irank==N ) THEN
          !
          !     STORE DIVISORS FOR THE TRIANGULAR SOLUTION
          !
          DO k = 1, N
            Div(k) = Diag(k)
          END DO
          GOTO 100
        ELSE
          !
          !     FOR RANK DEFICIENT PROBLEMS USE ADDITIONAL ORTHOGONAL
          !     TRANSFORMATIONS TO FURTHER REDUCE Q
          !
          IF( Irank/=0 ) CALL OHTROL(Q,N,Nrda,Diag,Irank,Div,Td)
          RETURN
        END IF
      ELSEIF( Iflag==1 ) THEN
        GOTO 100
      END IF
    END IF
  END IF
  !
  !     INVALID INPUT FOR LSSUDS
  Iflag = 2
  ERROR STOP 'LSSUDS : INVALID INPUT PARAMETERS.'
  RETURN
  !
  !
  100 CONTINUE
  IF( Irank>0 ) THEN
    !
    !     COPY CONSTANT VECTOR INTO S AFTER FIRST INTERCHANGING
    !     THE ELEMENTS ACCORDING TO THE PIVOTAL SEQUENCE
    !
    DO k = 1, N
      kp = Kpivot(k)
      X(k) = B(kp)
    END DO
    DO k = 1, N
      S(k) = X(k)
    END DO
    !
    irp = Irank + 1
    nu = 1
    IF( Mlso==0 ) nu = 0
    IF( Irank/=N ) THEN
      !
      !     FOR RANK DEFICIENT PROBLEMS WE MUST APPLY THE
      !     ORTHOGONAL TRANSFORMATION TO S
      !     WE ALSO CHECK TO SEE IF THE SYSTEM APPEARS TO BE INCONSISTENT
      !
      nmir = N - Irank
      ss = NORM2(S(1:N))**2
      DO l = 1, Irank
        k = irp - l
        gam = ((Td(k)*S(k))+DOT_PRODUCT(Q(irp:N,k),S(irp:N))) &
          /(Td(k)*Div(k))
        S(k) = S(k) + gam*Td(k)
        DO j = irp, N
          S(j) = S(j) + gam*Q(j,k)
        END DO
      END DO
      res = NORM2(S(irp:irp+nmir-1))**2
      IF( res>ss*(10._SP*MAX(10._SP**Isflg,10._SP*uro))**2 ) THEN
        !
        !     INCONSISTENT SYSTEM
        Iflag = 4
        nu = 0
      END IF
    END IF
    !
    !     APPLY FORWARD SUBSTITUTION TO SOLVE LOWER TRIANGULAR SYSTEM
    !
    S(1) = S(1)/Div(1)
    IF( Irank/=1 ) THEN
      DO k = 2, Irank
        S(k) = (S(k)-DOT_PRODUCT(Q(k,1:k-1),S(1:k-1)))/Div(k)
      END DO
    END IF
    !
    !     INITIALIZE X VECTOR AND THEN APPLY ORTHOGONAL TRANSFORMATION
    !
    DO k = 1, M
      X(k) = 0._SP
      IF( k<=Irank ) X(k) = S(k)
    END DO
    !
    DO jr = 1, Irank
      j = irp - jr
      gama = DOT_PRODUCT(Q(j,j:M),X(j:M))/(Diag(j)*Q(j,j))
      DO k = j, M
        X(k) = X(k) + gama*Q(j,k)
      END DO
    END DO
    !
    !     RESCALE ANSWERS AS DICTATED
    !
    DO k = 1, M
      X(k) = X(k)*Scales(k)
    END DO
    !
    IF( (nu==0) .OR. (M==Irank) ) RETURN
    !
    !     INITIALIZE U MATRIX AND THEN APPLY ORTHOGONAL TRANSFORMATION
    !
    l = M - Irank
    DO k = 1, l
      DO i = 1, M
        U(i,k) = 0._SP
        IF( i==Irank+k ) U(i,k) = 1._SP
      END DO
      !
      DO jr = 1, Irank
        j = irp - jr
        gama = DOT_PRODUCT(Q(j,j:M),U(j:M,k))/(Diag(j)*Q(j,j))
        DO i = j, M
          U(i,k) = U(i,k) + gama*Q(j,i)
        END DO
      END DO
    END DO
  ELSE
    !
    !     SPECIAL CASE FOR THE NULL MATRIX
    DO k = 1, M
      X(k) = 0._SP
      IF( Mlso/=0 ) THEN
        U(k,k) = 1._SP
        DO j = 1, M
          IF( j/=k ) U(j,k) = 0._SP
        END DO
      END IF
    END DO
    DO k = 1, N
      IF( B(k)>0. ) Iflag = 4
    END DO
    RETURN
  END IF
  !
END SUBROUTINE LSSUDS