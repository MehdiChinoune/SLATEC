!** FDJAC1
SUBROUTINE FDJAC1(FCN,N,X,Fvec,Fjac,Ldfjac,Iflag,Ml,Mu,Epsfcn,Wa1,Wa2)
  !>
  !***
  !  Subsidiary to SNSQ and SNSQE
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (FDJAC1-S, DFDJC1-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     This subroutine computes a forward-difference approximation
  !     to the N by N Jacobian matrix associated with a specified
  !     problem of N functions in N VARIABLES. If the Jacobian has
  !     a banded form, then function evaluations are saved by only
  !     approximating the nonzero terms.
  !
  !     The subroutine statement is
  !
  !       SUBROUTINE FDJAC1(FCN,N,X,FVEC,FJAC,LDFJAC,IFLAG,ML,MU,EPSFCN,
  !                         WA1,WA2)
  !
  !     where
  !
  !       FCN is the name of the user-supplied subroutine which
  !         calculates the functions. FCN must be declared
  !         in an external statement in the user calling
  !         program, and should be written as follows.
  !
  !         SUBROUTINE FCN(N,X,FVEC,IFLAG)
  !         INTEGER N,IFLAG
  !         REAL X(N),FVEC(N)
  !         ----------
  !         Calculate the functions at X and
  !         return this vector in FVEC.
  !         ----------
  !         RETURN
  !         END
  !
  !         The value of IFLAG should not be changed by FCN unless
  !         the user wants to terminate execution of FDJAC1.
  !         In this case set IFLAG to a negative integer.
  !
  !       N Is a positive integer input variable set to the number
  !         of functions and variables.
  !
  !       X is an input array of length N.
  !
  !       FVEC is an input array of length N which must contain the
  !         functions evaluated at X.
  !
  !       FJAC is an output N by N array which contains the
  !         approximation to the Jacobian matrix evaluated at X.
  !
  !       LDFJAC is a positive integer input variable not less than N
  !         which specifies the leading dimension of the array FJAC.
  !
  !       IFLAG is an integer variable which can be used to terminate
  !         the execution of FDJAC1. See description of FCN.
  !
  !       ML is a nonnegative integer input variable which specifies
  !         the number of subdiagonals within the band of the
  !         Jacobian matrix. If the Jacobian is not banded, set
  !         ML to at least N - 1.
  !
  !       EPSFCN is an input variable used in determining a suitable
  !         step length for the forward-difference approximation. This
  !         approximation assumes that the relative errors in the
  !         functions are of the order of EPSFCN. If EPSFCN is less
  !         than the machine precision, it is assumed that the relative
  !         errors in the functions are of the order of the machine
  !         precision.
  !
  !       MU is a nonnegative integer input variable which specifies
  !         the number of superdiagonals within the band of the
  !         Jacobian matrix. If the Jacobian is not banded, set
  !         MU to at least N - 1.
  !
  !       WA1 and WA2 are work arrays of length N. If ML + MU + 1 is at
  !         least N, then the Jacobian is considered dense, and WA2 is
  !         not referenced.
  !
  !***
  ! **See also:**  SNSQ, SNSQE
  !***
  ! **Routines called:**  R1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   800301  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   900328  Added TYPE section.  (WRB)

  INTEGER N, Ldfjac, Iflag, Ml, Mu
  REAL Epsfcn
  REAL X(*), Fvec(*), Fjac(Ldfjac,*), Wa1(*), Wa2(*)
  INTEGER i, j, k, msum
  REAL eps, epsmch, h, temp
  REAL, PARAMETER :: zero = 0.0E0
  !* FIRST EXECUTABLE STATEMENT  FDJAC1
  epsmch = R1MACH(4)
  !
  eps = SQRT(MAX(Epsfcn,epsmch))
  msum = Ml + Mu + 1
  IF ( msum<N ) THEN
    !
    !        COMPUTATION OF BANDED APPROXIMATE JACOBIAN.
    !
    DO k = 1, msum
      DO j = k, N, msum
        Wa2(j) = X(j)
        h = eps*ABS(Wa2(j))
        IF ( h==zero ) h = eps
        X(j) = Wa2(j) + h
      END DO
      CALL FCN(N,X,Wa1,Iflag)
      IF ( Iflag<0 ) EXIT
      DO j = k, N, msum
        X(j) = Wa2(j)
        h = eps*ABS(Wa2(j))
        IF ( h==zero ) h = eps
        DO i = 1, N
          Fjac(i,j) = zero
          IF ( i>=j-Mu.AND.i<=j+Ml ) Fjac(i,j) = (Wa1(i)-Fvec(i))/h
        END DO
      END DO
    END DO
  ELSE
    !
    !        COMPUTATION OF DENSE APPROXIMATE JACOBIAN.
    !
    DO j = 1, N
      temp = X(j)
      h = eps*ABS(temp)
      IF ( h==zero ) h = eps
      X(j) = temp + h
      CALL FCN(N,X,Wa1,Iflag)
      IF ( Iflag<0 ) EXIT
      X(j) = temp
      DO i = 1, N
        Fjac(i,j) = (Wa1(i)-Fvec(i))/h
      END DO
    END DO
  END IF
  !
  !     LAST CARD OF SUBROUTINE FDJAC1.
  !
END SUBROUTINE FDJAC1
