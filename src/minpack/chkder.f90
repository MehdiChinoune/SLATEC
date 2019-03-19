!** CHKDER
SUBROUTINE CHKDER(M,N,X,Fvec,Fjac,Ldfjac,Xp,Fvecp,Mode,Err)
  IMPLICIT NONE
  !>
  !***
  !  Check the gradients of M nonlinear functions in N
  !            variables, evaluated at a point X, for consistency
  !            with the functions themselves.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  F3, G4C
  !***
  ! **Type:**      SINGLE PRECISION (CHKDER-S, DCKDER-D)
  !***
  ! **Keywords:**  GRADIENTS, JACOBIAN, MINPACK, NONLINEAR
  !***
  ! **Author:**  Hiebert, K. L. (SNLA)
  !***
  ! **Description:**
  !
  !   This subroutine is a companion routine to SNLS1,SNLS1E,SNSQ,and
  !   SNSQE which may be used to check the calculation of the Jacobian.
  !
  !     SUBROUTINE CHKDER
  !
  !     This subroutine checks the gradients of M nonlinear functions
  !     in N variables, evaluated at a point X, for consistency with
  !     the functions themselves. The user must call CKDER twice,
  !     first with MODE = 1 and then with MODE = 2.
  !
  !     MODE = 1. On input, X must contain the point of evaluation.
  !               On output, XP is set to a neighboring point.
  !
  !     MODE = 2. On input, FVEC must contain the functions and the
  !                         rows of FJAC must contain the gradients
  !                         of the respective functions each evaluated
  !                         at X, and FVECP must contain the functions
  !                         evaluated at XP.
  !               On output, ERR contains measures of correctness of
  !                          the respective gradients.
  !
  !     The subroutine does not perform reliably if cancellation or
  !     rounding errors cause a severe loss of significance in the
  !     evaluation of a function. Therefore, none of the components
  !     of X should be unusually small (in particular, zero) or any
  !     other value which may cause loss of significance.
  !
  !     The SUBROUTINE statement is
  !
  !       SUBROUTINE CHKDER(M,N,X,FVEC,FJAC,LDFJAC,XP,FVECP,MODE,ERR)
  !
  !     where
  !
  !       M is a positive integer input variable set to the number
  !         of functions.
  !
  !       N is a positive integer input variable set to the number
  !         of variables.
  !
  !       X is an input array of length N.
  !
  !       FVEC is an array of length M. On input when MODE = 2,
  !         FVEC must contain the functions evaluated at X.
  !
  !       FJAC is an M by N array. On input when MODE = 2,
  !         the rows of FJAC must contain the gradients of
  !         the respective functions evaluated at X.
  !
  !       LDFJAC is a positive integer input parameter not less than M
  !         which specifies the leading dimension of the array FJAC.
  !
  !       XP is an array of length N. On output when MODE = 1,
  !         XP is set to a neighboring point of X.
  !
  !       FVECP is an array of length M. On input when MODE = 2,
  !         FVECP must contain the functions evaluated at XP.
  !
  !       MODE is an integer input variable set to 1 on the first call
  !         and 2 on the second. Other values of MODE are equivalent
  !         to MODE = 1.
  !
  !       ERR is an array of length M. On output when MODE = 2,
  !         ERR contains measures of correctness of the respective
  !         gradients. If there is no severe loss of significance,
  !         then if ERR(I) is 1.0 the I-th gradient is correct,
  !         while if ERR(I) is 0.0 the I-th gradient is incorrect.
  !         For values of ERR between 0.0 and 1.0, the categorization
  !         is less certain. In general, a value of ERR(I) greater
  !         than 0.5 indicates that the I-th gradient is probably
  !         correct, while a value of ERR(I) less than 0.5 indicates
  !         that the I-th gradient is probably incorrect.
  !
  !***
  ! **References:**  M. J. D. Powell, A hybrid method for nonlinear equa-
  !                 tions. In Numerical Methods for Nonlinear Algebraic
  !                 Equations, P. Rabinowitz, Editor.  Gordon and Breach,
  !                 1988.
  !***
  ! **Routines called:**  R1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   800301  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  
  INTEGER M, N, Ldfjac, Mode
  REAL X(*), Fvec(*), Fjac(Ldfjac,*), Xp(*), Fvecp(*), Err(*)
  INTEGER i, j
  REAL eps, epsf, epslog, epsmch, factor, one, temp, zero
  REAL R1MACH
  SAVE factor, one, zero
  !
  DATA factor, one, zero/1.0E2, 1.0E0, 0.0E0/
  !* FIRST EXECUTABLE STATEMENT  CHKDER
  epsmch = R1MACH(4)
  !
  eps = SQRT(epsmch)
  !
  IF ( Mode==2 ) THEN
    !
    !        MODE = 2.
    !
    epsf = factor*epsmch
    epslog = LOG10(eps)
    DO i = 1, M
      Err(i) = zero
    ENDDO
    DO j = 1, N
      temp = ABS(X(j))
      IF ( temp==zero ) temp = one
      DO i = 1, M
        Err(i) = Err(i) + temp*Fjac(i,j)
      ENDDO
    ENDDO
    DO i = 1, M
      temp = one
      IF ( Fvec(i)/=zero.AND.Fvecp(i)/=zero.AND.ABS(Fvecp(i)-Fvec(i))&
        >=epsf*ABS(Fvec(i)) )&
        temp = eps*ABS((Fvecp(i)-Fvec(i))/eps-Err(i))&
        /(ABS(Fvec(i))+ABS(Fvecp(i)))
      Err(i) = one
      IF ( temp>epsmch.AND.temp<eps ) Err(i) = (LOG10(temp)-epslog)/epslog
      IF ( temp>=eps ) Err(i) = zero
    ENDDO
  ELSE
    !
    !        MODE = 1.
    !
    DO j = 1, N
      temp = eps*ABS(X(j))
      IF ( temp==zero ) temp = eps
      Xp(j) = X(j) + temp
    ENDDO
  ENDIF
  !
  !
  !     LAST CARD OF SUBROUTINE CHKDER.
  !
END SUBROUTINE CHKDER
