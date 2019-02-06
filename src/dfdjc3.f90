!*==DFDJC3.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DFDJC3
      SUBROUTINE DFDJC3(FCN,M,N,X,Fvec,Fjac,Ldfjac,Iflag,Epsfcn,Wa)
      IMPLICIT NONE
!*--DFDJC35
!***BEGIN PROLOGUE  DFDJC3
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DNLS1 and DNLS1E
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (FDJAC3-S, DFDJC3-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!  **** Double Precision version of FDJAC3 ****
!
!     This subroutine computes a forward-difference approximation
!     to the M by N Jacobian matrix associated with a specified
!     problem of M functions in N variables.
!
!     The subroutine statement is
!
!       SUBROUTINE DFDJC3(FCN,M,N,X,FVEC,FJAC,LDFJAC,IFLAG,EPSFCN,WA)
!
!     where
!
!       FCN is the name of the user-supplied subroutine which
!         calculates the functions. FCN must be declared
!         in an external statement in the user calling
!         program, and should be written as follows.
!
!         SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
!         INTEGER LDFJAC,M,N,IFLAG
!         DOUBLE PRECISION X(N),FVEC(M),FJAC(LDFJAC,N)
!         ----------
!         When IFLAG.EQ.1 calculate the functions at X and
!         return this vector in FVEC.
!         ----------
!         RETURN
!         END
!
!         The value of IFLAG should not be changed by FCN unless
!         the user wants to terminate execution of DFDJC3.
!         In this case set IFLAG to a negative integer.
!
!       M is a positive integer input variable set to the number
!         of functions.
!
!       N is a positive integer input variable set to the number
!         of variables. N must not exceed M.
!
!       X is an input array of length N.
!
!       FVEC is an input array of length M which must contain the
!         functions evaluated at X.
!
!       FJAC is an output M by N array which contains the
!         approximation to the Jacobian matrix evaluated at X.
!
!       LDFJAC is a positive integer input variable not less than M
!         which specifies the leading dimension of the array FJAC.
!
!       IFLAG is an integer variable which can be used to terminate
!         THE EXECUTION OF DFDJC3. See description of FCN.
!
!       EPSFCN is an input variable used in determining a suitable
!         step length for the forward-difference approximation. This
!         approximation assumes that the relative errors in the
!         functions are of the order of EPSFCN. If EPSFCN is less
!         than the machine precision, it is assumed that the relative
!         errors in the functions are of the order of the machine
!         precision.
!
!       WA is a work array of length M.
!
!***SEE ALSO  DNLS1, DNLS1E
!***ROUTINES CALLED  D1MACH
!***REVISION HISTORY  (YYMMDD)
!   800301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  DFDJC3
      INTEGER M , N , Ldfjac , Iflag
      DOUBLE PRECISION Epsfcn
      DOUBLE PRECISION X(*) , Fvec(*) , Fjac(Ldfjac,*) , Wa(*)
      INTEGER i , j
      DOUBLE PRECISION eps , epsmch , h , temp , zero
      DOUBLE PRECISION D1MACH
      SAVE zero
      DATA zero/0.0D0/
!***FIRST EXECUTABLE STATEMENT  DFDJC3
      epsmch = D1MACH(4)
!
      eps = SQRT(MAX(Epsfcn,epsmch))
!      SET IFLAG=1 TO INDICATE THAT FUNCTION VALUES
!           ARE TO BE RETURNED BY FCN.
      Iflag = 1
      DO j = 1 , N
        temp = X(j)
        h = eps*ABS(temp)
        IF ( h==zero ) h = eps
        X(j) = temp + h
        CALL FCN(Iflag,M,N,X,Wa,Fjac,Ldfjac)
        IF ( Iflag<0 ) EXIT
        X(j) = temp
        DO i = 1 , M
          Fjac(i,j) = (Wa(i)-Fvec(i))/h
        ENDDO
      ENDDO
!
!     LAST CARD OF SUBROUTINE DFDJC3.
!
      END SUBROUTINE DFDJC3
