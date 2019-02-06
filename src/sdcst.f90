!*==SDCST.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK SDCST
      SUBROUTINE SDCST(Maxord,Mint,Iswflg,El,Tq)
      IMPLICIT NONE
!*--SDCST5
!***BEGIN PROLOGUE  SDCST
!***SUBSIDIARY
!***PURPOSE  SDCST sets coefficients used by the core integrator SDSTP.
!***LIBRARY   SLATEC (SDRIVE)
!***TYPE      SINGLE PRECISION (SDCST-S, DDCST-D, CDCST-C)
!***AUTHOR  Kahaner, D. K., (NIST)
!             National Institute of Standards and Technology
!             Gaithersburg, MD  20899
!           Sutherland, C. D., (LANL)
!             Mail Stop D466
!             Los Alamos National Laboratory
!             Los Alamos, NM  87545
!***DESCRIPTION
!
!  SDCST is called by SDNTL.  The array EL determines the basic method.
!  The array TQ is involved in adjusting the step size in relation
!  to truncation error.  EL and TQ depend upon MINT, and are calculated
!  for orders 1 to MAXORD(.LE. 12).  For each order NQ, the coefficients
!  EL are calculated from the generating polynomial:
!    L(T) = EL(1,NQ) + EL(2,NQ)*T + ... + EL(NQ+1,NQ)*T**NQ.
!  For the implicit Adams methods, L(T) is given by
!    dL/dT = (1+T)*(2+T)* ... *(NQ-1+T)/K,   L(-1) = 0,
!    where      K = factorial(NQ-1).
!  For the Gear methods,
!    L(T) = (1+T)*(2+T)* ... *(NQ+T)/K,
!    where      K = factorial(NQ)*(1 + 1/2 + ... + 1/NQ).
!  For each order NQ, there are three components of TQ.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   900329  Initial submission to SLATEC.
!***END PROLOGUE  SDCST
      REAL El(13,12) , factrl(12) , gamma(14) , sum , Tq(3,12)
      INTEGER i , Iswflg , j , Maxord , Mint , mxrd
!***FIRST EXECUTABLE STATEMENT  SDCST
      factrl(1) = 1.E0
      DO i = 2 , Maxord
        factrl(i) = i*factrl(i-1)
      ENDDO
!                                             Compute Adams coefficients
      IF ( Mint==1 ) THEN
        gamma(1) = 1.E0
        DO i = 1 , Maxord + 1
          sum = 0.E0
          DO j = 1 , i
            sum = sum - gamma(j)/(i-j+2)
          ENDDO
          gamma(i+1) = sum
        ENDDO
        El(1,1) = 1.E0
        El(2,1) = 1.E0
        El(2,2) = 1.E0
        El(3,2) = 1.E0
        DO j = 3 , Maxord
          El(2,j) = factrl(j-1)
          DO i = 3 , j
            El(i,j) = (j-1)*El(i,j-1) + El(i-1,j-1)
          ENDDO
          El(j+1,j) = 1.E0
        ENDDO
        DO j = 2 , Maxord
          El(1,j) = El(1,j-1) + gamma(j)
          El(2,j) = 1.E0
          DO i = 3 , j + 1
            El(i,j) = El(i,j)/((i-1)*factrl(j-1))
          ENDDO
        ENDDO
        DO j = 1 , Maxord
          Tq(1,j) = -1.E0/(factrl(j)*gamma(j))
          Tq(2,j) = -1.E0/gamma(j+1)
          Tq(3,j) = -1.E0/gamma(j+2)
        ENDDO
!                                              Compute Gear coefficients
      ELSEIF ( Mint==2 ) THEN
        El(1,1) = 1.E0
        El(2,1) = 1.E0
        DO j = 2 , Maxord
          El(1,j) = factrl(j)
          DO i = 2 , j
            El(i,j) = j*El(i,j-1) + El(i-1,j-1)
          ENDDO
          El(j+1,j) = 1.E0
        ENDDO
        sum = 1.E0
        DO j = 2 , Maxord
          sum = sum + 1.E0/j
          DO i = 1 , j + 1
            El(i,j) = El(i,j)/(factrl(j)*sum)
          ENDDO
        ENDDO
        DO j = 1 , Maxord
          IF ( j>1 ) Tq(1,j) = 1.E0/factrl(j-1)
          Tq(2,j) = (j+1)/El(1,j)
          Tq(3,j) = (j+2)/El(1,j)
        ENDDO
      ENDIF
!                          Compute constants used in the stiffness test.
!                          These are the ratio of TQ(2,NQ) for the Gear
!                          methods to those for the Adams methods.
      IF ( Iswflg==3 ) THEN
        mxrd = MIN(Maxord,5)
        IF ( Mint==2 ) THEN
          gamma(1) = 1.E0
          DO i = 1 , mxrd
            sum = 0.E0
            DO j = 1 , i
              sum = sum - gamma(j)/(i-j+2)
            ENDDO
            gamma(i+1) = sum
          ENDDO
        ENDIF
        sum = 1.E0
        DO i = 2 , mxrd
          sum = sum + 1.E0/i
          El(1+i,1) = -(i+1)*sum*gamma(i+1)
        ENDDO
      ENDIF
      END SUBROUTINE SDCST
