!** SDCST
SUBROUTINE SDCST(Maxord,Mint,Iswflg,El,Tq)
  !>
  !***
  !  SDCST sets coefficients used by the core integrator SDSTP.
  !***
  ! **Library:**   SLATEC (SDRIVE)
  !***
  ! **Type:**      SINGLE PRECISION (SDCST-S, DDCST-D, CDCST-C)
  !***
  ! **Author:**  Kahaner, D. K., (NIST)
  !             National Institute of Standards and Technology
  !             Gaithersburg, MD  20899
  !           Sutherland, C. D., (LANL)
  !             Mail Stop D466
  !             Los Alamos National Laboratory
  !             Los Alamos, NM  87545
  !***
  ! **Description:**
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
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   900329  Initial submission to SLATEC.

  REAL El(13,12), factrl(12), gama(14), summ, Tq(3,12)
  INTEGER i, Iswflg, j, Maxord, Mint, mxrd
  !* FIRST EXECUTABLE STATEMENT  SDCST
  factrl(1) = 1.E0
  DO i = 2, Maxord
    factrl(i) = i*factrl(i-1)
  END DO
  !                                             Compute Adams coefficients
  IF ( Mint==1 ) THEN
    gama(1) = 1.E0
    DO i = 1, Maxord + 1
      summ = 0.E0
      DO j = 1, i
        summ = summ - gama(j)/(i-j+2)
      END DO
      gama(i+1) = summ
    END DO
    El(1,1) = 1.E0
    El(2,1) = 1.E0
    El(2,2) = 1.E0
    El(3,2) = 1.E0
    DO j = 3, Maxord
      El(2,j) = factrl(j-1)
      DO i = 3, j
        El(i,j) = (j-1)*El(i,j-1) + El(i-1,j-1)
      END DO
      El(j+1,j) = 1.E0
    END DO
    DO j = 2, Maxord
      El(1,j) = El(1,j-1) + gama(j)
      El(2,j) = 1.E0
      DO i = 3, j + 1
        El(i,j) = El(i,j)/((i-1)*factrl(j-1))
      END DO
    END DO
    DO j = 1, Maxord
      Tq(1,j) = -1.E0/(factrl(j)*gama(j))
      Tq(2,j) = -1.E0/gama(j+1)
      Tq(3,j) = -1.E0/gama(j+2)
    END DO
    !                                              Compute Gear coefficients
  ELSEIF ( Mint==2 ) THEN
    El(1,1) = 1.E0
    El(2,1) = 1.E0
    DO j = 2, Maxord
      El(1,j) = factrl(j)
      DO i = 2, j
        El(i,j) = j*El(i,j-1) + El(i-1,j-1)
      END DO
      El(j+1,j) = 1.E0
    END DO
    summ = 1.E0
    DO j = 2, Maxord
      summ = summ + 1.E0/j
      DO i = 1, j + 1
        El(i,j) = El(i,j)/(factrl(j)*summ)
      END DO
    END DO
    DO j = 1, Maxord
      IF ( j>1 ) Tq(1,j) = 1.E0/factrl(j-1)
      Tq(2,j) = (j+1)/El(1,j)
      Tq(3,j) = (j+2)/El(1,j)
    END DO
  END IF
  !                          Compute constants used in the stiffness test.
  !                          These are the ratio of TQ(2,NQ) for the Gear
  !                          methods to those for the Adams methods.
  IF ( Iswflg==3 ) THEN
    mxrd = MIN(Maxord,5)
    IF ( Mint==2 ) THEN
      gama(1) = 1.E0
      DO i = 1, mxrd
        summ = 0.E0
        DO j = 1, i
          summ = summ - gama(j)/(i-j+2)
        END DO
        gama(i+1) = summ
      END DO
    END IF
    summ = 1.E0
    DO i = 2, mxrd
      summ = summ + 1.E0/i
      El(1+i,1) = -(i+1)*summ*gama(i+1)
    END DO
  END IF
END SUBROUTINE SDCST
