*DECK QXGBUN
      SUBROUTINE QXGBUN (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  QXGBUN
C***PURPOSE
C***LIBRARY   SLATEC
C***KEYWORDS  QUICK CHECK
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                        F I S H P A K                          *
C     *                                                               *
C     *                                                               *
C     *     A PACKAGE OF FORTRAN SUBPROGRAMS FOR THE SOLUTION OF      *
C     *                                                               *
C     *      SEPARABLE ELLIPTIC PARTIAL DIFFERENTIAL EQUATIONS        *
C     *                                                               *
C     *                  (VERSION  3 , JUNE 1979)                     *
C     *                                                               *
C     *                             BY                                *
C     *                                                               *
C     *        JOHN ADAMS, PAUL SWARZTRAUBER AND ROLAND SWEET         *
C     *                                                               *
C     *                             OF                                *
C     *                                                               *
C     *         THE NATIONAL CENTER FOR ATMOSPHERIC RESEARCH          *
C     *                                                               *
C     *                BOULDER, COLORADO  (80307)  U.S.A.             *
C     *                                                               *
C     *                   WHICH IS SPONSORED BY                       *
C     *                                                               *
C     *              THE NATIONAL SCIENCE FOUNDATION                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     PROGRAM TO ILLUSTRATE THE USE OF SUBROUTINE GENBUN
C
C***ROUTINES CALLED  GENBUN, PIMACH
C***REVISION HISTORY  (YYMMDD)
C   750701  DATE WRITTEN
C   890718  Changed computation of PI to use PIMACH.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   891009  Removed unreferenced variable.  (WRB)
C   891009  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   901010  Added PASS/FAIL message and cleaned up FORMATs.  (RWC)
C***END PROLOGUE  QXGBUN
      DIMENSION F(25,130), A(20), B(20), C(20), W(1200), X(20), Y(120)
C***FIRST EXECUTABLE STATEMENT  QXGBUN
C
C     FROM DIMENSION STATEMENT WE GET VALUE OF IDIMY.  ALSO NOTE THAT
C     W(.) IS DIMENSIONED 6*N + 5*M.
C
      ERMAX=1.E-2
      IDIMY = 25
      MPEROD = 1
      M = 20
      DELTAX = 1.0E0/M
      NPEROD = 0
      N = 120
      PI = PIMACH(DUM)
      DELTAY = 2.0E0*PI/N
C
C     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING
C     COEFFICIENTS AND RIGHT SIDE OF EQUATION.
C
      DO 100 I=1,M
         X(I) = (I-1)*DELTAX
  100    CONTINUE
      DO 105 J=1,N
         Y(J) = -PI + (J-1)*DELTAY
  105    CONTINUE
C
C     GENERATE COEFFICIENTS.
C
      S = (DELTAY/DELTAX)**2
      T = S*DELTAX
      A(1) = 0.
      B(1) = -2.0E0*S
      C(1) = 2.0E0*S
      DO 110 I=2,M
         A(I) = (1.+X(I))**2*S + (1.+X(I))*T
         C(I) = (1.+X(I))**2*S - (1.+X(I))*T
         B(I) = -2.0E0*(1.0E0+X(I))**2*S
  110    CONTINUE
      C(M) = 0.
C
C     GENERATE RIGHT SIDE OF EQUATION FOR I = 1 SHOWING INTRODUCTION OF
C     BOUNDARY DATA.
C
      DYSQ = DELTAY**2
      DO 115 J=1,N
         F(1,J) = DYSQ*(11. + 8./DELTAX)*SIN(Y(J))
  115    CONTINUE
C
C     GENERATE RIGHT SIDE.
C
      MM1 = M-1
      DO 125 I=2,MM1
         DO 120 J=1,N
            F(I,J) = DYSQ*3.*(1.+X(I))**4*SIN(Y(J))
  120       CONTINUE
  125    CONTINUE
C
C     GENERATE RIGHT SIDE FOR I = M SHOWING INTRODUCTION OF
C     BOUNDARY DATA.
C
      DO 130 J=1,N
         F(M,J) = DYSQ*(3.*(1.+X(M))**4 - 16.*((1.+X(M))/DELTAX)**2
     +                + 16.*(1.+X(M))/DELTAX)*SIN(Y(J))
  130    CONTINUE
      CALL GENBUN(NPEROD,N,MPEROD,M,A,B,C,IDIMY,F,IERROR,W)
C
C     COMPUTE DISCRETIZATION ERROR.  THE EXACT SOLUTION IS
C                   U(X,Y) = (1+X)**4*SIN(Y)
C
      ERR = 0.
      DO 140 I=1,M
         DO 135 J=1,N
            Z = ABS(F(I,J)-(1.+X(I))**4*SIN(Y(J)))
            IF (Z .GT. ERR) ERR = Z
  135       CONTINUE
  140    CONTINUE
C
      IPASS = 1
      IF (ERR.GT.ERMAX) IPASS = 0
      IF (KPRINT.EQ.0) RETURN
      IF (KPRINT.GE.2 .OR. IPASS.EQ.0) THEN
         WRITE (LUN,1001) IERROR, ERR, INT(W(1))
         IF (IPASS.EQ.1) THEN
            WRITE (LUN, 1002)
         ELSE
            WRITE (LUN, 1003)
         ENDIF
      ENDIF
      RETURN
C
 1001 FORMAT ('1',20X,'SUBROUTINE GENBUN EXAMPLE'///
     1        10X,'THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS'//
     2        32X,'IERROR = 0'/
     3        18X,'DISCRETIZATION ERROR = 7.94113E-03'/
     4        12X,'REQUIRED LENGTH OF W ARRAY = 740'//
     5        10X,'THE OUTPUT FROM YOUR COMPUTER IS'//
     6        32X,'IERROR =',I2/
     7        18X,'DISCRETIZATION ERROR =',1PE12.5/
     8        12X,'REQUIRED LENGTH OF W ARRAY =',I4)
 1002 FORMAT (60X,'PASS'/)
 1003 FORMAT (60X,'FAIL'/)
      END
