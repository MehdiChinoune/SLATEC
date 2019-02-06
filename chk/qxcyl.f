*DECK QXCYL
      SUBROUTINE QXCYL (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  QXCYL
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
C
C          PROGRAM TO ILLUSTRATE THE USE OF SUBROUTINE HWSCYL TO SOLVE
C     THE EQUATION
C
C     (1/R)(D/DR)(R*(DU/DR)) + (D/DZ)(DU/DZ)
C
C     = (2*R*Z)**2*(4*Z**2 + 3*R**2)
C
C     ON THE RECTANGLE 0 .LT. R .LT. 1, 0 .LT. Z .LT. 1 WITH THE
C     BOUNDARY CONDITIONS
C
C     U(0,Z) UNSPECIFIED
C                                            0 .LE. Z .LE. 1
C     (DU/DR)(1,Z) = 4*Z**4
C
C     AND
C
C     (DU/DZ)(R,0) = 0
C                                            0 .LE. R .LE. 1
C     (DU/DZ)(R,1) = 4*R**4 .
C
C          THE R-INTERVAL WILL BE DIVIDED INTO 50 PANELS AND THE
C     Z-INTERVAL WILL BE DIVIDED INTO 100 PANELS.
C
C***ROUTINES CALLED  HWSCYL
C***REVISION HISTORY  (YYMMDD)
C   800103  DATE WRITTEN
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   901010  Added PASS/FAIL message and cleaned up FORMATs.  (RWC)
C   930415  Test modified to use a 64 by 128 grid.  (WRB)
C***END PROLOGUE  QXCYL
      DIMENSION F(65,129), BDA(129), BDB(129), BDC(65), BDD(65),
     1   W(1400), R(65), Z(129)
C***FIRST EXECUTABLE STATEMENT  QXCYL
      IF (KPRINT .GE. 2) WRITE (LUN, 1000)
C
C     FROM DIMENSION STATEMENT WE GET VALUE OF IDIMF.
C
      IDIMF = 65
      ERMAX = 1.0E-3
      A = 0.0
      B = 1.0
      M = 64
      MBDCND = 6
      C = 0.0
      D = 1.0
      N = 128
      NBDCND = 3
      ELMBDA = 0.0
C
C     AUXILIARY QUANTITIES.
C
      MP1 = M+1
      NP1 = N+1
C
C     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING
C     BOUNDARY DATA AND THE RIGHT SIDE OF THE POISSON EQUATION.
C
      DO 101 I=1,MP1
         R(I) = (I-1)/64.0E0
  101 CONTINUE
      DO 102 J=1,NP1
         Z(J) = (J-1)/128.0E0
  102 CONTINUE
C
C     GENERATE BOUNDARY DATA.
C
      DO 103 J=1,NP1
         BDB(J) = 4.0*Z(J)**4
  103 CONTINUE
      DO 104 I=1,MP1
         BDC(I) = 0.0
         BDD(I) = 4.0*R(I)**4
  104 CONTINUE
C
C     BDA IS A DUMMY VARIABLE.
C
C     GENERATE RIGHT SIDE OF EQUATION.
C
      DO 106 I=1,MP1
         DO 105 J=1,NP1
            F(I,J) = 4.0*R(I)**2*Z(J)**2*(4.0*Z(J)**2+3.0*R(I)**2)
  105    CONTINUE
  106 CONTINUE
      CALL HWSCYL (A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,ELMBDA,F,
     1             IDIMF,PERTRB,IERROR,W)
C
C     COMPUTE DISCRETIZATION ERROR BY MINIMIZING OVER ALL A THE FUNCTION
C     NORM(F(I,J) - A - U(R(I),Z(J))).  THE EXACT SOLUTION IS
C     U(R,Z) = (R*Z)**4 + ARBITRARY CONSTANT.
C
      X = 0.0
      DO 108 I=1,MP1
         DO 107 J=1,NP1
            X = X+F(I,J)-(R(I)*Z(J))**4
  107    CONTINUE
  108 CONTINUE
      X = X/(NP1*MP1)
      DO 110 I=1,MP1
         DO 109 J=1,NP1
            F(I,J) = F(I,J)-X
  109    CONTINUE
  110 CONTINUE
      ERR = 0.0
      DO 112 I=1,MP1
         DO 111 J=1,NP1
            X = ABS(F(I,J)-(R(I)*Z(J))**4)
            IF (X .GT. ERR) ERR = X
  111    CONTINUE
  112 CONTINUE
C
      IPASS = 1
      IF (ERR .GT. ERMAX) IPASS = 0
      IF (KPRINT.GE.3 .OR. (KPRINT.GE.2.AND.IPASS.EQ.0))
     +   WRITE (LUN,1001) IERROR,PERTRB,ERR,INT(W(1))
C
C     Print PASS/FAIL message.
C
      IF (IPASS.EQ.1 .AND. KPRINT.GE.2) WRITE (LUN, 1002)
      IF (IPASS.EQ.0 .AND. KPRINT.GE.1) WRITE (LUN, 1003)
      RETURN
C
 1000 FORMAT ('1', 20X, 'SUBROUTINE HWSCYL EXAMPLE' //)
 1001 FORMAT (10X, 'THE OUTPUT FROM YOUR COMPUTER IS' //
     +        32X, 'IERROR =', I2 /
     +        32X, 'PERTRB =', E12.5 /
     +        18X, 'DISCRETIZATION ERROR =', 1PE12.5 /
     +        12X, 'REQUIRED LENGTH OF W ARRAY = ', I4)
 1002 FORMAT (25X, 'HWSCYL TEST PASSED' /)
 1003 FORMAT (25X, 'HWSCYL TEST FAILED' /)
      END
