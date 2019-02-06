*DECK QXCRT
      SUBROUTINE QXCRT (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  QXCRT
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
C          PROGRAM TO ILLUSTRATE THE USE OF SUBROUTINE HWSCRT TO SOLVE
C     THE EQUATION
C
C     (D/DX)(DU/DX) + (D/DY)(DU/DY) - 4*U
C
C     = (2 - (4 + PI**2/4)*X**2)*COS((Y+1)*PI/2)
C
C     WITH THE BOUNDARY CONDITIONS
C     ON THE RECTANGLE 0 .LT. X .LT. 2, -1 .LT. Y .LT. 3 WITH THE
C
C     U(0,Y) = 0
C                                          -1 .LE. Y .LE. 3
C     (DU/DX)(2,Y) = 4*COS((Y+1)*PI/2)
C
C     AND WITH U PERIODIC IN Y.
C          THE X-INTERVAL WILL BE DIVIDED INTO 40 PANELS AND THE
C     Y-INTERVAL WILL BE DIVIDED INTO 80 PANELS.
C
C***ROUTINES CALLED  HWSCRT, PIMACH
C***REVISION HISTORY  (YYMMDD)
C   800103  DATE WRITTEN
C   890718  Changed computation of PI to use PIMACH.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   901010  Added PASS/FAIL message and cleaned up FORMATs.  (RWC)
C***END PROLOGUE  QXCRT
      DIMENSION F(45,82), BDB(81), W(1200), X(41), Y(81)
C***FIRST EXECUTABLE STATEMENT  QXCRT
C
C     FROM DIMENSION STATEMENT WE GET VALUE OF IDIMF.  ALSO NOTE THAT W
C     IS DIMENSIONED 6*(N+1) + 8*(M+1).
C
      IDIMF = 45
      ERMAX=1.E-3
      A = 0.
      B = 2.
      M = 40
      MBDCND = 2
      C = -1.
      D = 3.
      N = 80
      NBDCND = 0
      ELMBDA = -4.
C
C     AUXILIARY QUANTITIES.
C
      PI = PIMACH(DUM)
      PIBY2 = PI/2.
      PISQ = PI**2
      MP1 = M+1
      NP1 = N+1
C
C     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING
C     BOUNDARY DATA AND THE RIGHT SIDE OF THE HELMHOLTZ EQUATION.
C
      DO 101 I=1,MP1
         X(I) = (I-1)/20.0E0
  101 CONTINUE
      DO 102 J=1,NP1
         Y(J) = -1.0E0+(J-1)/20.0E0
  102 CONTINUE
C
C     GENERATE BOUNDARY DATA.
C
      DO 103 J=1,NP1
         BDB(J) = 4.*COS((Y(J)+1.)*PIBY2)
  103 CONTINUE
C
C     BDA, BDC, AND BDD ARE DUMMY VARIABLES.
C
      DO 104 J=1,NP1
         F(1,J) = 0.
  104 CONTINUE
C
C     GENERATE RIGHT SIDE OF EQUATION.
C
      DO 106 I=2,MP1
         DO 105 J=1,NP1
            F(I,J) = (2.-(4.+PISQ/4.)*X(I)**2)*COS((Y(J)+1.)*PIBY2)
  105    CONTINUE
  106 CONTINUE
      CALL HWSCRT(A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,ELMBDA,F,
     1             IDIMF,PERTRB,IERROR,W)
C
C     COMPUTE DISCRETIZATION ERROR.  THE EXACT SOLUTION IS
C                U(X,Y) = X**2*COS((Y+1)*PIBY2)
C
      ERR = 0.
      DO 108 I=1,MP1
         DO 107 J=1,NP1
            Z = ABS(F(I,J)-X(I)**2*COS((Y(J)+1.)*PIBY2))
            IF (Z .GT. ERR) ERR = Z
  107    CONTINUE
  108 CONTINUE
C
      IPASS = 1
      IF (ERR.GT.ERMAX) IPASS = 0
      IF (KPRINT.EQ.0) RETURN
      IF (KPRINT.GE.2 .OR. IPASS.EQ.0) THEN
         WRITE (LUN,1001) IERROR,ERR,INT(W(1))
         IF (IPASS.EQ.1) THEN
            WRITE (LUN, 1002)
         ELSE
            WRITE (LUN, 1003)
         ENDIF
      ENDIF
      RETURN
C
 1001 FORMAT ('1',20X,'SUBROUTINE HWSCRT EXAMPLE'///
     1        10X,'THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS'//
     2        32X,'IERROR = 0'/
     3        18X,'DISCRETIZATION ERROR = 5.36508E-04'/
     4        12X,'REQUIRED LENGTH OF W ARRAY = 880'//
     5        10X,'THE OUTPUT FROM YOUR COMPUTER IS'//
     6        32X,'IERROR =',I2/
     7        18X,'DISCRETIZATION ERROR =',1PE12.5/
     8        12X,'REQUIRED LENGTH OF W ARRAY =',I4)
 1002 FORMAT (60X,'PASS'/)
 1003 FORMAT (60X,'FAIL'/)
      END
