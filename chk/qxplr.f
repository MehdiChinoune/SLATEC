*DECK QXPLR
      SUBROUTINE QXPLR (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  QXPLR
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
C          PROGRAM TO ILLUSTRATE THE USE OF SUBROUTINE HWSPLR TO SOLVE
C     THE EQUATION
C
C     (1/R)(D/DR)(R*(DU/DR)) + (1/R**2)(D/DTHETA)(DU/DTHETA) = 16*R**2
C
C     ON THE QUARTER-DISK 0 .LT. R .LT. 1, 0 .LT. THETA .LT. PI/2 WITH
C     WITH THE BOUNDARY CONDITIONS
C
C     U(1,THETA) = 1 - COS(4*THETA), 0 .LE. THETA .LE. 1
C
C     AND
C
C     (DU/DTHETA)(R,0) = (DU/DTHETA)(R,PI/2) = 0,  0 .LE. R .LE. 1.
C
C     (NOTE THAT THE SOLUTION U IS UNSPECIFIED AT R = 0.)
C          THE R-INTERVAL WILL BE DIVIDED INTO 50 PANELS AND THE
C     THETA-INTERVAL WILL BE DIVIDED INTO 48 PANELS.
C
C***ROUTINES CALLED  HWSPLR, PIMACH
C***REVISION HISTORY  (YYMMDD)
C   800103  DATE WRITTEN
C   890718  Changed computation of PI to use PIMACH.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   901010  Added PASS/FAIL message and cleaned up FORMATs.  (RWC)
C***END PROLOGUE  QXPLR
      DIMENSION F(100,50), BDC(51), BDD(51), W(1200), R(51), THETA(49)
C***FIRST EXECUTABLE STATEMENT  QXPLR
C
C     FROM DIMENSION STATEMENT WE GET VALUE OF IDIMF.  ALSO NOTE THAT W
C     IS DIMENSIONED 6*(N+1) + 8*(M+1).
C
      IDIMF = 100
      ERMAX=1.E-3
      A = 0.
      B = 1.
      M = 50
      MBDCND = 5
      C = 0.
      PI = PIMACH(DUM)
      D = PI/2.
      N = 48
      NBDCND = 3
      ELMBDA = 0.
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
         R(I) = (I-1)/50.0E0
  101 CONTINUE
      DO 102 J=1,NP1
         THETA(J) = (J-1)*PI/96.0E0
  102 CONTINUE
C
C     GENERATE BOUNDARY DATA.
C
      DO 103 I=1,MP1
         BDC(I) = 0.
         BDD(I) = 0.
  103 CONTINUE
C
C     BDA AND BDB ARE DUMMY VARIABLES.
C
      DO 104 J=1,NP1
         F(MP1,J) = 1.-COS(4.*THETA(J))
  104 CONTINUE
C
C     GENERATE RIGHT SIDE OF EQUATION.
C
      DO 106 I=1,M
         DO 105 J=1,NP1
            F(I,J) = 16.*R(I)**2
  105    CONTINUE
  106 CONTINUE
      CALL HWSPLR(A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,ELMBDA,F,
     1             IDIMF,PERTRB,IERROR,W)
C
C     COMPUTE DISCRETIZATION ERROR.  THE EXACT SOLUTION IS
C                U(R,THETA) = R**4*(1 - COS(4*THETA))
C
      ERR = 0.
      DO 108 I=1,MP1
         DO 107 J=1,NP1
            Z = ABS(F(I,J)-R(I)**4*(1.-COS(4.*THETA(J))))
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
 1001 FORMAT ('1',20X,'SUBROUTINE HWSPLR EXAMPLE'///
     1        10X,'THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS'//
     2        32X,'IERROR = 0'/
     3        18X,'DISCRETIZATION ERROR = 6.19134E-04'/
     4        12X,'REQUIRED LENGTH OF W ARRAY = 882'//
     5        10X,'THE OUTPUT FROM YOUR COMPUTER IS'//
     6        32X,'IERROR =',I2/
     7        18X,'DISCRETIZATION ERROR =',1PE12.5/
     8        12X,'REQUIRED LENGTH OF W ARRAY =',I4)
 1002 FORMAT (60X,'PASS'/)
 1003 FORMAT (60X,'FAIL'/)
      END
