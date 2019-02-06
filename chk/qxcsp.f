*DECK QXCSP
      SUBROUTINE QXCSP (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  QXCSP
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
C     PROGRAM TO ILLUSTRATE THE USE OF HWSCSP
C
C***ROUTINES CALLED  HWSCSP, PIMACH
C***REVISION HISTORY  (YYMMDD)
C   800103  DATE WRITTEN
C   890718  Changed computation of PI to use PIMACH.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   901010  Added PASS/FAIL message and cleaned up FORMATs.  (RWC)
C***END PROLOGUE  QXCSP
      DIMENSION F(48,33), BDTF(33), W(1200), R(33), THETA(48)
C***FIRST EXECUTABLE STATEMENT  QXCSP
C
C     THE VALUE OF IDIMF IS THE FIRST DIMENSION OF F.  SINCE M=36, N=32,
C     L=N THEREFORE K=5 AND W IS DIMENSIONED 2*(L+1)*(K-1) + 6*(M+N)
C     + MAX(4*N,6*M) + 14 = 902.
C
      ERMAX=1.E-3
      PI = PIMACH(DUM)
      INTL = 0
      TS = 0.
      TF = PI/2.
      M = 36
      MBDCND = 6
      RS = 0.
      RF = 1.
      N = 32
      NBDCND = 5
      ELMBDA = 0.
      IDIMF = 48
C
C     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING THE
C     BOUNDARY DATA AND THE RIGHT SIDE OF THE EQUATION.
C
      MP1 = M+1
      DTHETA = TF/M
      DO 101 I=1,MP1
         THETA(I) = (I-1)*DTHETA
  101 CONTINUE
      NP1 = N+1
      DR = 1.0E0/N
      DO 102 J=1,NP1
         R(J) = (J-1)*DR
  102 CONTINUE
C
C     GENERATE NORMAL DERIVATIVE DATA AT EQUATOR
C
      DO 103 J=1,NP1
         BDTF(J) = 0.
  103 CONTINUE
C
C     COMPUTE BOUNDARY DATA ON THE SURFACE OF THE SPHERE
C
      DO 104 I=1,MP1
         F(I,N+1) = COS(THETA(I))**4
  104 CONTINUE
C
C     COMPUTE RIGHT SIDE OF EQUATION
C
      DO 106 I=1,MP1
         CI4 = 12.0E0*COS(THETA(I))**2
         DO 105 J=1,N
            F(I,J) = CI4*R(J)**2
  105    CONTINUE
  106 CONTINUE
C
      CALL HWSCSP (INTL,TS,TF,M,MBDCND,BDTS,BDTF,RS,RF,N,NBDCND,BDRS,
     1             BDRF,ELMBDA,F,IDIMF,PERTRB,IERROR,W)
C
C     COMPUTE DISCRETIZATION ERROR
C
      ERR = 0.
      DO 108 I=1,MP1
         CI4 = COS(THETA(I))**4
         DO 107 J=1,N
            Z = ABS(F(I,J)-CI4*R(J)**4)
            IF (Z .GT. ERR) ERR = Z
  107    CONTINUE
  108 CONTINUE
C
      IPASS = 1
      IF (ERR.GT.ERMAX) IPASS = 0
      IF (KPRINT.NE.0) THEN
         IF (KPRINT.GE.2 .OR. IPASS.EQ.0) THEN
            WRITE (LUN,1001) IERROR,ERR,INT(W(1))
            IF (IPASS.EQ.1) THEN
               WRITE (LUN, 1003)
            ELSE
               WRITE (LUN, 1004)
            ENDIF
         ENDIF
      ENDIF
C
C     THE FOLLOWING PROGRAM ILLUSTRATES THE USE OF HWSCSP TO SOLVE
C     A THREE DIMENSIONAL PROBLEM WHICH HAS LONGITUDINAL DEPENDENCE
C
      MBDCND = 2
      NBDCND = 1
      DPHI = PI/72.
      ELMBDA = -2.0E0*(1.0E0-COS(DPHI))/DPHI**2
C
C     COMPUTE BOUNDARY DATA ON THE SURFACE OF THE SPHERE
C
      DO 109 I=1,MP1
         F(I,N+1) = SIN(THETA(I))
  109 CONTINUE
C
C     COMPUTE RIGHT SIDE OF THE EQUATION
C
      DO 111 J=1,N
         DO 110 I=1,MP1
            F(I,J) = 0.
  110    CONTINUE
  111 CONTINUE
C
      CALL HWSCSP (INTL,TS,TF,M,MBDCND,BDTS,BDTF,RS,RF,N,NBDCND,BDRS,
     1             BDRF,ELMBDA,F,IDIMF,PERTRB,IERROR,W)
C
C     COMPUTE DISCRETIZATION ERROR   (FOURIER COEFFICIENTS)
C
      ERR = 0.
      DO 113 I=1,MP1
         SI = SIN(THETA(I))
         DO 112 J=1,NP1
            Z = ABS(F(I,J)-R(J)*SI)
            IF (Z .GT. ERR) ERR = Z
  112    CONTINUE
  113 CONTINUE
C
      IF (ERR.GT.ERMAX) IPASS = 0
      IF (KPRINT.EQ.0) RETURN
      IF (KPRINT.GE.2 .OR. IPASS.EQ.0) THEN
         WRITE (LUN,1002) IERROR,ERR,INT(W(1))
         IF (IPASS.EQ.1) THEN
            WRITE (LUN, 1003)
         ELSE
            WRITE (LUN, 1004)
         ENDIF
      ENDIF
      RETURN
C
 1001 FORMAT ('1',20X,'SUBROUTINE HWSCSP EXAMPLE 1'///
     1        10X,'THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS'//
     2        32X,'IERROR = 0'/
     3        18X,'DISCRETIZATION ERROR = 7.99842E-04'/
     4        12X,'REQUIRED LENGTH OF W ARRAY = 775'//
     5        10X,'THE OUTPUT FROM YOUR COMPUTER IS'//
     6        32X,'IERROR =',I2/
     7        18X,'DISCRETIZATION ERROR =',1PE12.5/
     8        12X,'REQUIRED LENGTH OF W ARRAY =',I4)
 1002 FORMAT ('1',20X,'SUBROUTINE HWSCSP EXAMPLE 2'///
     1        10X,'THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS'//
     2        32X,'IERROR = 0'/
     3        18X,'DISCRETIZATION ERROR = 5.86824E-05'/
     4        12X,'REQUIRED LENGTH OF W ARRAY = 775'//
     5        10X,'THE OUTPUT FROM YOUR COMPUTER IS'//
     6        32X,'IERROR =',I2/
     7        18X,'DISCRETIZATION ERROR =',1PE12.5/
     8        12X,'REQUIRED LENGTH OF W ARRAY =',I4)
 1003 FORMAT (60X,'PASS'/)
 1004 FORMAT (60X,'FAIL'/)
      END
