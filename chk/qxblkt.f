*DECK QXBLKT
      SUBROUTINE QXBLKT (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  QXBLKT
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
C     PROGRAM TO ILLUSTRATE THE USE OF BLKTRI
C
C***ROUTINES CALLED  BLKTRI
C***REVISION HISTORY  (YYMMDD)
C   800103  DATE WRITTEN
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   901010  Added PASS/FAIL message and cleaned up FORMATs.  (RWC)
C***END PROLOGUE  QXBLKT
      DIMENSION Y(75,105), AM(75), BM(75), CM(75), AN(105), BN(105),
     1          CN(105), W(1952), S(75), T(105)
C***FIRST EXECUTABLE STATEMENT  QXBLKT
      ERMAX=1.E-3
      IFLG = 0
      NP = 1
      N = 63
      MP = 1
      M = 50
      IDIMY = 75
C
C     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING THE
C     COEFFICIENTS AND THE ARRAY Y.
C
      DELTAS = 1.0E0/(M+1)
      DO 101 I=1,M
         S(I) = I * DELTAS
  101 CONTINUE
      DELTAT = 1.0E0/(N+1)
      DO 102 J=1,N
         T(J) = J*DELTAT
  102 CONTINUE
C
C     COMPUTE THE COEFFICIENTS AM, BM AND CM CORRESPONDING TO THE S
C     DIRECTION.
C
      HDS = DELTAS/2.
      TDS = DELTAS+DELTAS
      DO 103 I=1,M
         TEMP1 = 1./(S(I)*TDS)
         TEMP2 = 1./((S(I)-HDS)*TDS)
         TEMP3 = 1./((S(I)+HDS)*TDS)
         AM(I) = TEMP1*TEMP2
         CM(I) = TEMP1*TEMP3
         BM(I) = -(AM(I)+CM(I))
  103 CONTINUE
C
C     COMPUTE THE COEFFICIENTS AN, BN AND CN CORRESPONDING TO THE T
C     DIRECTION.
C
      HDT = DELTAT/2.
      TDT = DELTAT+DELTAT
      DO 104 J=1,N
         TEMP1 = 1./(T(J)*TDT)
         TEMP2 = 1./((T(J)-HDT)*TDT)
         TEMP3 = 1./((T(J)+HDT)*TDT)
         AN(J) = TEMP1*TEMP2
         CN(J) = TEMP1*TEMP3
         BN(J) = -(AN(J)+CN(J))
  104 CONTINUE
C
C     COMPUTE RIGHT SIDE OF EQUATION
C
      DO 106 J=1,N
         DO 105 I=1,M
            Y(I,J) = 3.75*S(I)*T(J)*(S(I)**4.+T(J)**4.)
  105    CONTINUE
  106 CONTINUE
C
C     INCLUDE NONHOMOGENEOUS BOUNDARY INTO RIGHT SIDE. NOTE THAT THE
C     CORNER AT J=N,I=M INCLUDES CONTRIBUTIONS FROM BOTH BOUNDARIES.
C
      DO 107 J=1,N
         Y(M,J) = Y(M,J)-CM(M)*T(J)**5.
  107 CONTINUE
      DO 108 I=1,M
         Y(I,N) = Y(I,N)-CN(N)*S(I)**5.
  108 CONTINUE
C
  109 CALL BLKTRI (IFLG,NP,N,AN,BN,CN,MP,M,AM,BM,CM,IDIMY,Y,IERROR,W)
      IFLG = IFLG+1
      IF (IFLG-1) 109,109,110
C
C     COMPUTE DISCRETIZATION ERROR
C
  110 ERR = 0.
      DO 112 J=1,N
         DO 111 I=1,M
            Z = ABS(Y(I,J)-(S(I)*T(J))**5.)
            IF (Z .GT. ERR) ERR = Z
  111    CONTINUE
  112 CONTINUE
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
 1001 FORMAT ('1',20X,'SUBROUTINE BLKTRI EXAMPLE'///
     1        10X,'THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS'//
     2        32X,'IERROR = 0'/
     3        18X,'DISCRETIZATION ERROR = 1.6478E-05'/
     4        12X,'REQUIRED LENGTH OF W ARRAY = 823'//
     5        10X,'THE OUTPUT FROM YOUR COMPUTER IS'//
     6        32X,'IERROR =',I2/
     7        18X,'DISCRETIZATION ERROR =',1PE12.5/
     8        12X,'REQUIRED LENGTH OF W ARRAY =', I4)
 1002 FORMAT (60X,'PASS'/)
 1003 FORMAT (60X,'FAIL'/)
      END
