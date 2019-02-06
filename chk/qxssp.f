*DECK QXSSP
      SUBROUTINE QXSSP (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  QXSSP
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
C
C     PROGRAM TO ILLUSTRATE THE USE OF HWSSSP
C
C***ROUTINES CALLED  HWSSSP, PIMACH
C***REVISION HISTORY  (YYMMDD)
C   800103  DATE WRITTEN
C   890718  Changed computation of PI to use PIMACH.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   901010  Added PASS/FAIL message and cleaned up FORMATs.  (RWC)
C***END PROLOGUE  QXSSP
      DIMENSION F(19,73), BDTF(73), SINT(19), SINP(73), W(1200)
C***FIRST EXECUTABLE STATEMENT  QXSSP
C
C     THE VALUE OF IDIMF IS THE FIRST DIMENSION OF F.  W IS
C     DIMENSIONED 11*(M+1)+6*(N+1)=647 SINCE M=18 AND N=72.
C
      PI = PIMACH(DUM)
      ERMAX=5.E-3
      TS = 0.0
      TF = PI/2.
      M = 18
      MBDCND = 6
      PS = 0.0
      PF = PI+PI
      N = 72
      NBDCND = 0
      ELMBDA = 0.
      IDIMF = 19
C
C     GENERATE SINES FOR USE IN SUBSEQUENT COMPUTATIONS
C
      DTHETA = TF/M
      MP1 = M+1
      DO 101 I=1,MP1
         SINT(I) = SIN((I-1)*DTHETA)
  101 CONTINUE
      DPHI = (PI+PI)/N
      NP1 = N+1
      DO 102 J=1,NP1
         SINP(J) = SIN((J-1)*DPHI)
  102 CONTINUE
C
C     COMPUTE RIGHT SIDE OF EQUATION AND STORE IN F
C
      DO 104 J=1,NP1
         DO 103 I=1,MP1
            F(I,J) = 2.-6.*(SINT(I)*SINP(J))**2
  103    CONTINUE
  104 CONTINUE
C
C     STORE DERIVATIVE DATA AT THE EQUATOR
C
      DO 105 J=1,NP1
         BDTF(J) = 0.
  105 CONTINUE
C
      CALL HWSSSP(TS,TF,M,MBDCND,BDTS,BDTF,PS,PF,N,NBDCND,BDPS,BDPF,
     1             ELMBDA,F,IDIMF,PERTRB,IERROR,W)
C
C     COMPUTE DISCRETIZATION ERROR. SINCE PROBLEM IS SINGULAR, THE
C     SOLUTION MUST BE NORMALIZED.
C
      ERR = 0.0
      DO 107 J=1,NP1
         DO 106 I=1,MP1
            Z = ABS(F(I,J)-(SINT(I)*SINP(J))**2-F(1,1))
            IF (Z .GT. ERR) ERR = Z
  106    CONTINUE
  107 CONTINUE
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
 1001 FORMAT ('1',20X,'SUBROUTINE HWSSSP EXAMPLE'///
     1        10X,'THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS'//
     2        32X,'IERROR = 0'/
     3        18X,'DISCRETIZATION ERROR = 3.38107E-03'/
     4        12X,'REQUIRED LENGTH OF W ARRAY = 600'//
     5        10X,'THE OUTPUT FROM YOUR COMPUTER IS'//
     6        32X,'IERROR =',I2/
     7        18X,'DISCRETIZATION ERROR =',1PE12.5 /
     8        12X,'REQUIRED LENGTH OF W ARRAY =',I4)
 1002 FORMAT (60X,'PASS'/)
 1003 FORMAT (60X,'FAIL'/)
      END
