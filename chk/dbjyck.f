*DECK DBJYCK
      SUBROUTINE DBJYCK (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  DBJYCK
C***PURPOSE  Quick check for DBESJ and DBESY.
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (BJYCK-S, DBJYCK-D)
C***KEYWORDS  QUICK CHECK
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C   DBJYCK is a quick check routine for DBESJ and DBESY.  The main loops
C   evaluate the Wronskian and test the error.  Underflow and overflow
C   diagnostics are checked in addition to illegal arguments.
C
C***ROUTINES CALLED  D1MACH, DBESJ, DBESY, NUMXER, XERCLR, XGETF, XSETF
C***REVISION HISTORY  (YYMMDD)
C   750101  DATE WRITTEN
C   890911  Removed unnecessary intrinsics.  (WRB)
C   891004  Removed unreachable code.  (WRB)
C   891004  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   901205  Changed usage of D1MACH(3) to D1MACH(4).  (RWC)
C   910121  Editorial Changes.  (RWC)
C   910501  Added TYPE record.  (WRB)
C   910708  Code revised to test error returns for all values of
C           KPRINT.  (WRB)
C   910801  Editorial changes, some restructing and modifications to
C           obtain more information when there is failure of the
C           Wronskian.  (WRB)
C***END PROLOGUE  DBJYCK
      INTEGER I, IPASS, IX, K, KONTRL, LUN, M, N, NERR, NU, NY
      DOUBLE PRECISION ALP, DEL, ER, FNU, FNUP, RHPI, RX, TOL, X
      DOUBLE PRECISION FN(3), W(5), XX(5), Y(5)
      DOUBLE PRECISION D1MACH
      LOGICAL FATAL
C***FIRST EXECUTABLE STATEMENT  DBJYCK
      IF (KPRINT.GE.2) WRITE (LUN,90000)
C
      IPASS = 1
      RHPI  = 0.5D0/ATAN(1.0D0)
      XX(1) = 0.49D0
      XX(2) = 1.3D0
      XX(3) = 5.3D0
      XX(4) = 13.3D0
      XX(5) = 21.3D0
      FN(1) = 0.095D0
      FN(2) = 0.70D0
      FN(3) = 0.0D0
      TOL = MAX(500.0D0*D1MACH(4), 7.1D-12)
      DO 50 M=1,3
         DO 40 N=1,4
            DO 30 NU=1,4
               FNU = FN(M) + 12*(NU-1)
               DO 20 IX=1,5
                  IF (IX.LT.2 .AND. NU.GT.3) GO TO 20
                  X = XX(IX)
                  RX = RHPI/X
                  CALL DBESJ(X, FNU, N, Y, NY)
                  IF (NY.NE.0) GO TO 20
                  CALL DBESY(X, FNU, N, W)
                  FNUP = FNU + N
                  CALL DBESJ(X,FNUP,1,Y(N+1),NY)
                  IF (NY.NE.0) GO TO 20
                  CALL DBESY(X,FNUP,1,W(N+1))
                  DO 10 I=1,N
                     ER = Y(I+1)*W(I) - W(I+1)*Y(I) - RX
                     ER = ABS(ER)/RX
                     IF (ER.GT.TOL) THEN
                        IPASS = 0
                        IF (KPRINT.GE.2) WRITE (LUN,90010) M,N,NU,IX,I,
     *                     X,ER,TOL,Y(I),Y(I+1),W(I),W(I+1)
                     ENDIF
   10             CONTINUE
   20          CONTINUE
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE
C
C     Check small values of X and order
C
      N = 2
      FNU = 1.0D0
      X = D1MACH(4)/5.0D0
      RX = RHPI/X
      DO 60 I=1,3
         CALL DBESJ(X, FNU, N, Y, NY)
         CALL DBESY(X, FNU, N, W)
         ER = Y(2)*W(1) - W(2)*Y(1) - RX
         ER = ABS(ER)/RX
         IF (ER.GT.TOL) THEN
            IPASS = 0
            IF (KPRINT.GE.2) WRITE (LUN,90020) I,FNU,X,ER,TOL,
     +         Y(I),Y(I+1),W(I),W(I+1)
            GO TO 600
         ENDIF
         FNU = D1MACH(4)/100.0D0
         X = XX(2*I-1)
         RX = RHPI/X
   60 CONTINUE
C
C     Check large values of X and order
C
  600 DO 76 K=1,2
         DEL = 30*(K-1)
         FNU = 70.0D0+DEL
         DO 75 N=1,2
            X = 50.0D0 + DEL
            DO 70 I=1,5
               RX = RHPI/X
               CALL DBESJ(X, FNU, N, Y, NY)
               IF (NY.NE.0) GO TO 70
               CALL DBESY(X, FNU, N, W)
               IF (N.EQ.1) THEN
                  FNUP = FNU + 1.0D0
                  CALL DBESJ(X,FNUP,1,Y(2),NY)
                  IF (NY.NE.0) GO TO 70
                  CALL DBESY(X,FNUP,1,W(2))
               ENDIF
               ER = Y(2)*W(1) - Y(1)*W(2) - RX
               ER = ABS(ER)/RX
               IF (ER.GT.TOL) THEN
                  IPASS = 0
                  IF (KPRINT.GE.2) WRITE (LUN,90030) K,N,I,X,ER,TOL,
     *               Y(1),Y(2),W(1),W(2)
                  GO TO 800
               ENDIF
               X = X + 10.0D0
   70       CONTINUE
   75    CONTINUE
   76 CONTINUE
C
C     Check underflow flags
C
  800 X = D1MACH(1)*10.0D0
      ALP = 12.3D0
      N = 3
      CALL DBESJ(X, ALP, N, Y, NY)
      IF (NY.NE.3) THEN
         IPASS = 0
         IF (KPRINT.GE.2) WRITE (LUN,90040)
      ENDIF
C
C     Trigger 7 error conditions
C
      CALL XGETF (KONTRL)
      IF (KPRINT .LE. 2) THEN
         CALL XSETF (0)
      ELSE
         CALL XSETF (1)
      ENDIF
      FATAL = .FALSE.
      CALL XERCLR
C
      IF (KPRINT .GE. 3) WRITE (LUN,90050)
      XX(1) = 1.0D0
      XX(2) = 1.0D0
      XX(3) = 1.0D0
C
C     Illegal arguments
C
      DO 80 I=1,3
         XX(I) = -XX(I)
         N = INT(XX(3))
         CALL DBESJ(XX(1), XX(2), N, Y, NY)
         IF (NUMXER(NERR) .NE. 2) THEN
            IPASS = 0
            FATAL = .TRUE.
         ENDIF
         CALL XERCLR
         CALL DBESY(XX(1), XX(2), N, W)
         IF (NUMXER(NERR) .NE. 2) THEN
            IPASS = 0
            FATAL = .TRUE.
         ENDIF
         CALL XERCLR
         XX(I) = -XX(I)
   80 CONTINUE
C
C     Trigger overflow
C
      X = D1MACH(1)*10.0D0
      N = 3
      ALP = 2.3D0
      CALL DBESY(X, ALP, N, W)
      IF (NUMXER(NERR) .NE. 6) THEN
         IPASS = 0
         FATAL = .TRUE.
      ENDIF
      CALL XERCLR
      CALL XSETF (KONTRL)
      IF (FATAL) THEN
         IF (KPRINT .GE. 2) THEN
            WRITE (LUN, 90070)
         ENDIF
      ELSE
         IF (KPRINT .GE. 3) THEN
            WRITE (LUN, 90080)
         ENDIF
      ENDIF
C
      IF (IPASS.EQ.1 .AND. KPRINT.GE.2) WRITE (LUN,90100)
      IF (IPASS.EQ.0 .AND. KPRINT.GE.1) WRITE (LUN,90110)
      RETURN
90000 FORMAT (/ ' QUICK CHECKS FOR DBESJ AND DBESY' //)
90010 FORMAT (/ ' ERROR IN QUICK CHECK OF WRONSKIAN', 1P /
     +        ' M = ', I1,', N = ', I1, ', NU = ', I1, ', IX = ', I1,
     +        ', I = ', I1, /
     +        ' X = ', E14.7, ', ER   = ', E14.7, ', TOL = ', E14.7 /
     +        ' Y(I) = ', E14.7, ', Y(I+1) = ', E14.7 /
     +        ' W(I) = ', E14.7, ', W(I+1) = ', E14.7)
90020 FORMAT (/ ' ERROR IN QUICK CHECK OF SMALL X AND ORDER', 1P /
     +        ' I = ', I1,',  FNU = ', E14.7 /
     +        ' X = ', E14.7, ', ER = ', E14.7, ', TOL = ', E14.7 /
     +        ' Y(1) = ', E14.7, ', Y(2) = ', E14.7 /
     +        ' W(1) = ', E14.7, ', W(2) = ', E14.7)
90030 FORMAT (/ ' ERROR IN QUICK CHECK OF LARGE X AND ORDER', 1P /
     +        ' K = ', I1,', N = ', I1, ', I = ', I1 /
     +        ' X = ', E14.7, ', ER = ', E14.7, ', TOL = ', E14.7 /
     +        ' Y(1) = ', E14.7, ', Y(2) = ', E14.7 /
     +        ' W(1) = ', E14.7, ', W(2) = ', E14.7)
90040 FORMAT (/ ' ERROR IN DBESJ UNDERFLOW TEST' /)
90050 FORMAT (// ' TRIGGER 7 ERROR CONDITIONS' //)
90070 FORMAT (/ ' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
90080 FORMAT (/ ' ALL INCORRECT ARGUMENT TESTS PASSED')
90100 FORMAT (/' *********DBESJ AND DBESY PASSED ALL TESTS*********')
90110 FORMAT (/' *********DBESJ OR DBESY FAILED SOME TESTS*********')
      END
