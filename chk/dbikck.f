*DECK DBIKCK
      SUBROUTINE DBIKCK (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  DBIKCK
C***PURPOSE  Quick check for DBESI and DBESK.
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (BIKCK-S, DBIKCK-D)
C***KEYWORDS  QUICK CHECK
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C   DBIKCK is a quick check routine for DBESI and DBESK.  The main loops
C   evaluate the Wronskian and test the error.  Underflow and overflow
C   diagnostics are checked in addition to illegal arguments.
C
C***ROUTINES CALLED  D1MACH, DBESI, DBESK, NUMXER, XERCLR, XGETF, XSETF
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
C***END PROLOGUE  DBIKCK
      INTEGER I, IPASS, IX, K, KODE, KONTRL, LUN, M, N, NERR, NU, NW, NY
      DOUBLE PRECISION ALP, DEL, ER, FNU, FNUP, RX, TOL, X
      DOUBLE PRECISION FN(3), W(5), XX(5), Y(5)
      DOUBLE PRECISION D1MACH
      LOGICAL FATAL
C***FIRST EXECUTABLE STATEMENT  DBIKCK
      IF (KPRINT .GE. 2) WRITE (LUN,90000)
C
      IPASS = 1
      XX(1) = 0.49D0
      XX(2) = 1.3D0
      XX(3) = 5.3D0
      XX(4) = 13.3D0
      XX(5) = 21.3D0
      FN(1) = 0.095D0
      FN(2) = 0.70D0
      FN(3) = 0.0D0
      TOL = MAX(500.0D0*D1MACH(4), 7.1D-12)
      DO 60 KODE=1,2
         DO 50 M=1,3
            DO 40 N=1,4
               DO 30 NU=1,4
                  FNU = FN(M) + 12*(NU-1)
                  DO 20 IX=1,5
                     IF (IX.LT.2 .AND. NU.GT.3) GO TO 20
                     X = XX(IX)
                     RX = 1.0D0/X
                     CALL DBESI(X, FNU, KODE, N, Y, NY)
                     IF (NY.NE.0) GO TO 20
                     CALL DBESK(X, FNU, KODE, N, W, NW)
                     IF (NW.NE.0) GO TO 20
                     FNUP = FNU + N
                     CALL DBESI(X,FNUP,KODE,1,Y(N+1),NY)
                     IF (NY.NE.0) GO TO 20
                     CALL DBESK(X,FNUP,KODE,1,W(N+1),NW)
                     IF (NW.NE.0) GO TO 20
                     DO 10 I=1,N
                        ER = Y(I+1)*W(I) + W(I+1)*Y(I) - RX
                        ER = ABS(ER)*X
                        IF (ER.GT.TOL) THEN
                           IPASS = 0
                           IF (KPRINT.GE.2) WRITE (LUN,90010) KODE,M,N,
     *                        NU,IX,I,X,ER,TOL,
     *                        Y(I),Y(I+1),W(I),W(I+1)
                        ENDIF
   10                CONTINUE
   20             CONTINUE
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE
C
C     Check small values of X and order
C
      N = 2
      FNU = 1.0D0
      X = D1MACH(4)
      DO 80 I=1,3
         DO 70 KODE=1,2
            CALL DBESI(X, FNU, KODE, N, Y, NY)
            CALL DBESK(X, FNU, KODE, N, W, NW)
            ER = Y(2)*W(1) + W(2)*Y(1) - 1.0D0/X
            ER = ABS(ER)*X
            IF (ER.GT.TOL) THEN
               IPASS = 0
               IF (KPRINT.GE.2) WRITE (LUN,90020) I,KODE,FNU,X,ER,TOL,
     +            Y(1),Y(2),W(1),W(2)
               GO TO 700
            ENDIF
   70    CONTINUE
C
  700    FNU = D1MACH(4)/100.0D0
         X = XX(2*I-1)
   80 CONTINUE
C
C     Check large values of X and order
C
      KODE = 2
      DO 76 K=1,2
         DEL = 30*(K-1)
         FNU = 45.0D0+DEL
         DO 75 N=1,2
            X = 20.0D0 + DEL
            DO 71 I=1,5
               RX = 1.0D0/X
               CALL DBESI(X, FNU, KODE, N, Y, NY)
               IF (NY.NE.0) GO TO 71
               CALL DBESK(X, FNU, KODE, N, W, NW)
               IF (NW.NE.0) GO TO 71
               IF (N.EQ.1) THEN
                  FNUP = FNU + 1.0D0
                  CALL DBESI(X,FNUP,KODE,1,Y(2),NY)
                  IF (NY.NE.0) GO TO 71
                  CALL DBESK(X,FNUP,KODE,1,W(2),NW)
                  IF (NW.NE.0) GO TO 71
               ENDIF
               ER = Y(2)*W(1) + Y(1)*W(2) - RX
               ER = ABS(ER)*X
               IF (ER.GT.TOL) THEN
                  IPASS = 0
                  IF (KPRINT.GE.2) WRITE (LUN,90030) K,N,I,FNUP,X,
     +               ER,TOL,Y(1),Y(2),W(1),W(2)
                  GO TO 760
               ENDIF
               X = X + 10.0D0
   71       CONTINUE
   75    CONTINUE
   76 CONTINUE
C
C     Check underflow flags
C
  760 X = D1MACH(1)*10.0D0
      ALP = 12.3D0
      N = 3
      CALL DBESI(X, ALP, 1, N, Y, NY)
      IF (NY.NE.3) THEN
         IPASS = 0
         IF (KPRINT.GE.2) WRITE (LUN,90040)
      ENDIF
C
      X = LOG(D1MACH(2)/10.0D0) + 20.0D0
      ALP = 1.3D0
      N = 3
      CALL DBESK(X, ALP, 1, N, W, NW)
      IF (NW.NE.3) THEN
         IPASS = 0
         IF (KPRINT.GE.2) WRITE (LUN,90050)
      ENDIF
C
C     Trigger 10 error conditions
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
      IF (KPRINT .GE. 3) WRITE (LUN,90060)
      XX(1) = 1.0D0
      XX(2) = 1.0D0
      XX(3) = 1.0D0
      XX(4) = 1.0D0
C
C     Illegal arguments
C
      DO 90 I=1,4
         XX(I) = -XX(I)
         K = INT(XX(3))
         N = INT(XX(4))
         CALL DBESI(XX(1), XX(2), K, N, Y, NY)
         IF (NUMXER(NERR) .NE. 2) THEN
            IPASS = 0
            FATAL = .TRUE.
         ENDIF
         CALL XERCLR
         CALL DBESK(XX(1), XX(2), K, N, W, NW)
         IF (NUMXER(NERR) .NE. 2) THEN
            IPASS = 0
            FATAL = .TRUE.
         ENDIF
         CALL XERCLR
         XX(I) = -XX(I)
   90 CONTINUE
C
C     Trigger overflow
C
      X = LOG(D1MACH(2)/10.0D0) + 20.0D0
      N = 3
      ALP = 2.3D0
      CALL DBESI(X, ALP, 1, N, Y, NY)
      IF (NUMXER(NERR) .NE. 6) THEN
         IPASS = 0
         FATAL = .TRUE.
      ENDIF
      CALL XERCLR
C
      X = D1MACH(1)*10.0D0
      CALL DBESK(X, ALP, 1, N, W, NW)
      IF (NUMXER(NERR) .NE. 6) THEN
         IPASS = 0
         FATAL = .TRUE.
      ENDIF
      CALL XERCLR
C
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
90000 FORMAT (/ ' QUICK CHECKS FOR DBESI AND DBESK' //)
90010 FORMAT (/ ' ERROR IN QUICK CHECK OF WRONSKIAN', 1P /
     +        ' KODE = ', I1,', M = ', I1, ', N = ', I1, ', NU = ', I1,
     +        ', IX = ', I1, ', I = ', I1 /
     +        ' X = ', E14.7, ', ER   = ', E14.7, ', TOL = ', E14.7 /
     +        ' Y(I) = ', E14.7, ', Y(I+1) = ', E14.7 /
     +        ' W(I) = ', E14.7, ', W(I+1) = ', E14.7)
90020 FORMAT (/ ' ERROR IN QUICK CHECK OF SMALL X AND ORDER', 1P /
     +        ' I = ', I1,', KODE = ', I1, ', FNU = ', E14.7 /
     +        ' X = ', E14.7, ', ER = ', E14.7, ', TOL = ', E14.7 /
     +        ' Y(1) = ', E14.7, ', Y(2) = ', E14.7 /
     +        ' W(1) = ', E14.7, ', W(2) = ', E14.7)
90030 FORMAT (/ ' ERROR IN QUICK CHECK OF LARGE X AND ORDER', 1P /
     +        ' K = ', I1,', N = ', I1, ', I = ', I1,
     +        ', FNUP = ', E14.7 /
     +        ' X = ', E14.7, ', ER = ', E14.7, ', TOL = ', E14.7 /
     +        ' Y(1) = ', E14.7, ', Y(2) = ', E14.7 /
     +        ' W(1) = ', E14.7, ', W(2) = ', E14.7)
90040 FORMAT (/ ' ERROR IN DBESI UNDERFLOW TEST' /)
90050 FORMAT (/ ' ERROR IN DBESK UNDERFLOW TEST' /)
90060 FORMAT (// ' TRIGGER 10 ERROR CONDITIONS' //)
90070 FORMAT (/ ' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
90080 FORMAT (/ ' ALL INCORRECT ARGUMENT TESTS PASSED')
90100 FORMAT (/' *********DBESI AND DBESK PASSED ALL TESTS***********')
90110 FORMAT (/' *********DBESI OR DBESK FAILED SOME TESTS***********')
      END
