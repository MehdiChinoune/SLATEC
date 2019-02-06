*DECK DBOCQX
      SUBROUTINE DBOCQX (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  DBOCQX
C***PURPOSE  Quick check for DBOCLS.
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (SBOCQX-S, DBOCQX-D)
C***KEYWORDS  QUICK CHECK
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     MINIMAL TEST DRIVER FOR DBOCLS, BOUNDED CONSTRAINED LEAST
C     SQUARES SOLVER.  DELIVERS THE VALUE IPASS=1 IF 8 TESTS WERE
C     PASSED.  DELIVER THE VALUE IPASS=0 IF ANY ONE OF THEM FAILED.
C
C     RUN FOUR BOUNDED LEAST SQUARES PROBLEMS THAT COME FROM THE
C     DIPLOME WORK OF P. ZIMMERMANN.
C
C***ROUTINES CALLED  D1MACH, DBOCLS, DBOLS, DCOPY, DNRM2
C***REVISION HISTORY  (YYMMDD)
C   850310  DATE WRITTEN
C   890618  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   901010  Added PASS/FAIL message.  (RWC)
C***END PROLOGUE  DBOCQX
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION
     *  D(6,5),W(11,11),BL(5,2),BU(5,2),X(30),RW(55),XTRUE(9)
      DOUBLE PRECISION C(5,5)
      DOUBLE PRECISION BL1(10),BU1(10)
      INTEGER IND(10),IW(20),IOPT(40)
      DOUBLE PRECISION RHS(6,2)
      CHARACTER*4 MSG
C
      DATA ((C(I,J),I=1,5),J=1,5)/1.D0,10.D0,4.D0,8.D0,1.D0,1.D0,10.D0,
     +     2.D0,-1.D0,1.D0,1.D0,-3.D0,-3.D0,2.D0,1.D0,1.D0,5.D0,5.D0,
     +     5.D0,1.D0,1.D0,4.D0,-1.D0,-3.D0,1.D0/
      DATA ((D(I,J),I=1,6),J=1,5)/-74.D0,14.D0,66.D0,-12.D0,3.D0,4.D0,
     +     80.D0,-69.D0,-72.D0,66.D0,8.D0,-12.D0,18.D0,21.D0,-5.D0,
     +     -30.D0,-7.D0,4.D0,-11.D0,28.D0,7.D0,-23.D0,-4.D0,4.D0,-4.D0,
     +     0.D0,1.D0,3.D0,1.D0,0.D0/
      DATA ((BL(I,J),I=1,5),J=1,2)/1.D0,0.D0,-1.D0,1.D0,-4.D0,-1.D0,
     +     0.D0,-3.D0,1.D0,-6.D0/
      DATA ((BU(I,J),I=1,5),J=1,2)/3.D0,2.D0,1.D0,3.D0,-2.D0,3.D0,4.D0,
     +     1.D0,5.D0,-2.D0/
      DATA ((RHS(I,J),I=1,6),J=1,2)/51.D0,-61.D0,-56.D0,69.D0,10.D0,
     +     -12.D0,-5.D0,-9.D0,708.D0,4165.D0,-13266.D0,8409.D0/
      DATA (XTRUE(J),J=1,9)/1.D0,2.D0,-1.D0,3.D0,-4.D0,1.D0,32.D0,30.D0,
     +     31.D0/
C***FIRST EXECUTABLE STATEMENT  DBOCQX
      MDW = 11
      MROWS = 6
      NCOLS = 5
      MCON = 4
      IOPT(1) = 99
      IPASS = 1
      ITEST = 0
C
      IF (KPRINT.GE.2) WRITE (LUN, 99998)
C
      DO 50 IB = 1,2
          DO 40 IRHS = 1,2
C
C           TRANSFER DATA TO WORKING ARRAY W(*,*).
C
              DO 10 J = 1,NCOLS
                  CALL DCOPY(MROWS,D(1,J),1,W(1,J),1)
   10         CONTINUE
C
              CALL DCOPY(MROWS,RHS(1,IRHS),1,W(1,NCOLS+1),1)
C
C             SET BOUND INDICATOR FLAGS.
C
              DO 20 J = 1,NCOLS
                  IND(J) = 3
   20         CONTINUE
C
              CALL DBOLS(W,MDW,MROWS,NCOLS,BL(1,IB),BU(1,IB),IND,IOPT,X,
     *                   RNORM,MODE,RW,IW)
              DO 30 J = 1,NCOLS
                  X(J) = X(J) - XTRUE(J)
   30         CONTINUE
C
              SR = DNRM2(NCOLS,X,1)
              MPASS = 1
              IF (SR.GT.10.D2*SQRT(D1MACH(4))) MPASS = 0
              IPASS = IPASS*MPASS
              IF (KPRINT.GE.2) THEN
                 MSG = 'PASS'
                 IF (MPASS.EQ.0) MSG = 'FAIL'
                 ITEST = ITEST + 1
                 WRITE (LUN, 99999) ITEST, IB, IRHS, SR, MSG
              ENDIF
   40     CONTINUE
   50 CONTINUE
C
C     RUN STOER'S PROBLEM FROM 1971 SIAM J. N. ANAL. PAPER.
C
      DO 90 IB = 1,2
         DO 80 IRHS = 1,2
            CALL DCOPY(11*10,0.D0,0,W,1)
            CALL DCOPY(NCOLS,BL(1,IB),1,BL1,1)
            CALL DCOPY(NCOLS,BU(1,IB),1,BU1,1)
            IND(NCOLS+1) = 2
            IND(NCOLS+2) = 1
            IND(NCOLS+3) = 2
            IND(NCOLS+4) = 3
            BU1(NCOLS+1) = 5.
            BL1(NCOLS+2) = 20.
            BU1(NCOLS+3) = 30.
            BL1(NCOLS+4) = 11.
            BU1(NCOLS+4) = 40.
            DO 60 J = 1,NCOLS
               CALL DCOPY(MCON,C(1,J),1,W(1,J),1)
               CALL DCOPY(MROWS,D(1,J),1,W(MCON+1,J),1)
   60       CONTINUE
C
            CALL DCOPY(MROWS,RHS(1,IRHS),1,W(MCON+1,NCOLS+1),1)
C
C           CHECK LENGTHS OF REQD. ARRAYS.
C
            IOPT(01) = 2
            IOPT(02) = 11
            IOPT(03) = 11
            IOPT(04) = 10
            IOPT(05) = 30
            IOPT(06) = 55
            IOPT(07) = 20
            IOPT(08) = 40
            IOPT(09) = 99
            CALL DBOCLS(W,MDW,MCON,MROWS,NCOLS,BL1,BU1,IND,IOPT,X,
     *                  RNORMC,RNORM,MODE,RW,IW)
            DO 70 J = 1,NCOLS + MCON
               X(J) = X(J) - XTRUE(J)
   70       CONTINUE
C
            SR = DNRM2(NCOLS+MCON,X,1)
            MPASS = 1
            IF (SR.GT.10.D2*SQRT(D1MACH(4))) MPASS = 0
            IPASS = IPASS*MPASS
            IF (KPRINT.GE.2) THEN
               MSG = 'PASS'
               IF (MPASS.EQ.0) MSG = 'FAIL'
               ITEST = ITEST + 1
               WRITE (LUN, 99999) ITEST, IB, IRHS, SR, MSG
            ENDIF
   80    CONTINUE
   90 CONTINUE
C
C     HERE THE VALUE OF IPASS=1 SAYS THAT DBOCLS HAS PASSED ITS TESTS.
C          THE VALUE OF IPASS=0 SAYS THAT DBOCLS HAS NOT PASSED.
C
      IF(KPRINT.GE.3)
     *WRITE(LUN,'('' IPASS VALUE. (A 1 IS GOOD, 0 IS BAD.)'',I4)')IPASS
      IF(KPRINT.GE.2.AND.IPASS.EQ.0) WRITE(LUN,10789)
      RETURN
C
10789 FORMAT (' ERROR IN DBOCLS OR DBOLS')
99998 FORMAT (' TEST   IB IRHS             SR')
99999 FORMAT (3I5, 1P,E20.6, ' TEST ', A, 'ED.')
      END
