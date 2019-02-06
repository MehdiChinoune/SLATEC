*DECK DQCK
      SUBROUTINE DQCK (LUN, KPRINT, NERR)
C***BEGIN PROLOGUE  DQCK
C***PURPOSE  Quick check for DPOFS AND DNBFS.
C***LIBRARY   SLATEC
C***KEYWORDS  QUICK CHECK
C***AUTHOR  Voorhees, E. A., (LANL)
C***DESCRIPTION
C
C    QUICK CHECK SUBROUTINE DQCK TESTS THE EXECUTION OF THE
C    SLATEC SUBROUTINES DPOFS AND DNBFS.
C    A TITLE LINE AND A SUMMARY LINE ARE ALWAYS OUTPUTTED.
C
C    THE SUMMARY LINE GIVES A COUNT OF THE NUMBER OF
C    PROBLEMS ENCOUNTERED IN THE TEST IF ANY EXIST.  DQCK
C    CHECKS COMPUTED VS. EXACT SOLUTIONS TO AGREE TO
C    WITHIN 0.8 TIMES THE WORD LENGTH OF THE COMPUTER
C    (1.6 IF DOUBLE PRECISION) FOR CASE 1.  DQCK ALSO
C    TESTS ERROR HANDLING BY THE SUBROUTINE (CALLS TO
C    XERMSG (DQCK SETS IFLAG/KONTRL TO 0))
C    USING A SINGULAR MATRIX FOR CASE 2.  EACH EXECUTION
C    PROBLEM DETECTED BY DQCK RESULTS IN AN ADDITIONAL
C    EXPLANATORY LINE OF OUTPUT.
C
C    DQCK REQUIRES NO INPUT ARGUMENTS.
C    ON RETURN, NERR (INTEGER TYPE) CONTAINS THE TOTAL COUNT
C    OF ALL PROBLEMS DETECTED BY DQCK.
C
C***ROUTINES CALLED  D1MACH, DNBFS, DPOFS
C***REVISION HISTORY  (YYMMDD)
C   801002  DATE WRITTEN
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   901009  Restructured using IF-THEN-ELSE-ENDIF, cleaned up FORMATs,
C           including removing an illegal character from column 1, and
C           editorial changes.  (RWC)
C***END PROLOGUE  DQCK
      DOUBLE PRECISION A(4,4),AT(5,4),ABE(5,7),ABET(5,7),B(4),BT(4),
     1 C(4),WORK(35),SIGN,D1MACH
      REAL R,DELX,DELMAX
      CHARACTER*4 LIST(2)
      INTEGER LDA,N,ML,MU,IND,IWORK(4),NERR,I,J,J1,J2,JD,MLP,K,KCASE,
     1 KPROG
      DATA A/5.0D0,4.0D0,1.0D0,1.0D0,4.0D0,5.0D0,1.0D0,1.0D0,
     1 1.0D0,1.0D0,4.0D0,2.0D0,1.0D0,1.0D0,2.0D0,4.0D0/
      DATA LIST/'POFS', 'NBFS'/
C***FIRST EXECUTABLE STATEMENT  DQCK
      IF (KPRINT.GE.3) WRITE (LUN,800)
      LDA = 5
      N = 4
      ML = 2
      MU = 1
      JD = 2*ML+MU+1
      NERR = 0
      R = D1MACH(4)**0.8E0
C
C     COMPUTE C VECTOR.
C
      SIGN = 1.0D0
      DO 10 I=1,N
         C(I) = SIGN/I
         SIGN = -SIGN
   10 CONTINUE
C
C     CASE 1 FOR WELL-CONDITIONED MATRIX, CASE 2 FOR SINGULAR MATRIX.
C
      DO 170 KCASE=1,2
         DO 140 KPROG=1,2
C           SET VECTOR B TO ZERO.
            DO 11 I=1,N
               B(I) = 0.0D0
   11       CONTINUE
C
C           FORM VECTOR B FOR NON-BANDED.
C
            IF (KPROG.EQ.1) THEN
               DO 13 I=1,N
                  DO 12 J=1,N
                     B(I) = B(I)+A(I,J)*C(J)
   12             CONTINUE
   13          CONTINUE
            ELSE
C
C              FORM ABE(NB ARRAY) FROM MATRIX A
C              AND FORM VECTOR B FOR BANDED.
C
               DO 30 J=1,JD
                  DO 20 I=1,N
                     ABE(I,J) = 0.0D0
   20             CONTINUE
   30          CONTINUE
C
               MLP = ML+1
               DO 50 I=1,N
                  J1 = MAX(1,I-ML)
                  J2 = MIN(N,I+MU)
                  DO 40 J=J1,J2
                     K = J-I+MLP
                     ABE(I,K) = A(I,J)
                     B(I) = B(I)+(A(I,J)*C(J))
   40             CONTINUE
   50          CONTINUE
            ENDIF
C
C           FORM BT FROM B, AT FROM A, AND ABET FROM ABE.
C
            DO 60 I=1,N
               BT(I) = B(I)
               DO 58 J=1,N
                  AT(I,J) = A(I,J)
   58          CONTINUE
   60       CONTINUE
C
            DO 80 J=1,JD
               DO 70 I=1,N
                  ABET(I,J) = ABE(I,J)
   70          CONTINUE
   80       CONTINUE
C
C           MAKE AT AND ABET SINGULAR FOR CASE  =  2
C
            IF (KCASE.EQ.2) THEN
               DO 88 J=1,N
                  AT(1,J) = 0.0D0
   88          CONTINUE
C
               DO 90 J=1,JD
                  ABET(1,J) = 0.0D0
   90          CONTINUE
            ENDIF
C
C           SOLVE FOR X
C
            IF (KPROG.EQ.1) CALL DPOFS (AT,LDA,N,BT,1,IND,WORK)
            IF (KPROG.EQ.2) CALL DNBFS (ABET,LDA,N,ML,MU,BT,1,IND,WORK,
     *         IWORK)
C
C           COMPARE EXACT AND COMPUTED SOLUTIONS FOR CASE 1
C
            IF (KCASE.EQ.1) THEN
               DELMAX = 0.0E0
               DO 110 I=1,N
                  DELX = ABS(BT(I)-C(I))
                  DELMAX = MAX(DELMAX,DELX)
  110          CONTINUE
C
               IF (R.LE.DELMAX) THEN
                  NERR = NERR+1
                  WRITE (LUN,801) LIST(KPROG),KCASE,DELMAX
               ENDIF
            ELSE
C
C              CHECK CONTROL FOR SINGULAR MATRIX FOR CASE 2
C
               IF (IND.NE.-4) THEN
                  NERR = NERR+1
                  WRITE (LUN,802) LIST(KPROG),KCASE,IND
               ENDIF
            ENDIF
  140    CONTINUE
  170 CONTINUE
C
C     SUMMARY PRINT
C
      IF (NERR.NE.0) WRITE (LUN,803) NERR
      IF (KPRINT.GE.2 .AND. NERR.EQ.0) WRITE (LUN,804)
      RETURN
C
  800 FORMAT (/' *    DQCK - QUICK CHECK FOR  DPOFS AND DNBFS'/)
  801 FORMAT ('   PROBLEM WITH D', A, ', CASE ', I1,
     1   '.  MAX ABS ERROR OF', E11.4/)
  802 FORMAT ('   PROBLEM WITH D', A, ', CASE ', I1, '.  IND = ', I2,
     1   ' INSTEAD OF -4'/)
  803 FORMAT (/' **** DQCK DETECTED A TOTAL OF ', I2,' PROBLEMS. ****'/)
  804 FORMAT ('     DQCK DETECTED NO PROBLEMS.'/)
      END
