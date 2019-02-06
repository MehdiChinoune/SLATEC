*DECK ISRTQC
      SUBROUTINE ISRTQC (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  ISRTQC
C***SUBSIDIARY
C***PURPOSE  Quick check for SLATEC routines ISORT, IPSORT, IPPERM
C***LIBRARY   SLATEC
C***CATEGORY  N6A
C***TYPE      INTEGER (SSRTQC-S, DSRTQC-D, ISRTQC-I, HSRTQC-H)
C***KEYWORDS  IPPERM, IPSORT, ISORT, QUICK CHECK
C***AUTHOR  Boisvert, Ronald, (NIST)
C***REFERENCES  (NONE)
C***ROUTINES CALLED  IPPERM, IPSORT, ISORT
C***REVISION HISTORY  (YYMMDD)
C   890620  DATE WRITTEN
C   901005  Included test of IPPERM.  (MAM)
C   920511  Added error message tests.  (MAM)
C***END PROLOGUE  ISRTQC
C
      INTEGER N, NTEST
      PARAMETER (N=9,NTEST=4)
C
      LOGICAL FAIL
      INTEGER X(N,NTEST), XS(N,NTEST), Y(N), YC(N)
      INTEGER IX(N,NTEST), IY(N), KFLAG(NTEST), KPRINT, LUN, IPASS, J,
     +        I, KABS, IER, NERR, NUMXER, NN, KKFLAG
C
C     ---------
C     TEST DATA
C     ---------
C
C         X   = TEST VECTOR
C         XS  = TEST VECTOR IN SORTED ORDER
C         IX  = PERMUTATION VECTOR, I.E.  X(IX(J)) = XS(J)
C
      DATA KFLAG(1)       /  2 /
      DATA (X(I,1),I=1,N) / 36, 54, -1, 29,  1, 80, 98, 99, 55/
      DATA (IX(I,1),I=1,N)/  3,  5,  4,  1,  2,  9,  6,  7,  8/
      DATA (XS(I,1),I=1,N)/ -1,  1, 29, 36, 54, 55, 80, 98, 99/
C
      DATA KFLAG(2)       / -1 /
      DATA (X(I,2),I=1,N) /  1,  2,  3,  4,  5,  6,  7,  8,  9/
      DATA (IX(I,2),I=1,N)/  9,  8,  7,  6,  5,  4,  3,  2,  1/
      DATA (XS(I,2),I=1,N)/  9,  8,  7,  6,  5,  4,  3,  2,  1/
C
      DATA KFLAG(3)       / -2 /
      DATA (X(I,3),I=1,N) / -9, -8, -7, -6, -5, -4, -3, -2, -1/
      DATA (IX(I,3),I=1,N)/  9,  8,  7,  6,  5,  4,  3,  2,  1/
      DATA (XS(I,3),I=1,N)/ -1, -2, -3, -4, -5, -6, -7, -8, -9/
C
      DATA KFLAG(4)       /  1 /
      DATA (X(I,4),I=1,N) / 36, 54, -1, 29,  1, 80, 98, 99, 55/
      DATA (IX(I,4),I=1,N)/  3,  5,  4,  1,  2,  9,  6,  7,  8/
      DATA (XS(I,4),I=1,N)/ -1,  1, 29, 36, 54, 55, 80, 98, 99/
C
C***FIRST EXECUTABLE STATEMENT  ISRTQC
      IF ( KPRINT .GE. 2 ) THEN
         WRITE (LUN,2001) '================='
         WRITE (LUN,2002) 'OUTPUT FROM ISRTQC'
         WRITE (LUN,2002) '================='
      ENDIF
      IPASS = 1
C
C     -------------------------------------------------------------
C                          CHECK ISORT
C     -------------------------------------------------------------
C
      DO 200 J=1,NTEST
C
C        ... SETUP PROBLEM
C
         DO 110 I=1,N
            Y(I) = X(I,J)
            YC(I) = X(I,J)
  110    CONTINUE
C
C        ... CALL ROUTINE TO BE TESTED
C
         CALL ISORT(Y,YC,N,KFLAG(J))
C
C        ... EVALUATE RESULTS
C
         KABS = ABS(KFLAG(J))
         FAIL = .FALSE.
         DO 120 I=1,N
            FAIL = FAIL .OR. (Y(I).NE.XS(I,J))
     +                  .OR. ((KABS.EQ.1).AND.(YC(I).NE.X(I,J)))
     +                  .OR. ((KABS.EQ.2).AND.(YC(I).NE.XS(I,J)))
  120    CONTINUE
C
C        ... PRODUCE REQUIRED OUTPUT
C
         IF (FAIL) THEN
             IPASS = 0
             IF (KPRINT .GT. 0) WRITE(LUN,2001) 'ISORT FAILED TEST ',J
         ELSE
             IF (KPRINT .GE. 2) WRITE(LUN,2001) 'ISORT PASSED TEST ',J
         ENDIF
         IF ((FAIL .AND. (KPRINT .GE. 2)) .OR. (KPRINT .GE. 3)) THEN
            WRITE(LUN,2001) '------------------------'
            WRITE(LUN,2002) 'DETAILS OF ISORT TEST ',J
            WRITE(LUN,2002) '------------------------'
            WRITE(LUN,2002) '1ST ARGUMENT (VECTOR TO BE SORTED)'
            WRITE(LUN,2003) '             INPUT = ',(X(I,J),I=1,N)
            WRITE(LUN,2003) '   COMPUTED OUTPUT = ',(Y(I),I=1,N)
            WRITE(LUN,2003) '    CORRECT OUTPUT = ',(XS(I,J),I=1,N)
            WRITE(LUN,2002) '2ND ARGUMENT (VECTOR CARRIED ALONG)'
            WRITE(LUN,2003) '             INPUT = ',(X(I,J),I=1,N)
            WRITE(LUN,2003) '   COMPUTED OUTPUT = ',(YC(I),I=1,N)
            IF (KABS .EQ. 1) THEN
               WRITE(LUN,2003) '    CORRECT OUTPUT = ',(X(I,J),I=1,N)
            ELSE
               WRITE(LUN,2003) '    CORRECT OUTPUT = ',(XS(I,J),I=1,N)
            ENDIF
            WRITE(LUN,2002) '3RD ARGUMENT (VECTOR LENGTH)'
            WRITE(LUN,2004) '             INPUT = ',N
            WRITE(LUN,2002) '4TH ARGUMENT (TYPE OF SORT)'
            WRITE(LUN,2004) '             INPUT = ',KFLAG(J)
         ENDIF
  200 CONTINUE
C
C     -------------------------------------------------------------
C                            CHECK IPSORT
C     -------------------------------------------------------------
C
      DO 300 J=1,NTEST
C
C        ... SETUP PROBLEM
C
         DO 210 I=1,N
            Y(I) = X(I,J)
  210    CONTINUE
C
C        ... CALL ROUTINE TO BE TESTED
C
         CALL IPSORT(Y,N,IY,KFLAG(J),IER)
C
C        ... EVALUATE RESULTS
C
         KABS = ABS(KFLAG(J))
         FAIL = .FALSE. .OR. (IER.GT.0)
         DO 220 I=1,N
            FAIL = FAIL .OR. (IY(I).NE.IX(I,J))
     +                  .OR. ((KABS.EQ.1).AND.(Y(I).NE.X(I,J)))
     +                  .OR. ((KABS.EQ.2).AND.(Y(I).NE.XS(I,J)))
  220    CONTINUE
C
C        ... PRODUCE REQUIRED OUTPUT
C
         IF (FAIL) THEN
            IPASS = 0
            IF (KPRINT .GT. 0) WRITE(LUN,2001) 'IPSORT FAILED TEST ',J
         ELSE
            IF (KPRINT .GE. 2) WRITE(LUN,2001) 'IPSORT PASSED TEST ',J
         ENDIF
         IF ((FAIL .AND. (KPRINT .GE. 2)) .OR. (KPRINT .GE. 3)) THEN
            WRITE(LUN,2001) '-------------------------'
            WRITE(LUN,2002) 'DETAILS OF IPSORT TEST ',J
            WRITE(LUN,2002) '-------------------------'
            WRITE(LUN,2002) '1ST ARGUMENT (VECTOR TO BE SORTED)'
            WRITE(LUN,2003) '             INPUT = ',(X(I,J),I=1,N)
            WRITE(LUN,2003) '   COMPUTED OUTPUT = ',(Y(I),I=1,N)
            IF (KABS .EQ. 1) THEN
               WRITE(LUN,2003) '    CORRECT OUTPUT = ',(X(I,J),I=1,N)
            ELSE
               WRITE(LUN,2003) '    CORRECT OUTPUT = ',(XS(I,J),I=1,N)
            ENDIF
            WRITE(LUN,2002) '2ND ARGUMENT (VECTOR LENGTH)'
            WRITE(LUN,2004) '             INPUT = ',N
            WRITE(LUN,2002) '3RD ARGUMENT (PERMUTATION VECTOR)'
            WRITE(LUN,2004) '   COMPUTED OUTPUT = ',(IY(I),I=1,N)
            WRITE(LUN,2004) '    CORRECT OUTPUT = ',(IX(I,J),I=1,N)
            WRITE(LUN,2002) '4TH ARGUMENT (TYPE OF SORT)'
            WRITE(LUN,2004) '             INPUT = ',KFLAG(J)
         ENDIF
C
  300 CONTINUE
C
C     ... TEST ERROR MESSAGES
C
      IF(KPRINT.LE.2)THEN
         CALL XSETF(0)
      ELSE
         CALL XSETF(-1)
      ENDIF
C
      NN=-1
      KKFLAG=1
      IF(KPRINT.GE.3)WRITE(LUN,*)
      CALL XERCLR
      CALL IPSORT(Y,NN,IY,KKFLAG,IER)
      IF(NUMXER(NERR).NE.IER)IPASS=0
C
      NN=1
      KKFLAG=0
      IF(KPRINT.GE.3)WRITE(LUN,*)
      CALL XERCLR
      CALL IPSORT(Y,NN,IY,KKFLAG,IER)
      IF(NUMXER(NERR).NE.IER)IPASS=0
C
      IF((KPRINT.GE.2).AND.(IPASS.EQ.1))THEN
         WRITE(LUN,*)
         WRITE(LUN,*)' IPSORT PASSED ERROR MESSAGE TESTS'
      ELSE IF((KPRINT.GE.1).AND.(IPASS.EQ.0))THEN
         WRITE(LUN,*)
         WRITE(LUN,*)' IPSORT FAILED ERROR MESSAGE TESTS'
      ENDIF
C
C     -------------------------------------------------------------
C                            CHECK IPPERM
C     -------------------------------------------------------------
C
      DO 400 J=1,NTEST
C
C        ... SETUP PROBLEM
C
         KABS = ABS(KFLAG(J))
         DO 310 I=1,N
            Y(I) = X(I,J)
            IF(KABS.EQ.1)THEN
               IY(I) = I
            ELSE
               IY(I) = IX(I,J)
            ENDIF
  310    CONTINUE
C
C        ... CALL ROUTINE TO BE TESTED
C
         CALL IPPERM(Y,N,IY,IER)
C
C        ... EVALUATE RESULTS
C
         FAIL = .FALSE. .OR. (IER.GT.0)
         DO 320 I=1,N
            FAIL = FAIL .OR. ((KABS.EQ.1).AND.(IY(I).NE.I))
     +                  .OR. ((KABS.EQ.2).AND.(IY(I).NE.IX(I,J)))
     +                  .OR. ((KABS.EQ.1).AND.(Y(I).NE.X(I,J)))
     +                  .OR. ((KABS.EQ.2).AND.(Y(I).NE.XS(I,J)))
  320    CONTINUE
C
C        ... PRODUCE REQUIRED OUTPUT
C
         IF (FAIL) THEN
            IPASS = 0
            IF (KPRINT.GT.0) WRITE(LUN,2001)'IPPERM FAILED TEST ',J
         ELSE
            IF (KPRINT.GE.2) WRITE(LUN,2001)'IPPERM PASSED TEST ',J
         ENDIF
         IF ((FAIL .AND. (KPRINT.GE.2)) .OR. (KPRINT.GE.3)) THEN
            WRITE(LUN,2001)'------------------------'
            WRITE(LUN,2002)'DETAILS OF IPPERM TEST',J
            WRITE(LUN,2002)'------------------------'
            WRITE(LUN,2002)'1ST ARGUMENT (VECTOR TO BE PERMUTED)'
            WRITE(LUN,2003)'             INPUT =',(X(I,J),I=1,N)
            WRITE(LUN,2003)'   COMPUTED OUTPUT =',(Y(I),I=1,N)
            IF(KABS.EQ.1)THEN
               WRITE(LUN,2003)'    CORRECT OUTPUT =',(X(I,J),I=1,N)
            ELSE
               WRITE(LUN,2003)'    CORRECT OUTPUT =',(XS(I,J),I=1,N)
            ENDIF
            WRITE(LUN,2002)'2ND ARGUMENT (VECTOR LENGTH)'
            WRITE(LUN,2004)'             INPUT =',N
            WRITE(LUN,2002)'3RD ARGUMENT (PERMUTATION VECTOR)'
            WRITE(LUN,2004)'             INPUT =',(IY(I),I=1,N)
            WRITE(LUN,2002)'4TH ARGUMENT (ERROR FLAG)'
            WRITE(LUN,2004)'             OUTPUT =',IER
         ENDIF
C
  400 CONTINUE
C
C     ... TEST ERROR MESSAGES
C
      IF(KPRINT.LE.2)THEN
         CALL XSETF(0)
      ELSE
         CALL XSETF(-1)
      ENDIF
C
      NN=-1
      IF(KPRINT.GE.3)WRITE(LUN,*)
      CALL XERCLR
      CALL IPPERM(Y,NN,IY,IER)
      IF(NUMXER(NERR).NE.IER)IPASS=0
C
      NN=1
      IY(1)=5
      IF(KPRINT.GE.3)WRITE(LUN,*)
      CALL XERCLR
      CALL IPPERM(Y,NN,IY,IER)
      IF(NUMXER(NERR).NE.IER)IPASS=0
C
      IF((KPRINT.GE.2).AND.(IPASS.EQ.1))THEN
         WRITE(LUN,*)
         WRITE(LUN,*)' IPPERM PASSED ERROR MESSAGE TESTS'
      ELSE IF((KPRINT.GE.1).AND.(IPASS.EQ.0))THEN
         WRITE(LUN,*)
         WRITE(LUN,*)' IPPERM FAILED ERROR MESSAGE TESTS'
      ENDIF
C
      RETURN
C
 2001 FORMAT(/ 1X,A,I2)
 2002 FORMAT(1X,A,I2)
 2003 FORMAT(1X,A,9I4)
 2004 FORMAT(1X,A,9I4)
      END
