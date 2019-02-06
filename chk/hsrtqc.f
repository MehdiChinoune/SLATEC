*DECK HSRTQC
      SUBROUTINE HSRTQC (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  HSRTQC
C***SUBSIDIARY
C***PURPOSE  Quick check for SLATEC routine HPSORT, HPPERM
C***LIBRARY   SLATEC
C***CATEGORY  N6A
C***TYPE      CHARACTER (SSRTQC-S, DSRTQC-D, ISRTQC-I, HSRTQC-H)
C***KEYWORDS  HPPERM, HPSORT, QUICK CHECK
C***AUTHOR  Boisvert, Ronald, (NIST)
C***REFERENCES  (NONE)
C***ROUTINES CALLED  HPPERM, HPSORT
C***REVISION HISTORY  (YYMMDD)
C   890620  DATE WRITTEN
C   901005  Included test of HPPERM.  (MAM)
C   920511  Added error message tests.  (MAM)
C***END PROLOGUE  HSRTQC
C
      INTEGER N, NTEST
      PARAMETER (N=9,NTEST=4)
C
      LOGICAL FAIL
      CHARACTER*1 SHORT
      CHARACTER*2 X(N,NTEST), XS(N,NTEST), Y(N), WORK(N)
      INTEGER IX(N,NTEST), IY(N), KFLAG(NTEST), KPRINT, LUN, IPASS, J,
     +        I, KABS, IER, NERR, NUMXER, NN, KKFLAG, STRBEG, STREND
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
      DATA (X(I,1),I=1,N) /'AC','AZ','AD','AA','AB','ZZ','ZA','ZX','ZY'/
      DATA (IX(I,1),I=1,N)/  4,   5,   1,   3,   2,   7,   8,   9,   6/
      DATA (XS(I,1),I=1,N)/'AA','AB','AC','AD','AZ','ZA','ZX','ZY','ZZ'/
C
      DATA KFLAG(2)       / -1 /
      DATA (X(I,2),I=1,N) /'AA','BB','CC','DD','EE','FF','GG','HH','II'/
      DATA (IX(I,2),I=1,N)/  9,   8,   7,   6,   5,   4,   3,   2,   1 /
      DATA (XS(I,2),I=1,N)/'II','HH','GG','FF','EE','DD','CC','BB','AA'/
C
      DATA KFLAG(3)       / -2 /
      DATA (X(I,3),I=1,N) /'AA','BB','CC','DD','EE','FF','GG','HH','II'/
      DATA (IX(I,3),I=1,N)/  9,   8,   7,   6,   5,   4,   3,   2,   1/
      DATA (XS(I,3),I=1,N)/'II','HH','GG','FF','EE','DD','CC','BB','AA'/
C
      DATA KFLAG(4)       /  1 /
      DATA (X(I,4),I=1,N) /'AC','AZ','AD','AA','AB','ZZ','ZA','ZX','ZY'/
      DATA (IX(I,4),I=1,N)/  4,   5,   1,   3,   2,   7,   8,   9,   6/
      DATA (XS(I,4),I=1,N)/'AA','AB','AC','AD','AZ','ZA','ZX','ZY','ZZ'/
C
C***FIRST EXECUTABLE STATEMENT  HSRTQC
      IF ( KPRINT .GE. 2 ) THEN
         WRITE (LUN,2001) '================='
         WRITE (LUN,2002) 'OUTPUT FROM HSRTQC'
         WRITE (LUN,2002) '================='
      ENDIF
      IPASS = 1
C
C     -------------------------------------------------------------
C                            CHECK HPSORT
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
         CALL HPSORT(Y,N,1,2,IY,KFLAG(J),WORK,IER)
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
            IF (KPRINT .GT. 0) WRITE(LUN,2001) 'HPSORT FAILED TEST ',J
         ELSE
            IF (KPRINT .GE. 2) WRITE(LUN,2001) 'HPSORT PASSED TEST ',J
         ENDIF
         IF ((FAIL .AND. (KPRINT .GE. 2)) .OR. (KPRINT .GE. 3)) THEN
            WRITE(LUN,2001) '-------------------------'
            WRITE(LUN,2002) 'DETAILS OF HPSORT TEST ',J
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
      STRBEG=1
      STREND=2
      KKFLAG=1
      IF(KPRINT.GE.3)WRITE(LUN,*)
      CALL XERCLR
      CALL HPSORT(Y,NN,STRBEG,STREND,IY,KKFLAG,WORK,IER)
      IF(NUMXER(NERR).NE.IER)IPASS=0
C
      NN=1
      STRBEG=1
      STREND=2
      KKFLAG=0
      IF(KPRINT.GE.3)WRITE(LUN,*)
      CALL XERCLR
      CALL HPSORT(Y,NN,STRBEG,STREND,IY,KKFLAG,WORK,IER)
      IF(NUMXER(NERR).NE.IER)IPASS=0
C
      NN=1
      STRBEG=1
      STREND=2
      KKFLAG=1
      IF(KPRINT.GE.3)WRITE(LUN,*)
      CALL XERCLR
      CALL HPSORT(Y,NN,STRBEG,STREND,IY,KKFLAG,SHORT,IER)
      IF(NUMXER(NERR).NE.IER)IPASS=0
C
      NN=1
      STRBEG=2
      STREND=1
      KKFLAG=1
      IF(KPRINT.GE.3)WRITE(LUN,*)
      CALL XERCLR
      CALL HPSORT(Y,NN,STRBEG,STREND,IY,KKFLAG,WORK,IER)
      IF(NUMXER(NERR).NE.IER)IPASS=0
C
      NN=1
      STRBEG=-1
      STREND=2
      KKFLAG=1
      IF(KPRINT.GE.3)WRITE(LUN,*)
      CALL XERCLR
      CALL HPSORT(Y,NN,STRBEG,STREND,IY,KKFLAG,WORK,IER)
      IF(NUMXER(NERR).NE.IER)IPASS=0
C
      NN=1
      STRBEG=1
      STREND=3
      KKFLAG=1
      IF(KPRINT.GE.3)WRITE(LUN,*)
      CALL XERCLR
      CALL HPSORT(Y,NN,STRBEG,STREND,IY,KKFLAG,WORK,IER)
      IF(NUMXER(NERR).NE.IER)IPASS=0
C
      IF((KPRINT.GE.2).AND.(IPASS.EQ.1))THEN
         WRITE(LUN,*)
         WRITE(LUN,*)' HPSORT PASSED ERROR MESSAGE TESTS'
      ELSE IF((KPRINT.GE.1).AND.(IPASS.EQ.0))THEN
         WRITE(LUN,*)
         WRITE(LUN,*)' HPSORT FAILED ERROR MESSAGE TESTS'
      ENDIF
C
C     -------------------------------------------------------------
C                            CHECK HPPERM
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
         CALL HPPERM(Y,N,IY,WORK,IER)
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
            IF (KPRINT.GT.0) WRITE(LUN,2001)'HPPERM FAILED TEST ',J
         ELSE
            IF (KPRINT.GE.2) WRITE(LUN,2001)'HPPERM PASSED TEST ',J
         ENDIF
         IF ((FAIL .AND. (KPRINT.GE.2)) .OR. (KPRINT.GE.3)) THEN
            WRITE(LUN,2001)'------------------------'
            WRITE(LUN,2002)'DETAILS OF HPPERM TEST',J
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
      CALL HPPERM(Y,NN,IY,WORK,IER)
      IF(NUMXER(NERR).NE.IER)IPASS=0
C
      NN=1
      IY(1)=5
      IF(KPRINT.GE.3)WRITE(LUN,*)
      CALL XERCLR
      CALL HPPERM(Y,NN,IY,WORK,IER)
      IF(NUMXER(NERR).NE.IER)IPASS=0
C
      IF((KPRINT.GE.2).AND.(IPASS.EQ.1))THEN
         WRITE(LUN,*)
         WRITE(LUN,*)' HPPERM PASSED ERROR MESSAGE TESTS'
      ELSE IF((KPRINT.GE.1).AND.(IPASS.EQ.0))THEN
         WRITE(LUN,*)
         WRITE(LUN,*)' HPPERM FAILED ERROR MESSAGE TESTS'
      ENDIF
C
      RETURN
C
 2001 FORMAT(/ 1X,A,I2)
 2002 FORMAT(1X,A,I2)
 2003 FORMAT(1X,A,9(2X,A2))
 2004 FORMAT(1X,A,9I4)
      END
