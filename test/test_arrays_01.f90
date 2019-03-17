MODULE TEST54_MOD
  IMPLICIT NONE

CONTAINS
  !DECK ISRTQC
  SUBROUTINE ISRTQC(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !***BEGIN PROLOGUE  ISRTQC
    !***SUBSIDIARY
    !***PURPOSE  Quick check for SLATEC routines ISORT, IPSORT, IPPERM
    !***LIBRARY   SLATEC
    !***CATEGORY  N6A
    !***TYPE      INTEGER (SSRTQC-S, DSRTQC-D, ISRTQC-I, HSRTQC-H)
    !***KEYWORDS  IPPERM, IPSORT, ISORT, QUICK CHECK
    !***AUTHOR  Boisvert, Ronald, (NIST)
    !***REFERENCES  (NONE)
    !***ROUTINES CALLED  IPPERM, IPSORT, ISORT
    !***REVISION HISTORY  (YYMMDD)
    !   890620  DATE WRITTEN
    !   901005  Included test of IPPERM.  (MAM)
    !   920511  Added error message tests.  (MAM)
    !***END PROLOGUE  ISRTQC
    !
    INTEGER N, NTEST
    PARAMETER (N=9,NTEST=4)
    !
    LOGICAL fail
    INTEGER x(N,NTEST), xs(N,NTEST), y(N), yc(N)
    INTEGER ix(N,NTEST), iy(N), kflag(NTEST), Kprint, Lun, Ipass, j, &
      i, kabs, ier, nerr, NUMXER, nn, kkflag
    !
    !     ---------
    !     TEST DATA
    !     ---------
    !
    !         X   = TEST VECTOR
    !         XS  = TEST VECTOR IN SORTED ORDER
    !         IX  = PERMUTATION VECTOR, I.E.  X(IX(J)) = XS(J)
    !
    DATA kflag(1)/2/
    DATA (x(i,1),i=1,N)/36, 54, -1, 29, 1, 80, 98, 99, 55/
    DATA (ix(i,1),i=1,N)/3, 5, 4, 1, 2, 9, 6, 7, 8/
    DATA (xs(i,1),i=1,N)/ - 1, 1, 29, 36, 54, 55, 80, 98, 99/
    !
    DATA kflag(2)/ - 1/
    DATA (x(i,2),i=1,N)/1, 2, 3, 4, 5, 6, 7, 8, 9/
    DATA (ix(i,2),i=1,N)/9, 8, 7, 6, 5, 4, 3, 2, 1/
    DATA (xs(i,2),i=1,N)/9, 8, 7, 6, 5, 4, 3, 2, 1/
    !
    DATA kflag(3)/ - 2/
    DATA (x(i,3),i=1,N)/ - 9, -8, -7, -6, -5, -4, -3, -2, -1/
    DATA (ix(i,3),i=1,N)/9, 8, 7, 6, 5, 4, 3, 2, 1/
    DATA (xs(i,3),i=1,N)/ - 1, -2, -3, -4, -5, -6, -7, -8, -9/
    !
    DATA kflag(4)/1/
    DATA (x(i,4),i=1,N)/36, 54, -1, 29, 1, 80, 98, 99, 55/
    DATA (ix(i,4),i=1,N)/3, 5, 4, 1, 2, 9, 6, 7, 8/
    DATA (xs(i,4),i=1,N)/ - 1, 1, 29, 36, 54, 55, 80, 98, 99/
    !
    !***FIRST EXECUTABLE STATEMENT  ISRTQC
    IF ( Kprint>=2 ) THEN
      WRITE (Lun,99001) '================='
      WRITE (Lun,99002) 'OUTPUT FROM ISRTQC'
      WRITE (Lun,99002) '================='
    ENDIF
    Ipass = 1
    !
    !     -------------------------------------------------------------
    !                          CHECK ISORT
    !     -------------------------------------------------------------
    !
    DO j = 1, NTEST
      !
      !        ... SETUP PROBLEM
      !
      DO i = 1, N
        y(i) = x(i,j)
        yc(i) = x(i,j)
      ENDDO
      !
      !        ... CALL ROUTINE TO BE TESTED
      !
      CALL ISORT(y,yc,N,kflag(j))
      !
      !        ... EVALUATE RESULTS
      !
      kabs = ABS(kflag(j))
      fail = .FALSE.
      DO i = 1, N
        fail = fail .OR. (y(i)/=xs(i,j)) .OR. ((kabs==1).AND.(yc(i)/=x(i,j)))&
          .OR. ((kabs==2).AND.(yc(i)/=xs(i,j)))
      ENDDO
      !
      !        ... PRODUCE REQUIRED OUTPUT
      !
      IF ( fail ) THEN
        Ipass = 0
        IF ( Kprint>0 ) WRITE (Lun,99001) 'ISORT FAILED TEST ', j
      ELSE
        IF ( Kprint>=2 ) WRITE (Lun,99001) 'ISORT PASSED TEST ', j
      ENDIF
      IF ( (fail.AND.(Kprint>=2)).OR.(Kprint>=3) ) THEN
        WRITE (Lun,99001) '------------------------'
        WRITE (Lun,99002) 'DETAILS OF ISORT TEST ', j
        WRITE (Lun,99002) '------------------------'
        WRITE (Lun,99002) '1ST ARGUMENT (VECTOR TO BE SORTED)'
        WRITE (Lun,99003) '             INPUT = ', (x(i,j),i=1,N)
        WRITE (Lun,99003) '   COMPUTED OUTPUT = ', (y(i),i=1,N)
        WRITE (Lun,99003) '    CORRECT OUTPUT = ', (xs(i,j),i=1,N)
        WRITE (Lun,99002) '2ND ARGUMENT (VECTOR CARRIED ALONG)'
        WRITE (Lun,99003) '             INPUT = ', (x(i,j),i=1,N)
        WRITE (Lun,99003) '   COMPUTED OUTPUT = ', (yc(i),i=1,N)
        IF ( kabs==1 ) THEN
          WRITE (Lun,99003) '    CORRECT OUTPUT = ', (x(i,j),i=1,N)
        ELSE
          WRITE (Lun,99003) '    CORRECT OUTPUT = ', (xs(i,j),i=1,N)
        ENDIF
        WRITE (Lun,99002) '3RD ARGUMENT (VECTOR LENGTH)'
        WRITE (Lun,99004) '             INPUT = ', N
        WRITE (Lun,99002) '4TH ARGUMENT (TYPE OF SORT)'
        WRITE (Lun,99004) '             INPUT = ', kflag(j)
      ENDIF
    ENDDO
    !
    !     -------------------------------------------------------------
    !                            CHECK IPSORT
    !     -------------------------------------------------------------
    !
    DO j = 1, NTEST
      !
      !        ... SETUP PROBLEM
      !
      DO i = 1, N
        y(i) = x(i,j)
      ENDDO
      !
      !        ... CALL ROUTINE TO BE TESTED
      !
      CALL IPSORT(y,N,iy,kflag(j),ier)
      !
      !        ... EVALUATE RESULTS
      !
      kabs = ABS(kflag(j))
      fail = .FALSE. .OR. (ier>0)
      DO i = 1, N
        fail = fail .OR. (iy(i)/=ix(i,j)) .OR. ((kabs==1).AND.(y(i)/=x(i,j)))&
          .OR. ((kabs==2).AND.(y(i)/=xs(i,j)))
      ENDDO
      !
      !        ... PRODUCE REQUIRED OUTPUT
      !
      IF ( fail ) THEN
        Ipass = 0
        IF ( Kprint>0 ) WRITE (Lun,99001) 'IPSORT FAILED TEST ', j
      ELSE
        IF ( Kprint>=2 ) WRITE (Lun,99001) 'IPSORT PASSED TEST ', j
      ENDIF
      IF ( (fail.AND.(Kprint>=2)).OR.(Kprint>=3) ) THEN
        WRITE (Lun,99001) '-------------------------'
        WRITE (Lun,99002) 'DETAILS OF IPSORT TEST ', j
        WRITE (Lun,99002) '-------------------------'
        WRITE (Lun,99002) '1ST ARGUMENT (VECTOR TO BE SORTED)'
        WRITE (Lun,99003) '             INPUT = ', (x(i,j),i=1,N)
        WRITE (Lun,99003) '   COMPUTED OUTPUT = ', (y(i),i=1,N)
        IF ( kabs==1 ) THEN
          WRITE (Lun,99003) '    CORRECT OUTPUT = ', (x(i,j),i=1,N)
        ELSE
          WRITE (Lun,99003) '    CORRECT OUTPUT = ', (xs(i,j),i=1,N)
        ENDIF
        WRITE (Lun,99002) '2ND ARGUMENT (VECTOR LENGTH)'
        WRITE (Lun,99004) '             INPUT = ', N
        WRITE (Lun,99002) '3RD ARGUMENT (PERMUTATION VECTOR)'
        WRITE (Lun,99004) '   COMPUTED OUTPUT = ', (iy(i),i=1,N)
        WRITE (Lun,99004) '    CORRECT OUTPUT = ', (ix(i,j),i=1,N)
        WRITE (Lun,99002) '4TH ARGUMENT (TYPE OF SORT)'
        WRITE (Lun,99004) '             INPUT = ', kflag(j)
      ENDIF
      !
    ENDDO
    !
    !     ... TEST ERROR MESSAGES
    !
    IF ( Kprint<=2 ) THEN
      CALL XSETF(0)
    ELSE
      CALL XSETF(-1)
    ENDIF
    !
    nn = -1
    kkflag = 1
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL IPSORT(y,nn,iy,kkflag,ier)
    IF ( NUMXER(nerr)/=ier ) Ipass = 0
    !
    nn = 1
    kkflag = 0
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL IPSORT(y,nn,iy,kkflag,ier)
    IF ( NUMXER(nerr)/=ier ) Ipass = 0
    !
    IF ( (Kprint>=2).AND.(Ipass==1) ) THEN
      WRITE (Lun,*)
      WRITE (Lun,*) ' IPSORT PASSED ERROR MESSAGE TESTS'
    ELSEIF ( (Kprint>=1).AND.(Ipass==0) ) THEN
      WRITE (Lun,*)
      WRITE (Lun,*) ' IPSORT FAILED ERROR MESSAGE TESTS'
    ENDIF
    !
    !     -------------------------------------------------------------
    !                            CHECK IPPERM
    !     -------------------------------------------------------------
    !
    DO j = 1, NTEST
      !
      !        ... SETUP PROBLEM
      !
      kabs = ABS(kflag(j))
      DO i = 1, N
        y(i) = x(i,j)
        IF ( kabs==1 ) THEN
          iy(i) = i
        ELSE
          iy(i) = ix(i,j)
        ENDIF
      ENDDO
      !
      !        ... CALL ROUTINE TO BE TESTED
      !
      CALL IPPERM(y,N,iy,ier)
      !
      !        ... EVALUATE RESULTS
      !
      fail = .FALSE. .OR. (ier>0)
      DO i = 1, N
        fail = fail .OR. ((kabs==1).AND.(iy(i)/=i)) .OR. &
          ((kabs==2).AND.(iy(i)/=ix(i,j))) .OR. &
          ((kabs==1).AND.(y(i)/=x(i,j))) .OR. &
          ((kabs==2).AND.(y(i)/=xs(i,j)))
      ENDDO
      !
      !        ... PRODUCE REQUIRED OUTPUT
      !
      IF ( fail ) THEN
        Ipass = 0
        IF ( Kprint>0 ) WRITE (Lun,99001) 'IPPERM FAILED TEST ', j
      ELSE
        IF ( Kprint>=2 ) WRITE (Lun,99001) 'IPPERM PASSED TEST ', j
      ENDIF
      IF ( (fail.AND.(Kprint>=2)).OR.(Kprint>=3) ) THEN
        WRITE (Lun,99001) '------------------------'
        WRITE (Lun,99002) 'DETAILS OF IPPERM TEST', j
        WRITE (Lun,99002) '------------------------'
        WRITE (Lun,99002) '1ST ARGUMENT (VECTOR TO BE PERMUTED)'
        WRITE (Lun,99003) '             INPUT =', (x(i,j),i=1,N)
        WRITE (Lun,99003) '   COMPUTED OUTPUT =', (y(i),i=1,N)
        IF ( kabs==1 ) THEN
          WRITE (Lun,99003) '    CORRECT OUTPUT =', (x(i,j),i=1,N)
        ELSE
          WRITE (Lun,99003) '    CORRECT OUTPUT =', (xs(i,j),i=1,N)
        ENDIF
        WRITE (Lun,99002) '2ND ARGUMENT (VECTOR LENGTH)'
        WRITE (Lun,99004) '             INPUT =', N
        WRITE (Lun,99002) '3RD ARGUMENT (PERMUTATION VECTOR)'
        WRITE (Lun,99004) '             INPUT =', (iy(i),i=1,N)
        WRITE (Lun,99002) '4TH ARGUMENT (ERROR FLAG)'
        WRITE (Lun,99004) '             OUTPUT =', ier
      ENDIF
      !
    ENDDO
    !
    !     ... TEST ERROR MESSAGES
    !
    IF ( Kprint<=2 ) THEN
      CALL XSETF(0)
    ELSE
      CALL XSETF(-1)
    ENDIF
    !
    nn = -1
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL IPPERM(y,nn,iy,ier)
    IF ( NUMXER(nerr)/=ier ) Ipass = 0
    !
    nn = 1
    iy(1) = 5
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL IPPERM(y,nn,iy,ier)
    IF ( NUMXER(nerr)/=ier ) Ipass = 0
    !
    IF ( (Kprint>=2).AND.(Ipass==1) ) THEN
      WRITE (Lun,*)
      WRITE (Lun,*) ' IPPERM PASSED ERROR MESSAGE TESTS'
    ELSEIF ( (Kprint>=1).AND.(Ipass==0) ) THEN
      WRITE (Lun,*)
      WRITE (Lun,*) ' IPPERM FAILED ERROR MESSAGE TESTS'
    ENDIF
    !
    RETURN
    !
    99001 FORMAT (/1X,A,I2)
    99002 FORMAT (1X,A,I2)
    99003 FORMAT (1X,A,9I4)
    99004 FORMAT (1X,A,9I4)
  END SUBROUTINE ISRTQC
  !DECK HSRTQC
  SUBROUTINE HSRTQC(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !***BEGIN PROLOGUE  HSRTQC
    !***SUBSIDIARY
    !***PURPOSE  Quick check for SLATEC routine HPSORT, HPPERM
    !***LIBRARY   SLATEC
    !***CATEGORY  N6A
    !***TYPE      CHARACTER (SSRTQC-S, DSRTQC-D, ISRTQC-I, HSRTQC-H)
    !***KEYWORDS  HPPERM, HPSORT, QUICK CHECK
    !***AUTHOR  Boisvert, Ronald, (NIST)
    !***REFERENCES  (NONE)
    !***ROUTINES CALLED  HPPERM, HPSORT
    !***REVISION HISTORY  (YYMMDD)
    !   890620  DATE WRITTEN
    !   901005  Included test of HPPERM.  (MAM)
    !   920511  Added error message tests.  (MAM)
    !***END PROLOGUE  HSRTQC
    !
    INTEGER N, NTEST
    PARAMETER (N=9,NTEST=4)
    !
    LOGICAL fail
    CHARACTER :: short
    CHARACTER(2) :: x(N,NTEST), xs(N,NTEST), y(N), work(N)
    INTEGER ix(N,NTEST), iy(N), kflag(NTEST), Kprint, Lun, Ipass, j, &
      i, kabs, ier, nerr, NUMXER, nn, kkflag, strbeg, strend
    !
    !     ---------
    !     TEST DATA
    !     ---------
    !
    !         X   = TEST VECTOR
    !         XS  = TEST VECTOR IN SORTED ORDER
    !         IX  = PERMUTATION VECTOR, I.E.  X(IX(J)) = XS(J)
    !
    DATA kflag(1)/2/
    DATA (x(i,1),i=1,N)/'AC', 'AZ', 'AD', 'AA', 'AB', 'ZZ', 'ZA', &
      'ZX', 'ZY'/
    DATA (ix(i,1),i=1,N)/4, 5, 1, 3, 2, 7, 8, 9, 6/
    DATA (xs(i,1),i=1,N)/'AA', 'AB', 'AC', 'AD', 'AZ', 'ZA', 'ZX', &
      'ZY', 'ZZ'/
    !
    DATA kflag(2)/ - 1/
    DATA (x(i,2),i=1,N)/'AA', 'BB', 'CC', 'DD', 'EE', 'FF', 'GG', &
      'HH', 'II'/
    DATA (ix(i,2),i=1,N)/9, 8, 7, 6, 5, 4, 3, 2, 1/
    DATA (xs(i,2),i=1,N)/'II', 'HH', 'GG', 'FF', 'EE', 'DD', 'CC', &
      'BB', 'AA'/
    !
    DATA kflag(3)/ - 2/
    DATA (x(i,3),i=1,N)/'AA', 'BB', 'CC', 'DD', 'EE', 'FF', 'GG', &
      'HH', 'II'/
    DATA (ix(i,3),i=1,N)/9, 8, 7, 6, 5, 4, 3, 2, 1/
    DATA (xs(i,3),i=1,N)/'II', 'HH', 'GG', 'FF', 'EE', 'DD', 'CC', &
      'BB', 'AA'/
    !
    DATA kflag(4)/1/
    DATA (x(i,4),i=1,N)/'AC', 'AZ', 'AD', 'AA', 'AB', 'ZZ', 'ZA', &
      'ZX', 'ZY'/
    DATA (ix(i,4),i=1,N)/4, 5, 1, 3, 2, 7, 8, 9, 6/
    DATA (xs(i,4),i=1,N)/'AA', 'AB', 'AC', 'AD', 'AZ', 'ZA', 'ZX', &
      'ZY', 'ZZ'/
    !
    !***FIRST EXECUTABLE STATEMENT  HSRTQC
    IF ( Kprint>=2 ) THEN
      WRITE (Lun,99001) '================='
      WRITE (Lun,99002) 'OUTPUT FROM HSRTQC'
      WRITE (Lun,99002) '================='
    ENDIF
    Ipass = 1
    !
    !     -------------------------------------------------------------
    !                            CHECK HPSORT
    !     -------------------------------------------------------------
    !
    DO j = 1, NTEST
      !
      !        ... SETUP PROBLEM
      !
      DO i = 1, N
        y(i) = x(i,j)
      ENDDO
      !
      !        ... CALL ROUTINE TO BE TESTED
      !
      CALL HPSORT(y,N,1,2,iy,kflag(j),work,ier)
      !
      !        ... EVALUATE RESULTS
      !
      kabs = ABS(kflag(j))
      fail = .FALSE. .OR. (ier>0)
      DO i = 1, N
        fail = fail .OR. (iy(i)/=ix(i,j)) .OR. ((kabs==1).AND.(y(i)/=x(i,j)))&
          .OR. ((kabs==2).AND.(y(i)/=xs(i,j)))
      ENDDO
      !
      !        ... PRODUCE REQUIRED OUTPUT
      !
      IF ( fail ) THEN
        Ipass = 0
        IF ( Kprint>0 ) WRITE (Lun,99001) 'HPSORT FAILED TEST ', j
      ELSE
        IF ( Kprint>=2 ) WRITE (Lun,99001) 'HPSORT PASSED TEST ', j
      ENDIF
      IF ( (fail.AND.(Kprint>=2)).OR.(Kprint>=3) ) THEN
        WRITE (Lun,99001) '-------------------------'
        WRITE (Lun,99002) 'DETAILS OF HPSORT TEST ', j
        WRITE (Lun,99002) '-------------------------'
        WRITE (Lun,99002) '1ST ARGUMENT (VECTOR TO BE SORTED)'
        WRITE (Lun,99003) '             INPUT = ', (x(i,j),i=1,N)
        WRITE (Lun,99003) '   COMPUTED OUTPUT = ', (y(i),i=1,N)
        IF ( kabs==1 ) THEN
          WRITE (Lun,99003) '    CORRECT OUTPUT = ', (x(i,j),i=1,N)
        ELSE
          WRITE (Lun,99003) '    CORRECT OUTPUT = ', (xs(i,j),i=1,N)
        ENDIF
        WRITE (Lun,99002) '2ND ARGUMENT (VECTOR LENGTH)'
        WRITE (Lun,99004) '             INPUT = ', N
        WRITE (Lun,99002) '3RD ARGUMENT (PERMUTATION VECTOR)'
        WRITE (Lun,99004) '   COMPUTED OUTPUT = ', (iy(i),i=1,N)
        WRITE (Lun,99004) '    CORRECT OUTPUT = ', (ix(i,j),i=1,N)
        WRITE (Lun,99002) '4TH ARGUMENT (TYPE OF SORT)'
        WRITE (Lun,99004) '             INPUT = ', kflag(j)
      ENDIF
      !
    ENDDO
    !
    !     ... TEST ERROR MESSAGES
    !
    IF ( Kprint<=2 ) THEN
      CALL XSETF(0)
    ELSE
      CALL XSETF(-1)
    ENDIF
    !
    nn = -1
    strbeg = 1
    strend = 2
    kkflag = 1
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL HPSORT(y,nn,strbeg,strend,iy,kkflag,work,ier)
    IF ( NUMXER(nerr)/=ier ) Ipass = 0
    !
    nn = 1
    strbeg = 1
    strend = 2
    kkflag = 0
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL HPSORT(y,nn,strbeg,strend,iy,kkflag,work,ier)
    IF ( NUMXER(nerr)/=ier ) Ipass = 0
    !
    nn = 1
    strbeg = 1
    strend = 2
    kkflag = 1
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL HPSORT(y,nn,strbeg,strend,iy,kkflag,short,ier)
    IF ( NUMXER(nerr)/=ier ) Ipass = 0
    !
    nn = 1
    strbeg = 2
    strend = 1
    kkflag = 1
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL HPSORT(y,nn,strbeg,strend,iy,kkflag,work,ier)
    IF ( NUMXER(nerr)/=ier ) Ipass = 0
    !
    nn = 1
    strbeg = -1
    strend = 2
    kkflag = 1
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL HPSORT(y,nn,strbeg,strend,iy,kkflag,work,ier)
    IF ( NUMXER(nerr)/=ier ) Ipass = 0
    !
    nn = 1
    strbeg = 1
    strend = 3
    kkflag = 1
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL HPSORT(y,nn,strbeg,strend,iy,kkflag,work,ier)
    IF ( NUMXER(nerr)/=ier ) Ipass = 0
    !
    IF ( (Kprint>=2).AND.(Ipass==1) ) THEN
      WRITE (Lun,*)
      WRITE (Lun,*) ' HPSORT PASSED ERROR MESSAGE TESTS'
    ELSEIF ( (Kprint>=1).AND.(Ipass==0) ) THEN
      WRITE (Lun,*)
      WRITE (Lun,*) ' HPSORT FAILED ERROR MESSAGE TESTS'
    ENDIF
    !
    !     -------------------------------------------------------------
    !                            CHECK HPPERM
    !     -------------------------------------------------------------
    !
    DO j = 1, NTEST
      !
      !        ... SETUP PROBLEM
      !
      kabs = ABS(kflag(j))
      DO i = 1, N
        y(i) = x(i,j)
        IF ( kabs==1 ) THEN
          iy(i) = i
        ELSE
          iy(i) = ix(i,j)
        ENDIF
      ENDDO
      !
      !        ... CALL ROUTINE TO BE TESTED
      !
      CALL HPPERM(y,N,iy,work,ier)
      !
      !        ... EVALUATE RESULTS
      !
      fail = .FALSE. .OR. (ier>0)
      DO i = 1, N
        fail = fail .OR. ((kabs==1).AND.(iy(i)/=i)) .OR.&
          ((kabs==2).AND.(iy(i)/=ix(i,j))) .OR.&
          ((kabs==1).AND.(y(i)/=x(i,j))) .OR.&
          ((kabs==2).AND.(y(i)/=xs(i,j)))
      ENDDO
      !
      !        ... PRODUCE REQUIRED OUTPUT
      !
      IF ( fail ) THEN
        Ipass = 0
        IF ( Kprint>0 ) WRITE (Lun,99001) 'HPPERM FAILED TEST ', j
      ELSE
        IF ( Kprint>=2 ) WRITE (Lun,99001) 'HPPERM PASSED TEST ', j
      ENDIF
      IF ( (fail.AND.(Kprint>=2)).OR.(Kprint>=3) ) THEN
        WRITE (Lun,99001) '------------------------'
        WRITE (Lun,99002) 'DETAILS OF HPPERM TEST', j
        WRITE (Lun,99002) '------------------------'
        WRITE (Lun,99002) '1ST ARGUMENT (VECTOR TO BE PERMUTED)'
        WRITE (Lun,99003) '             INPUT =', (x(i,j),i=1,N)
        WRITE (Lun,99003) '   COMPUTED OUTPUT =', (y(i),i=1,N)
        IF ( kabs==1 ) THEN
          WRITE (Lun,99003) '    CORRECT OUTPUT =', (x(i,j),i=1,N)
        ELSE
          WRITE (Lun,99003) '    CORRECT OUTPUT =', (xs(i,j),i=1,N)
        ENDIF
        WRITE (Lun,99002) '2ND ARGUMENT (VECTOR LENGTH)'
        WRITE (Lun,99004) '             INPUT =', N
        WRITE (Lun,99002) '3RD ARGUMENT (PERMUTATION VECTOR)'
        WRITE (Lun,99004) '             INPUT =', (iy(i),i=1,N)
        WRITE (Lun,99002) '4TH ARGUMENT (ERROR FLAG)'
        WRITE (Lun,99004) '             OUTPUT =', ier
      ENDIF
      !
    ENDDO
    !
    !     ... TEST ERROR MESSAGES
    !
    IF ( Kprint<=2 ) THEN
      CALL XSETF(0)
    ELSE
      CALL XSETF(-1)
    ENDIF
    !
    nn = -1
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL HPPERM(y,nn,iy,work,ier)
    IF ( NUMXER(nerr)/=ier ) Ipass = 0
    !
    nn = 1
    iy(1) = 5
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL HPPERM(y,nn,iy,work,ier)
    IF ( NUMXER(nerr)/=ier ) Ipass = 0
    !
    IF ( (Kprint>=2).AND.(Ipass==1) ) THEN
      WRITE (Lun,*)
      WRITE (Lun,*) ' HPPERM PASSED ERROR MESSAGE TESTS'
    ELSEIF ( (Kprint>=1).AND.(Ipass==0) ) THEN
      WRITE (Lun,*)
      WRITE (Lun,*) ' HPPERM FAILED ERROR MESSAGE TESTS'
    ENDIF
    !
    RETURN
    !
    99001 FORMAT (/1X,A,I2)
    99002 FORMAT (1X,A,I2)
    99003 FORMAT (1X,A,9(2X,A2))
    99004 FORMAT (1X,A,9I4)
  END SUBROUTINE HSRTQC
  !DECK SSRTQC
  SUBROUTINE SSRTQC(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !***BEGIN PROLOGUE  SSRTQC
    !***SUBSIDIARY
    !***PURPOSE  Quick check for SLATEC routines SSORT, SPSORT, SPPERM
    !***LIBRARY   SLATEC
    !***CATEGORY  N6A
    !***TYPE      SINGLE PRECISION (SSRTQC-S, DSRTQC-D, ISRTQC-I, HSRTQC-H)
    !***KEYWORDS  QUICK CHECK, SPPERM, SPSORT, SSORT
    !***AUTHOR  Boisvert, Ronald, (NIST)
    !***REFERENCES  (NONE)
    !***ROUTINES CALLED  SPPERM, SPSORT, SSORT
    !***REVISION HISTORY  (YYMMDD)
    !   890620  DATE WRITTEN
    !   901005  Included test of SPPERM.  (MAM)
    !   920511  Added error message tests.  (MAM)
    !***END PROLOGUE  SSRTQC
    !
    INTEGER N, NTEST
    PARAMETER (N=9,NTEST=4)
    !
    LOGICAL fail
    REAL x(N,NTEST), xs(N,NTEST), y(N), yc(N)
    INTEGER ix(N,NTEST), iy(N), kflag(NTEST), Kprint, Lun, Ipass, j, &
      i, kabs, ier, nerr, NUMXER, nn, kkflag
    !
    !     ---------
    !     TEST DATA
    !     ---------
    !
    !         X   = TEST VECTOR
    !         XS  = TEST VECTOR IN SORTED ORDER
    !         IX  = PERMUTATION VECTOR, I.E.  X(IX(J)) = XS(J)
    !
    DATA kflag(1)/2/
    DATA (x(i,1),i=1,N)/36., 54., -1., 29., 1., 80., 98., 99., 55./
    DATA (ix(i,1),i=1,N)/3, 5, 4, 1, 2, 9, 6, 7, 8/
    DATA (xs(i,1),i=1,N)/ - 1., 1., 29., 36., 54., 55., 80., 98., 99./
    !
    DATA kflag(2)/ - 1/
    DATA (x(i,2),i=1,N)/1., 2., 3., 4., 5., 6., 7., 8., 9./
    DATA (ix(i,2),i=1,N)/9, 8, 7, 6, 5, 4, 3, 2, 1/
    DATA (xs(i,2),i=1,N)/9., 8., 7., 6., 5., 4., 3., 2., 1./
    !
    DATA kflag(3)/ - 2/
    DATA (x(i,3),i=1,N)/ - 9., -8., -7., -6., -5., -4., -3., -2., -1./
    DATA (ix(i,3),i=1,N)/9, 8, 7, 6, 5, 4, 3, 2, 1/
    DATA (xs(i,3),i=1,N)/ - 1., -2., -3., -4., -5., -6., -7., -8., &
      -9./
    !
    DATA kflag(4)/1/
    DATA (x(i,4),i=1,N)/36., 54., -1., 29., 1., 80., 98., 99., 55./
    DATA (ix(i,4),i=1,N)/3, 5, 4, 1, 2, 9, 6, 7, 8/
    DATA (xs(i,4),i=1,N)/ - 1., 1., 29., 36., 54., 55., 80., 98., 99./
    !
    !***FIRST EXECUTABLE STATEMENT  SSRTQC
    IF ( Kprint>=2 ) THEN
      WRITE (Lun,99001) '================='
      WRITE (Lun,99002) 'OUTPUT FROM SSRTQC'
      WRITE (Lun,99002) '================='
    ENDIF
    Ipass = 1
    !
    !     -------------------------------------------------------------
    !                          CHECK SSORT
    !     -------------------------------------------------------------
    !
    DO j = 1, NTEST
      !
      !        ... SETUP PROBLEM
      !
      DO i = 1, N
        y(i) = x(i,j)
        yc(i) = x(i,j)
      ENDDO
      !
      !        ... CALL ROUTINE TO BE TESTED
      !
      CALL SSORT(y,yc,N,kflag(j))
      !
      !        ... EVALUATE RESULTS
      !
      kabs = ABS(kflag(j))
      fail = .FALSE.
      DO i = 1, N
        fail = fail .OR. (y(i)/=xs(i,j)) .OR. ((kabs==1).AND.(yc(i)/=x(i,j)))&
          .OR. ((kabs==2).AND.(yc(i)/=xs(i,j)))
      ENDDO
      !
      !        ... PRODUCE REQUIRED OUTPUT
      !
      IF ( fail ) THEN
        Ipass = 0
        IF ( Kprint>0 ) WRITE (Lun,99001) 'SSORT FAILED TEST ', j
      ELSE
        IF ( Kprint>=2 ) WRITE (Lun,99001) 'SSORT PASSED TEST ', j
      ENDIF
      IF ( (fail.AND.(Kprint>=2)).OR.(Kprint>=3) ) THEN
        WRITE (Lun,99001) '------------------------'
        WRITE (Lun,99002) 'DETAILS OF SSORT TEST ', j
        WRITE (Lun,99002) '------------------------'
        WRITE (Lun,99002) '1ST ARGUMENT (VECTOR TO BE SORTED)'
        WRITE (Lun,99003) '             INPUT = ', (x(i,j),i=1,N)
        WRITE (Lun,99003) '   COMPUTED OUTPUT = ', (y(i),i=1,N)
        WRITE (Lun,99003) '    CORRECT OUTPUT = ', (xs(i,j),i=1,N)
        WRITE (Lun,99002) '2ND ARGUMENT (VECTOR CARRIED ALONG)'
        WRITE (Lun,99003) '             INPUT = ', (x(i,j),i=1,N)
        WRITE (Lun,99003) '   COMPUTED OUTPUT = ', (yc(i),i=1,N)
        IF ( kabs==1 ) THEN
          WRITE (Lun,99003) '    CORRECT OUTPUT = ', (x(i,j),i=1,N)
        ELSE
          WRITE (Lun,99003) '    CORRECT OUTPUT = ', (xs(i,j),i=1,N)
        ENDIF
        WRITE (Lun,99002) '3RD ARGUMENT (VECTOR LENGTH)'
        WRITE (Lun,99004) '             INPUT = ', N
        WRITE (Lun,99002) '4TH ARGUMENT (TYPE OF SORT)'
        WRITE (Lun,99004) '             INPUT = ', kflag(j)
      ENDIF
    ENDDO
    !
    !     -------------------------------------------------------------
    !                            CHECK SPSORT
    !     -------------------------------------------------------------
    !
    DO j = 1, NTEST
      !
      !        ... SETUP PROBLEM
      !
      DO i = 1, N
        y(i) = x(i,j)
      ENDDO
      !
      !        ... CALL ROUTINE TO BE TESTED
      !
      CALL SPSORT(y,N,iy,kflag(j),ier)
      !
      !        ... EVALUATE RESULTS
      !
      kabs = ABS(kflag(j))
      fail = .FALSE. .OR. (ier>0)
      DO i = 1, N
        fail = fail .OR. (iy(i)/=ix(i,j)) .OR. ((kabs==1).AND.(y(i)/=x(i,j)))&
          .OR. ((kabs==2).AND.(y(i)/=xs(i,j)))
      ENDDO
      !
      !        ... PRODUCE REQUIRED OUTPUT
      !
      IF ( fail ) THEN
        Ipass = 0
        IF ( Kprint>0 ) WRITE (Lun,99001) 'SPSORT FAILED TEST ', j
      ELSE
        IF ( Kprint>=2 ) WRITE (Lun,99001) 'SPSORT PASSED TEST ', j
      ENDIF
      IF ( (fail.AND.(Kprint>=2)).OR.(Kprint>=3) ) THEN
        WRITE (Lun,99001) '-------------------------'
        WRITE (Lun,99002) 'DETAILS OF SPSORT TEST ', j
        WRITE (Lun,99002) '-------------------------'
        WRITE (Lun,99002) '1ST ARGUMENT (VECTOR TO BE SORTED)'
        WRITE (Lun,99003) '             INPUT = ', (x(i,j),i=1,N)
        WRITE (Lun,99003) '   COMPUTED OUTPUT = ', (y(i),i=1,N)
        IF ( kabs==1 ) THEN
          WRITE (Lun,99003) '    CORRECT OUTPUT = ', (x(i,j),i=1,N)
        ELSE
          WRITE (Lun,99003) '    CORRECT OUTPUT = ', (xs(i,j),i=1,N)
        ENDIF
        WRITE (Lun,99002) '2ND ARGUMENT (VECTOR LENGTH)'
        WRITE (Lun,99004) '             INPUT = ', N
        WRITE (Lun,99002) '3RD ARGUMENT (PERMUTATION VECTOR)'
        WRITE (Lun,99004) '   COMPUTED OUTPUT = ', (iy(i),i=1,N)
        WRITE (Lun,99004) '    CORRECT OUTPUT = ', (ix(i,j),i=1,N)
        WRITE (Lun,99002) '4TH ARGUMENT (TYPE OF SORT)'
        WRITE (Lun,99004) '             INPUT = ', kflag(j)
      ENDIF
      !
    ENDDO
    !
    !     ... TEST ERROR MESSAGES
    !
    IF ( Kprint<=2 ) THEN
      CALL XSETF(0)
    ELSE
      CALL XSETF(-1)
    ENDIF
    !
    nn = -1
    kkflag = 1
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL SPSORT(y,nn,iy,kkflag,ier)
    IF ( NUMXER(nerr)/=ier ) Ipass = 0
    !
    nn = 1
    kkflag = 0
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL SPSORT(y,nn,iy,kkflag,ier)
    IF ( NUMXER(nerr)/=ier ) Ipass = 0
    !
    IF ( (Kprint>=2).AND.(Ipass==1) ) THEN
      WRITE (Lun,*)
      WRITE (Lun,*) ' SPSORT PASSED ERROR MESSAGE TESTS'
    ELSEIF ( (Kprint>=1).AND.(Ipass==0) ) THEN
      WRITE (Lun,*)
      WRITE (Lun,*) ' SPSORT FAILED ERROR MESSAGE TESTS'
    ENDIF
    !
    !     -------------------------------------------------------------
    !                            CHECK SPPERM
    !     -------------------------------------------------------------
    !
    DO j = 1, NTEST
      !
      !        ... SETUP PROBLEM
      !
      kabs = ABS(kflag(j))
      DO i = 1, N
        y(i) = x(i,j)
        IF ( kabs==1 ) THEN
          iy(i) = i
        ELSE
          iy(i) = ix(i,j)
        ENDIF
      ENDDO
      !
      !        ... CALL ROUTINE TO BE TESTED
      !
      CALL SPPERM(y,N,iy,ier)
      !
      !        ... EVALUATE RESULTS
      !
      fail = .FALSE. .OR. (ier>0)
      DO i = 1, N
        fail = fail .OR. ((kabs==1).AND.(iy(i)/=i)) .OR. &
          ((kabs==2).AND.(iy(i)/=ix(i,j))) .OR. &
          ((kabs==1).AND.(y(i)/=x(i,j))) .OR. &
          ((kabs==2).AND.(y(i)/=xs(i,j)))
      ENDDO
      !
      !        ... PRODUCE REQUIRED OUTPUT
      !
      IF ( fail ) THEN
        Ipass = 0
        IF ( Kprint>0 ) WRITE (Lun,99001) 'SPPERM FAILED TEST ', j
      ELSE
        IF ( Kprint>=2 ) WRITE (Lun,99001) 'SPPERM PASSED TEST ', j
      ENDIF
      IF ( (fail.AND.(Kprint>=2)).OR.(Kprint>=3) ) THEN
        WRITE (Lun,99001) '------------------------'
        WRITE (Lun,99002) 'DETAILS OF SPPERM TEST', j
        WRITE (Lun,99002) '------------------------'
        WRITE (Lun,99002) '1ST ARGUMENT (VECTOR TO BE PERMUTED)'
        WRITE (Lun,99003) '             INPUT =', (x(i,j),i=1,N)
        WRITE (Lun,99003) '   COMPUTED OUTPUT =', (y(i),i=1,N)
        IF ( kabs==1 ) THEN
          WRITE (Lun,99003) '    CORRECT OUTPUT =', (x(i,j),i=1,N)
        ELSE
          WRITE (Lun,99003) '    CORRECT OUTPUT =', (xs(i,j),i=1,N)
        ENDIF
        WRITE (Lun,99002) '2ND ARGUMENT (VECTOR LENGTH)'
        WRITE (Lun,99004) '             INPUT =', N
        WRITE (Lun,99002) '3RD ARGUMENT (PERMUTATION VECTOR)'
        WRITE (Lun,99004) '             INPUT =', (iy(i),i=1,N)
        WRITE (Lun,99002) '4TH ARGUMENT (ERROR FLAG)'
        WRITE (Lun,99004) '             OUTPUT =', ier
      ENDIF
      !
    ENDDO
    !
    !     ... TEST ERROR MESSAGES
    !
    IF ( Kprint<=2 ) THEN
      CALL XSETF(0)
    ELSE
      CALL XSETF(-1)
    ENDIF
    !
    nn = -1
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL SPPERM(y,nn,iy,ier)
    IF ( NUMXER(nerr)/=ier ) Ipass = 0
    !
    nn = 1
    iy(1) = 5
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL SPPERM(y,nn,iy,ier)
    IF ( NUMXER(nerr)/=ier ) Ipass = 0
    !
    IF ( (Kprint>=2).AND.(Ipass==1) ) THEN
      WRITE (Lun,*)
      WRITE (Lun,*) ' SPPERM PASSED ERROR MESSAGE TESTS'
    ELSEIF ( (Kprint>=1).AND.(Ipass==0) ) THEN
      WRITE (Lun,*)
      WRITE (Lun,*) ' SPPERM FAILED ERROR MESSAGE TESTS'
    ENDIF
    !
    RETURN
    !
    99001 FORMAT (/1X,A,I2)
    99002 FORMAT (1X,A,I2)
    99003 FORMAT (1X,A,9F4.0)
    99004 FORMAT (1X,A,9I4)
  END SUBROUTINE SSRTQC
  !DECK DSRTQC
  SUBROUTINE DSRTQC(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !***BEGIN PROLOGUE  DSRTQC
    !***SUBSIDIARY
    !***PURPOSE  Quick check for SLATEC routines DSORT, DPSORT, DPPERM
    !***LIBRARY   SLATEC
    !***CATEGORY  N6A
    !***TYPE      DOUBLE PRECISION (SSRTQC-S, DSRTQC-D, ISRTQC-I, HSRTQC-H)
    !***KEYWORDS  DPPERM, DPSORT, DSORT, QUICK CHECK
    !***AUTHOR  Boisvert, Ronald, (NIST)
    !***REFERENCES  (NONE)
    !***ROUTINES CALLED  DPPERM, DPSORT, DSORT
    !***REVISION HISTORY  (YYMMDD)
    !   890620  DATE WRITTEN
    !   901005  Included test of DPPERM.  (MAM)
    !   920511  Added error message tests.  (MAM)
    !***END PROLOGUE  DSRTQC
    !
    INTEGER N, NTEST
    PARAMETER (N=9,NTEST=4)
    !
    LOGICAL fail
    REAL(8) :: x(N,NTEST), xs(N,NTEST), y(N), yc(N)
    INTEGER ix(N,NTEST), iy(N), kflag(NTEST), Kprint, Lun, Ipass, j, &
      i, kabs, ier, nerr, NUMXER, nn, kkflag
    !
    !     ---------
    !     TEST DATA
    !     ---------
    !
    !         X   = TEST VECTOR
    !         XS  = TEST VECTOR IN SORTED ORDER
    !         IX  = PERMUTATION VECTOR, I.E.  X(IX(J)) = XS(J)
    !
    DATA kflag(1)/2/
    DATA (x(i,1),i=1,N)/36D0, 54D0, -1D0, 29D0, 1D0, 80D0, 98D0, 99D0, &
      55D0/
    DATA (ix(i,1),i=1,N)/3, 5, 4, 1, 2, 9, 6, 7, 8/
    DATA (xs(i,1),i=1,N)/ - 1D0, 1D0, 29D0, 36D0, 54D0, 55D0, 80D0, &
      98D0, 99D0/
    !
    DATA kflag(2)/ - 1/
    DATA (x(i,2),i=1,N)/1D0, 2D0, 3D0, 4D0, 5D0, 6D0, 7D0, 8D0, 9D0/
    DATA (ix(i,2),i=1,N)/9, 8, 7, 6, 5, 4, 3, 2, 1/
    DATA (xs(i,2),i=1,N)/9D0, 8D0, 7D0, 6D0, 5D0, 4D0, 3D0, 2D0, 1D0/
    !
    DATA kflag(3)/ - 2/
    DATA (x(i,3),i=1,N)/ - 9D0, -8D0, -7D0, -6D0, -5D0, -4D0, -3D0, &
      -2D0, -1D0/
    DATA (ix(i,3),i=1,N)/9, 8, 7, 6, 5, 4, 3, 2, 1/
    DATA (xs(i,3),i=1,N)/ - 1D0, -2D0, -3D0, -4D0, -5D0, -6D0, -7D0, &
      -8D0, -9D0/
    !
    DATA kflag(4)/1/
    DATA (x(i,4),i=1,N)/36D0, 54D0, -1D0, 29D0, 1D0, 80D0, 98D0, 99D0, &
      55D0/
    DATA (ix(i,4),i=1,N)/3, 5, 4, 1, 2, 9, 6, 7, 8/
    DATA (xs(i,4),i=1,N)/ - 1D0, 1D0, 29D0, 36D0, 54D0, 55D0, 80D0, &
      98D0, 99D0/
    !
    !***FIRST EXECUTABLE STATEMENT  DSRTQC
    IF ( Kprint>=2 ) THEN
      WRITE (Lun,99001) '================='
      WRITE (Lun,99002) 'OUTPUT FROM DSRTQC'
      WRITE (Lun,99002) '================='
    ENDIF
    Ipass = 1
    !
    !     -------------------------------------------------------------
    !                          CHECK DSORT
    !     -------------------------------------------------------------
    !
    DO j = 1, NTEST
      !
      !        ... SETUP PROBLEM
      !
      DO i = 1, N
        y(i) = x(i,j)
        yc(i) = x(i,j)
      ENDDO
      !
      !        ... CALL ROUTINE TO BE TESTED
      !
      CALL DSORT(y,yc,N,kflag(j))
      !
      !        ... EVALUATE RESULTS
      !
      kabs = ABS(kflag(j))
      fail = .FALSE.
      DO i = 1, N
        fail = fail .OR. (y(i)/=xs(i,j)) .OR. ((kabs==1).AND.(yc(i)/=x(i,j)))&
          .OR. ((kabs==2).AND.(yc(i)/=xs(i,j)))
      ENDDO
      !
      !        ... PRODUCE REQUIRED OUTPUT
      !
      IF ( fail ) THEN
        Ipass = 0
        IF ( Kprint>0 ) WRITE (Lun,99001) 'DSORT FAILED TEST ', j
      ELSE
        IF ( Kprint>=2 ) WRITE (Lun,99001) 'DSORT PASSED TEST ', j
      ENDIF
      IF ( (fail.AND.(Kprint>=2)).OR.(Kprint>=3) ) THEN
        WRITE (Lun,99001) '------------------------'
        WRITE (Lun,99002) 'DETAILS OF DSORT TEST ', j
        WRITE (Lun,99002) '------------------------'
        WRITE (Lun,99002) '1ST ARGUMENT (VECTOR TO BE SORTED)'
        WRITE (Lun,99003) '             INPUT = ', (x(i,j),i=1,N)
        WRITE (Lun,99003) '   COMPUTED OUTPUT = ', (y(i),i=1,N)
        WRITE (Lun,99003) '    CORRECT OUTPUT = ', (xs(i,j),i=1,N)
        WRITE (Lun,99002) '2ND ARGUMENT (VECTOR CARRIED ALONG)'
        WRITE (Lun,99003) '             INPUT = ', (x(i,j),i=1,N)
        WRITE (Lun,99003) '   COMPUTED OUTPUT = ', (yc(i),i=1,N)
        IF ( kabs==1 ) THEN
          WRITE (Lun,99003) '    CORRECT OUTPUT = ', (x(i,j),i=1,N)
        ELSE
          WRITE (Lun,99003) '    CORRECT OUTPUT = ', (xs(i,j),i=1,N)
        ENDIF
        WRITE (Lun,99002) '3RD ARGUMENT (VECTOR LENGTH)'
        WRITE (Lun,99004) '             INPUT = ', N
        WRITE (Lun,99002) '4TH ARGUMENT (TYPE OF SORT)'
        WRITE (Lun,99004) '             INPUT = ', kflag(j)
      ENDIF
    ENDDO
    !
    !     -------------------------------------------------------------
    !                            CHECK DPSORT
    !     -------------------------------------------------------------
    !
    DO j = 1, NTEST
      !
      !        ... SETUP PROBLEM
      !
      DO i = 1, N
        y(i) = x(i,j)
      ENDDO
      !
      !        ... CALL ROUTINE TO BE TESTED
      !
      CALL DPSORT(y,N,iy,kflag(j),ier)
      !
      !        ... EVALUATE RESULTS
      !
      kabs = ABS(kflag(j))
      fail = .FALSE. .OR. (ier>0)
      DO i = 1, N
        fail = fail .OR. (iy(i)/=ix(i,j)) .OR. ((kabs==1).AND.(y(i)/=x(i,j)))&
          .OR. ((kabs==2).AND.(y(i)/=xs(i,j)))
      ENDDO
      !
      !        ... PRODUCE REQUIRED OUTPUT
      !
      IF ( fail ) THEN
        Ipass = 0
        IF ( Kprint>0 ) WRITE (Lun,99001) 'DPSORT FAILED TEST ', j
      ELSE
        IF ( Kprint>=2 ) WRITE (Lun,99001) 'DPSORT PASSED TEST ', j
      ENDIF
      IF ( (fail.AND.(Kprint>=2)).OR.(Kprint>=3) ) THEN
        WRITE (Lun,99001) '-------------------------'
        WRITE (Lun,99002) 'DETAILS OF DPSORT TEST ', j
        WRITE (Lun,99002) '-------------------------'
        WRITE (Lun,99002) '1ST ARGUMENT (VECTOR TO BE SORTED)'
        WRITE (Lun,99003) '             INPUT = ', (x(i,j),i=1,N)
        WRITE (Lun,99003) '   COMPUTED OUTPUT = ', (y(i),i=1,N)
        IF ( kabs==1 ) THEN
          WRITE (Lun,99003) '    CORRECT OUTPUT = ', (x(i,j),i=1,N)
        ELSE
          WRITE (Lun,99003) '    CORRECT OUTPUT = ', (xs(i,j),i=1,N)
        ENDIF
        WRITE (Lun,99002) '2ND ARGUMENT (VECTOR LENGTH)'
        WRITE (Lun,99004) '             INPUT = ', N
        WRITE (Lun,99002) '3RD ARGUMENT (PERMUTATION VECTOR)'
        WRITE (Lun,99004) '   COMPUTED OUTPUT = ', (iy(i),i=1,N)
        WRITE (Lun,99004) '    CORRECT OUTPUT = ', (ix(i,j),i=1,N)
        WRITE (Lun,99002) '4TH ARGUMENT (TYPE OF SORT)'
        WRITE (Lun,99004) '             INPUT = ', kflag(j)
      ENDIF
      !
    ENDDO
    !
    !     ... TEST ERROR MESSAGES
    !
    IF ( Kprint<=2 ) THEN
      CALL XSETF(0)
    ELSE
      CALL XSETF(-1)
    ENDIF
    !
    nn = -1
    kkflag = 1
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL DPSORT(y,nn,iy,kkflag,ier)
    IF ( NUMXER(nerr)/=ier ) Ipass = 0
    !
    nn = 1
    kkflag = 0
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL DPSORT(y,nn,iy,kkflag,ier)
    IF ( NUMXER(nerr)/=ier ) Ipass = 0
    !
    IF ( (Kprint>=2).AND.(Ipass==1) ) THEN
      WRITE (Lun,*)
      WRITE (Lun,*) ' DPSORT PASSED ERROR MESSAGE TESTS'
    ELSEIF ( (Kprint>=1).AND.(Ipass==0) ) THEN
      WRITE (Lun,*)
      WRITE (Lun,*) ' DPSORT FAILED ERROR MESSAGE TESTS'
    ENDIF
    !
    !     -------------------------------------------------------------
    !                            CHECK DPPERM
    !     -------------------------------------------------------------
    !
    DO j = 1, NTEST
      !
      !        ... SETUP PROBLEM
      !
      kabs = ABS(kflag(j))
      DO i = 1, N
        y(i) = x(i,j)
        IF ( kabs==1 ) THEN
          iy(i) = i
        ELSE
          iy(i) = ix(i,j)
        ENDIF
      ENDDO
      !
      !        ... CALL ROUTINE TO BE TESTED
      !
      CALL DPPERM(y,N,iy,ier)
      !
      !        ... EVALUATE RESULTS
      !
      fail = .FALSE. .OR. (ier>0)
      DO i = 1, N
        fail = fail .OR. ((kabs==1).AND.(iy(i)/=i)) .OR. &
          ((kabs==2).AND.(iy(i)/=ix(i,j))) .OR. &
          ((kabs==1).AND.(y(i)/=x(i,j))) .OR. &
          ((kabs==2).AND.(y(i)/=xs(i,j)))
      ENDDO
      !
      !        ... PRODUCE REQUIRED OUTPUT
      !
      IF ( fail ) THEN
        Ipass = 0
        IF ( Kprint>0 ) WRITE (Lun,99001) 'DPPERM FAILED TEST ', j
      ELSE
        IF ( Kprint>=2 ) WRITE (Lun,99001) 'DPPERM PASSED TEST ', j
      ENDIF
      IF ( (fail.AND.(Kprint>=2)).OR.(Kprint>=3) ) THEN
        WRITE (Lun,99001) '------------------------'
        WRITE (Lun,99002) 'DETAILS OF DPPERM TEST', j
        WRITE (Lun,99002) '------------------------'
        WRITE (Lun,99002) '1ST ARGUMENT (VECTOR TO BE PERMUTED)'
        WRITE (Lun,99003) '             INPUT =', (x(i,j),i=1,N)
        WRITE (Lun,99003) '   COMPUTED OUTPUT =', (y(i),i=1,N)
        IF ( kabs==1 ) THEN
          WRITE (Lun,99003) '    CORRECT OUTPUT =', (x(i,j),i=1,N)
        ELSE
          WRITE (Lun,99003) '    CORRECT OUTPUT =', (xs(i,j),i=1,N)
        ENDIF
        WRITE (Lun,99002) '2ND ARGUMENT (VECTOR LENGTH)'
        WRITE (Lun,99004) '             INPUT =', N
        WRITE (Lun,99002) '3RD ARGUMENT (PERMUTATION VECTOR)'
        WRITE (Lun,99004) '             INPUT =', (iy(i),i=1,N)
        WRITE (Lun,99002) '4TH ARGUMENT (ERROR FLAG)'
        WRITE (Lun,99004) '             OUTPUT =', ier
      ENDIF
      !
    ENDDO
    !
    !     ... TEST ERROR MESSAGES
    !
    IF ( Kprint<=2 ) THEN
      CALL XSETF(0)
    ELSE
      CALL XSETF(-1)
    ENDIF
    !
    nn = -1
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL DPPERM(y,nn,iy,ier)
    IF ( NUMXER(nerr)/=ier ) Ipass = 0
    !
    nn = 1
    iy(1) = 5
    IF ( Kprint>=3 ) WRITE (Lun,*)
    CALL XERCLR
    CALL DPPERM(y,nn,iy,ier)
    IF ( NUMXER(nerr)/=ier ) Ipass = 0
    !
    IF ( (Kprint>=2).AND.(Ipass==1) ) THEN
      WRITE (Lun,*)
      WRITE (Lun,*) ' DPPERM PASSED ERROR MESSAGE TESTS'
    ELSEIF ( (Kprint>=1).AND.(Ipass==0) ) THEN
      WRITE (Lun,*)
      WRITE (Lun,*) ' DPPERM FAILED ERROR MESSAGE TESTS'
    ENDIF
    !
    RETURN
    !
    99001 FORMAT (/1X,A,I2)
    99002 FORMAT (1X,A,I2)
    99003 FORMAT (1X,A,9F4.0)
    99004 FORMAT (1X,A,9I4)
  END SUBROUTINE DSRTQC
END MODULE TEST54_MOD
!DECK TEST54
PROGRAM TEST54
  USE TEST54_MOD
  IMPLICIT NONE
  !***BEGIN PROLOGUE  TEST54
  !***PURPOSE  Driver for testing SLATEC subprograms
  !***LIBRARY   SLATEC
  !***CATEGORY  N6A
  !***TYPE      ALL (TEST54-A)
  !***KEYWORDS  QUICK CHECK DRIVER
  !***AUTHOR  SLATEC Common Mathematical Library Committee
  !***DESCRIPTION
  !
  ! *Usage:
  !     One input data record is required
  !         READ (LIN, '(I1)') KPRINT
  !
  ! *Arguments:
  !     KPRINT = 0  Quick checks - No printing.
  !                 Driver       - Short pass or fail message printed.
  !              1  Quick checks - No message printed for passed tests,
  !                                short message printed for failed tests.
  !                 Driver       - Short pass or fail message printed.
  !              2  Quick checks - Print short message for passed tests,
  !                                fuller information for failed tests.
  !                 Driver       - Pass or fail message printed.
  !              3  Quick checks - Print complete quick check results.
  !                 Driver       - Pass or fail message printed.
  !
  ! *Description:
  !     Driver for testing SLATEC subprograms
  !        ISORT   SSORT   DSORT
  !        IPSORT  SPSORT  DPSORT  HPSORT
  !        IPPERM  SPPERM  DPPERM  HPPERM
  !
  !***REFERENCES  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***ROUTINES CALLED  DSRTQC, HSRTQC, I1MACH, ISRTQC, SSRTQC, XERMAX,
  !                    XSETF, XSETUN
  !***REVISION HISTORY  (YYMMDD)
  !   890620  DATE WRITTEN
  !***END PROLOGUE  TEST54
  INTEGER lun, I1MACH, lin, nfail, kprint, ipass
  !***FIRST EXECUTABLE STATEMENT  TEST54
  lun = I1MACH(2)
  lin = I1MACH(1)
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  CALL GET_ARGUMENT(kprint)
  CALL XERMAX(1000)
  CALL XSETUN(lun)
  IF ( kprint<=1 ) THEN
    CALL XSETF(0)
  ELSE
    CALL XSETF(1)
  ENDIF
  !
  !     Test ISORT, IPSORT and IPPERM
  !
  CALL ISRTQC(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test SSORT, SPSORT and SPPERM
  !
  CALL SSRTQC(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test DSORT, DPSORT and DPPERM
  !
  CALL DSRTQC(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test HPSORT and HPPERM
  !
  CALL HSRTQC(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST54 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST54 *************')
  ENDIF
  STOP
END PROGRAM TEST54
