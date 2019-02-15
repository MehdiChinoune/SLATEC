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
