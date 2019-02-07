!*==HSRTQC.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK HSRTQC
SUBROUTINE HSRTQC(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--HSRTQC5
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
  INTEGER N , NTEST
  PARAMETER (N=9,NTEST=4)
  !
  LOGICAL fail
  CHARACTER :: short
  CHARACTER(2) :: x(N,NTEST) , xs(N,NTEST) , y(N) , work(N)
  INTEGER ix(N,NTEST) , iy(N) , kflag(NTEST) , Kprint , Lun , Ipass , j ,&
    i , kabs , ier , nerr , NUMXER , nn , kkflag , strbeg , strend
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
  DATA (x(i,1),i=1,N)/'AC' , 'AZ' , 'AD' , 'AA' , 'AB' , 'ZZ' , 'ZA' ,&
    'ZX' , 'ZY'/
  DATA (ix(i,1),i=1,N)/4 , 5 , 1 , 3 , 2 , 7 , 8 , 9 , 6/
  DATA (xs(i,1),i=1,N)/'AA' , 'AB' , 'AC' , 'AD' , 'AZ' , 'ZA' , 'ZX' ,&
    'ZY' , 'ZZ'/
  !
  DATA kflag(2)/ - 1/
  DATA (x(i,2),i=1,N)/'AA' , 'BB' , 'CC' , 'DD' , 'EE' , 'FF' , 'GG' ,&
    'HH' , 'II'/
  DATA (ix(i,2),i=1,N)/9 , 8 , 7 , 6 , 5 , 4 , 3 , 2 , 1/
  DATA (xs(i,2),i=1,N)/'II' , 'HH' , 'GG' , 'FF' , 'EE' , 'DD' , 'CC' ,&
    'BB' , 'AA'/
  !
  DATA kflag(3)/ - 2/
  DATA (x(i,3),i=1,N)/'AA' , 'BB' , 'CC' , 'DD' , 'EE' , 'FF' , 'GG' ,&
    'HH' , 'II'/
  DATA (ix(i,3),i=1,N)/9 , 8 , 7 , 6 , 5 , 4 , 3 , 2 , 1/
  DATA (xs(i,3),i=1,N)/'II' , 'HH' , 'GG' , 'FF' , 'EE' , 'DD' , 'CC' ,&
    'BB' , 'AA'/
  !
  DATA kflag(4)/1/
  DATA (x(i,4),i=1,N)/'AC' , 'AZ' , 'AD' , 'AA' , 'AB' , 'ZZ' , 'ZA' ,&
    'ZX' , 'ZY'/
  DATA (ix(i,4),i=1,N)/4 , 5 , 1 , 3 , 2 , 7 , 8 , 9 , 6/
  DATA (xs(i,4),i=1,N)/'AA' , 'AB' , 'AC' , 'AD' , 'AZ' , 'ZA' , 'ZX' ,&
    'ZY' , 'ZZ'/
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
  DO j = 1 , NTEST
    !
    !        ... SETUP PROBLEM
    !
    DO i = 1 , N
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
    DO i = 1 , N
      fail = fail .OR. (iy(i)/=ix(i,j)) .OR. ((kabs==1).AND.(y(i)/=x(i,j)))&
        .OR. ((kabs==2).AND.(y(i)/=xs(i,j)))
    ENDDO
    !
    !        ... PRODUCE REQUIRED OUTPUT
    !
    IF ( fail ) THEN
      Ipass = 0
      IF ( Kprint>0 ) WRITE (Lun,99001) 'HPSORT FAILED TEST ' , j
    ELSE
      IF ( Kprint>=2 ) WRITE (Lun,99001) 'HPSORT PASSED TEST ' , j
    ENDIF
    IF ( (fail.AND.(Kprint>=2)).OR.(Kprint>=3) ) THEN
      WRITE (Lun,99001) '-------------------------'
      WRITE (Lun,99002) 'DETAILS OF HPSORT TEST ' , j
      WRITE (Lun,99002) '-------------------------'
      WRITE (Lun,99002) '1ST ARGUMENT (VECTOR TO BE SORTED)'
      WRITE (Lun,99003) '             INPUT = ' , (x(i,j),i=1,N)
      WRITE (Lun,99003) '   COMPUTED OUTPUT = ' , (y(i),i=1,N)
      IF ( kabs==1 ) THEN
        WRITE (Lun,99003) '    CORRECT OUTPUT = ' , (x(i,j),i=1,N)
      ELSE
        WRITE (Lun,99003) '    CORRECT OUTPUT = ' , (xs(i,j),i=1,N)
      ENDIF
      WRITE (Lun,99002) '2ND ARGUMENT (VECTOR LENGTH)'
      WRITE (Lun,99004) '             INPUT = ' , N
      WRITE (Lun,99002) '3RD ARGUMENT (PERMUTATION VECTOR)'
      WRITE (Lun,99004) '   COMPUTED OUTPUT = ' , (iy(i),i=1,N)
      WRITE (Lun,99004) '    CORRECT OUTPUT = ' , (ix(i,j),i=1,N)
      WRITE (Lun,99002) '4TH ARGUMENT (TYPE OF SORT)'
      WRITE (Lun,99004) '             INPUT = ' , kflag(j)
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
  DO j = 1 , NTEST
    !
    !        ... SETUP PROBLEM
    !
    kabs = ABS(kflag(j))
    DO i = 1 , N
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
    DO i = 1 , N
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
      IF ( Kprint>0 ) WRITE (Lun,99001) 'HPPERM FAILED TEST ' , j
    ELSE
      IF ( Kprint>=2 ) WRITE (Lun,99001) 'HPPERM PASSED TEST ' , j
    ENDIF
    IF ( (fail.AND.(Kprint>=2)).OR.(Kprint>=3) ) THEN
      WRITE (Lun,99001) '------------------------'
      WRITE (Lun,99002) 'DETAILS OF HPPERM TEST' , j
      WRITE (Lun,99002) '------------------------'
      WRITE (Lun,99002) '1ST ARGUMENT (VECTOR TO BE PERMUTED)'
      WRITE (Lun,99003) '             INPUT =' , (x(i,j),i=1,N)
      WRITE (Lun,99003) '   COMPUTED OUTPUT =' , (y(i),i=1,N)
      IF ( kabs==1 ) THEN
        WRITE (Lun,99003) '    CORRECT OUTPUT =' , (x(i,j),i=1,N)
      ELSE
        WRITE (Lun,99003) '    CORRECT OUTPUT =' , (xs(i,j),i=1,N)
      ENDIF
      WRITE (Lun,99002) '2ND ARGUMENT (VECTOR LENGTH)'
      WRITE (Lun,99004) '             INPUT =' , N
      WRITE (Lun,99002) '3RD ARGUMENT (PERMUTATION VECTOR)'
      WRITE (Lun,99004) '             INPUT =' , (iy(i),i=1,N)
      WRITE (Lun,99002) '4TH ARGUMENT (ERROR FLAG)'
      WRITE (Lun,99004) '             OUTPUT =' , ier
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
