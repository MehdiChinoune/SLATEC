!*==DPCHQ4.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DPCHQ4
SUBROUTINE DPCHQ4(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--DPCHQ45
  !***BEGIN PROLOGUE  DPCHQ4
  !***PURPOSE  Test the PCHIP monotonicity checker DPCHCM.
  !***LIBRARY   SLATEC (PCHIP)
  !***TYPE      DOUBLE PRECISION (PCHQK4-S, DPCHQ4-D)
  !***KEYWORDS  PCHIP MONOTONICITY CHECKER QUICK CHECK
  !***AUTHOR  Fritsch, F. N., (LLNL)
  !***DESCRIPTION
  !
  !             DPCHIP QUICK CHECK NUMBER 4
  !
  !     TESTS THE MONOTONICITY CHECKER:  DPCHCM.
  ! *Usage:
  !
  !        INTEGER  LUN, KPRINT, IPASS
  !
  !        CALL DPCHQ4 (LUN, KPRINT, IPASS)
  !
  ! *Arguments:
  !
  !     LUN   :IN  is the unit number to which output is to be written.
  !
  !     KPRINT:IN  controls the amount of output, as specified in the
  !                SLATEC Guidelines.
  !
  !     IPASS:OUT  will contain a pass/fail flag.  IPASS=1 is good.
  !                IPASS=0 indicates one or more tests failed.
  !
  ! *Description:
  !
  !   This routine tests a constructed data set with three different
  !   INCFD settings and compares with the expected results.  It then
  !   runs a special test to check for bug in overall monotonicity found
  !   in DPCHMC.  Finally, it reverses the data and repeats all tests.
  !
  !***ROUTINES CALLED  DPCHCM
  !***REVISION HISTORY  (YYMMDD)
  !   890208  DATE WRITTEN
  !   890306  Changed LOUT to LUN and added it to call list.  (FNF)
  !   890316  Removed DATA statements to suit new quick check standards.
  !   890410  Changed PCHMC to PCHCM.
  !   890410  Added a SLATEC 4.0 format prologue.
  !   900314  Changed name from PCHQK3 to PCHQK4 and improved some output
  !           formats.
  !   900315  Revised prologue and improved some output formats.  (FNF)
  !   900320  Converted to double precision.
  !   900321  Removed IFAIL from call sequence for SLATEC standards and
  !           made miscellaneous cosmetic changes.  (FNF)
  !   900322  Added declarations so all variables are declared.  (FNF)
  !   910708  Minor modifications in use of KPRINT.  (WRB)
  !   930317  Improved output formats.  (FNF)
  !***END PROLOGUE  DPCHQ4
  !
  !*Internal Notes:
  !
  !     Data set-up is done via assignment statements to avoid modifying
  !     DATA-loaded arrays, as required by the 1989 SLATEC Guidelines.
  !     Run with KPRINT=3 to display the data.
  !**End
  !
  !  Declare arguments.
  !
  INTEGER Lun , Kprint , Ipass
  !
  !  DECLARE VARIABLES.
  !
  INTEGER MAXN , MAXN2 , MAXN3 , NB
  PARAMETER (MAXN=16,MAXN2=8,MAXN3=6,NB=7)
  INTEGER i , ierr , ifail , incfd , ismex1(MAXN) , ismex2(MAXN2) , &
    ismex3(MAXN3) , ismexb(NB) , ismon(MAXN) , k , n , ns(3)
  DOUBLE PRECISION d(MAXN) , db(NB) , f(MAXN) , fb(NB) , x(MAXN)
  LOGICAL skip
  !
  !  DEFINE EXPECTED RESULTS.
  !
  DATA ismex1/1 , 1 , -1 , 1 , 1 , -1 , 1 , 1 , -1 , 1 , 1 , -1 , 1 , 1 , &
    -1 , 2/
  DATA ismex2/1 , 2 , 2 , 1 , 2 , 2 , 1 , 2/
  DATA ismex3/1 , 1 , 1 , 1 , 1 , 1/
  DATA ismexb/1 , 3 , 1 , -1 , -3 , -1 , 2/
  !
  !  DEFINE TEST DATA.
  !
  DATA ns/16 , 8 , 6/
  !
  !***FIRST EXECUTABLE STATEMENT  DPCHQ4
  IF ( Kprint>=3 ) WRITE (Lun,99001)
  !
  !  FORMATS.
  !
  99001 FORMAT ('1'//10X,'TEST DPCHIP MONOTONICITY CHECKER')
  IF ( Kprint>=2 ) WRITE (Lun,99002)
  99002 FORMAT (//10X,'DPCHQ4 RESULTS'/10X,'--------------')
  !
  !       Define X, F, D.
  DO i = 1 , MAXN
    x(i) = i
    d(i) = 0.D0
  ENDDO
  DO i = 2 , MAXN , 3
    d(i) = 2.D0
  ENDDO
  DO i = 1 , 3
    f(i) = x(i)
    f(i+3) = f(i) + 1.D0
    f(i+6) = f(i+3) + 1.D0
    f(i+9) = f(i+6) + 1.D0
    f(i+12) = f(i+9) + 1.D0
  ENDDO
  f(16) = 6.D0
  !       Define FB, DB.
  fb(1) = 0.D0
  fb(2) = 2.D0
  fb(3) = 3.D0
  fb(4) = 5.D0
  db(1) = 1.D0
  db(2) = 3.D0
  db(3) = 3.D0
  db(4) = 0.D0
  DO i = 1 , 3
    fb(NB-i+1) = fb(i)
    db(NB-i+1) = -db(i)
  ENDDO
  !
  !  INITIALIZE.
  !
  ifail = 0
  !
  IF ( Kprint>=3 ) THEN
    WRITE (Lun,99003)
    99003   FORMAT (//5X,'DATA:'//9X,'I',4X,'X',5X,'F',5X,'D',5X,'FB',4X,'DB')
    DO i = 1 , NB
      WRITE (Lun,99010) i , x(i) , f(i) , d(i) , fb(i) , db(i)
    ENDDO
    DO i = NB + 1 , MAXN
      WRITE (Lun,99010) i , x(i) , f(i) , d(i)
    ENDDO
  ENDIF
  !
  !  TRANSFER POINT FOR SECOND SET OF TESTS.
  !
  !
  !  Loop over a series of values of INCFD.
  !
  100  DO incfd = 1 , 3
  n = ns(incfd)
  skip = .FALSE.
  !        -------------------------------------------------
  CALL DPCHCM(n,x,f,d,incfd,skip,ismon,ierr)
  !        -------------------------------------------------
  IF ( Kprint>=3 ) WRITE (Lun,99004) incfd , ierr , (ismon(i),i=1,n)
  99004   FORMAT (/4X,'INCFD =',I2,':  IERR =',I3/15X,'ISMON =',16I3)
  IF ( ierr/=0 ) THEN
    ifail = ifail + 1
    IF ( Kprint>=3 ) WRITE (Lun,99011)
  ELSE
    DO i = 1 , n
      IF ( incfd==1 ) THEN
        IF ( ismon(i)/=ismex1(i) ) THEN
          ifail = ifail + 1
          IF ( Kprint>=3 ) WRITE (Lun,99012) (ismex1(k),k=1,n)
          EXIT
        ENDIF
      ELSEIF ( incfd==2 ) THEN
        IF ( ismon(i)/=ismex2(i) ) THEN
          ifail = ifail + 1
          IF ( Kprint>=3 ) WRITE (Lun,99012) (ismex2(k),k=1,n)
          EXIT
        ENDIF
      ELSEIF ( ismon(i)/=ismex3(i) ) THEN
        ifail = ifail + 1
        IF ( Kprint>=3 ) WRITE (Lun,99012) (ismex3(k),k=1,n)
        EXIT
      ENDIF
    ENDDO
  ENDIF
ENDDO
!
!  Test for -1,3,1 bug.
!
skip = .FALSE.
!     ------------------------------------------------
CALL DPCHCM(NB,x,fb,db,1,skip,ismon,ierr)
!     ------------------------------------------------
IF ( Kprint>=3 ) WRITE (Lun,99005) ierr , (ismon(i),i=1,NB)
99005 FORMAT (/4X,' Bug test:  IERR =',I3/15X,'ISMON =',7I3)
IF ( ierr/=0 ) THEN
  ifail = ifail + 1
  IF ( Kprint>=3 ) WRITE (Lun,99011)
ELSE
  DO i = 1 , NB
    IF ( ismon(i)/=ismexb(i) ) THEN
      ifail = ifail + 1
      IF ( Kprint>=3 ) WRITE (Lun,99012) (ismexb(k),k=1,NB)
      EXIT
    ENDIF
  ENDDO
ENDIF
!
IF ( f(1)<0. ) THEN
  !
  !  PRINT SUMMARY AND TERMINATE.
  !
  IF ( (Kprint>=2).AND.(ifail/=0) ) WRITE (Lun,99006) ifail
  99006   FORMAT (/' *** TROUBLE ***',I5,' MONOTONICITY TESTS FAILED.')
  !
  IF ( ifail==0 ) THEN
    Ipass = 1
    IF ( Kprint>=2 ) WRITE (Lun,99007)
    99007     FORMAT (/' ------------ DPCHIP PASSED  ALL MONOTONICITY TESTS',&
      ' ------------')
  ELSE
    Ipass = 0
    IF ( Kprint>=1 ) WRITE (Lun,99008)
    99008     FORMAT (/' ************ DPCHIP FAILED SOME MONOTONICITY TESTS',&
      ' ************')
  ENDIF
  !
  RETURN
ELSE
  !
  !  Change sign and do again.
  !
  IF ( Kprint>=3 ) WRITE (Lun,99009)
  99009   FORMAT (/4X,'Changing sign of data.....')
  DO i = 1 , MAXN
    f(i) = -f(i)
    d(i) = -d(i)
    IF ( ismex1(i)/=2 ) ismex1(i) = -ismex1(i)
  ENDDO
  DO i = 1 , MAXN2
    IF ( ismex2(i)/=2 ) ismex2(i) = -ismex2(i)
  ENDDO
  DO i = 1 , MAXN3
    IF ( ismex3(i)/=2 ) ismex3(i) = -ismex3(i)
  ENDDO
  DO i = 1 , NB
    fb(i) = -fb(i)
    db(i) = -db(i)
    IF ( ismexb(i)/=2 ) ismexb(i) = -ismexb(i)
  ENDDO
  GOTO 100
ENDIF
99010 FORMAT (5X,I5,5F6.1)
99011 FORMAT (' *** Failed -- bad IERR value.')
99012 FORMAT (' *** Failed -- expect:',16I3)
!------------- LAST LINE OF DPCHQ4 FOLLOWS -----------------------------
END SUBROUTINE DPCHQ4
