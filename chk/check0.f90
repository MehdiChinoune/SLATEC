!*==CHECK0.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CHECK0
SUBROUTINE CHECK0(Sfac,Dfac,Kprint)
  IMPLICIT NONE
  !*--CHECK05
  !*** Start of declarations inserted by SPAG
  INTEGER i , ICAse , INCx , INCy , jump , k , Kprint , MODe , N , NPRint
  REAL sa , sb , sc , Sfac , ss , zero
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  CHECK0
  !***PURPOSE  (UNKNOWN)
  !***LIBRARY   SLATEC
  !***AUTHOR  Lawson, C. L., (JPL)
  !           Hanson, R. J., (SNLA)
  !           Wisniewski, J. A., (SNLA)
  !***DESCRIPTION
  !
  !     THIS SUBROUTINE TESTS SUBPROGRAMS 12-13 AND 16-17.
  !     THESE SUBPROGRAMS HAVE NO ARRAY ARGUMENTS.
  !
  !     C. L. LAWSON, JPL, 1975 MAR 07, MAY 28
  !     R. J. HANSON, J. A. WISNIEWSKI, SANDIA LABS, APRIL 25,1977.
  !
  !***ROUTINES CALLED  DROTG, DROTMG, DTEST, SROTG, SROTMG, STEST
  !***COMMON BLOCKS    COMBLA
  !***REVISION HISTORY  (YYMMDD)
  !   750307  DATE WRITTEN
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   890911  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  CHECK0
  COMMON /COMBLA/ NPRint , ICAse , N , INCx , INCy , MODe , PASs
  LOGICAL PASs
  REAL strue(9) , stemp(9)
  DOUBLE PRECISION dc , ds , da1(8) , db1(8) , dc1(8) , ds1(8)
  DOUBLE PRECISION da , datrue(8) , dbtrue(8) , dzero , Dfac , db
  DOUBLE PRECISION dab(4,9) , dtemp(9) , dtrue(9,9) , d12
  DATA zero , dzero/0. , 0.D0/
  DATA da1/.3D0 , .4D0 , -.3D0 , -.4D0 , -.3D0 , 0.D0 , 0.D0 , 1.D0/
  DATA db1/.4D0 , .3D0 , .4D0 , .3D0 , -.4D0 , 0.D0 , 1.D0 , 0.D0/
  DATA dc1/.6D0 , .8D0 , -.6D0 , .8D0 , .6D0 , 1.D0 , 0.D0 , 1.D0/
  DATA ds1/.8D0 , .6D0 , .8D0 , -.6D0 , .8D0 , 0.D0 , 1.D0 , 0.D0/
  DATA datrue/.5D0 , .5D0 , .5D0 , -.5D0 , -.5D0 , 0.D0 , 1.D0 , 1.D0/
  DATA dbtrue/0.D0 , .6D0 , 0.D0 , -.6D0 , 0.D0 , 0.D0 , 1.D0 , 0.D0/
  !                                              INPUT FOR MODIFIED GIVENS
  DATA dab/.1D0 , .3D0 , 1.2D0 , .2D0 , .7D0 , .2D0 , .6D0 , 4.2D0 , 0.D0 , &
    0.D0 , 0.D0 , 0.D0 , 4.D0 , -1.D0 , 2.D0 , 4.D0 , 6.D-10 , 2.D-2 , &
    1.D5 , 10.D0 , 4.D10 , 2.D-2 , 1.D-5 , 10.D0 , 2.D-10 , 4.D-2 , &
    1.D5 , 10.D0 , 2.D10 , 4.D-2 , 1.D-5 , 10.D0 , 4.D0 , -2.D0 , 8.D0 , &
    4.D0/
  !                                       TRUE RESULTS FOR MODIFIED GIVENS
  DATA dtrue/0.D0 , 0.D0 , 1.3D0 , .2D0 , 0.D0 , 0.D0 , 0.D0 , .5D0 , 0.D0 , &
    0.D0 , 0.D0 , 4.5D0 , 4.2D0 , 1.D0 , .5D0 , 0.D0 , 0.D0 , 0.D0 , &
    0.D0 , 0.D0 , 0.D0 , 0.D0 , -2.D0 , 0.D0 , 0.D0 , 0.D0 , 0.D0 , &
    0.D0 , 0.D0 , 0.D0 , 4.D0 , -1.D0 , 0.D0 , 0.D0 , 0.D0 , 0.D0 , &
    0.D0 , 15.D-3 , 0.D0 , 10.D0 , -1.D0 , 0.D0 , -1.D-4 , 0.D0 , 1.D0 , &
    0.D0 , 0.D0 , 6144.D-5 , 10.D0 , -1.D0 , 4096.D0 , -1.D6 , 0.D0 , &
    1.D0 , 0.D0 , 0.D0 , 15.D0 , 10.D0 , -1.D0 , 5.D-5 , 0.D0 , 1.D0 , &
    0.D0 , 0.D0 , 0.D0 , 15.D0 , 10.D0 , -1.D0 , 5.D5 , -4096.D0 , 1.D0 , &
    4096.D-6 , 0.D0 , 0.D0 , 7.D0 , 4.D0 , 0.D0 , 0.D0 , -.5D0 , -.25D0 , &
    0.D0/
  !                   4096 = 2 ** 12
  DATA d12/4096.D0/
  !***FIRST EXECUTABLE STATEMENT  CHECK0
  !
  !                   COMPUTE TRUE VALUES WHICH CANNOT BE PRESTORED
  !                   IN DECIMAL NOTATION.
  dtrue(1,1) = 12.D0/130.D0
  dtrue(2,1) = 36.D0/130.D0
  dtrue(7,1) = -1.D0/6.D0
  dtrue(1,2) = 14.D0/75.D0
  dtrue(2,2) = 49.D0/75.D0
  dtrue(9,2) = 1.D0/7.D0
  dtrue(1,5) = 45.D-11*(d12*d12)
  dtrue(3,5) = 4.D5/(3.D0*d12)
  dtrue(6,5) = 1.D0/d12
  dtrue(8,5) = 1.D4/(3.D0*d12)
  dtrue(1,6) = 4.D10/(1.5D0*d12*d12)
  dtrue(2,6) = 2.D-2/1.5D0
  dtrue(8,6) = 5.D-7*d12
  dtrue(1,7) = 4.D0/150.D0
  dtrue(2,7) = (2.D-10/1.5D0)*(d12*d12)
  dtrue(7,7) = -dtrue(6,5)
  dtrue(9,7) = 1.D4/d12
  dtrue(1,8) = dtrue(1,7)
  dtrue(2,8) = 2.D10/(1.5D0*d12*d12)
  dtrue(1,9) = 32.D0/7.D0
  dtrue(2,9) = -16.D0/7.D0
  dbtrue(1) = 1.D0/.6D0
  dbtrue(3) = -1.D0/.6D0
  dbtrue(5) = 1.D0/.6D0
  !
  jump = ICAse - 11
  DO k = 1 , 9
    !                        SET N=K FOR IDENTIFICATION IN OUTPUT IF ANY.
    N = k
    !                             BRANCH TO SELECT SUBPROGRAM TO BE TESTED.
    !
    SELECT CASE (jump)
      CASE (2)
        !                                                             13. DROTG
        IF ( k>8 ) EXIT
        da = da1(k)
        db = db1(k)
        CALL DROTG(da,db,dc,ds)
        CALL DTEST(1,da,datrue(k),datrue(k),Dfac,Kprint)
        CALL DTEST(1,db,dbtrue(k),dbtrue(k),Dfac,Kprint)
        CALL DTEST(1,dc,dc1(k),dc1(k),Dfac,Kprint)
        CALL DTEST(1,ds,ds1(k),ds1(k),Dfac,Kprint)
      CASE (3,4)
        GOTO 100
      CASE (5)
        !                                                             16. SROTMG
        DO i = 1 , 4
          stemp(i) = dab(i,k)
          stemp(i+4) = zero
        ENDDO
        stemp(9) = zero
        CALL SROTMG(stemp(1),stemp(2),stemp(3),stemp(4),stemp(5))
        !
        DO i = 1 , 9
          strue(i) = dtrue(i,k)
        ENDDO
        CALL STEST(9,stemp,strue,strue,Sfac,Kprint)
      CASE (6)
        !                                                             17. DROTMG
        DO i = 1 , 4
          dtemp(i) = dab(i,k)
          dtemp(i+4) = dzero
        ENDDO
        dtemp(9) = dzero
        CALL DROTMG(dtemp(1),dtemp(2),dtemp(3),dtemp(4),dtemp(5))
        CALL DTEST(9,dtemp,dtrue(1,k),dtrue(1,k),Dfac,Kprint)
      CASE DEFAULT
        !                                                             12. SROTG
        IF ( k>8 ) EXIT
        sa = da1(k)
        sb = db1(k)
        CALL SROTG(sa,sb,sc,ss)
        CALL STEST(1,sa,REAL(datrue(k)),REAL(datrue(k)),Sfac,Kprint)
        CALL STEST(1,sb,REAL(dbtrue(k)),REAL(dbtrue(k)),Sfac,Kprint)
        CALL STEST(1,sc,REAL(dc1(k)),REAL(dc1(k)),Sfac,Kprint)
        CALL STEST(1,ss,REAL(ds1(k)),REAL(ds1(k)),Sfac,Kprint)
    END SELECT
  ENDDO
  RETURN
  !                     THE FOLLOWING STOP SHOULD NEVER BE REACHED.
  100  STOP
END SUBROUTINE CHECK0
