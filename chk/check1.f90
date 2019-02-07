!*==CHECK1.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CHECK1
SUBROUTINE CHECK1(Sfac,Dfac,Kprint)
  IMPLICIT NONE
  !*--CHECK15
  !*** Start of declarations inserted by SPAG
  INTEGER i , ICAMAX , ICAse , IDAMAX , INCx , INCy , ISAMAX , jump , &
    Kprint , len , MODe , N , np1 , NPRint
  REAL sa , SASUM , SCASUM , SCNRM2 , Sfac , SNRM2 , stemp
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  CHECK1
  !***PURPOSE  (UNKNOWN)
  !***LIBRARY   SLATEC
  !***AUTHOR  Lawson, C. L., (JPL)
  !***DESCRIPTION
  !
  !     THIS SUBPROGRAM TESTS THE INCREMENTING AND ACCURACY OF THE LINEAR
  !     ALGEBRA SUBPROGRAMS 26 - 38 (SNRM2 TO ICAMAX). STORED RESULTS ARE
  !     COMPARED WITH THE RESULT RETURNED BY THE SUBPROGRAM.
  !
  !     THESE SUBPROGRAMS REQUIRE A SINGLE VECTOR ARGUMENT.
  !
  !     ICASE            DESIGNATES WHICH SUBPROGRAM TO TEST.
  !                      26 .LE. ICASE .LE. 38
  !     C. L. LAWSON, JPL, 1974 DEC 10, MAY 28
  !
  !***ROUTINES CALLED  CSCAL, CSSCAL, DASUM, DNRM2, DSCAL, DTEST, ICAMAX,
  !                    IDAMAX, ISAMAX, ITEST, SASUM, SCASUM, SCNRM2,
  !                    SNRM2, SSCAL, STEST
  !***COMMON BLOCKS    COMBLA
  !***REVISION HISTORY  (YYMMDD)
  !   741210  DATE WRITTEN
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   890911  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  CHECK1
  COMMON /COMBLA/ NPRint , ICAse , N , INCx , INCy , MODe , PASs
  LOGICAL PASs
  INTEGER itrue2(5) , itrue3(5)
  DOUBLE PRECISION da , dx(8)
  DOUBLE PRECISION dv(8,5,2)
  DOUBLE PRECISION Dfac
  DOUBLE PRECISION DNRM2 , DASUM
  DOUBLE PRECISION dtrue1(5) , dtrue3(5) , dtrue5(8,5,2)
  REAL strue2(5) , strue4(5) , strue(8) , sx(8)
  COMPLEX ca , cv(8,5,2) , ctrue5(8,5,2) , ctrue6(8,5,2) , cx(8)
  !
  DATA sa , da , ca/.3 , .3D0 , (.4,-.7)/
  DATA dv/.1D0 , 2.D0 , 2.D0 , 2.D0 , 2.D0 , 2.D0 , 2.D0 , 2.D0 , .3D0 , &
    3.D0 , 3.D0 , 3.D0 , 3.D0 , 3.D0 , 3.D0 , 3.D0 , .3D0 , -.4D0 , &
    4.D0 , 4.D0 , 4.D0 , 4.D0 , 4.D0 , 4.D0 , .2D0 , -.6D0 , .3D0 , &
    5.D0 , 5.D0 , 5.D0 , 5.D0 , 5.D0 , .1D0 , -.3D0 , .5D0 , -.1D0 , &
    6.D0 , 6.D0 , 6.D0 , 6.D0 , .1D0 , 8.D0 , 8.D0 , 8.D0 , 8.D0 , 8.D0 , &
    8.D0 , 8.D0 , .3D0 , 9.D0 , 9.D0 , 9.D0 , 9.D0 , 9.D0 , 9.D0 , 9.D0 , &
    .3D0 , 2.D0 , -.4D0 , 2.D0 , 2.D0 , 2.D0 , 2.D0 , 2.D0 , .2D0 , &
    3.D0 , -.6D0 , 5.D0 , .3D0 , 2.D0 , 2.D0 , 2.D0 , .1D0 , 4.D0 , &
    -.3D0 , 6.D0 , -.5D0 , 7.D0 , -.1D0 , 3.D0/
  !     COMPLEX TEST VECTORS
  DATA cv/(.1,.1) , (1.,2.) , (1.,2.) , (1.,2.) , (1.,2.) , (1.,2.) , &
    (1.,2.) , (1.,2.) , (.3,-.4) , (3.,4.) , (3.,4.) , (3.,4.) , (3.,4.)&
    , (3.,4.) , (3.,4.) , (3.,4.) , (.1,-.3) , (.5,-.1) , (5.,6.) , &
    (5.,6.) , (5.,6.) , (5.,6.) , (5.,6.) , (5.,6.) , (.1,.1) , (-.6,.1)&
    , (.1,-.3) , (7.,8.) , (7.,8.) , (7.,8.) , (7.,8.) , (7.,8.) , &
    (.3,.1) , (.1,.4) , (.4,.1) , (.1,.2) , (2.,3.) , (2.,3.) , (2.,3.) , &
    (2.,3.) , (.1,.1) , (4.,5.) , (4.,5.) , (4.,5.) , (4.,5.) , (4.,5.) , &
    (4.,5.) , (4.,5.) , (.3,-.4) , (6.,7.) , (6.,7.) , (6.,7.) , (6.,7.)&
    , (6.,7.) , (6.,7.) , (6.,7.) , (.1,-.3) , (8.,9.) , (.5,-.1) , &
    (2.,5.) , (2.,5.) , (2.,5.) , (2.,5.) , (2.,5.) , (.1,.1) , (3.,6.) , &
    (-.6,.1) , (4.,7.) , (.1,-.3) , (7.,2.) , (7.,2.) , (7.,2.) , (.3,.1)&
    , (5.,8.) , (.1,.4) , (6.,9.) , (.4,.1) , (8.,3.) , (.1,.2) , (9.,4.)&
    /
  !
  DATA strue2/.0 , .5 , .6 , .7 , .7/
  DATA strue4/.0 , .7 , 1. , 1.3 , 1.7/
  DATA dtrue1/.0D0 , .3D0 , .5D0 , .7D0 , .6D0/
  DATA dtrue3/.0D0 , .3D0 , .7D0 , 1.1D0 , 1.D0/
  DATA dtrue5/.10D0 , 2.D0 , 2.D0 , 2.D0 , 2.D0 , 2.D0 , 2.D0 , 2.D0 , &
    .09D0 , 3.D0 , 3.D0 , 3.D0 , 3.D0 , 3.D0 , 3.D0 , 3.D0 , .09D0 , &
    -.12D0 , 4.D0 , 4.D0 , 4.D0 , 4.D0 , 4.D0 , 4.D0 , .06D0 , -.18D0 , &
    .09D0 , 5.D0 , 5.D0 , 5.D0 , 5.D0 , 5.D0 , .03D0 , -.09D0 , .15D0 , &
    -.03D0 , 6.D0 , 6.D0 , 6.D0 , 6.D0 , .10D0 , 8.D0 , 8.D0 , 8.D0 , &
    8.D0 , 8.D0 , 8.D0 , 8.D0 , .09D0 , 9.D0 , 9.D0 , 9.D0 , 9.D0 , &
    9.D0 , 9.D0 , 9.D0 , .09D0 , 2.D0 , -.12D0 , 2.D0 , 2.D0 , 2.D0 , &
    2.D0 , 2.D0 , .06D0 , 3.D0 , -.18D0 , 5.D0 , .09D0 , 2.D0 , 2.D0 , &
    2.D0 , .03D0 , 4.D0 , -.09D0 , 6.D0 , -.15D0 , 7.D0 , -.03D0 , 3.D0/
  !
  DATA ctrue5/(.1,.1) , (1.,2.) , (1.,2.) , (1.,2.) , (1.,2.) , (1.,2.) , &
    (1.,2.) , (1.,2.) , (-.16,-.37) , (3.,4.) , (3.,4.) , (3.,4.) , &
    (3.,4.) , (3.,4.) , (3.,4.) , (3.,4.) , (-.17,-.19) , (.13,-.39) , &
    (5.,6.) , (5.,6.) , (5.,6.) , (5.,6.) , (5.,6.) , (5.,6.) , &
    (.11,-.03) , (-.17,.46) , (-.17,-.19) , (7.,8.) , (7.,8.) , (7.,8.) , &
    (7.,8.) , (7.,8.) , (.19,-.17) , (.32,.09) , (.23,-.24) , (.18,.01) , &
    (2.,3.) , (2.,3.) , (2.,3.) , (2.,3.) , (.1,.1) , (4.,5.) , (4.,5.) , &
    (4.,5.) , (4.,5.) , (4.,5.) , (4.,5.) , (4.,5.) , (-.16,-.37) , &
    (6.,7.) , (6.,7.) , (6.,7.) , (6.,7.) , (6.,7.) , (6.,7.) , (6.,7.) , &
    (-.17,-.19) , (8.,9.) , (.13,-.39) , (2.,5.) , (2.,5.) , (2.,5.) , &
    (2.,5.) , (2.,5.) , (.11,-.03) , (3.,6.) , (-.17,.46) , (4.,7.) , &
    (-.17,-.19) , (7.,2.) , (7.,2.) , (7.,2.) , (.19,-.17) , (5.,8.) , &
    (.32,.09) , (6.,9.) , (.23,-.24) , (8.,3.) , (.18,.01) , (9.,4.)/
  !
  DATA ctrue6/(.1,.1) , (1.,2.) , (1.,2.) , (1.,2.) , (1.,2.) , (1.,2.) , &
    (1.,2.) , (1.,2.) , (.09,-.12) , (3.,4.) , (3.,4.) , (3.,4.) , &
    (3.,4.) , (3.,4.) , (3.,4.) , (3.,4.) , (.03,-.09) , (.15,-.03) , &
    (5.,6.) , (5.,6.) , (5.,6.) , (5.,6.) , (5.,6.) , (5.,6.) , (.03,.03)&
    , (-.18,.03) , (.03,-.09) , (7.,8.) , (7.,8.) , (7.,8.) , (7.,8.) , &
    (7.,8.) , (.09,.03) , (.03,.12) , (.12,.03) , (.03,.06) , (2.,3.) , &
    (2.,3.) , (2.,3.) , (2.,3.) , (.1,.1) , (4.,5.) , (4.,5.) , (4.,5.) , &
    (4.,5.) , (4.,5.) , (4.,5.) , (4.,5.) , (.09,-.12) , (6.,7.) , &
    (6.,7.) , (6.,7.) , (6.,7.) , (6.,7.) , (6.,7.) , (6.,7.) , &
    (.03,-.09) , (8.,9.) , (.15,-.03) , (2.,5.) , (2.,5.) , (2.,5.) , &
    (2.,5.) , (2.,5.) , (.03,.03) , (3.,6.) , (-.18,.03) , (4.,7.) , &
    (.03,-.09) , (7.,2.) , (7.,2.) , (7.,2.) , (.09,.03) , (5.,8.) , &
    (.03,.12) , (6.,9.) , (.12,.03) , (8.,3.) , (.03,.06) , (9.,4.)/
  !
  !
  DATA itrue2/0 , 1 , 2 , 2 , 3/
  DATA itrue3/0 , 1 , 2 , 2 , 2/
  !***FIRST EXECUTABLE STATEMENT  CHECK1
  jump = ICAse - 25
  DO INCx = 1 , 2
    DO np1 = 1 , 5
      N = np1 - 1
      len = 2*MAX(N,1)
      !                                                  SET VECTOR ARGUMENTS.
      DO i = 1 , len
        sx(i) = dv(i,np1,INCx)
        dx(i) = dv(i,np1,INCx)
        cx(i) = cv(i,np1,INCx)
      ENDDO
      !
      !                        BRANCH TO INVOKE SUBPROGRAM TO BE TESTED.
      !
      SELECT CASE (jump)
        CASE (2)
          !                                                             27. DNRM2
          CALL DTEST(1,DNRM2(N,dx,INCx),dtrue1(np1),dtrue1(np1),Dfac,Kprint)
        CASE (3)
          !                                                             28. SCNRM2
          CALL STEST(1,SCNRM2(N,cx,INCx),strue2(np1),strue2(np1),Sfac,Kprint)
        CASE (4)
          !                                                             29. SASUM
          stemp = dtrue3(np1)
          CALL STEST(1,SASUM(N,sx,INCx),stemp,stemp,Sfac,Kprint)
        CASE (5)
          !                                                             30. DASUM
          CALL DTEST(1,DASUM(N,dx,INCx),dtrue3(np1),dtrue3(np1),Dfac,Kprint)
        CASE (6)
          !                                                             31. SCASUM
          CALL STEST(1,SCASUM(N,cx,INCx),strue4(np1),strue4(np1),Sfac,Kprint)
        CASE (7)
          !                                                             32. SSCALE
          CALL SSCAL(N,sa,sx,INCx)
          DO i = 1 , len
            strue(i) = dtrue5(i,np1,INCx)
          ENDDO
          CALL STEST(len,sx,strue,strue,Sfac,Kprint)
        CASE (8)
          !                                                             33. DSCALE
          CALL DSCAL(N,da,dx,INCx)
          CALL DTEST(len,dx,dtrue5(1,np1,INCx),dtrue5(1,np1,INCx),Dfac,Kprint)
        CASE (9)
          !                                                             34. CSCALE
          CALL CSCAL(N,ca,cx,INCx)
          CALL STEST(2*len,cx,ctrue5(1,np1,INCx),ctrue5(1,np1,INCx),Sfac,&
            Kprint)
        CASE (10)
          !                                                             35. CSSCAL
          CALL CSSCAL(N,sa,cx,INCx)
          CALL STEST(2*len,cx,ctrue6(1,np1,INCx),ctrue6(1,np1,INCx),Sfac,&
            Kprint)
        CASE (11)
          !                                                             36. ISAMAX
          CALL ITEST(1,ISAMAX(N,sx,INCx),itrue2(np1),Kprint)
        CASE (12)
          !                                                             37. IDAMAX
          CALL ITEST(1,IDAMAX(N,dx,INCx),itrue2(np1),Kprint)
        CASE (13)
          !                                                             38. ICAMAX
          CALL ITEST(1,ICAMAX(N,cx,INCx),itrue3(np1),Kprint)
        CASE DEFAULT
          !                                                             26. SNRM2
          stemp = dtrue1(np1)
          CALL STEST(1,SNRM2(N,sx,INCx),stemp,stemp,Sfac,Kprint)
      END SELECT
      !
    ENDDO
  ENDDO
END SUBROUTINE CHECK1
