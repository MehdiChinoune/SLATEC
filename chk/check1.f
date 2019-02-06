*DECK CHECK1
      SUBROUTINE CHECK1 (SFAC, DFAC, KPRINT)
C***BEGIN PROLOGUE  CHECK1
C***PURPOSE  (UNKNOWN)
C***LIBRARY   SLATEC
C***AUTHOR  Lawson, C. L., (JPL)
C***DESCRIPTION
C
C     THIS SUBPROGRAM TESTS THE INCREMENTING AND ACCURACY OF THE LINEAR
C     ALGEBRA SUBPROGRAMS 26 - 38 (SNRM2 TO ICAMAX). STORED RESULTS ARE
C     COMPARED WITH THE RESULT RETURNED BY THE SUBPROGRAM.
C
C     THESE SUBPROGRAMS REQUIRE A SINGLE VECTOR ARGUMENT.
C
C     ICASE            DESIGNATES WHICH SUBPROGRAM TO TEST.
C                      26 .LE. ICASE .LE. 38
C     C. L. LAWSON, JPL, 1974 DEC 10, MAY 28
C
C***ROUTINES CALLED  CSCAL, CSSCAL, DASUM, DNRM2, DSCAL, DTEST, ICAMAX,
C                    IDAMAX, ISAMAX, ITEST, SASUM, SCASUM, SCNRM2,
C                    SNRM2, SSCAL, STEST
C***COMMON BLOCKS    COMBLA
C***REVISION HISTORY  (YYMMDD)
C   741210  DATE WRITTEN
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  CHECK1
      COMMON /COMBLA/ NPRINT, ICASE, N, INCX, INCY, MODE, PASS
      LOGICAL          PASS
      INTEGER          ITRUE2(5),ITRUE3(5)
      DOUBLE PRECISION DA,DX(8)
      DOUBLE PRECISION DV(8,5,2)
      DOUBLE PRECISION DFAC
      DOUBLE PRECISION DNRM2,DASUM
      DOUBLE PRECISION DTRUE1(5),DTRUE3(5),DTRUE5(8,5,2)
      REAL             STRUE2(5),STRUE4(5),STRUE(8),SX(8)
      COMPLEX          CA,CV(8,5,2),CTRUE5(8,5,2),CTRUE6(8,5,2),CX(8)
C
      DATA SA, DA, CA        / .3, .3D0, (.4,-.7)    /
      DATA DV/.1D0,2.D0,2.D0,2.D0,2.D0,2.D0,2.D0,2.D0,
     1        .3D0,3.D0,3.D0,3.D0,3.D0,3.D0,3.D0,3.D0,
     2        .3D0,-.4D0,4.D0,4.D0,4.D0,4.D0,4.D0,4.D0,
     3        .2D0,-.6D0,.3D0,5.D0,5.D0,5.D0,5.D0,5.D0,
     4        .1D0,-.3D0,.5D0,-.1D0,6.D0,6.D0,6.D0,6.D0,
     5        .1D0,8.D0,8.D0,8.D0,8.D0,8.D0,8.D0,8.D0,
     6        .3D0,9.D0,9.D0,9.D0,9.D0,9.D0,9.D0,9.D0,
     7        .3D0,2.D0,-.4D0,2.D0,2.D0,2.D0,2.D0,2.D0,
     8        .2D0,3.D0,-.6D0,5.D0,.3D0,2.D0,2.D0,2.D0,
     9         .1D0,4.D0,-.3D0,6.D0,-.5D0,7.D0,-.1D0,              3.D0/
C     COMPLEX TEST VECTORS
      DATA CV/
     1(.1,.1),(1.,2.),(1.,2.),(1.,2.),(1.,2.),(1.,2.),(1.,2.),(1.,2.),
     2(.3,-.4),(3.,4.),(3.,4.),(3.,4.),(3.,4.),(3.,4.),(3.,4.),(3.,4.),
     3(.1,-.3),(.5,-.1),(5.,6.),(5.,6.),(5.,6.),(5.,6.),(5.,6.),(5.,6.),
     4(.1,.1),(-.6,.1),(.1,-.3),(7.,8.),(7.,8.),(7.,8.),(7.,8.),(7.,8.),
     5(.3,.1),(.1,.4),(.4,.1),(.1,.2),(2.,3.),(2.,3.),(2.,3.),(2.,3.),
     6(.1,.1),(4.,5.),(4.,5.),(4.,5.),(4.,5.),(4.,5.),(4.,5.),(4.,5.),
     7(.3,-.4),(6.,7.),(6.,7.),(6.,7.),(6.,7.),(6.,7.),(6.,7.),(6.,7.),
     8(.1,-.3),(8.,9.),(.5,-.1),(2.,5.),(2.,5.),(2.,5.),(2.,5.),(2.,5.),
     9(.1,.1),(3.,6.),(-.6,.1),(4.,7.),(.1,-.3),(7.,2.),(7.,2.),(7.,2.),
     T(.3,.1),(5.,8.),(.1,.4),(6.,9.),(.4,.1),(8.,3.),(.1,.2),(9.,4.) /
C
      DATA STRUE2/.0,.5,.6,.7,.7/
      DATA STRUE4/.0,.7,1.,1.3,1.7/
      DATA DTRUE1/.0D0,.3D0,.5D0,.7D0,.6D0/
      DATA DTRUE3/.0D0,.3D0,.7D0,1.1D0,1.D0/
      DATA DTRUE5/.10D0,2.D0,2.D0,2.D0,2.D0,2.D0,2.D0,2.D0,
     1            .09D0,3.D0,3.D0,3.D0,3.D0,3.D0,3.D0,3.D0,
     2            .09D0,-.12D0,4.D0,4.D0,4.D0,4.D0,4.D0,4.D0,
     3            .06D0,-.18D0,.09D0,5.D0,5.D0,5.D0,5.D0,5.D0,
     4            .03D0,-.09D0,.15D0,-.03D0,6.D0,6.D0,6.D0,6.D0,
     5            .10D0,8.D0,8.D0,8.D0,8.D0,8.D0,8.D0,8.D0,
     6            .09D0,9.D0,9.D0,9.D0,9.D0,9.D0,9.D0,9.D0,
     7            .09D0,2.D0,-.12D0,2.D0,2.D0,2.D0,2.D0,2.D0,
     8            .06D0,3.D0,-.18D0,5.D0,.09D0,2.D0,2.D0,2.D0,
     9            .03D0,4.D0, -.09D0,6.D0, -.15D0,7.D0, -.03D0,  3.D0/
C
      DATA CTRUE5/
     A(.1,.1),(1.,2.),(1.,2.),(1.,2.),(1.,2.),(1.,2.),(1.,2.),(1.,2.),
     B(-.16,-.37),(3.,4.),(3.,4.),(3.,4.),(3.,4.),(3.,4.),(3.,4.),
     C                                                         (3.,4.),
     D(-.17,-.19),(.13,-.39),(5.,6.),(5.,6.),(5.,6.),(5.,6.),(5.,6.),
     E                                                         (5.,6.),
     F(.11,-.03),(-.17,.46),(-.17,-.19),(7.,8.),(7.,8.),(7.,8.),(7.,8.),
     G                                                         (7.,8.),
     H(.19,-.17),(.32,.09),(.23,-.24),(.18,.01),(2.,3.),(2.,3.),(2.,3.),
     I                                                         (2.,3.),
     J(.1,.1),(4.,5.),(4.,5.),(4.,5.),(4.,5.),(4.,5.),(4.,5.),(4.,5.),
     K(-.16,-.37),(6.,7.),(6.,7.),(6.,7.),(6.,7.),(6.,7.),(6.,7.),
     L                                                         (6.,7.),
     M(-.17,-.19),(8.,9.),(.13,-.39),(2.,5.),(2.,5.),(2.,5.),(2.,5.),
     N                                                         (2.,5.),
     O(.11,-.03),(3.,6.),(-.17,.46),(4.,7.),(-.17,-.19),(7.,2.),(7.,2.),
     P                                                         (7.,2.),
     Q(.19,-.17),(5.,8.),(.32,.09),(6.,9.),(.23,-.24),(8.,3.),(.18,.01),
     R                                                         (9.,4.) /
C
      DATA CTRUE6/
     A(.1,.1),(1.,2.),(1.,2.),(1.,2.),(1.,2.),(1.,2.),(1.,2.),(1.,2.),
     B(.09,-.12),(3.,4.),(3.,4.),(3.,4.),(3.,4.),(3.,4.),(3.,4.),
     C                                                         (3.,4.),
     D(.03,-.09),(.15,-.03),(5.,6.),(5.,6.),(5.,6.),(5.,6.),(5.,6.),
     E                                                         (5.,6.),
     F(.03,.03),(-.18,.03),(.03,-.09),(7.,8.),(7.,8.),(7.,8.),(7.,8.),
     G                                                         (7.,8.),
     H(.09,.03),(.03,.12),(.12,.03),(.03,.06),(2.,3.),(2.,3.),(2.,3.),
     I                                                         (2.,3.),
     J(.1,.1),(4.,5.),(4.,5.),(4.,5.),(4.,5.),(4.,5.),(4.,5.),(4.,5.),
     K(.09,-.12),(6.,7.),(6.,7.),(6.,7.),(6.,7.),(6.,7.),(6.,7.),
     L                                                         (6.,7.),
     M(.03,-.09),(8.,9.),(.15,-.03),(2.,5.),(2.,5.),(2.,5.),(2.,5.),
     N                                                         (2.,5.),
     O(.03,.03),(3.,6.),(-.18,.03),(4.,7.),(.03,-.09),(7.,2.),(7.,2.),
     P                                                         (7.,2.),
     Q(.09,.03),(5.,8.),(.03,.12),(6.,9.),(.12,.03),(8.,3.),(.03,.06),
     R                                                         (9.,4.) /
C
C
      DATA ITRUE2/ 0, 1, 2, 2, 3/
      DATA ITRUE3/ 0, 1, 2, 2, 2/
C***FIRST EXECUTABLE STATEMENT  CHECK1
      JUMP=ICASE-25
         DO 520 INCX=1,2
            DO 500 NP1=1,5
            N=NP1-1
            LEN= 2*MAX(N,1)
C                                                  SET VECTOR ARGUMENTS.
                    DO 22 I = 1, LEN
                    SX(I) = DV(I,NP1,INCX)
                    DX(I) = DV(I,NP1,INCX)
   22               CX(I) = CV(I,NP1,INCX)
C
C                        BRANCH TO INVOKE SUBPROGRAM TO BE TESTED.
C
               GO TO (260,270,280,290,300,310,320,
     *                330,340,350,360,370,380),JUMP
C                                                             26. SNRM2
  260       STEMP = DTRUE1(NP1)
            CALL STEST(1,SNRM2(N,SX,INCX),STEMP,STEMP,SFAC,KPRINT)
            GO TO 500
C                                                             27. DNRM2
  270       CALL DTEST(1,DNRM2(N,DX,INCX),DTRUE1(NP1),DTRUE1(NP1),DFAC,
     1                 KPRINT)
            GO TO 500
C                                                             28. SCNRM2
  280       CALL STEST(1,SCNRM2(N,CX,INCX),STRUE2(NP1),STRUE2(NP1),
     1                 SFAC,KPRINT)
            GO TO 500
C                                                             29. SASUM
  290       STEMP = DTRUE3(NP1)
            CALL STEST(1,SASUM(N,SX,INCX),STEMP,STEMP,SFAC,KPRINT)
            GO TO 500
C                                                             30. DASUM
  300       CALL DTEST(1,DASUM(N,DX,INCX),DTRUE3(NP1),DTRUE3(NP1),DFAC,
     1                 KPRINT)
            GO TO 500
C                                                             31. SCASUM
  310       CALL STEST(1,SCASUM(N,CX,INCX),STRUE4(NP1),STRUE4(NP1),SFAC,
     1                 KPRINT)
            GO TO 500
C                                                             32. SSCALE
  320       CALL SSCAL(N,SA,SX,INCX)
               DO 322 I = 1, LEN
  322          STRUE(I) = DTRUE5(I,NP1,INCX)
            CALL STEST(LEN,SX,STRUE,STRUE,SFAC,KPRINT)
            GO TO 500
C                                                             33. DSCALE
  330       CALL DSCAL(N,DA,DX,INCX)
           CALL DTEST(LEN,DX,DTRUE5(1,NP1,INCX),DTRUE5(1,NP1,INCX),
     1                 DFAC,KPRINT)
            GO TO 500
C                                                             34. CSCALE
  340       CALL CSCAL(N,CA,CX,INCX)
        CALL STEST(2*LEN,CX,CTRUE5(1,NP1,INCX),CTRUE5(1,NP1,INCX),
     1                 SFAC,KPRINT)
            GO TO 500
C                                                             35. CSSCAL
  350       CALL CSSCAL(N,SA,CX,INCX)
         CALL STEST(2*LEN,CX,CTRUE6(1,NP1,INCX),CTRUE6(1,NP1,INCX),
     1                 SFAC,KPRINT)
            GO TO 500
C                                                             36. ISAMAX
  360       CALL ITEST(1,ISAMAX(N,SX,INCX),ITRUE2(NP1),KPRINT)
            GO TO 500
C                                                             37. IDAMAX
  370       CALL ITEST(1,IDAMAX(N,DX,INCX),ITRUE2(NP1),KPRINT)
            GO TO 500
C                                                             38. ICAMAX
  380       CALL ITEST(1,ICAMAX(N,CX,INCX),ITRUE3(NP1),KPRINT)
C
  500       CONTINUE
  520    CONTINUE
      RETURN
      END
