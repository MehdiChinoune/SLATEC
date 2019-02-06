*DECK CHECK0
      SUBROUTINE CHECK0 (SFAC, DFAC, KPRINT)
C***BEGIN PROLOGUE  CHECK0
C***PURPOSE  (UNKNOWN)
C***LIBRARY   SLATEC
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Wisniewski, J. A., (SNLA)
C***DESCRIPTION
C
C     THIS SUBROUTINE TESTS SUBPROGRAMS 12-13 AND 16-17.
C     THESE SUBPROGRAMS HAVE NO ARRAY ARGUMENTS.
C
C     C. L. LAWSON, JPL, 1975 MAR 07, MAY 28
C     R. J. HANSON, J. A. WISNIEWSKI, SANDIA LABS, APRIL 25,1977.
C
C***ROUTINES CALLED  DROTG, DROTMG, DTEST, SROTG, SROTMG, STEST
C***COMMON BLOCKS    COMBLA
C***REVISION HISTORY  (YYMMDD)
C   750307  DATE WRITTEN
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  CHECK0
      COMMON /COMBLA/ NPRINT, ICASE, N, INCX, INCY, MODE, PASS
      LOGICAL          PASS
      REAL             STRUE(9),STEMP(9)
      DOUBLE PRECISION DC,DS,DA1(8),DB1(8),DC1(8),DS1(8)
      DOUBLE PRECISION DA,DATRUE(8),DBTRUE(8),DZERO,DFAC,DB
      DOUBLE PRECISION DAB(4,9),DTEMP(9),DTRUE(9,9),D12
      DATA ZERO, DZERO / 0., 0.D0 /
      DATA DA1/ .3D0,  .4D0, -.3D0, -.4D0, -.3D0,  0.D0,  0.D0,  1.D0/
      DATA DB1/ .4D0,  .3D0,  .4D0,  .3D0, -.4D0,  0.D0,  1.D0,  0.D0/
      DATA DC1/ .6D0,  .8D0, -.6D0,  .8D0,  .6D0,  1.D0,  0.D0,  1.D0/
      DATA DS1/ .8D0,  .6D0,  .8D0, -.6D0,  .8D0,  0.D0,  1.D0,  0.D0/
      DATA DATRUE/ .5D0,  .5D0,  .5D0, -.5D0, -.5D0, 0.D0, 1.D0, 1.D0/
      DATA DBTRUE/ 0.D0,  .6D0,  0.D0, -.6D0,  0.D0, 0.D0, 1.D0, 0.D0/
C                                              INPUT FOR MODIFIED GIVENS
      DATA DAB/ .1D0,.3D0,1.2D0,.2D0,
     A          .7D0, .2D0, .6D0, 4.2D0,
     B          0.D0,0.D0,0.D0,0.D0,
     C          4.D0, -1.D0, 2.D0, 4.D0,
     D          6.D-10, 2.D-2, 1.D5, 10.D0,
     E          4.D10, 2.D-2, 1.D-5, 10.D0,
     F          2.D-10, 4.D-2, 1.D5, 10.D0,
     G          2.D10, 4.D-2, 1.D-5, 10.D0,
     H          4.D0, -2.D0, 8.D0, 4.D0    /
C                                       TRUE RESULTS FOR MODIFIED GIVENS
      DATA DTRUE/0.D0,0.D0, 1.3D0, .2D0, 0.D0,0.D0,0.D0, .5D0, 0.D0,
     A           0.D0,0.D0, 4.5D0, 4.2D0, 1.D0, .5D0, 0.D0,0.D0,0.D0,
     B           0.D0,0.D0,0.D0,0.D0, -2.D0, 0.D0,0.D0,0.D0,0.D0,
     C           0.D0,0.D0,0.D0, 4.D0, -1.D0, 0.D0,0.D0,0.D0,0.D0,
     D           0.D0, 15.D-3, 0.D0, 10.D0, -1.D0, 0.D0, -1.D-4,
     E           0.D0, 1.D0,
     F           0.D0,0.D0, 6144.D-5, 10.D0, -1.D0, 4096.D0, -1.D6,
     G           0.D0, 1.D0,
     H           0.D0,0.D0,15.D0,10.D0,-1.D0, 5.D-5, 0.D0,1.D0,0.D0,
     I           0.D0,0.D0, 15.D0, 10.D0, -1. D0, 5.D5, -4096.D0,
     J           1.D0, 4096.D-6,
     K           0.D0,0.D0, 7.D0, 4.D0, 0.D0,0.D0, -.5D0, -.25D0, 0.D0/
C                   4096 = 2 ** 12
      DATA D12  /4096.D0/
C***FIRST EXECUTABLE STATEMENT  CHECK0
C
C                   COMPUTE TRUE VALUES WHICH CANNOT BE PRESTORED
C                   IN DECIMAL NOTATION.
      DTRUE(1,1) = 12.D0 / 130.D0
      DTRUE(2,1) = 36.D0 / 130.D0
      DTRUE(7,1) = -1.D0 / 6.D0
      DTRUE(1,2) = 14.D0 / 75.D0
      DTRUE(2,2) = 49.D0 / 75.D0
      DTRUE(9,2) = 1.D0 / 7.D0
      DTRUE(1,5) = 45.D-11 * (D12 * D12)
      DTRUE(3,5) = 4.D5 / (3.D0 * D12)
      DTRUE(6,5) = 1.D0 / D12
      DTRUE(8,5) = 1.D4 / (3.D0 * D12)
      DTRUE(1,6) = 4.D10 / (1.5D0 * D12 * D12)
      DTRUE(2,6) = 2.D-2 / 1.5D0
      DTRUE(8,6) = 5.D-7 * D12
      DTRUE(1,7) = 4.D0 / 150.D0
      DTRUE(2,7) = (2.D-10 / 1.5D0) * (D12 * D12)
      DTRUE(7,7) = -DTRUE(6,5)
      DTRUE(9,7) = 1.D4 / D12
      DTRUE(1,8) = DTRUE(1,7)
      DTRUE(2,8) = 2.D10 / (1.5D0 * D12 * D12)
      DTRUE(1,9) = 32.D0 / 7.D0
      DTRUE(2,9) = -16.D0 / 7.D0
      DBTRUE(1) = 1.D0/.6D0
      DBTRUE(3) = -1.D0/.6D0
      DBTRUE(5) = 1.D0/.6D0
C
      JUMP= ICASE-11
          DO 500 K = 1, 9
C                        SET N=K FOR IDENTIFICATION IN OUTPUT IF ANY.
          N=K
C                             BRANCH TO SELECT SUBPROGRAM TO BE TESTED.
C
          GO TO (120,130,999,999,160,170), JUMP
C                                                             12. SROTG
  120 IF(K.GT.8) GO TO 600
          SA = DA1(K)
          SB = DB1(K)
          CALL SROTG(SA,SB,SC,SS)
          CALL STEST(1,SA,REAL(DATRUE(K)),REAL(DATRUE(K)),SFAC,KPRINT)
          CALL STEST(1,SB,REAL(DBTRUE(K)),REAL(DBTRUE(K)),SFAC,KPRINT)
          CALL STEST(1,SC,REAL(DC1(K)),REAL(DC1(K)),SFAC,KPRINT)
          CALL STEST(1,SS,REAL(DS1(K)),REAL(DS1(K)),SFAC,KPRINT)
          GO TO 500
C                                                             13. DROTG
  130 IF(K.GT.8) GO TO 600
          DA = DA1(K)
          DB = DB1(K)
          CALL DROTG(DA,DB,DC,DS)
          CALL DTEST(1,DA,DATRUE(K),DATRUE(K),DFAC,KPRINT)
          CALL DTEST(1,DB,DBTRUE(K),DBTRUE(K),DFAC,KPRINT)
          CALL DTEST(1,DC,DC1(K),DC1(K),DFAC,KPRINT)
          CALL DTEST(1,DS,DS1(K),DS1(K),DFAC,KPRINT)
          GO TO 500
C                                                             16. SROTMG
  160     CONTINUE
               DO 162 I = 1, 4
               STEMP(I) = DAB(I,K)
               STEMP(I+4) = ZERO
  162          CONTINUE
           STEMP(9) = ZERO
           CALL SROTMG(STEMP(1),STEMP(2),STEMP(3),STEMP(4),STEMP(5))
C
               DO 166 I = 1, 9
  166          STRUE(I) = DTRUE(I,K)
          CALL STEST(9,STEMP,STRUE,STRUE,SFAC,KPRINT)
          GO TO 500
C                                                             17. DROTMG
  170     CONTINUE
               DO 172 I = 1, 4
               DTEMP(I) = DAB(I,K)
               DTEMP(I+4) = DZERO
  172          CONTINUE
          DTEMP(9) = DZERO
          CALL DROTMG(DTEMP(1),DTEMP(2),DTEMP(3),DTEMP(4),DTEMP(5))
          CALL DTEST(9,DTEMP,DTRUE(1,K),DTRUE(1,K),DFAC,KPRINT)
  500     CONTINUE
  600 RETURN
C                     THE FOLLOWING STOP SHOULD NEVER BE REACHED.
  999 STOP
      END
