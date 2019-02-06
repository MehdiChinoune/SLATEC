!*==SPLPQX.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK SPLPQX
      SUBROUTINE SPLPQX(Lun,Kprint,Ipass)
      IMPLICIT NONE
!*--SPLPQX5
!*** Start of declarations inserted by SPAG
      REAL bl , bu , d , dattrv , duals , prgopt , primal , USRMAT , work , zero
      INTEGER i , ibasis , ic , icnt , ind , info , Ipass , isoln , iv , ivv , 
     &        iwork , j , kk , kount , Kprint , liw , Lun , lw , mm , mrelas
      INTEGER nvars
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  SPLPQX
!***PURPOSE  Quick check for SPLP.
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (SPLPQX-S, DPLPQX-D)
!***AUTHOR  (UNKNOWN)
!***ROUTINES CALLED  PASS, SCOPY, SPLP, USRMAT
!***REVISION HISTORY  (YYMMDD)
!   ??????  DATE WRITTEN
!   890911  Removed unnecessary intrinsics.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   901013  Added additional printout on failure.  (RWC)
!***END PROLOGUE  SPLPQX
      EXTERNAL USRMAT
      REAL costs(37)
      DIMENSION prgopt(50) , dattrv(210) , bl(60) , bu(60)
      DIMENSION ind(60) , primal(60) , duals(60) , ibasis(60)
      DIMENSION work(800) , iwork(900) , isoln(14)
      DIMENSION d(14,37)
!***FIRST EXECUTABLE STATEMENT  SPLPQX
      IF ( Kprint>=2 ) WRITE (Lun,99001)
99001 FORMAT ('1 SPLP QUICK CHECK')
      icnt = 1
      zero = 0.0
!
!     DEFINE WORKING ARRAY LENGTHS
!
      liw = 900
      lw = 800
      mrelas = 14
      nvars = 37
!
!     DEFINE THE ARRAY COSTS(*) FOR THE OBJECTIVE FUNCTION
!
      costs(1) = 1.030
      costs(2) = 0.985
      costs(3) = 0.997
      costs(4) = 1.036
      costs(5) = 1.005
      costs(6) = 0.980
      costs(7) = 1.004
      costs(8) = 0.993
      costs(9) = 1.018
      costs(10) = 0.947
      costs(11) = 0.910
      costs(12) = 1.028
      costs(13) = 0.957
      costs(14) = 1.025
      costs(15) = 1.036
      costs(16) = 1.060
      costs(17) = 0.954
      costs(18) = 0.891
      costs(19) = 0.921
      costs(20) = 1.040
      costs(21) = 0.912
      costs(22) = 0.926
      costs(23) = 1.000
      costs(24) = 0.000
      costs(25) = 0.000
      costs(26) = 0.000
      costs(27) = 0.000
      costs(28) = 0.000
      costs(29) = 0.000
      costs(30) = 0.000
      costs(31) = 0.000
      costs(32) = 0.000
      costs(33) = 0.000
      costs(34) = 0.000
      costs(35) = 0.000
      costs(36) = 0.000
      costs(37) = 0.000
!
!     PLACE THE NONZERO INFORMATION ABOUT THE MATRIX IN DATTRV(*)
!
      CALL SCOPY(14*37,zero,0,d,1)
      d(1,1) = 1.04000
      d(1,23) = 1.00000
      d(1,24) = -1.00000
      d(2,6) = 0.04125
      d(2,7) = 0.05250
      d(2,17) = 0.04875
      d(2,24) = 1.00000
      d(2,25) = -1.00000
      d(3,8) = 0.05625
      d(3,9) = 0.06875
      d(3,11) = 0.02250
      d(3,25) = 1.00000
      d(3,26) = -1.00000
      d(4,2) = 1.04000
      d(4,3) = 1.05375
      d(4,5) = 1.06125
      d(4,12) = 0.08000
      d(4,16) = 0.09375
      d(4,18) = 0.03750
      d(4,19) = 0.04625
      d(4,20) = 0.08125
      d(4,22) = 0.05250
      d(4,26) = 1.00000
      d(4,27) = -1.00000
      d(5,10) = 0.04375
      d(5,27) = 1.00000
      d(5,28) = -1.00000
      d(6,4) = 1.05875
      d(6,13) = 0.04500
      d(6,14) = 0.06375
      d(6,15) = 0.06625
      d(6,21) = 0.05000
      d(6,28) = 1.00000
      d(6,29) = -1.00000
      d(7,6) = 1.04125
      d(7,7) = 1.05250
      d(7,8) = 1.05625
      d(7,9) = 1.06875
      d(7,11) = 0.02250
      d(7,17) = 0.04875
      d(7,29) = 1.00000
      d(7,30) = -1.00000
      d(8,10) = 1.04375
      d(8,12) = 0.08000
      d(8,13) = 0.04500
      d(8,14) = 0.06375
      d(8,15) = 0.06625
      d(8,16) = 0.09375
      d(8,18) = 0.03750
      d(8,19) = 0.04625
      d(8,20) = 0.08125
      d(8,21) = 0.05000
      d(8,22) = 0.05250
      d(8,30) = 1.00000
      d(8,31) = -1.00000
      d(9,11) = 1.02250
      d(9,17) = 0.04875
      d(9,31) = 1.00000
      d(9,32) = -1.00000
      d(10,12) = 1.08000
      d(10,13) = 1.04500
      d(10,14) = 1.06375
      d(10,15) = 1.06625
      d(10,16) = 1.09375
      d(10,18) = 0.03750
      d(10,19) = 0.04625
      d(10,20) = 0.08125
      d(10,21) = 0.05000
      d(10,22) = 0.05250
      d(10,32) = 1.00000
      d(10,33) = -1.00000
      d(11,17) = 1.04875
      d(11,33) = 1.00000
      d(11,34) = -1.00000
      d(12,18) = 1.03750
      d(12,19) = 1.04625
      d(12,20) = 1.08125
      d(12,21) = 1.05000
      d(12,22) = 0.05250
      d(12,34) = 1.00000
      d(12,35) = -1.00000
      d(13,35) = 1.00000
      d(13,36) = -1.00000
      d(14,22) = 1.05250
      d(14,36) = 1.00000
      d(14,37) = -1.00000
      kount = 1
      DO mm = 1 , nvars
        dattrv(kount) = -mm
        DO kk = 1 , mrelas
          IF ( d(kk,mm)/=zero ) THEN
            kount = kount + 1
            dattrv(kount) = kk
            kount = kount + 1
            dattrv(kount) = d(kk,mm)
          ENDIF
        ENDDO
        kount = kount + 1
      ENDDO
      dattrv(kount) = zero
!
!     NON-NEGATIVITY CONSTRAINT
!
      DO ic = 1 , nvars
        bl(ic) = zero
        ind(ic) = 3
        bu(ic) = 10000000.000
      ENDDO
!
!     LE CONSTRAINTS
!
      DO iv = 1 , mrelas
        ivv = iv + nvars
        ind(ivv) = 3
        bl(ivv) = 100.00000
        bu(ivv) = 100000000.00000
      ENDDO
      prgopt(01) = 18
      prgopt(02) = 59
      prgopt(03) = 0
      prgopt(04) = 1
      prgopt(05) = 3
      prgopt(06) = 8
      prgopt(07) = 10
      prgopt(08) = 11
      prgopt(09) = 16
      prgopt(10) = 17
      prgopt(11) = 21
      prgopt(12) = 22
      prgopt(13) = 24
      prgopt(14) = 25
      prgopt(15) = 27
      prgopt(16) = 28
      prgopt(17) = 35
      prgopt(18) = 21
      prgopt(19) = 51
      prgopt(20) = 0
      prgopt(21) = 1
      CALL SPLP(USRMAT,mrelas,nvars,costs,prgopt,dattrv,bl,bu,ind,info,primal,
     &          duals,ibasis,work,lw,iwork,liw)
!
!     LOOK FOR THE KNOWN BASIS AT THE SOLN., NOW IS ISOLN(*).
!
      DO i = 1 , mrelas
        isoln(i) = prgopt(i+3)
      ENDDO
!
      Ipass = 1
      DO j = 1 , mrelas
        DO i = 1 , mrelas
          IF ( isoln(i)==ibasis(j) ) GOTO 100
        ENDDO
        Ipass = 0
        EXIT
 100  ENDDO
!
      IF ( Kprint>=2 ) WRITE (Lun,99002) (isoln(i),ibasis(i),i=1,mrelas)
!
99002 FORMAT (/'     ISOLN    IBASIS'/(2I10))
!
      IF ( Kprint>=2.OR.(Kprint==1.AND.Ipass/=1) ) CALL PASS(Lun,icnt,Ipass)
!
!     HERE IPASS=0 IF CODE FAILED QUICK CHECK;
!               =1 IF CODE PASSED QUICK CHECK.
!
      IF ( Kprint>=1.AND.Ipass/=1 ) WRITE (Lun,99003)
99003 FORMAT (/' ************ SPLP FAILED SOME TESTS ****************')
      IF ( Kprint>=2.AND.Ipass==1 ) WRITE (Lun,99004)
99004 FORMAT (/' ************ SPLP PASSED ALL TESTS *****************')
      RETURN
      END SUBROUTINE SPLPQX
