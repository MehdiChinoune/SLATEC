!*==DPLPQX.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DPLPQX
SUBROUTINE DPLPQX(Lun,Kprint,Ipass)
  !***BEGIN PROLOGUE  DPLPQX
  !***PURPOSE  Quick check for DSPLP.
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (SPLPQX-S, DPLPQX-D)
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  DCOPY, DSPLP, DUSRMT, PASS
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901013  Added additional printout on failure.  (RWC)
  !***END PROLOGUE  DPLPQX
  IMPLICIT NONE
  !*--DPLPQX17
  !*** Start of declarations inserted by SPAG
  DOUBLE PRECISION DUSRMT
  INTEGER i , ic , iv , ivv , j , kk , kount , Kprint , Lun , mm
  !*** End of declarations inserted by SPAG
  EXTERNAL DUSRMT
  INTEGER icnt , ind(60) , ibasis(60) , Ipass , iwork(900) , isoln(14)
  DOUBLE PRECISION costs(37)
  DOUBLE PRECISION prgopt(50) , dattrv(210) , bl(60) , bu(60)
  DOUBLE PRECISION primal(60) , duals(60)
  DOUBLE PRECISION work(800)
  DOUBLE PRECISION d(14,37)
  DOUBLE PRECISION zero
  INTEGER mrelas , nvars , info , lw , liw
  !***FIRST EXECUTABLE STATEMENT  DPLPQX
  IF ( Kprint>=2 ) WRITE (Lun,99001)
  99001 FORMAT ('1 DSPLP QUICK CHECK')
  icnt = 1
  zero = 0.0D0
  Ipass = 0
  !     DEFINE WORKING ARRAY LENGTHS
  liw = 900
  lw = 800
  mrelas = 14
  nvars = 37
  !     DEFINE THE ARRAY COSTS(*) FOR THE OBJECTIVE FUNCTION
  costs(1) = 1.030D0
  costs(2) = 0.985D0
  costs(3) = 0.997D0
  costs(4) = 1.036D0
  costs(5) = 1.005D0
  costs(6) = 0.980D0
  costs(7) = 1.004D0
  costs(8) = 0.993D0
  costs(9) = 1.018D0
  costs(10) = 0.947D0
  costs(11) = 0.910D0
  costs(12) = 1.028D0
  costs(13) = 0.957D0
  costs(14) = 1.025D0
  costs(15) = 1.036D0
  costs(16) = 1.060D0
  costs(17) = 0.954D0
  costs(18) = 0.891D0
  costs(19) = 0.921D0
  costs(20) = 1.040D0
  costs(21) = 0.912D0
  costs(22) = 0.926D0
  costs(23) = 1.000D0
  costs(24) = 0.000D0
  costs(25) = 0.000D0
  costs(26) = 0.000D0
  costs(27) = 0.000D0
  costs(28) = 0.000D0
  costs(29) = 0.000D0
  costs(30) = 0.000D0
  costs(31) = 0.000D0
  costs(32) = 0.000D0
  costs(33) = 0.000D0
  costs(34) = 0.000D0
  costs(35) = 0.000D0
  costs(36) = 0.000D0
  costs(37) = 0.000D0
  !     PLACE THE NONZERO INFORMATION ABOUT THE MATRIX IN DATTRV(*)
  CALL DCOPY(14*37,zero,0,d,1)
  d(1,1) = 1.04000D0
  d(1,23) = 1.00000D0
  d(1,24) = -1.00000D0
  d(2,6) = 0.04125D0
  d(2,7) = 0.05250D0
  d(2,17) = 0.04875D0
  d(2,24) = 1.00000D0
  d(2,25) = -1.00000D0
  d(3,8) = 0.05625D0
  d(3,9) = 0.06875D0
  d(3,11) = 0.02250D0
  d(3,25) = 1.00000D0
  d(3,26) = -1.00000D0
  d(4,2) = 1.04000D0
  d(4,3) = 1.05375D0
  d(4,5) = 1.06125D0
  d(4,12) = 0.08000D0
  d(4,16) = 0.09375D0
  d(4,18) = 0.03750D0
  d(4,19) = 0.04625D0
  d(4,20) = 0.08125D0
  d(4,22) = 0.05250D0
  d(4,26) = 1.00000D0
  d(4,27) = -1.00000D0
  d(5,10) = 0.04375D0
  d(5,27) = 1.00000D0
  d(5,28) = -1.00000D0
  d(6,4) = 1.05875D0
  d(6,13) = 0.04500D0
  d(6,14) = 0.06375D0
  d(6,15) = 0.06625D0
  d(6,21) = 0.05000D0
  d(6,28) = 1.00000D0
  d(6,29) = -1.00000D0
  d(7,6) = 1.04125D0
  d(7,7) = 1.05250D0
  d(7,8) = 1.05625D0
  d(7,9) = 1.06875D0
  d(7,11) = 0.02250D0
  d(7,17) = 0.04875D0
  d(7,29) = 1.00000D0
  d(7,30) = -1.00000D0
  d(8,10) = 1.04375D0
  d(8,12) = 0.08000D0
  d(8,13) = 0.04500D0
  d(8,14) = 0.06375D0
  d(8,15) = 0.06625D0
  d(8,16) = 0.09375D0
  d(8,18) = 0.03750D0
  d(8,19) = 0.04625D0
  d(8,20) = 0.08125D0
  d(8,21) = 0.05000D0
  d(8,22) = 0.05250D0
  d(8,30) = 1.00000D0
  d(8,31) = -1.00000D0
  d(9,11) = 1.02250D0
  d(9,17) = 0.04875D0
  d(9,31) = 1.00000D0
  d(9,32) = -1.00000D0
  d(10,12) = 1.08000D0
  d(10,13) = 1.04500D0
  d(10,14) = 1.06375D0
  d(10,15) = 1.06625D0
  d(10,16) = 1.09375D0
  d(10,18) = 0.03750D0
  d(10,19) = 0.04625D0
  d(10,20) = 0.08125D0
  d(10,21) = 0.05000D0
  d(10,22) = 0.05250D0
  d(10,32) = 1.00000D0
  d(10,33) = -1.00000D0
  d(11,17) = 1.04875D0
  d(11,33) = 1.00000D0
  d(11,34) = -1.00000D0
  d(12,18) = 1.03750D0
  d(12,19) = 1.04625D0
  d(12,20) = 1.08125D0
  d(12,21) = 1.05000D0
  d(12,22) = 0.05250D0
  d(12,34) = 1.00000D0
  d(12,35) = -1.00000D0
  d(13,35) = 1.00000D0
  d(13,36) = -1.00000D0
  d(14,22) = 1.05250D0
  d(14,36) = 1.00000D0
  d(14,37) = -1.00000D0
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
  !     NON-NEGATIVITY CONSTRAINT
  DO ic = 1 , nvars
    bl(ic) = zero
    ind(ic) = 3
    bu(ic) = 10000000.000D0
  ENDDO
  !     LE CONSTRAINTS
  DO iv = 1 , mrelas
    ivv = iv + nvars
    ind(ivv) = 3
    bl(ivv) = 100.00000D0
    bu(ivv) = 100000000.00000D0
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
  CALL DSPLP(DUSRMT,mrelas,nvars,costs,prgopt,dattrv,bl,bu,ind,info,primal,&
    duals,ibasis,work,lw,iwork,liw)
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
    99003 FORMAT (/' ************ DSPLP FAILED SOME TESTS ***************')
    IF ( Kprint>=2.AND.Ipass==1 ) WRITE (Lun,99004)
    99004 FORMAT (/' ************ DSPLP PASSED ALL TESTS ****************')
    RETURN
END SUBROUTINE DPLPQX
