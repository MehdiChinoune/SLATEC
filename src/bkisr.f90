!*==BKISR.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK BKISR
SUBROUTINE BKISR(X,N,Sum,Ierr)
  IMPLICIT NONE
  !*--BKISR5
  !***BEGIN PROLOGUE  BKISR
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to BSKIN
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (BKISR-S, DBKISR-D)
  !***AUTHOR  Amos, D. E., (SNLA)
  !***DESCRIPTION
  !
  !     BKISR computes repeated integrals of the K0 Bessel function
  !     by the series for N=0,1, and 2.
  !
  !***SEE ALSO  BSKIN
  !***ROUTINES CALLED  PSIXN, R1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   820601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !***END PROLOGUE  BKISR
  INTEGER i , Ierr , k , kk , kkn , k1 , N , np
  REAL ak , atol , bk , c , fk , fn , hx , hxs , pol , pr , Sum , tkp , &
    tol , trm , X , xln
  REAL PSIXN , R1MACH
  DIMENSION c(2)
  SAVE c
  !
  DATA c(1) , c(2)/1.57079632679489662E+00 , 1.0E0/
  !***FIRST EXECUTABLE STATEMENT  BKISR
  Ierr = 0
  tol = MAX(R1MACH(4),1.0E-18)
  IF ( X>=tol ) THEN
    pr = 1.0E0
    pol = 0.0E0
    IF ( N/=0 ) THEN
      DO i = 1 , N
        pol = -pol*X + c(i)
        pr = pr*X/i
      ENDDO
    ENDIF
    hx = X*0.5E0
    hxs = hx*hx
    xln = LOG(hx)
    np = N + 1
    tkp = 3.0E0
    fk = 2.0E0
    fn = N
    bk = 4.0E0
    ak = 2.0E0/((fn+1.0E0)*(fn+2.0E0))
    Sum = ak*(PSIXN(N+3)-PSIXN(3)+PSIXN(2)-xln)
    atol = Sum*tol*0.75E0
    DO k = 2 , 20
      ak = ak*(hxs/bk)*((tkp+1.0E0)/(tkp+fn+1.0E0))*(tkp/(tkp+fn))
      k1 = k + 1
      kk = k1 + k
      kkn = kk + N
      trm = (PSIXN(k1)+PSIXN(kkn)-PSIXN(kk)-xln)*ak
      Sum = Sum + trm
      IF ( ABS(trm)<=atol ) GOTO 100
      tkp = tkp + 2.0E0
      bk = bk + tkp
      fk = fk + 1.0E0
    ENDDO
    Ierr = 2
    GOTO 99999
    !-----------------------------------------------------------------------
    !     SMALL X CASE, X.LT.WORD TOLERANCE
    !-----------------------------------------------------------------------
  ELSEIF ( N>0 ) THEN
    Sum = c(N)
    RETURN
  ELSE
    hx = X*0.5E0
    Sum = PSIXN(1) - LOG(hx)
    RETURN
  ENDIF
  100  Sum = (Sum*hxs+PSIXN(np)-xln)*pr
  IF ( N==1 ) Sum = -Sum
  Sum = pol + Sum
  RETURN
  99999 END SUBROUTINE BKISR
