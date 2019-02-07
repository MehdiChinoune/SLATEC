!*==COMPB.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK COMPB
SUBROUTINE COMPB(N,Ierror,An,Bn,Cn,B,Ah,Bh)
  IMPLICIT NONE
  !*--COMPB5
  !*** Start of declarations inserted by SPAG
  REAL Ah , An , arg , B , Bh , Bn , bnorm , Cn , CNV , d1 , d2 , d3 , EPS , &
    R1MACH
  INTEGER i , i2 , i4 , ib , Ierror , if , ifd , IK , ipl , ir , j , j1 , &
    j2 , jf , js , K , kdo , l , l1 , l2
  INTEGER lh , ls , N , n2m2 , nb , NCMplx , NM , nmp , NPP
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  COMPB
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to BLKTRI
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (COMPB-S, CCMPB-C)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !     COMPB computes the roots of the B polynomials using subroutine
  !     TEVLS which is a modification the EISPACK program TQLRAT.
  !     IERROR is set to 4 if either TEVLS fails or if A(J+1)*C(J) is
  !     less than zero for some J.  AH,BH are temporary work arrays.
  !
  !***SEE ALSO  BLKTRI
  !***ROUTINES CALLED  INDXB, PPADD, R1MACH, TEVLS
  !***COMMON BLOCKS    CBLKT
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  COMPB
  !
  DIMENSION An(*) , Bn(*) , Cn(*) , B(*) , Ah(*) , Bh(*)
  COMMON /CBLKT / NPP , K , EPS , CNV , NM , NCMplx , IK
  !***FIRST EXECUTABLE STATEMENT  COMPB
  EPS = R1MACH(4)
  bnorm = ABS(Bn(1))
  DO j = 2 , NM
    bnorm = MAX(bnorm,ABS(Bn(j)))
    arg = An(j)*Cn(j-1)
    IF ( arg<0 ) GOTO 200
    B(j) = SIGN(SQRT(arg),An(j))
  ENDDO
  CNV = EPS*bnorm
  if = 2**K
  kdo = K - 1
  DO l = 1 , kdo
    ir = l - 1
    i2 = 2**ir
    i4 = i2 + i2
    ipl = i4 - 1
    ifd = if - i4
    DO i = i4 , ifd , i4
      CALL INDXB(i,l,ib,nb)
      IF ( nb<=0 ) EXIT
      js = i - ipl
      jf = js + nb - 1
      ls = 0
      DO j = js , jf
        ls = ls + 1
        Bh(ls) = Bn(j)
        Ah(ls) = B(j)
      ENDDO
      CALL TEVLS(nb,Bh,Ah,Ierror)
      IF ( Ierror/=0 ) GOTO 100
      lh = ib - 1
      DO j = 1 , nb
        lh = lh + 1
        B(lh) = -Bh(j)
      ENDDO
    ENDDO
  ENDDO
  DO j = 1 , NM
    B(j) = -Bn(j)
  ENDDO
  IF ( NPP==0 ) THEN
    nmp = NM + 1
    nb = NM + nmp
    DO j = 1 , nb
      l1 = MOD(j-1,nmp) + 1
      l2 = MOD(j+NM-1,nmp) + 1
      arg = An(l1)*Cn(l2)
      IF ( arg<0 ) GOTO 200
      Bh(j) = SIGN(SQRT(arg),-An(l1))
      Ah(j) = -Bn(l1)
    ENDDO
    CALL TEVLS(nb,Ah,Bh,Ierror)
    IF ( Ierror/=0 ) GOTO 100
    CALL INDXB(if,K-1,j2,lh)
    CALL INDXB(if/2,K-1,j1,lh)
    j2 = j2 + 1
    lh = j2
    n2m2 = j2 + NM + NM - 2
    DO
      d1 = ABS(B(j1)-B(j2-1))
      d2 = ABS(B(j1)-B(j2))
      d3 = ABS(B(j1)-B(j2+1))
      IF ( (d2<d1).AND.(d2<d3) ) THEN
        j2 = j2 + 1
        j1 = j1 + 1
        IF ( j2>n2m2 ) EXIT
      ELSE
        B(lh) = B(j2)
        j2 = j2 + 1
        lh = lh + 1
        IF ( j2>n2m2 ) EXIT
      ENDIF
    ENDDO
    B(lh) = B(n2m2+1)
    CALL INDXB(if,K-1,j1,j2)
    j2 = j1 + nmp + nmp
    CALL PPADD(NM+1,Ierror,An,Cn,B(j1),B(j1),B(j2))
  ENDIF
  RETURN
  100  Ierror = 4
  RETURN
  200  Ierror = 5
END SUBROUTINE COMPB
