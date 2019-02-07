!*==ZBINU.f90  processed by SPAG 6.72Dc at 11:02 on  6 Feb 2019
!DECK ZBINU
SUBROUTINE ZBINU(Zr,Zi,Fnu,Kode,N,Cyr,Cyi,Nz,Rl,Fnul,Tol,Elim,Alim)
  IMPLICIT NONE
  !*--ZBINU5
  !***BEGIN PROLOGUE  ZBINU
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK and ZBIRY
  !***LIBRARY   SLATEC
  !***TYPE      ALL (CBINU-A, ZBINU-A)
  !***AUTHOR  Amos, D. E., (SNL)
  !***DESCRIPTION
  !
  !     ZBINU COMPUTES THE I FUNCTION IN THE RIGHT HALF Z PLANE
  !
  !***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBIRY
  !***ROUTINES CALLED  ZABS, ZASYI, ZBUNI, ZMLRI, ZSERI, ZUOIK, ZWRSK
  !***REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  ZBINU
  DOUBLE PRECISION Alim , az , cwi , cwr , Cyi , Cyr , dfnu , Elim , Fnu , &
    Fnul , Rl , Tol , zeroi , zeror , Zi , Zr , ZABS
  INTEGER i , inw , Kode , N , nlast , nn , nui , nw , Nz
  DIMENSION Cyr(N) , Cyi(N) , cwr(2) , cwi(2)
  EXTERNAL ZABS
  DATA zeror , zeroi/0.0D0 , 0.0D0/
  !***FIRST EXECUTABLE STATEMENT  ZBINU
  Nz = 0
  az = ZABS(Zr,Zi)
  nn = N
  dfnu = Fnu + (N-1)
  IF ( az>2.0D0 ) THEN
    IF ( az*az*0.25D0>dfnu+1.0D0 ) GOTO 100
  ENDIF
  !-----------------------------------------------------------------------
  !     POWER SERIES
  !-----------------------------------------------------------------------
  CALL ZSERI(Zr,Zi,Fnu,Kode,nn,Cyr,Cyi,nw,Tol,Elim,Alim)
  inw = ABS(nw)
  Nz = Nz + inw
  nn = nn - inw
  IF ( nn==0 ) RETURN
  IF ( nw>=0 ) GOTO 600
  dfnu = Fnu + (nn-1)
  100  IF ( az>=Rl ) THEN
  IF ( dfnu>1.0D0 ) THEN
    IF ( az+az<dfnu*dfnu ) GOTO 200
  ENDIF
  !-----------------------------------------------------------------------
  !     ASYMPTOTIC EXPANSION FOR LARGE Z
  !-----------------------------------------------------------------------
  CALL ZASYI(Zr,Zi,Fnu,Kode,nn,Cyr,Cyi,nw,Rl,Tol,Elim,Alim)
  IF ( nw>=0 ) GOTO 600
  GOTO 700
ELSEIF ( dfnu<=1.0D0 ) THEN
  GOTO 400
ENDIF
!-----------------------------------------------------------------------
!     OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM
!-----------------------------------------------------------------------
200  CALL ZUOIK(Zr,Zi,Fnu,Kode,1,nn,Cyr,Cyi,nw,Tol,Elim,Alim)
IF ( nw<0 ) GOTO 700
Nz = Nz + nw
nn = nn - nw
IF ( nn==0 ) RETURN
dfnu = Fnu + (nn-1)
IF ( dfnu>Fnul ) GOTO 500
IF ( az>Fnul ) GOTO 500
300  IF ( az>Rl ) THEN
!-----------------------------------------------------------------------
!     MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     OVERFLOW TEST ON K FUNCTIONS USED IN WRONSKIAN
!-----------------------------------------------------------------------
CALL ZUOIK(Zr,Zi,Fnu,Kode,2,2,cwr,cwi,nw,Tol,Elim,Alim)
IF ( nw>=0 ) THEN
  IF ( nw>0 ) GOTO 700
  CALL ZWRSK(Zr,Zi,Fnu,Kode,nn,Cyr,Cyi,nw,cwr,cwi,Tol,Elim,Alim)
  IF ( nw>=0 ) GOTO 600
  GOTO 700
ELSE
  Nz = nn
  DO i = 1 , nn
    Cyr(i) = zeror
    Cyi(i) = zeroi
  ENDDO
  RETURN
ENDIF
ENDIF
!-----------------------------------------------------------------------
!     MILLER ALGORITHM NORMALIZED BY THE SERIES
!-----------------------------------------------------------------------
400  CALL ZMLRI(Zr,Zi,Fnu,Kode,nn,Cyr,Cyi,nw,Tol)
IF ( nw>=0 ) GOTO 600
GOTO 700
!-----------------------------------------------------------------------
!     INCREMENT FNU+NN-1 UP TO FNUL, COMPUTE AND RECUR BACKWARD
!-----------------------------------------------------------------------
500  nui = Fnul - dfnu + 1
nui = MAX(nui,0)
CALL ZBUNI(Zr,Zi,Fnu,Kode,nn,Cyr,Cyi,nw,nui,nlast,Fnul,Tol,Elim,Alim)
IF ( nw<0 ) GOTO 700
Nz = Nz + nw
IF ( nlast/=0 ) THEN
nn = nlast
GOTO 300
ENDIF
600  RETURN
700  Nz = -1
IF ( nw==(-2) ) Nz = -2
END SUBROUTINE ZBINU
