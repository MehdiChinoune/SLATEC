!** CBINU
SUBROUTINE CBINU(Z,Fnu,Kode,N,Cy,Nz,Rl,Fnul,Tol,Elim,Alim)
  !> Subsidiary to CAIRY, CBESH, CBESI, CBESJ, CBESK and CBIRY
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CBINU-A, ZBINU-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     CBINU COMPUTES THE I FUNCTION IN THE RIGHT HALF Z PLANE
  !
  !***
  ! **See also:**  CAIRY, CBESH, CBESI, CBESJ, CBESK, CBIRY
  !***
  ! **Routines called:**  CASYI, CBUNI, CMLRI, CSERI, CUOIK, CWRSK

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)

  INTEGER :: i, inw, Kode, N, nlast, nn, nui, nw, Nz
  COMPLEX(SP) :: cw(2), Cy(N), Z
  REAL(SP) :: Alim, az, dfnu, Elim, Fnu, Fnul, Rl, Tol
  COMPLEX(SP), PARAMETER :: czero = (0._SP,0._SP)
  !* FIRST EXECUTABLE STATEMENT  CBINU
  Nz = 0
  az = ABS(Z)
  nn = N
  dfnu = Fnu + (N-1)
  IF( az>2._SP ) THEN
    IF( az*az*0.25E0>dfnu+1._SP ) GOTO 100
  END IF
  !-----------------------------------------------------------------------
  !     POWER SERIES
  !-----------------------------------------------------------------------
  CALL CSERI(Z,Fnu,Kode,nn,Cy,nw,Tol,Elim,Alim)
  inw = ABS(nw)
  Nz = Nz + inw
  nn = nn - inw
  IF( nn==0 ) RETURN
  IF( nw>=0 ) GOTO 600
  dfnu = Fnu + (nn-1)
  100 CONTINUE
  IF( az>=Rl ) THEN
    IF( dfnu>1._SP ) THEN
      IF( az+az<dfnu*dfnu ) GOTO 200
    END IF
    !-----------------------------------------------------------------------
    !     ASYMPTOTIC EXPANSION FOR LARGE Z
    !-----------------------------------------------------------------------
    CALL CASYI(Z,Fnu,Kode,nn,Cy,nw,Rl,Tol,Elim,Alim)
    IF( nw>=0 ) GOTO 600
    GOTO 700
  ELSEIF( dfnu<=1._SP ) THEN
    GOTO 400
  END IF
  !-----------------------------------------------------------------------
  !     OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM
  !-----------------------------------------------------------------------
  200  CALL CUOIK(Z,Fnu,Kode,1,nn,Cy,nw,Tol,Elim,Alim)
  IF( nw<0 ) GOTO 700
  Nz = Nz + nw
  nn = nn - nw
  IF( nn==0 ) RETURN
  dfnu = Fnu + (nn-1)
  IF( dfnu>Fnul ) GOTO 500
  IF( az>Fnul ) GOTO 500
  300 CONTINUE
  IF( az>Rl ) THEN
    !-----------------------------------------------------------------------
    !     MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    !     OVERFLOW TEST ON K FUNCTIONS USED IN WRONSKIAN
    !-----------------------------------------------------------------------
    CALL CUOIK(Z,Fnu,Kode,2,2,cw,nw,Tol,Elim,Alim)
    IF( nw>=0 ) THEN
      IF( nw>0 ) GOTO 700
      CALL CWRSK(Z,Fnu,Kode,nn,Cy,nw,cw,Tol,Elim,Alim)
      IF( nw>=0 ) GOTO 600
      GOTO 700
    ELSE
      Nz = nn
      DO i = 1, nn
        Cy(i) = czero
      END DO
      RETURN
    END IF
  END IF
  !-----------------------------------------------------------------------
  !     MILLER ALGORITHM NORMALIZED BY THE SERIES
  !-----------------------------------------------------------------------
  400  CALL CMLRI(Z,Fnu,Kode,nn,Cy,nw,Tol)
  IF( nw>=0 ) GOTO 600
  GOTO 700
  !-----------------------------------------------------------------------
  !     INCREMENT FNU+NN-1 UP TO FNUL, COMPUTE AND RECUR BACKWARD
  !-----------------------------------------------------------------------
  500  nui = INT( Fnul - dfnu ) + 1
  nui = MAX(nui,0)
  CALL CBUNI(Z,Fnu,Kode,nn,Cy,nw,nui,nlast,Fnul,Tol,Elim,Alim)
  IF( nw<0 ) GOTO 700
  Nz = Nz + nw
  IF( nlast/=0 ) THEN
    nn = nlast
    GOTO 300
  END IF
  600  RETURN
  700  Nz = -1
  IF( nw==(-2) ) Nz = -2
END SUBROUTINE CBINU
