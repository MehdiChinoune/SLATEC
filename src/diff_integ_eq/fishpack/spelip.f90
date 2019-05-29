!** SPELIP
SUBROUTINE SPELIP(Intl,Iorder,A,B,M,Mbdcnd,Bda,Alpha,Bdb,Beta,C,D,N,&
    Nbdcnd,Bdc,Gama,Bdd,Xnu,COFX,COFY,An,Bn,Cn,Dn,Un,Zn,Am,&
    Bm,Cm,Dm,Um,Zm,Grhs,Usol,Idmn,W,Pertrb,Ierror)
  !>
  !  Subsidiary to SEPELI
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (SPELIP-S)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     SPELIP sets up vectors and arrays for input to BLKTRI
  !     and computes a second order solution in USOL.  A return jump to
  !     SEPELI occurs if IORDER=2.  If IORDER=4 a fourth order
  !     solution is generated in USOL.
  !
  !***
  ! **See also:**  SEPELI
  !***
  ! **Routines called:**  BLKTRI, CHKSNG, DEFER, MINSOL, ORTHOG, TRISP
  !***
  ! COMMON BLOCKS    SPLPCM

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  USE SPLPCM, ONLY : L, AIT, BIT, CIT, DIT, DLX, DLX4, DLY, DLY4, IS, JS, K, KSWx, &
    KSWy, MIT, MS, NIT, NS, TDLx3, TDLy3
  INTEGER :: Idmn, Ierror, Intl, Iorder, M, Mbdcnd, N, Nbdcnd
  REAL :: A, Alpha, B, Beta, C, D, Gama, Pertrb, Xnu
  REAL :: Am(M), An(N), Bda(N+1), Bdb(N+1), Bdc(M+1), Bdd(M+1), Bm(M), Bn(N), Cm(M), &
    Cn(N), Dm(M), Dn(N), Grhs(Idmn,N+1), Um(M), Un(N), Usol(Idmn,N+1), W(:), &
    Zm(M), Zn(N)
  EXTERNAL :: COFX, COFY
  INTEGER :: i, i1, iord, j, mp, np
  REAL :: ai, ax1, axi, bi, bxi, ci, cxi, cxm, dj, dy1, dyj, ej, eyj, fj, fyj, &
    fyn, xi, yj
  LOGICAL :: singlr
  !* FIRST EXECUTABLE STATEMENT  SPELIP
  KSWx = Mbdcnd + 1
  KSWy = Nbdcnd + 1
  K = M + 1
  L = N + 1
  AIT = A
  BIT = B
  CIT = C
  DIT = D
  !
  !     SET RIGHT HAND SIDE VALUES FROM GRHS IN USOL ON THE INTERIOR
  !     AND NON-SPECIFIED BOUNDARIES.
  !
  DO i = 2, M
    DO j = 2, N
      Usol(i,j) = Grhs(i,j)
    END DO
  END DO
  IF ( KSWx/=2.AND.KSWx/=3 ) THEN
    DO j = 2, N
      Usol(1,j) = Grhs(1,j)
    END DO
  END IF
  IF ( KSWx/=2.AND.KSWx/=5 ) THEN
    DO j = 2, N
      Usol(K,j) = Grhs(K,j)
    END DO
  END IF
  IF ( KSWy/=2.AND.KSWy/=3 ) THEN
    DO i = 2, M
      Usol(i,1) = Grhs(i,1)
    END DO
  END IF
  IF ( KSWy/=2.AND.KSWy/=5 ) THEN
    DO i = 2, M
      Usol(i,L) = Grhs(i,L)
    END DO
  END IF
  IF ( KSWx/=2.AND.KSWx/=3.AND.KSWy/=2.AND.KSWy/=3 ) Usol(1,1) = Grhs(1,1)
  IF ( KSWx/=2.AND.KSWx/=5.AND.KSWy/=2.AND.KSWy/=3 ) Usol(K,1) = Grhs(K,1)
  IF ( KSWx/=2.AND.KSWx/=3.AND.KSWy/=2.AND.KSWy/=5 ) Usol(1,L) = Grhs(1,L)
  IF ( KSWx/=2.AND.KSWx/=5.AND.KSWy/=2.AND.KSWy/=5 ) Usol(K,L) = Grhs(K,L)
  i1 = 1
  !
  !     SET SWITCHES FOR PERIODIC OR NON-PERIODIC BOUNDARIES
  !
  mp = 1
  np = 1
  IF ( KSWx==1 ) mp = 0
  IF ( KSWy==1 ) np = 0
  !
  !     SET DLX,DLY AND SIZE OF BLOCK TRI-DIAGONAL SYSTEM GENERATED
  !     IN NINT,MINT
  !
  DLX = (BIT-AIT)/M
  MIT = K - 1
  IF ( KSWx==2 ) MIT = K - 2
  IF ( KSWx==4 ) MIT = K
  DLY = (DIT-CIT)/N
  NIT = L - 1
  IF ( KSWy==2 ) NIT = L - 2
  IF ( KSWy==4 ) NIT = L
  TDLx3 = 2.0*DLX**3
  DLX4 = DLX**4
  TDLy3 = 2.0*DLY**3
  DLY4 = DLY**4
  !
  !     SET SUBSCRIPT LIMITS FOR PORTION OF ARRAY TO INPUT TO BLKTRI
  !
  IS = 1
  JS = 1
  IF ( KSWx==2.OR.KSWx==3 ) IS = 2
  IF ( KSWy==2.OR.KSWy==3 ) JS = 2
  NS = NIT + JS - 1
  MS = MIT + IS - 1
  !
  !     SET X - DIRECTION
  !
  DO i = 1, MIT
    xi = AIT + (IS+i-2)*DLX
    CALL COFX(xi,ai,bi,ci)
    axi = (ai/DLX-0.5*bi)/DLX
    bxi = -2.*ai/DLX**2 + ci
    cxi = (ai/DLX+0.5*bi)/DLX
    Am(i) = axi
    Bm(i) = bxi
    Cm(i) = cxi
  END DO
  !
  !     SET Y DIRECTION
  !
  DO j = 1, NIT
    yj = CIT + (JS+j-2)*DLY
    CALL COFY(yj,dj,ej,fj)
    dyj = (dj/DLY-0.5*ej)/DLY
    eyj = (-2.*dj/DLY**2+fj)
    fyj = (dj/DLY+0.5*ej)/DLY
    An(j) = dyj
    Bn(j) = eyj
    Cn(j) = fyj
  END DO
  !
  !     ADJUST EDGES IN X DIRECTION UNLESS PERIODIC
  !
  ax1 = Am(1)
  cxm = Cm(MIT)
  SELECT CASE (KSWx)
    CASE (1)
    CASE (3)
      !
      !     DIRICHLET-MIXED IN X DIRECTION
      !
      Am(1) = 0.0
      Am(MIT) = Am(MIT) + cxm
      Bm(MIT) = Bm(MIT) - 2.*Beta*DLX*cxm
      Cm(MIT) = 0.0
    CASE (4)
      !
      !     MIXED - MIXED IN X DIRECTION
      !
      Am(1) = 0.0
      Bm(1) = Bm(1) + 2.*DLX*Alpha*ax1
      Cm(1) = Cm(1) + ax1
      Am(MIT) = Am(MIT) + cxm
      Bm(MIT) = Bm(MIT) - 2.*DLX*Beta*cxm
      Cm(MIT) = 0.0
    CASE (5)
      !
      !     MIXED-DIRICHLET IN X DIRECTION
      !
      Am(1) = 0.0
      Bm(1) = Bm(1) + 2.*Alpha*DLX*ax1
      Cm(1) = Cm(1) + ax1
      Cm(MIT) = 0.0
    CASE DEFAULT
      !
      !     DIRICHLET-DIRICHLET IN X DIRECTION
      !
      Am(1) = 0.0
      Cm(MIT) = 0.0
  END SELECT
  !
  !     ADJUST IN Y DIRECTION UNLESS PERIODIC
  !
  dy1 = An(1)
  fyn = Cn(NIT)
  SELECT CASE (KSWy)
    CASE (1)
    CASE (3)
      !
      !     DIRICHLET-MIXED IN Y DIRECTION
      !
      An(1) = 0.0
      An(NIT) = An(NIT) + fyn
      Bn(NIT) = Bn(NIT) - 2.*DLY*Xnu*fyn
      Cn(NIT) = 0.0
    CASE (4)
      !
      !     MIXED - MIXED DIRECTION IN Y DIRECTION
      !
      An(1) = 0.0
      Bn(1) = Bn(1) + 2.*DLY*Gama*dy1
      Cn(1) = Cn(1) + dy1
      An(NIT) = An(NIT) + fyn
      Bn(NIT) = Bn(NIT) - 2.0*DLY*Xnu*fyn
      Cn(NIT) = 0.0
    CASE (5)
      !
      !     MIXED-DIRICHLET IN Y DIRECTION
      !
      An(1) = 0.0
      Bn(1) = Bn(1) + 2.*DLY*Gama*dy1
      Cn(1) = Cn(1) + dy1
      Cn(NIT) = 0.0
    CASE DEFAULT
      !
      !     DIRICHLET-DIRICHLET IN Y DIRECTION
      !
      An(1) = 0.0
      Cn(NIT) = 0.0
  END SELECT
  IF ( KSWx/=1 ) THEN
    !
    !     ADJUST USOL ALONG X EDGE
    !
    DO j = JS, NS
      IF ( KSWx/=2.AND.KSWx/=3 ) THEN
        Usol(IS,j) = Usol(IS,j) + 2.0*DLX*ax1*Bda(j)
      ELSE
        Usol(IS,j) = Usol(IS,j) - ax1*Usol(1,j)
      END IF
      IF ( KSWx/=2.AND.KSWx/=5 ) THEN
        Usol(MS,j) = Usol(MS,j) - 2.0*DLX*cxm*Bdb(j)
      ELSE
        Usol(MS,j) = Usol(MS,j) - cxm*Usol(K,j)
      END IF
    END DO
  END IF
  IF ( KSWy/=1 ) THEN
    !
    !     ADJUST USOL ALONG Y EDGE
    !
    DO i = IS, MS
      IF ( KSWy/=2.AND.KSWy/=3 ) THEN
        Usol(i,JS) = Usol(i,JS) + 2.0*DLY*dy1*Bdc(i)
      ELSE
        Usol(i,JS) = Usol(i,JS) - dy1*Usol(i,1)
      END IF
      IF ( KSWy/=2.AND.KSWy/=5 ) THEN
        Usol(i,NS) = Usol(i,NS) - 2.0*DLY*fyn*Bdd(i)
      ELSE
        Usol(i,NS) = Usol(i,NS) - fyn*Usol(i,L)
      END IF
    END DO
  END IF
  !
  !     SAVE ADJUSTED EDGES IN GRHS IF IORDER=4
  !
  IF ( Iorder==4 ) THEN
    DO j = JS, NS
      Grhs(IS,j) = Usol(IS,j)
      Grhs(MS,j) = Usol(MS,j)
    END DO
    DO i = IS, MS
      Grhs(i,JS) = Usol(i,JS)
      Grhs(i,NS) = Usol(i,NS)
    END DO
  END IF
  iord = Iorder
  Pertrb = 0.0
  !
  !     CHECK IF OPERATOR IS SINGULAR
  !
  CALL CHKSNG(Mbdcnd,Nbdcnd,Alpha,Beta,Gama,Xnu,COFX,COFY,singlr)
  !
  !     COMPUTE NON-ZERO EIGENVECTOR IN NULL SPACE OF TRANSPOSE
  !     IF SINGULAR
  !
  IF ( singlr ) CALL TRISP(MIT,Am,Bm,Cm,Dm,Um,Zm)
  IF ( singlr ) CALL TRISP(NIT,An,Bn,Cn,Dn,Un,Zn)
  !
  !     MAKE INITIALIZATION CALL TO BLKTRI
  !
  IF ( Intl==0 ) CALL BLKTRI(Intl,np,NIT,An,Bn,Cn,mp,MIT,Am,Bm,Cm,Idmn,&
    Usol(IS,JS),Ierror,W)
  IF ( Ierror/=0 ) RETURN
  !
  !     ADJUST RIGHT HAND SIDE IF NECESSARY
  !
  100 CONTINUE
  IF ( singlr ) CALL ORTHOG(Usol,Idmn,Zn,Zm,Pertrb)
  !
  !     COMPUTE SOLUTION
  !
  CALL BLKTRI(i1,np,NIT,An,Bn,Cn,mp,MIT,Am,Bm,Cm,Idmn,Usol(IS,JS),Ierror,W)
  IF ( Ierror/=0 ) RETURN
  !
  !     SET PERIODIC BOUNDARIES IF NECESSARY
  !
  IF ( KSWx==1 ) THEN
    DO j = 1, L
      Usol(K,j) = Usol(1,j)
    END DO
  END IF
  IF ( KSWy==1 ) THEN
    DO i = 1, K
      Usol(i,L) = Usol(i,1)
    END DO
  END IF
  !
  !     MINIMIZE SOLUTION WITH RESPECT TO WEIGHTED LEAST SQUARES
  !     NORM IF OPERATOR IS SINGULAR
  !
  IF ( singlr ) CALL MINSOL(Usol,Idmn,Zn,Zm)
  !
  !     RETURN IF DEFERRED CORRECTIONS AND A FOURTH ORDER SOLUTION ARE
  !     NOT FLAGGED
  !
  IF ( iord==2 ) RETURN
  iord = 2
  !
  !     COMPUTE NEW RIGHT HAND SIDE FOR FOURTH ORDER SOLUTION
  !
  CALL DEFER(COFX,COFY,Idmn,Usol,Grhs)
  GOTO 100
END SUBROUTINE SPELIP
