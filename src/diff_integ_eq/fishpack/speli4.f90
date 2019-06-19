!** SPELI4
SUBROUTINE SPELI4(Iorder,A,B,M,Mbdcnd,Bda,Alpha,Bdb,Beta,C,D,N,Nbdcnd,Bdc,&
    Bdd,COFX,An,Bn,Cn,Dn,Un,Zn,Am,Bm,Cm,Dm,Um,Zm,Grhs,Usol,Idmn,W,Pertrb,Ierror)
  !> Subsidiary to SEPX4
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (SPELI4-S)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     SPELI4 sets up vectors and arrays for input to BLKTRI
  !     and computes a second order solution in USOL.  A return jump to
  !     SEPX4 occurs if IORDER=2.  If IORDER=4 a fourth order
  !     solution is generated in USOL.
  !
  !***
  ! **See also:**  SEPX4
  !***
  ! **Routines called:**  CHKSN4, DEFE4, GENBUN, MINSO4, ORTHO4, TRIS4
  !***
  ! COMMON BLOCKS    SPL4

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891009  Removed unreferenced variable.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  USE SPL4, ONLY : l_com, ait_com, bit_com, cit_com, dit_com, dlx_com, dlx4_com, &
    dly_com, dly4_com, is_com, js_com, k_com, kswx_com, kswy_com, mit_com, ms_com, &
    nit_com, ns_com, tdlx3_com, tdly3_com
  INTERFACE
    SUBROUTINE COFX(X,A,B,C)
      IMPORT SP
      REAL(SP) :: X, A, B, C
    END SUBROUTINE COFX
  END INTERFACE
  INTEGER :: Idmn, Ierror, Iorder, M, Mbdcnd, N, Nbdcnd
  REAL(SP) :: A, Alpha, B, Beta, C, D, Pertrb
  REAL(SP) :: Am(M+1), An(N+1), Bda(N+1), Bdb(N+1), Bdc(M+1), Bdd(M+1), Bm(M+1), &
    Bn(N+1), Cm(M+1), Cn(N+1), Dm(M+1), Dn(N+1), Grhs(Idmn,N), Um(M+1), Un(N+1), &
    Usol(Idmn,N+1), W(:), Zm(M+1), Zn(N+1)
  INTEGER :: i, ieror, iord, j, mp, np
  REAL(SP) :: ai, ax1, axi, bi, bxi, ci, cxi, cxm, dy1, dyj, eyj, fyj, fyn, gama, xi, xnu
  LOGICAL :: singlr
  !* FIRST EXECUTABLE STATEMENT  SPELI4
  kswx_com = Mbdcnd + 1
  kswy_com = Nbdcnd + 1
  k_com = M + 1
  l_com = N + 1
  ait_com = A
  bit_com = B
  cit_com = C
  dit_com = D
  dly_com = (dit_com-cit_com)/N
  !
  !     SET RIGHT HAND SIDE VALUES FROM GRHS IN USOL ON THE INTERIOR
  !     AND NON-SPECIFIED BOUNDARIES.
  !
  DO i = 2, M
    DO j = 2, N
      Usol(i,j) = dly_com**2*Grhs(i,j)
    END DO
  END DO
  IF( kswx_com/=2 .AND. kswx_com/=3 ) THEN
    DO j = 2, N
      Usol(1,j) = dly_com**2*Grhs(1,j)
    END DO
  END IF
  IF( kswx_com/=2 .AND. kswx_com/=5 ) THEN
    DO j = 2, N
      Usol(k_com,j) = dly_com**2*Grhs(k_com,j)
    END DO
  END IF
  IF( kswy_com/=2 .AND. kswy_com/=3 ) THEN
    DO i = 2, M
      Usol(i,1) = dly_com**2*Grhs(i,1)
    END DO
  END IF
  IF( kswy_com/=2 .AND. kswy_com/=5 ) THEN
    DO i = 2, M
      Usol(i,l_com) = dly_com**2*Grhs(i,l_com)
    END DO
  END IF
  IF( kswx_com/=2 .AND. kswx_com/=3 .AND. kswy_com/=2 .AND. kswy_com/=3 ) Usol(1,1) = dly_com**2*Grhs(1,1)
  IF( kswx_com/=2 .AND. kswx_com/=5 .AND. kswy_com/=2 .AND. kswy_com/=3 ) Usol(k_com,1) = dly_com**2*Grhs(k_com,1)
  IF( kswx_com/=2 .AND. kswx_com/=3 .AND. kswy_com/=2 .AND. kswy_com/=5 ) Usol(1,l_com) = dly_com**2*Grhs(1,l_com)
  IF( kswx_com/=2 .AND. kswx_com/=5 .AND. kswy_com/=2 .AND. kswy_com/=5 ) Usol(k_com,l_com) = dly_com**2*Grhs(k_com,l_com)
  !
  !     SET SWITCHES FOR PERIODIC OR NON-PERIODIC BOUNDARIES
  !
  mp = 1
  IF( kswx_com==1 ) mp = 0
  np = Nbdcnd
  !
  !     SET DLX,DLY AND SIZE OF BLOCK TRI-DIAGONAL SYSTEM GENERATED
  !     IN NINT,MINT
  !
  dlx_com = (bit_com-ait_com)/M
  mit_com = k_com - 1
  IF( kswx_com==2 ) mit_com = k_com - 2
  IF( kswx_com==4 ) mit_com = k_com
  dly_com = (dit_com-cit_com)/N
  nit_com = l_com - 1
  IF( kswy_com==2 ) nit_com = l_com - 2
  IF( kswy_com==4 ) nit_com = l_com
  tdlx3_com = 2.0*dlx_com**3
  dlx4_com = dlx_com**4
  tdly3_com = 2.0*dly_com**3
  dly4_com = dly_com**4
  !
  !     SET SUBSCRIPT LIMITS FOR PORTION OF ARRAY TO INPUT TO BLKTRI
  !
  is_com = 1
  js_com = 1
  IF( kswx_com==2 .OR. kswx_com==3 ) is_com = 2
  IF( kswy_com==2 .OR. kswy_com==3 ) js_com = 2
  ns_com = nit_com + js_com - 1
  ms_com = mit_com + is_com - 1
  !
  !     SET X - DIRECTION
  !
  DO i = 1, mit_com
    xi = ait_com + (is_com+i-2)*dlx_com
    CALL COFX(xi,ai,bi,ci)
    axi = (ai/dlx_com-0.5*bi)/dlx_com
    bxi = -2.*ai/dlx_com**2 + ci
    cxi = (ai/dlx_com+0.5*bi)/dlx_com
    Am(i) = dly_com**2*axi
    Bm(i) = dly_com**2*bxi
    Cm(i) = dly_com**2*cxi
  END DO
  !
  !     SET Y DIRECTION
  !
  DO j = 1, nit_com
    dyj = 1.0
    eyj = -2.0
    fyj = 1.0
    An(j) = dyj
    Bn(j) = eyj
    Cn(j) = fyj
  END DO
  !
  !     ADJUST EDGES IN X DIRECTION UNLESS PERIODIC
  !
  ax1 = Am(1)
  cxm = Cm(mit_com)
  SELECT CASE (kswx_com)
    CASE (1)
    CASE (3)
      !
      !     DIRICHLET-MIXED IN X DIRECTION
      !
      Am(1) = 0.0
      Am(mit_com) = Am(mit_com) + cxm
      Bm(mit_com) = Bm(mit_com) - 2.*Beta*dlx_com*cxm
      Cm(mit_com) = 0.0
    CASE (4)
      !
      !     MIXED - MIXED IN X DIRECTION
      !
      Am(1) = 0.0
      Bm(1) = Bm(1) + 2.*dlx_com*Alpha*ax1
      Cm(1) = Cm(1) + ax1
      Am(mit_com) = Am(mit_com) + cxm
      Bm(mit_com) = Bm(mit_com) - 2.*dlx_com*Beta*cxm
      Cm(mit_com) = 0.0
    CASE (5)
      !
      !     MIXED-DIRICHLET IN X DIRECTION
      !
      Am(1) = 0.0
      Bm(1) = Bm(1) + 2.*Alpha*dlx_com*ax1
      Cm(1) = Cm(1) + ax1
      Cm(mit_com) = 0.0
    CASE DEFAULT
      !
      !     DIRICHLET-DIRICHLET IN X DIRECTION
      !
      Am(1) = 0.0
      Cm(mit_com) = 0.0
  END SELECT
  !
  !     ADJUST IN Y DIRECTION UNLESS PERIODIC
  !
  dy1 = An(1)
  fyn = Cn(nit_com)
  gama = 0.0
  xnu = 0.0
  SELECT CASE (kswy_com)
    CASE (1)
    CASE (3)
      !
      !     DIRICHLET-MIXED IN Y DIRECTION
      !
      An(1) = 0.0
      An(nit_com) = An(nit_com) + fyn
      Bn(nit_com) = Bn(nit_com) - 2.*dly_com*xnu*fyn
      Cn(nit_com) = 0.0
    CASE (4)
      !
      !     MIXED - MIXED DIRECTION IN Y DIRECTION
      !
      An(1) = 0.0
      Bn(1) = Bn(1) + 2.*dly_com*gama*dy1
      Cn(1) = Cn(1) + dy1
      An(nit_com) = An(nit_com) + fyn
      Bn(nit_com) = Bn(nit_com) - 2.0*dly_com*xnu*fyn
      Cn(nit_com) = 0.0
    CASE (5)
      !
      !     MIXED-DIRICHLET IN Y DIRECTION
      !
      An(1) = 0.0
      Bn(1) = Bn(1) + 2.*dly_com*gama*dy1
      Cn(1) = Cn(1) + dy1
      Cn(nit_com) = 0.0
    CASE DEFAULT
      !
      !     DIRICHLET-DIRICHLET IN Y DIRECTION
      !
      An(1) = 0.0
      Cn(nit_com) = 0.0
  END SELECT
  IF( kswx_com/=1 ) THEN
    !
    !     ADJUST USOL ALONG X EDGE
    !
    DO j = js_com, ns_com
      IF( kswx_com/=2 .AND. kswx_com/=3 ) THEN
        Usol(is_com,j) = Usol(is_com,j) + 2.0*dlx_com*ax1*Bda(j)
      ELSE
        Usol(is_com,j) = Usol(is_com,j) - ax1*Usol(1,j)
      END IF
      IF( kswx_com/=2 .AND. kswx_com/=5 ) THEN
        Usol(ms_com,j) = Usol(ms_com,j) - 2.0*dlx_com*cxm*Bdb(j)
      ELSE
        Usol(ms_com,j) = Usol(ms_com,j) - cxm*Usol(k_com,j)
      END IF
    END DO
  END IF
  IF( kswy_com/=1 ) THEN
    !
    !     ADJUST USOL ALONG Y EDGE
    !
    DO i = is_com, ms_com
      IF( kswy_com/=2 .AND. kswy_com/=3 ) THEN
        Usol(i,js_com) = Usol(i,js_com) + 2.0*dly_com*dy1*Bdc(i)
      ELSE
        Usol(i,js_com) = Usol(i,js_com) - dy1*Usol(i,1)
      END IF
      IF( kswy_com/=2 .AND. kswy_com/=5 ) THEN
        Usol(i,ns_com) = Usol(i,ns_com) - 2.0*dly_com*fyn*Bdd(i)
      ELSE
        Usol(i,ns_com) = Usol(i,ns_com) - fyn*Usol(i,l_com)
      END IF
    END DO
  END IF
  !
  !     SAVE ADJUSTED EDGES IN GRHS IF IORDER=4
  !
  IF( Iorder==4 ) THEN
    DO j = js_com, ns_com
      Grhs(is_com,j) = Usol(is_com,j)
      Grhs(ms_com,j) = Usol(ms_com,j)
    END DO
    DO i = is_com, ms_com
      Grhs(i,js_com) = Usol(i,js_com)
      Grhs(i,ns_com) = Usol(i,ns_com)
    END DO
  END IF
  iord = Iorder
  Pertrb = 0.0
  !
  !     CHECK IF OPERATOR IS SINGULAR
  !
  CALL CHKSN4(Mbdcnd,Nbdcnd,Alpha,Beta,COFX,singlr)
  !
  !     COMPUTE NON-ZERO EIGENVECTOR IN NULL SPACE OF TRANSPOSE
  !     IF SINGULAR
  !
  IF( singlr ) CALL TRIS4(mit_com,Am,Bm,Cm,Dm,Um,Zm)
  IF( singlr ) CALL TRIS4(nit_com,An,Bn,Cn,Dn,Un,Zn)
  DO
    !
    !     ADJUST RIGHT HAND SIDE IF NECESSARY
    !
    IF( singlr ) CALL ORTHO4(Usol,Idmn,Zn,Zm,Pertrb)
    !
    !     COMPUTE SOLUTION
    !
    !     SAVE ADJUSTED RIGHT HAND SIDE IN GRHS
    DO j = js_com, ns_com
      DO i = is_com, ms_com
        Grhs(i,j) = Usol(i,j)
      END DO
    END DO
    CALL GENBUN(np,nit_com,mp,mit_com,Am,Bm,Cm,Idmn,Usol(is_com,js_com),ieror,W)
    !     CHECK IF ERROR DETECTED IN POIS
    !     THIS CAN ONLY CORRESPOND TO IERROR=12
    IF( ieror==0 ) THEN
      IF( Ierror/=0 ) RETURN
      !
      !     SET PERIODIC BOUNDARIES IF NECESSARY
      !
      IF( kswx_com==1 ) THEN
        DO j = 1, l_com
          Usol(k_com,j) = Usol(1,j)
        END DO
      END IF
      IF( kswy_com==1 ) THEN
        DO i = 1, k_com
          Usol(i,l_com) = Usol(i,1)
        END DO
      END IF
      !
      !     MINIMIZE SOLUTION WITH RESPECT TO WEIGHTED LEAST SQUARES
      !     NORM IF OPERATOR IS SINGULAR
      !
      IF( singlr ) CALL MINSO4(Usol,Idmn,Zn,Zm)
      !
      !     RETURN IF DEFERRED CORRECTIONS AND A FOURTH ORDER SOLUTION ARE
      !     NOT FLAGGED
      !
      IF( iord==2 ) RETURN
      iord = 2
      !
      !     COMPUTE NEW RIGHT HAND SIDE FOR FOURTH ORDER SOLUTION
      !
      CALL DEFE4(COFX,Idmn,Usol,Grhs)
    ELSE
      !     SET ERROR FLAG IF IMPROPER COEFFICIENTS INPUT TO POIS
      Ierror = 12
      RETURN
    END IF
  END DO
END SUBROUTINE SPELI4
