!** SPELIP
SUBROUTINE SPELIP(Intl,Iorder,A,B,M,Mbdcnd,Bda,Alpha,Bdb,Beta,C,D,N,&
    Nbdcnd,Bdc,Gama,Bdd,Xnu,COFX,COFY,An,Bn,Cn,Dn,Un,Zn,Am,&
    Bm,Cm,Dm,Um,Zm,Grhs,Usol,Idmn,W,Pertrb,Ierror)
  !> Subsidiary to SEPELI
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
  USE SPLPCM, ONLY : l_com, ait_com, bit_com, cit_com, dit_com, dlx_com, dlx4_com, &
    dly_com, dly4_com, is_com, js_com, k_com, kswx_com, kswy_com, mit_com, ms_com, &
    nit_com, ns_com, tdlx3_com, tdly3_com
  INTERFACE
    SUBROUTINE COFX(X,A,B,C)
      IMPORT SP
      REAL(SP) :: X, A, B, C
    END SUBROUTINE COFX
    SUBROUTINE COFY(Y,D,E,F)
      IMPORT SP
      REAL(SP) :: Y, D, E, F
    END SUBROUTINE COFY
  END INTERFACE
  INTEGER :: Idmn, Ierror, Intl, Iorder, M, Mbdcnd, N, Nbdcnd
  REAL(SP) :: A, Alpha, B, Beta, C, D, Gama, Pertrb, Xnu
  REAL(SP) :: Am(M), An(N), Bda(N+1), Bdb(N+1), Bdc(M+1), Bdd(M+1), Bm(M), Bn(N), Cm(M), &
    Cn(N), Dm(M), Dn(N), Grhs(Idmn,N+1), Um(M), Un(N), Usol(Idmn,N+1), W(:), &
    Zm(M), Zn(N)
  INTEGER :: i, i1, iord, j, mp, np
  REAL(SP) :: ai, ax1, axi, bi, bxi, ci, cxi, cxm, dj, dy1, dyj, ej, eyj, fj, fyj, &
    fyn, xi, yj
  LOGICAL :: singlr
  !* FIRST EXECUTABLE STATEMENT  SPELIP
  kswx_com = Mbdcnd + 1
  kswy_com = Nbdcnd + 1
  k_com = M + 1
  l_com = N + 1
  ait_com = A
  bit_com = B
  cit_com = C
  dit_com = D
  !
  !     SET RIGHT HAND SIDE VALUES FROM GRHS IN USOL ON THE INTERIOR
  !     AND NON-SPECIFIED BOUNDARIES.
  !
  DO i = 2, M
    DO j = 2, N
      Usol(i,j) = Grhs(i,j)
    END DO
  END DO
  IF( kswx_com/=2 .AND. kswx_com/=3 ) THEN
    DO j = 2, N
      Usol(1,j) = Grhs(1,j)
    END DO
  END IF
  IF( kswx_com/=2 .AND. kswx_com/=5 ) THEN
    DO j = 2, N
      Usol(k_com,j) = Grhs(k_com,j)
    END DO
  END IF
  IF( kswy_com/=2 .AND. kswy_com/=3 ) THEN
    DO i = 2, M
      Usol(i,1) = Grhs(i,1)
    END DO
  END IF
  IF( kswy_com/=2 .AND. kswy_com/=5 ) THEN
    DO i = 2, M
      Usol(i,l_com) = Grhs(i,l_com)
    END DO
  END IF
  IF( kswx_com/=2 .AND. kswx_com/=3 .AND. kswy_com/=2 .AND. kswy_com/=3 ) Usol(1,1) = Grhs(1,1)
  IF( kswx_com/=2 .AND. kswx_com/=5 .AND. kswy_com/=2 .AND. kswy_com/=3 ) Usol(k_com,1) = Grhs(k_com,1)
  IF( kswx_com/=2 .AND. kswx_com/=3 .AND. kswy_com/=2 .AND. kswy_com/=5 ) Usol(1,l_com) = Grhs(1,l_com)
  IF( kswx_com/=2 .AND. kswx_com/=5 .AND. kswy_com/=2 .AND. kswy_com/=5 ) Usol(k_com,l_com) = Grhs(k_com,l_com)
  i1 = 1
  !
  !     SET SWITCHES FOR PERIODIC OR NON-PERIODIC BOUNDARIES
  !
  mp = 1
  np = 1
  IF( kswx_com==1 ) mp = 0
  IF( kswy_com==1 ) np = 0
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
    Am(i) = axi
    Bm(i) = bxi
    Cm(i) = cxi
  END DO
  !
  !     SET Y DIRECTION
  !
  DO j = 1, nit_com
    yj = cit_com + (js_com+j-2)*dly_com
    CALL COFY(yj,dj,ej,fj)
    dyj = (dj/dly_com-0.5*ej)/dly_com
    eyj = (-2.*dj/dly_com**2+fj)
    fyj = (dj/dly_com+0.5*ej)/dly_com
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
  SELECT CASE (kswy_com)
    CASE (1)
    CASE (3)
      !
      !     DIRICHLET-MIXED IN Y DIRECTION
      !
      An(1) = 0.0
      An(nit_com) = An(nit_com) + fyn
      Bn(nit_com) = Bn(nit_com) - 2.*dly_com*Xnu*fyn
      Cn(nit_com) = 0.0
    CASE (4)
      !
      !     MIXED - MIXED DIRECTION IN Y DIRECTION
      !
      An(1) = 0.0
      Bn(1) = Bn(1) + 2.*dly_com*Gama*dy1
      Cn(1) = Cn(1) + dy1
      An(nit_com) = An(nit_com) + fyn
      Bn(nit_com) = Bn(nit_com) - 2.0*dly_com*Xnu*fyn
      Cn(nit_com) = 0.0
    CASE (5)
      !
      !     MIXED-DIRICHLET IN Y DIRECTION
      !
      An(1) = 0.0
      Bn(1) = Bn(1) + 2.*dly_com*Gama*dy1
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
  CALL CHKSNG(Mbdcnd,Nbdcnd,Alpha,Beta,Gama,Xnu,COFX,COFY,singlr)
  !
  !     COMPUTE NON-ZERO EIGENVECTOR IN NULL SPACE OF TRANSPOSE
  !     IF SINGULAR
  !
  IF( singlr ) CALL TRISP(mit_com,Am,Bm,Cm,Dm,Um,Zm)
  IF( singlr ) CALL TRISP(nit_com,An,Bn,Cn,Dn,Un,Zn)
  !
  !     MAKE INITIALIZATION CALL TO BLKTRI
  !
  IF( Intl==0 ) CALL BLKTRI(Intl,np,nit_com,An,Bn,Cn,mp,mit_com,Am,Bm,Cm,Idmn,&
    Usol(is_com,js_com),Ierror,W)
  IF( Ierror/=0 ) RETURN
  !
  !     ADJUST RIGHT HAND SIDE IF NECESSARY
  !
  100 CONTINUE
  IF( singlr ) CALL ORTHOG(Usol,Idmn,Zn,Zm,Pertrb)
  !
  !     COMPUTE SOLUTION
  !
  CALL BLKTRI(i1,np,nit_com,An,Bn,Cn,mp,mit_com,Am,Bm,Cm,Idmn,Usol(is_com,js_com),Ierror,W)
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
  IF( singlr ) CALL MINSOL(Usol,Idmn,Zn,Zm)
  !
  !     RETURN IF DEFERRED CORRECTIONS AND A FOURTH ORDER SOLUTION ARE
  !     NOT FLAGGED
  !
  IF( iord==2 ) RETURN
  iord = 2
  !
  !     COMPUTE NEW RIGHT HAND SIDE FOR FOURTH ORDER SOLUTION
  !
  CALL DEFER(COFX,COFY,Idmn,Usol,Grhs)
  GOTO 100
END SUBROUTINE SPELIP
