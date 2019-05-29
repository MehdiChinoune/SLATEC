!** HWSSS1
SUBROUTINE HWSSS1(Ts,Tf,M,Mbdcnd,Bdts,Bdtf,Ps,Pf,N,Nbdcnd,Bdps,Bdpf,&
    Elmbda,F,Idimf,Pertrb,Am,Bm,Cm,Sn,Ss,Sint,D)
  !>
  !  Subsidiary to HWSSSP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (HWSSS1-S)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **See also:**  HWSSSP
  !***
  ! **Routines called:**  GENBUN

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891009  Removed unreferenced variables.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)

  INTEGER :: Idimf, M, Mbdcnd, N, Nbdcnd
  REAL :: Elmbda, Pertrb, Pf, Ps, Tf, Ts
  REAL :: Am(M+1), Bdpf(M+1), Bdps(M+1), Bdtf(N+1), Bdts(N+1), Bm(M+1), Cm(M+1), &
    D(:), F(Idimf,N+1), Sint(M+1), Sn(M+1), Ss(M+1)
  INTEGER :: i, ierror, ii, iid, inp, ising, isp, itf, itfm, its, itsp, j, jpf, &
    jpfm, jps, jpsp, mbr, mp1, munk, nbr, np1, nunk
  REAL :: at, cf, cnp, cp, csp, ct, den, dfn, dfs, dnn, dns, dphi, dphi2, dsn, &
    dss, dth, dth2, fim1, fjj, fm, fn, hdth, hld, hne, rtn, rts, summ, sum1, sum2, &
    t1, tdp, tdt, theta, wp, wpf, wps, wtf, wts, yhld
  !* FIRST EXECUTABLE STATEMENT  HWSSS1
  mp1 = M + 1
  np1 = N + 1
  fn = N
  fm = M
  dth = (Tf-Ts)/fm
  hdth = dth/2.
  tdt = dth + dth
  dphi = (Pf-Ps)/fn
  tdp = dphi + dphi
  dphi2 = dphi*dphi
  dth2 = dth*dth
  cp = 4./(fn*dth2)
  wp = fn*SIN(hdth)/4.
  DO i = 1, mp1
    fim1 = i - 1
    theta = fim1*dth + Ts
    Sint(i) = SIN(theta)
    IF ( Sint(i)/=0 ) THEN
      t1 = 1./(dth2*Sint(i))
      Am(i) = t1*SIN(theta-hdth)
      Cm(i) = t1*SIN(theta+hdth)
      Bm(i) = -Am(i) - Cm(i) + Elmbda
    END IF
  END DO
  inp = 0
  isp = 0
  !
  ! BOUNDARY CONDITION AT THETA=TS
  !
  mbr = Mbdcnd + 1
  SELECT CASE (mbr)
    CASE (2,3,8)
      at = Am(2)
      its = 2
    CASE (4,5,9)
      at = Am(1)
      its = 1
      Cm(1) = Am(1) + Cm(1)
    CASE (6,7,10)
      at = Am(2)
      inp = 1
      its = 2
    CASE DEFAULT
      its = 1
  END SELECT
  !
  ! BOUNDARY CONDITION THETA=TF
  !
  SELECT CASE (mbr)
    CASE (2,5,6)
      ct = Cm(M)
      itf = M
    CASE (3,4,7)
      ct = Cm(M+1)
      Am(M+1) = Am(M+1) + Cm(M+1)
      itf = M + 1
    CASE (8,9,10)
      itf = M
      isp = 1
      ct = Cm(M)
    CASE DEFAULT
      itf = M
  END SELECT
  !
  ! COMPUTE HOMOGENEOUS SOLUTION WITH SOLUTION AT POLE EQUAL TO ONE
  !
  itsp = its + 1
  itfm = itf - 1
  wts = Sint(its+1)*Am(its+1)/Cm(its)
  wtf = Sint(itf-1)*Cm(itf-1)/Am(itf)
  munk = itf - its + 1
  IF ( isp>0 ) THEN
    D(its) = Cm(its)/Bm(its)
    DO i = itsp, M
      D(i) = Cm(i)/(Bm(i)-Am(i)*D(i-1))
    END DO
    Ss(M) = -D(M)
    iid = M - its
    DO ii = 1, iid
      i = M - ii
      Ss(i) = -D(i)*Ss(i+1)
    END DO
    Ss(M+1) = 1.
  END IF
  IF ( inp>0 ) THEN
    Sn(1) = 1.
    D(itf) = Am(itf)/Bm(itf)
    iid = itf - 2
    DO ii = 1, iid
      i = itf - ii
      D(i) = Am(i)/(Bm(i)-Cm(i)*D(i+1))
    END DO
    Sn(2) = -D(2)
    DO i = 3, itf
      Sn(i) = -D(i)*Sn(i-1)
    END DO
  END IF
  !
  ! BOUNDARY CONDITIONS AT PHI=PS
  !
  nbr = Nbdcnd + 1
  wps = 1.
  wpf = 1.
  SELECT CASE (nbr)
    CASE (2,3)
      jps = 2
    CASE (4,5)
      jps = 1
      wps = .5
    CASE DEFAULT
      jps = 1
  END SELECT
  !
  ! BOUNDARY CONDITION AT PHI=PF
  !
  SELECT CASE (nbr)
    CASE (2,5)
      jpf = N
    CASE (3,4)
      wpf = .5
      jpf = N + 1
    CASE DEFAULT
      jpf = N
  END SELECT
  jpsp = jps + 1
  jpfm = jpf - 1
  nunk = jpf - jps + 1
  fjj = jpfm - jpsp + 1
  !
  ! SCALE COEFFICIENTS FOR SUBROUTINE GENBUN
  !
  DO i = its, itf
    cf = dphi2*Sint(i)*Sint(i)
    Am(i) = cf*Am(i)
    Bm(i) = cf*Bm(i)
    Cm(i) = cf*Cm(i)
  END DO
  Am(its) = 0.
  Cm(itf) = 0.
  ising = 0
  SELECT CASE (mbr)
    CASE (2,3,5,6,8)
    CASE DEFAULT
      SELECT CASE (nbr)
        CASE (2,3,5)
        CASE DEFAULT
          IF ( Elmbda>=0 ) THEN
            ising = 1
            summ = wts*wps + wts*wpf + wtf*wps + wtf*wpf
            IF ( inp>0 ) summ = summ + wp
            IF ( isp>0 ) summ = summ + wp
            sum1 = 0.
            DO i = itsp, itfm
              sum1 = sum1 + Sint(i)
            END DO
            summ = summ + fjj*(sum1+wts+wtf)
            summ = summ + (wps+wpf)*sum1
            hne = summ
          END IF
      END SELECT
  END SELECT
  SELECT CASE (mbr)
    CASE (1)
    CASE (2,3,8)
      DO j = jps, jpf
        F(2,j) = F(2,j) - at*F(1,j)
      END DO
    CASE (4,5,9)
      DO j = jps, jpf
        F(1,j) = F(1,j) + tdt*Bdts(j)*at
      END DO
    CASE DEFAULT
      IF ( Nbdcnd==3 ) THEN
        yhld = F(1,jps) - 4./(fn*dphi*dth2)*(Bdpf(2)-Bdps(2))
        DO j = 1, np1
          F(1,j) = yhld
        END DO
      END IF
  END SELECT
  SELECT CASE (mbr)
    CASE (1)
    CASE (2,5,6)
      DO j = jps, jpf
        F(M,j) = F(M,j) - ct*F(M+1,j)
      END DO
    CASE (3,4,7)
      DO j = jps, jpf
        F(M+1,j) = F(M+1,j) - tdt*Bdtf(j)*ct
      END DO
    CASE DEFAULT
      IF ( Nbdcnd==3 ) THEN
        yhld = F(M+1,jps) - 4./(fn*dphi*dth2)*(Bdpf(M)-Bdps(M))
        DO j = 1, np1
          F(M+1,j) = yhld
        END DO
      END IF
  END SELECT
  SELECT CASE (nbr)
    CASE (1)
    CASE (4,5)
      DO i = its, itf
        F(i,1) = F(i,1) + tdp*Bdps(i)/(dphi2*Sint(i)*Sint(i))
      END DO
    CASE DEFAULT
      DO i = its, itf
        F(i,2) = F(i,2) - F(i,1)/(dphi2*Sint(i)*Sint(i))
      END DO
  END SELECT
  SELECT CASE (nbr)
    CASE (1)
    CASE (3,4)
      DO i = its, itf
        F(i,N+1) = F(i,N+1) - tdp*Bdpf(i)/(dphi2*Sint(i)*Sint(i))
      END DO
    CASE DEFAULT
      DO i = its, itf
        F(i,N) = F(i,N) - F(i,N+1)/(dphi2*Sint(i)*Sint(i))
      END DO
  END SELECT
  Pertrb = 0.
  IF ( ising/=0 ) THEN
    summ = wts*wps*F(its,jps) + wts*wpf*F(its,jpf) + wtf*wps*F(itf,jps)&
      + wtf*wpf*F(itf,jpf)
    IF ( inp>0 ) summ = summ + wp*F(1,jps)
    IF ( isp>0 ) summ = summ + wp*F(M+1,jps)
    DO i = itsp, itfm
      sum1 = 0.
      DO j = jpsp, jpfm
        sum1 = sum1 + F(i,j)
      END DO
      summ = summ + Sint(i)*sum1
    END DO
    sum1 = 0.
    sum2 = 0.
    DO j = jpsp, jpfm
      sum1 = sum1 + F(its,j)
      sum2 = sum2 + F(itf,j)
    END DO
    summ = summ + wts*sum1 + wtf*sum2
    sum1 = 0.
    sum2 = 0.
    DO i = itsp, itfm
      sum1 = sum1 + Sint(i)*F(i,jps)
      sum2 = sum2 + Sint(i)*F(i,jpf)
    END DO
    summ = summ + wps*sum1 + wpf*sum2
    Pertrb = summ/hne
    DO j = 1, np1
      DO i = 1, mp1
        F(i,j) = F(i,j) - Pertrb
      END DO
    END DO
  END IF
  !
  ! SCALE RIGHT SIDE FOR SUBROUTINE GENBUN
  !
  DO i = its, itf
    cf = dphi2*Sint(i)*Sint(i)
    DO j = jps, jpf
      F(i,j) = cf*F(i,j)
    END DO
  END DO
  CALL GENBUN(Nbdcnd,nunk,1,munk,Am(its),Bm(its),Cm(its),Idimf,F(its,jps),&
    ierror,D)
  IF ( ising>0 ) THEN
    IF ( inp<=0 ) THEN
      IF ( isp>0 ) THEN
        DO j = 1, np1
          F(M+1,j) = 0.
        END DO
        GOTO 100
      END IF
    ELSEIF ( isp<=0 ) THEN
      DO j = 1, np1
        F(1,j) = 0.
      END DO
      GOTO 100
    END IF
  END IF
  IF ( inp>0 ) THEN
    summ = wps*F(its,jps) + wpf*F(its,jpf)
    DO j = jpsp, jpfm
      summ = summ + F(its,j)
    END DO
    dfn = cp*summ
    dnn = cp*((wps+wpf+fjj)*(Sn(2)-1.)) + Elmbda
    dsn = cp*(wps+wpf+fjj)*Sn(M)
    IF ( isp<=0 ) THEN
      cnp = (F(1,1)-dfn)/dnn
      DO i = its, itf
        hld = cnp*Sn(i)
        DO j = jps, jpf
          F(i,j) = F(i,j) + hld
        END DO
      END DO
      DO j = 1, np1
        F(1,j) = cnp
      END DO
      GOTO 100
    END IF
  ELSEIF ( isp<=0 ) THEN
    GOTO 100
  END IF
  summ = wps*F(itf,jps) + wpf*F(itf,jpf)
  DO j = jpsp, jpfm
    summ = summ + F(itf,j)
  END DO
  dfs = cp*summ
  dss = cp*((wps+wpf+fjj)*(Ss(M)-1.)) + Elmbda
  dns = cp*(wps+wpf+fjj)*Ss(2)
  IF ( inp<=0 ) THEN
    csp = (F(M+1,1)-dfs)/dss
    DO i = its, itf
      hld = csp*Ss(i)
      DO j = jps, jpf
        F(i,j) = F(i,j) + hld
      END DO
    END DO
    DO j = 1, np1
      F(M+1,j) = csp
    END DO
  ELSE
    rtn = F(1,1) - dfn
    rts = F(M+1,1) - dfs
    IF ( ising>0 ) THEN
      csp = 0.
      cnp = rtn/dnn
    ELSEIF ( ABS(dnn)<=ABS(dsn) ) THEN
      den = dns - dss*dnn/dsn
      rtn = rtn - rts*dnn/dsn
      csp = rtn/den
      cnp = (rts-dss*csp)/dsn
    ELSE
      den = dss - dns*dsn/dnn
      rts = rts - rtn*dsn/dnn
      csp = rts/den
      cnp = (rtn-csp*dns)/dnn
    END IF
    DO i = its, itf
      hld = cnp*Sn(i) + csp*Ss(i)
      DO j = jps, jpf
        F(i,j) = F(i,j) + hld
      END DO
    END DO
    DO j = 1, np1
      F(1,j) = cnp
      F(M+1,j) = csp
    END DO
  END IF
  100 CONTINUE
  IF ( Nbdcnd==0 ) THEN
    DO i = 1, mp1
      F(i,jpf+1) = F(i,jps)
    END DO
  END IF
END SUBROUTINE HWSSS1
