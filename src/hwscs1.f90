!*==HWSCS1.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK HWSCS1
SUBROUTINE HWSCS1(Intl,Ts,Tf,M,Mbdcnd,Bdts,Bdtf,Rs,Rf,N,Nbdcnd,Bdrs,Bdrf,&
    Elmbda,F,Idimf,Pertrb,W,S,An,Bn,Cn,R,Am,Bm,Cm,Sint,Bmh)
  IMPLICIT NONE
  !*--HWSCS16
  !*** Start of declarations inserted by SPAG
  REAL Am , An , ar , at , Bdrf , Bdrs , Bdtf , Bdts , Bm , Bmh , Bn , Cm , &
    Cn , cr , ct , czr , dr , dr2 , dth , Elmbda
  REAL F , hdr , hdth , hne , Pertrb , R , r2 , Rf , rf2 , Rs , rs2 , rsq , &
    S , sdts , Sint , sum , t1 , tdr , tdt , Tf
  REAL theta , Ts , W , wrf , wrs , wrz , wtf , wtnm , wts , xp , xps , &
    yhld , yph , yps
  INTEGER i , ictr , Idimf , ierror , iflg , Intl , ising , itf , itfm , &
    its , itsp , j , jrf , jrfm , jrs , jrsp , l , M , Mbdcnd , mp
  INTEGER mp1 , munk , N , Nbdcnd , np , np1 , nunk
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  HWSCS1
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to HWSCSP
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (HWSCS1-S)
  !***AUTHOR  (UNKNOWN)
  !***SEE ALSO  HWSCSP
  !***ROUTINES CALLED  BLKTRI
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891009  Removed unreferenced variables.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  HWSCS1
  DIMENSION F(Idimf,*) , Bdrs(*) , Bdrf(*) , Bdts(*) , Bdtf(*) , Am(*) , &
    Bm(*) , Cm(*) , An(*) , Bn(*) , Cn(*) , S(*) , R(*) , Sint(*) , &
    Bmh(*) , W(*)
  !***FIRST EXECUTABLE STATEMENT  HWSCS1
  mp1 = M + 1
  dth = (Tf-Ts)/M
  tdt = dth + dth
  hdth = dth/2.
  sdts = 1./(dth*dth)
  DO i = 1 , mp1
    theta = Ts + (i-1)*dth
    Sint(i) = SIN(theta)
    IF ( Sint(i)/=0 ) THEN
      t1 = sdts/Sint(i)
      Am(i) = t1*SIN(theta-hdth)
      Cm(i) = t1*SIN(theta+hdth)
      Bm(i) = -(Am(i)+Cm(i))
    ENDIF
  ENDDO
  np1 = N + 1
  dr = (Rf-Rs)/N
  hdr = dr/2.
  tdr = dr + dr
  dr2 = dr*dr
  czr = 6.*dth/(dr2*(COS(Ts)-COS(Tf)))
  DO j = 1 , np1
    R(j) = Rs + (j-1)*dr
    An(j) = (R(j)-hdr)**2/dr2
    Cn(j) = (R(j)+hdr)**2/dr2
    Bn(j) = -(An(j)+Cn(j))
  ENDDO
  mp = 1
  np = 1
  !
  ! BOUNDARY CONDITION AT PHI=PS
  !
  SELECT CASE (Mbdcnd)
    CASE (3,4,8)
      at = Am(1)
      its = 1
      Cm(1) = Cm(1) + Am(1)
    CASE (5,6,9)
      its = 1
      Bm(1) = -4.*sdts
      Cm(1) = -Bm(1)
    CASE DEFAULT
      at = Am(2)
      its = 2
  END SELECT
  !
  ! BOUNDARY CONDITION AT PHI=PF
  !
  SELECT CASE (Mbdcnd)
    CASE (2,3,6)
      ct = Cm(M+1)
      Am(M+1) = Am(M+1) + Cm(M+1)
      itf = M + 1
    CASE (7,8,9)
      itf = M + 1
      Am(M+1) = 4.*sdts
      Bm(M+1) = -Am(M+1)
    CASE DEFAULT
      ct = Cm(M)
      itf = M
  END SELECT
  wts = Sint(its+1)*Am(its+1)/Cm(its)
  wtf = Sint(itf-1)*Cm(itf-1)/Am(itf)
  itsp = its + 1
  itfm = itf - 1
  !
  ! BOUNDARY CONDITION AT R=RS
  !
  ictr = 0
  SELECT CASE (Nbdcnd)
    CASE (3,4)
      ar = An(1)
      jrs = 1
      Cn(1) = Cn(1) + An(1)
    CASE (5,6)
      jrs = 2
      ictr = 1
      S(N) = An(N)/Bn(N)
      DO j = 3 , N
        l = N - j + 2
        S(l) = An(l)/(Bn(l)-Cn(l)*S(l+1))
      ENDDO
      S(2) = -S(2)
      DO j = 3 , N
        S(j) = -S(j)*S(j-1)
      ENDDO
      wtnm = wts + wtf
      DO i = itsp , itfm
        wtnm = wtnm + Sint(i)
      ENDDO
      yps = czr*wtnm*(S(2)-1.)
    CASE DEFAULT
      ar = An(2)
      jrs = 2
  END SELECT
  !
  ! BOUNDARY CONDITION AT R=RF
  !
  SELECT CASE (Nbdcnd)
    CASE (2,3,6)
      cr = Cn(N+1)
      An(N+1) = An(N+1) + Cn(N+1)
      jrf = N + 1
    CASE DEFAULT
      cr = Cn(N)
      jrf = N
  END SELECT
  wrs = An(jrs+1)*R(jrs)**2/Cn(jrs)
  wrf = Cn(jrf-1)*R(jrf)**2/An(jrf)
  wrz = An(jrs)/czr
  jrsp = jrs + 1
  jrfm = jrf - 1
  munk = itf - its + 1
  nunk = jrf - jrs + 1
  DO i = its , itf
    Bmh(i) = Bm(i)
  ENDDO
  ising = 0
  SELECT CASE (Nbdcnd)
    CASE (1,2,4,5)
    CASE DEFAULT
      SELECT CASE (Mbdcnd)
        CASE (1,2,4,5,7)
        CASE DEFAULT
          IF ( Elmbda>=0 ) THEN
            ising = 1
            sum = wts*wrs + wts*wrf + wtf*wrs + wtf*wrf
            IF ( ictr/=0 ) sum = sum + wrz
            DO j = jrsp , jrfm
              r2 = R(j)**2
              DO i = itsp , itfm
                sum = sum + r2*Sint(i)
              ENDDO
            ENDDO
            DO j = jrsp , jrfm
              sum = sum + (wts+wtf)*R(j)**2
            ENDDO
            DO i = itsp , itfm
              sum = sum + (wrs+wrf)*Sint(i)
            ENDDO
            hne = sum
          ENDIF
      END SELECT
  END SELECT
  SELECT CASE (Mbdcnd)
    CASE (5,6,9)
    CASE DEFAULT
      Bm(its) = Bmh(its) + Elmbda/Sint(its)**2
  END SELECT
  SELECT CASE (Mbdcnd)
    CASE (7,8,9)
    CASE DEFAULT
      Bm(itf) = Bmh(itf) + Elmbda/Sint(itf)**2
  END SELECT
  DO i = itsp , itfm
    Bm(i) = Bmh(i) + Elmbda/Sint(i)**2
  ENDDO
  SELECT CASE (Mbdcnd)
    CASE (3,4,8)
      DO j = jrs , jrf
        F(1,j) = F(1,j) + tdt*Bdts(j)*at/R(j)**2
      ENDDO
    CASE (5,6,9)
    CASE DEFAULT
      DO j = jrs , jrf
        F(2,j) = F(2,j) - at*F(1,j)/R(j)**2
      ENDDO
  END SELECT
  SELECT CASE (Mbdcnd)
    CASE (2,3,6)
      DO j = jrs , jrf
        F(M+1,j) = F(M+1,j) - tdt*Bdtf(j)*ct/R(j)**2
      ENDDO
    CASE (7,8,9)
    CASE DEFAULT
      DO j = jrs , jrf
        F(M,j) = F(M,j) - ct*F(M+1,j)/R(j)**2
      ENDDO
  END SELECT
  SELECT CASE (Nbdcnd)
    CASE (1,2)
      rs2 = (Rs+dr)**2
      DO i = its , itf
        F(i,2) = F(i,2) - ar*F(i,1)/rs2
      ENDDO
    CASE (3,4)
      DO i = its , itf
        F(i,1) = F(i,1) + tdr*Bdrs(i)*ar/Rs**2
      ENDDO
    CASE DEFAULT
      IF ( Mbdcnd==3 ) THEN
        yhld = F(its,1) - czr/tdt*(SIN(Tf)*Bdtf(2)-SIN(Ts)*Bdts(2))
        DO i = 1 , mp1
          F(i,1) = yhld
        ENDDO
      ENDIF
  END SELECT
  SELECT CASE (Nbdcnd)
    CASE (2,3,6)
      DO i = its , itf
        F(i,N+1) = F(i,N+1) - tdr*Bdrf(i)*cr/Rf**2
      ENDDO
    CASE DEFAULT
      rf2 = (Rf-dr)**2
      DO i = its , itf
        F(i,N) = F(i,N) - cr*F(i,N+1)/rf2
      ENDDO
  END SELECT
  Pertrb = 0.
  IF ( ising/=0 ) THEN
    sum = wts*wrs*F(its,jrs) + wts*wrf*F(its,jrf) + wtf*wrs*F(itf,jrs)&
      + wtf*wrf*F(itf,jrf)
    IF ( ictr/=0 ) sum = sum + wrz*F(its,1)
    DO j = jrsp , jrfm
      r2 = R(j)**2
      DO i = itsp , itfm
        sum = sum + r2*Sint(i)*F(i,j)
      ENDDO
    ENDDO
    DO j = jrsp , jrfm
      sum = sum + R(j)**2*(wts*F(its,j)+wtf*F(itf,j))
    ENDDO
    DO i = itsp , itfm
      sum = sum + Sint(i)*(wrs*F(i,jrs)+wrf*F(i,jrf))
    ENDDO
    Pertrb = sum/hne
    DO j = 1 , np1
      DO i = 1 , mp1
        F(i,j) = F(i,j) - Pertrb
      ENDDO
    ENDDO
  ENDIF
  DO j = jrs , jrf
    rsq = R(j)**2
    DO i = its , itf
      F(i,j) = rsq*F(i,j)
    ENDDO
  ENDDO
  iflg = Intl
  DO
    CALL BLKTRI(iflg,np,nunk,An(jrs),Bn(jrs),Cn(jrs),mp,munk,Am(its),Bm(its)&
      ,Cm(its),Idimf,F(its,jrs),ierror,W)
    iflg = iflg + 1
    IF ( iflg/=1 ) THEN
      IF ( Nbdcnd==0 ) THEN
        DO i = 1 , mp1
          F(i,jrf+1) = F(i,jrs)
        ENDDO
      ENDIF
      IF ( Mbdcnd==0 ) THEN
        DO j = 1 , np1
          F(itf+1,j) = F(its,j)
        ENDDO
      ENDIF
      xp = 0.
      IF ( ictr/=0 ) THEN
        IF ( ising==0 ) THEN
          sum = wts*F(its,2) + wtf*F(itf,2)
          DO i = itsp , itfm
            sum = sum + Sint(i)*F(i,2)
          ENDDO
          yph = czr*sum
          xp = (F(its,1)-yph)/yps
          DO j = jrs , jrf
            xps = xp*S(j)
            DO i = its , itf
              F(i,j) = F(i,j) + xps
            ENDDO
          ENDDO
        ENDIF
        DO i = 1 , mp1
          F(i,1) = xp
        ENDDO
      ENDIF
      EXIT
    ENDIF
  ENDDO
END SUBROUTINE HWSCS1
