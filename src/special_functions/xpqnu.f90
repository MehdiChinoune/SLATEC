!** XPQNU
SUBROUTINE XPQNU(Nu1,Nu2,Mu,Theta,Id,Pqa,Ipqa,Ierror)
  !>
  !  To compute the values of Legendre functions for XLEGF.
  !            This subroutine calculates initial values of P or Q using
  !            power series, then performs forward nu-wise recurrence to
  !            obtain P(-MU,NU,X), Q(0,NU,X), or Q(1,NU,X). The nu-wise
  !            recurrence is stable for P for all mu and for Q for mu=0,1.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  C3A2, C9
  !***
  ! **Type:**      SINGLE PRECISION (XPQNU-S, DXPQNU-D)
  !***
  ! **Keywords:**  LEGENDRE FUNCTIONS
  !***
  ! **Author:**  Smith, John M., (NBS and George Mason University)
  !***
  ! **Routines called:**  XADD, XADJ, XPSI
  !***
  ! COMMON BLOCKS    XBLK1

  !* REVISION HISTORY  (YYMMDD)
  !   820728  DATE WRITTEN
  !   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
  !   901019  Revisions to prologue.  (DWL and WRB)
  !   901106  Changed all specific intrinsics to generic.  (WRB)
  !           Corrected order of sections in prologue and added TYPE
  !           section.  (WRB)
  !   920127  Revised PURPOSE section of prologue.  (DWL)
  USE XBLK ,ONLY: nbitsf_com
  INTEGER :: Id, Ierror, Mu, Ipqa(*)
  REAL :: Nu1, Nu2, Theta, Pqa(*)
  INTEGER :: i, ia, if, ipq, ipq1, ipq2, ipsik, ipsix, ix1, ixs, j, j0, k
  REAL :: a, nu, pq, r, w, x, x1, x2, xs, y, z, di, dmu, pq1, pq2, factmu, flok
  !
  !        J0, IPSIK, AND IPSIX ARE INITIALIZED IN THIS SUBROUTINE.
  !        J0 IS THE NUMBER OF TERMS USED IN SERIES EXPANSION
  !        IN SUBROUTINE XPQNU.
  !        IPSIK, IPSIX ARE VALUES OF K AND X RESPECTIVELY
  !        USED IN THE CALCULATION OF THE XPSI FUNCTION.
  !
  !* FIRST EXECUTABLE STATEMENT  XPQNU
  Ierror = 0
  j0 = nbitsf_com
  ipsik = 1 + (nbitsf_com/10)
  ipsix = 5*ipsik
  ipq = 0
  !        FIND NU IN INTERVAL [-.5,.5) IF ID=2  ( CALCULATION OF Q )
  nu = MOD(Nu1,1.)
  IF ( nu>=.5 ) nu = nu - 1.
  !        FIND NU IN INTERVAL (-1.5,-.5] IF ID=1,3, OR 4  ( CALC. OF P )
  IF ( Id/=2.AND.nu>-.5 ) nu = nu - 1.
  !        CALCULATE MU FACTORIAL
  k = Mu
  dmu = Mu
  IF ( Mu>0 ) THEN
    factmu = 1.
    if = 0
    DO i = 1, k
      factmu = factmu*i
      CALL XADJ(factmu,if,Ierror)
    END DO
    IF ( Ierror/=0 ) RETURN
  END IF
  IF ( k==0 ) factmu = 1.
  IF ( k==0 ) if = 0
  !
  !        X=COS(THETA)
  !        Y=SIN(THETA/2)**2=(1-X)/2=.5-.5*X
  !        R=TAN(THETA/2)=SQRT((1-X)/(1+X)
  !
  x = COS(Theta)
  y = SIN(Theta/2.)**2
  r = TAN(Theta/2.)
  !
  !        USE ASCENDING SERIES TO CALCULATE TWO VALUES OF P OR Q
  !        FOR USE AS STARTING VALUES IN RECURRENCE RELATION.
  !
  pq2 = 0.0
  DO j = 1, 2
    ipq1 = 0
    IF ( Id==2 ) THEN
      !
      !        Z=-LN(R)=.5*LN((1+X)/(1-X))
      !
      z = -LOG(r)
      w = XPSI(nu+1.,ipsik,ipsix)
      xs = 1./SIN(Theta)
      !
      !        SERIES SUMMATION FOR Q ( ID = 2 )
      !        Q(0,NU,X)=SUM(FROM 0 TO J0-1)((.5*LN((1+X)/(1-X))
      !    +XPSI(J+1,IPSIK,IPSIX)-XPSI(NU+1,IPSIK,IPSIX)))*A(J)*(.5-.5*X)**J
      !
      !        Q(1,NU,X)=-SQRT(1./(1.-X**2))+SQRT((1-X)/(1+X))
      !             *SUM(FROM 0 T0 J0-1)(-NU*(NU+1)/2*LN((1+X)/(1-X))
      !                 +(J-NU)*(J+NU+1)/(2*(J+1))+NU*(NU+1)*
      !     (XPSI(NU+1,IPSIK,IPSIX)-XPSI(J+1,IPSIK,IPSIX))*A(J)*(.5-.5*X)**J
      !
      !        NOTE, IN THIS LOOP K=J+1
      !
      pq = 0.
      ipq = 0
      ia = 0
      a = 1.
      DO k = 1, j0
        flok = k
        IF ( k/=1 ) THEN
          a = a*y*(flok-2.-nu)*(flok-1.+nu)/((flok-1.+dmu)*(flok-1.))
          CALL XADJ(a,ia,Ierror)
          IF ( Ierror/=0 ) RETURN
        END IF
        IF ( Mu>=1 ) THEN
          x1 = (nu*(nu+1.)*(z-w+XPSI(flok,ipsik,ipsix))+(nu-flok+1.)&
            *(nu+flok)/(2.*flok))*a
          ix1 = ia
          CALL XADD(pq,ipq,x1,ix1,pq,ipq,Ierror)
          IF ( Ierror/=0 ) RETURN
        ELSE
          x1 = (XPSI(flok,ipsik,ipsix)-w+z)*a
          ix1 = ia
          CALL XADD(pq,ipq,x1,ix1,pq,ipq,Ierror)
          IF ( Ierror/=0 ) RETURN
        END IF
      END DO
      IF ( Mu>=1 ) pq = -r*pq
      ixs = 0
      IF ( Mu>=1 ) CALL XADD(pq,ipq,-xs,ixs,pq,ipq,Ierror)
      IF ( Ierror/=0 ) RETURN
      IF ( j==2 ) Mu = -Mu
      IF ( j==2 ) dmu = -dmu
    ELSE
      !
      !        SERIES FOR P ( ID = 1, 3, OR 4 )
      !        P(-MU,NU,X)=1./FACTORIAL(MU)*SQRT(((1.-X)/(1.+X))**MU)
      !                *SUM(FROM 0 TO J0-1)A(J)*(.5-.5*X)**J
      !
      ipq = 0
      pq = 1.
      a = 1.
      ia = 0
      DO i = 2, j0
        di = i
        a = a*y*(di-2.-nu)*(di-1.+nu)/((di-1.+dmu)*(di-1.))
        CALL XADJ(a,ia,Ierror)
        IF ( Ierror/=0 ) RETURN
        IF ( a==0. ) EXIT
        CALL XADD(pq,ipq,a,ia,pq,ipq,Ierror)
        IF ( Ierror/=0 ) RETURN
      END DO
      IF ( Mu>0 ) THEN
        x2 = r
        x1 = pq
        k = Mu
        DO i = 1, k
          x1 = x1*x2
          CALL XADJ(x1,ipq,Ierror)
        END DO
        IF ( Ierror/=0 ) RETURN
        pq = x1/factmu
        ipq = ipq - if
        CALL XADJ(pq,ipq,Ierror)
        IF ( Ierror/=0 ) RETURN
      END IF
    END IF
    IF ( j==1 ) pq2 = pq
    IF ( j==1 ) ipq2 = ipq
    nu = nu + 1.
  END DO
  k = 0
  IF ( nu-1.5>=Nu1 ) THEN
    k = k + 1
    Pqa(k) = pq2
    Ipqa(k) = ipq2
    IF ( nu>Nu2+.5 ) RETURN
  END IF
  100  pq1 = pq
  ipq1 = ipq
  IF ( nu>=Nu1+.5 ) THEN
    k = k + 1
    Pqa(k) = pq
    Ipqa(k) = ipq
    IF ( nu>Nu2+.5 ) RETURN
  END IF
  !
  !        FORWARD NU-WISE RECURRENCE FOR F(MU,NU,X) FOR FIXED MU
  !        USING
  !        (NU+MU+1)*F(MU,NU,X)=(2.*NU+1)*F(MU,NU,X)-(NU-MU)*F(MU,NU-1,X)
  !        WHERE F(MU,NU,X) MAY BE P(-MU,NU,X) OR IF MU IS REPLACED
  !        BY -MU THEN F(MU,NU,X) MAY BE Q(MU,NU,X).
  !        NOTE, IN THIS LOOP, NU=NU+1
  !
  x1 = (2.*nu-1.)/(nu+dmu)*x*pq1
  x2 = (nu-1.-dmu)/(nu+dmu)*pq2
  CALL XADD(x1,ipq1,-x2,ipq2,pq,ipq,Ierror)
  IF ( Ierror/=0 ) RETURN
  CALL XADJ(pq,ipq,Ierror)
  IF ( Ierror/=0 ) RETURN
  nu = nu + 1.
  pq2 = pq1
  ipq2 = ipq1
  GOTO 100
  !
END SUBROUTINE XPQNU
