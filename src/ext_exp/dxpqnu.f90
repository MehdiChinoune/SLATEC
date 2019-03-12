!DECK DXPQNU
SUBROUTINE DXPQNU(Nu1,Nu2,Mu,Theta,Id,Pqa,Ipqa,Ierror)
  IMPLICIT NONE
  INTEGER i, ia, Id, Ierror, if, ipq, ipq1, ipq2, Ipqa, ipsik, &
    ipsix, ix1, ixs, j, j0, k, Mu, NBItsf
  !***BEGIN PROLOGUE  DXPQNU
  !***SUBSIDIARY
  !***PURPOSE  To compute the values of Legendre functions for DXLEGF.
  !            This subroutine calculates initial values of P or Q using
  !            power series, then performs forward nu-wise recurrence to
  !            obtain P(-MU,NU,X), Q(0,NU,X), or Q(1,NU,X). The nu-wise
  !            recurrence is stable for P for all mu and for Q for mu=0,1.
  !***LIBRARY   SLATEC
  !***CATEGORY  C3A2, C9
  !***TYPE      DOUBLE PRECISION (XPQNU-S, DXPQNU-D)
  !***KEYWORDS  LEGENDRE FUNCTIONS
  !***AUTHOR  Smith, John M., (NBS and George Mason University)
  !***ROUTINES CALLED  DXADD, DXADJ, DXPSI
  !***COMMON BLOCKS    DXBLK1
  !***REVISION HISTORY  (YYMMDD)
  !   820728  DATE WRITTEN
  !   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
  !   901019  Revisions to prologue.  (DWL and WRB)
  !   901106  Changed all specific intrinsics to generic.  (WRB)
  !           Corrected order of sections in prologue and added TYPE
  !           section.  (WRB)
  !   920127  Revised PURPOSE section of prologue.  (DWL)
  !***END PROLOGUE  DXPQNU
  REAL(8) :: a, nu, Nu1, Nu2, pq, Pqa, DXPSI, r, Theta, w, &
    x, x1, x2, xs, y, z
  REAL(8) :: di, dmu, pq1, pq2, factmu, flok
  DIMENSION Pqa(*), Ipqa(*)
  COMMON /DXBLK1/ NBItsf
  SAVE /DXBLK1/
  !
  !        J0, IPSIK, AND IPSIX ARE INITIALIZED IN THIS SUBROUTINE.
  !        J0 IS THE NUMBER OF TERMS USED IN SERIES EXPANSION
  !        IN SUBROUTINE DXPQNU.
  !        IPSIK, IPSIX ARE VALUES OF K AND X RESPECTIVELY
  !        USED IN THE CALCULATION OF THE DXPSI FUNCTION.
  !
  !***FIRST EXECUTABLE STATEMENT  DXPQNU
  Ierror = 0
  j0 = NBItsf
  ipsik = 1 + (NBItsf/10)
  ipsix = 5*ipsik
  ipq = 0
  !        FIND NU IN INTERVAL [-.5,.5) IF ID=2  ( CALCULATION OF Q )
  nu = MOD(Nu1,1.D0)
  IF ( nu>=.5D0 ) nu = nu - 1.D0
  !        FIND NU IN INTERVAL (-1.5,-.5] IF ID=1,3, OR 4  ( CALC. OF P )
  IF ( Id/=2.AND.nu>-.5D0 ) nu = nu - 1.D0
  !        CALCULATE MU FACTORIAL
  k = Mu
  dmu = Mu
  IF ( Mu>0 ) THEN
    factmu = 1.D0
    if = 0
    DO i = 1, k
      factmu = factmu*i
      CALL DXADJ(factmu,if,Ierror)
    ENDDO
    IF ( Ierror/=0 ) RETURN
  ENDIF
  IF ( k==0 ) factmu = 1.D0
  IF ( k==0 ) if = 0
  !
  !        X=COS(THETA)
  !        Y=SIN(THETA/2)**2=(1-X)/2=.5-.5*X
  !        R=TAN(THETA/2)=SQRT((1-X)/(1+X)
  !
  x = COS(Theta)
  y = SIN(Theta/2.D0)**2
  r = TAN(Theta/2.D0)
  !
  !        USE ASCENDING SERIES TO CALCULATE TWO VALUES OF P OR Q
  !        FOR USE AS STARTING VALUES IN RECURRENCE RELATION.
  !
  pq2 = 0.0D0
  DO j = 1, 2
    ipq1 = 0
    IF ( Id==2 ) THEN
      !
      !        Z=-LN(R)=.5*LN((1+X)/(1-X))
      !
      z = -LOG(r)
      w = DXPSI(nu+1.D0,ipsik,ipsix)
      xs = 1.D0/SIN(Theta)
      !
      !        SERIES SUMMATION FOR Q ( ID = 2 )
      !        Q(0,NU,X)=SUM(FROM 0 TO J0-1)((.5*LN((1+X)/(1-X))
      !    +DXPSI(J+1,IPSIK,IPSIX)-DXPSI(NU+1,IPSIK,IPSIX)))*A(J)*(.5-.5*X)**J
      !
      !        Q(1,NU,X)=-SQRT(1./(1.-X**2))+SQRT((1-X)/(1+X))
      !             *SUM(FROM 0 T0 J0-1)(-NU*(NU+1)/2*LN((1+X)/(1-X))
      !                 +(J-NU)*(J+NU+1)/(2*(J+1))+NU*(NU+1)*
      !     (DXPSI(NU+1,IPSIK,IPSIX)-DXPSI(J+1,IPSIK,IPSIX))*A(J)*(.5-.5*X)**J
      !
      !        NOTE, IN THIS LOOP K=J+1
      !
      pq = 0.D0
      ipq = 0
      ia = 0
      a = 1.D0
      DO k = 1, j0
        flok = k
        IF ( k/=1 ) THEN
          a = a*y*(flok-2.D0-nu)*(flok-1.D0+nu)&
            /((flok-1.D0+dmu)*(flok-1.D0))
          CALL DXADJ(a,ia,Ierror)
          IF ( Ierror/=0 ) RETURN
        ENDIF
        IF ( Mu>=1 ) THEN
          x1 = (nu*(nu+1.D0)*(z-w+DXPSI(flok,ipsik,ipsix))+(nu-flok+1.D0)&
            *(nu+flok)/(2.D0*flok))*a
          ix1 = ia
          CALL DXADD(pq,ipq,x1,ix1,pq,ipq,Ierror)
          IF ( Ierror/=0 ) RETURN
        ELSE
          x1 = (DXPSI(flok,ipsik,ipsix)-w+z)*a
          ix1 = ia
          CALL DXADD(pq,ipq,x1,ix1,pq,ipq,Ierror)
          IF ( Ierror/=0 ) RETURN
        ENDIF
      ENDDO
      IF ( Mu>=1 ) pq = -r*pq
      ixs = 0
      IF ( Mu>=1 ) CALL DXADD(pq,ipq,-xs,ixs,pq,ipq,Ierror)
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
      pq = 1.D0
      a = 1.D0
      ia = 0
      DO i = 2, j0
        di = i
        a = a*y*(di-2.D0-nu)*(di-1.D0+nu)/((di-1.D0+dmu)*(di-1.D0))
        CALL DXADJ(a,ia,Ierror)
        IF ( Ierror/=0 ) RETURN
        IF ( a==0.D0 ) EXIT
        CALL DXADD(pq,ipq,a,ia,pq,ipq,Ierror)
        IF ( Ierror/=0 ) RETURN
      ENDDO
      IF ( Mu>0 ) THEN
        x2 = r
        x1 = pq
        k = Mu
        DO i = 1, k
          x1 = x1*x2
          CALL DXADJ(x1,ipq,Ierror)
        ENDDO
        IF ( Ierror/=0 ) RETURN
        pq = x1/factmu
        ipq = ipq - if
        CALL DXADJ(pq,ipq,Ierror)
        IF ( Ierror/=0 ) RETURN
      ENDIF
    ENDIF
    IF ( j==1 ) pq2 = pq
    IF ( j==1 ) ipq2 = ipq
    nu = nu + 1.D0
  ENDDO
  k = 0
  IF ( nu-1.5D0>=Nu1 ) THEN
    k = k + 1
    Pqa(k) = pq2
    Ipqa(k) = ipq2
    IF ( nu>Nu2+.5D0 ) RETURN
  ENDIF
  100  pq1 = pq
  ipq1 = ipq
  IF ( nu>=Nu1+.5D0 ) THEN
    k = k + 1
    Pqa(k) = pq
    Ipqa(k) = ipq
    IF ( nu>Nu2+.5D0 ) RETURN
  ENDIF
  !
  !        FORWARD NU-WISE RECURRENCE FOR F(MU,NU,X) FOR FIXED MU
  !        USING
  !        (NU+MU+1)*F(MU,NU,X)=(2.*NU+1)*F(MU,NU,X)-(NU-MU)*F(MU,NU-1,X)
  !        WHERE F(MU,NU,X) MAY BE P(-MU,NU,X) OR IF MU IS REPLACED
  !        BY -MU THEN F(MU,NU,X) MAY BE Q(MU,NU,X).
  !        NOTE, IN THIS LOOP, NU=NU+1
  !
  x1 = (2.D0*nu-1.D0)/(nu+dmu)*x*pq1
  x2 = (nu-1.D0-dmu)/(nu+dmu)*pq2
  CALL DXADD(x1,ipq1,-x2,ipq2,pq,ipq,Ierror)
  IF ( Ierror/=0 ) RETURN
  CALL DXADJ(pq,ipq,Ierror)
  IF ( Ierror/=0 ) RETURN
  nu = nu + 1.D0
  pq2 = pq1
  ipq2 = ipq1
  GOTO 100
  !
END SUBROUTINE DXPQNU
