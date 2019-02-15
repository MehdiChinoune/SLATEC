!DECK DXQNU
SUBROUTINE DXQNU(Nu1,Nu2,Mu1,Theta,X,Sx,Id,Pqa,Ipqa,Ierror)
  IMPLICIT NONE
  INTEGER Id, Ierror, ipq, ipq1, ipq2, Ipqa, ipql1, ipql2, k, mu, &
    Mu1
  !***BEGIN PROLOGUE  DXQNU
  !***SUBSIDIARY
  !***PURPOSE  To compute the values of Legendre functions for DXLEGF.
  !            Method: backward nu-wise recurrence for Q(MU,NU,X) for
  !            fixed mu to obtain Q(MU1,NU1,X), Q(MU1,NU1+1,X), ...,
  !            Q(MU1,NU2,X).
  !***LIBRARY   SLATEC
  !***CATEGORY  C3A2, C9
  !***TYPE      DOUBLE PRECISION (XQNU-S, DXQNU-D)
  !***KEYWORDS  LEGENDRE FUNCTIONS
  !***AUTHOR  Smith, John M., (NBS and George Mason University)
  !***ROUTINES CALLED  DXADD, DXADJ, DXPQNU
  !***REVISION HISTORY  (YYMMDD)
  !   820728  DATE WRITTEN
  !   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
  !   901019  Revisions to prologue.  (DWL and WRB)
  !   901106  Corrected order of sections in prologue and added TYPE
  !           section.  (WRB)
  !   920127  Revised PURPOSE section of prologue.  (DWL)
  !***END PROLOGUE  DXQNU
  DIMENSION Pqa(*), Ipqa(*)
  REAL(8) :: dmu, nu, Nu1, Nu2, pq, Pqa, pq1, pq2, Sx, X, &
    x1, x2
  REAL(8) :: Theta, pql1, pql2
  !***FIRST EXECUTABLE STATEMENT  DXQNU
  Ierror = 0
  k = 0
  pq2 = 0.0D0
  ipq2 = 0
  pql2 = 0.0D0
  ipql2 = 0
  IF ( Mu1/=1 ) THEN
    mu = 0
    !
    !        CALL DXPQNU TO OBTAIN Q(0.,NU2,X) AND Q(0.,NU2-1,X)
    !
    CALL DXPQNU(Nu1,Nu2,mu,Theta,Id,Pqa,Ipqa,Ierror)
    IF ( Ierror/=0 ) RETURN
    IF ( Mu1==0 ) RETURN
    k = (Nu2-Nu1+1.5D0)
    pq2 = Pqa(k)
    ipq2 = Ipqa(k)
    pql2 = Pqa(k-1)
    ipql2 = Ipqa(k-1)
  ENDIF
  mu = 1
  !
  !        CALL DXPQNU TO OBTAIN Q(1.,NU2,X) AND Q(1.,NU2-1,X)
  !
  CALL DXPQNU(Nu1,Nu2,mu,Theta,Id,Pqa,Ipqa,Ierror)
  IF ( Ierror/=0 ) RETURN
  IF ( Mu1==1 ) RETURN
  nu = Nu2
  pq1 = Pqa(k)
  ipq1 = Ipqa(k)
  pql1 = Pqa(k-1)
  ipql1 = Ipqa(k-1)
  100  mu = 1
  dmu = 1.D0
  DO
    !
    !        FORWARD RECURRENCE IN MU TO OBTAIN Q(MU1,NU2,X) AND
    !              Q(MU1,NU2-1,X) USING
    !              Q(MU+1,NU,X)=-2.*MU*X*SQRT(1./(1.-X**2))*Q(MU,NU,X)
    !                   -(NU+MU)*(NU-MU+1.)*Q(MU-1,NU,X)
    !
    !              FIRST FOR NU=NU2
    !
    x1 = -2.D0*dmu*X*Sx*pq1
    x2 = (nu+dmu)*(nu-dmu+1.D0)*pq2
    CALL DXADD(x1,ipq1,-x2,ipq2,pq,ipq,Ierror)
    IF ( Ierror/=0 ) RETURN
    CALL DXADJ(pq,ipq,Ierror)
    IF ( Ierror/=0 ) RETURN
    pq2 = pq1
    ipq2 = ipq1
    pq1 = pq
    ipq1 = ipq
    mu = mu + 1
    dmu = dmu + 1.D0
    IF ( mu>=Mu1 ) THEN
      Pqa(k) = pq
      Ipqa(k) = ipq
      IF ( k==1 ) RETURN
      IF ( nu<Nu2 ) THEN
        !
        !         BACKWARD RECURRENCE IN NU TO OBTAIN
        !              Q(MU1,NU1,X),Q(MU1,NU1+1,X),....,Q(MU1,NU2,X)
        !              USING
        !              (NU-MU+1.)*Q(MU,NU+1,X)=
        !                       (2.*NU+1.)*X*Q(MU,NU,X)-(NU+MU)*Q(MU,NU-1,X)
        !
        pq1 = Pqa(k)
        ipq1 = Ipqa(k)
        pq2 = Pqa(k+1)
        ipq2 = Ipqa(k+1)
        DO
          IF ( nu<=Nu1 ) RETURN
          k = k - 1
          x1 = (2.D0*nu+1.D0)*X*pq1/(nu+dmu)
          x2 = -(nu-dmu+1.D0)*pq2/(nu+dmu)
          CALL DXADD(x1,ipq1,x2,ipq2,pq,ipq,Ierror)
          IF ( Ierror/=0 ) RETURN
          CALL DXADJ(pq,ipq,Ierror)
          IF ( Ierror/=0 ) RETURN
          pq2 = pq1
          ipq2 = ipq1
          pq1 = pq
          ipq1 = ipq
          Pqa(k) = pq
          Ipqa(k) = ipq
          nu = nu - 1.D0
        ENDDO
      ELSE
        !
        !              THEN FOR NU=NU2-1
        !
        nu = nu - 1.D0
        pq2 = pql2
        ipq2 = ipql2
        pq1 = pql1
        ipq1 = ipql1
        k = k - 1
        GOTO 100
      ENDIF
    ENDIF
  ENDDO
END SUBROUTINE DXQNU
