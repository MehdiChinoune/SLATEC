!DECK XQMU
SUBROUTINE XQMU(Nu1,Nu2,Mu1,Mu2,Theta,X,Sx,Id,Pqa,Ipqa,Ierror)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  XQMU
  !***SUBSIDIARY
  !***PURPOSE  To compute the values of Legendre functions for XLEGF.
  !            Method: forward mu-wise recurrence for Q(MU,NU,X) for fixed
  !            nu to obtain Q(MU1,NU,X), Q(MU1+1,NU,X), ..., Q(MU2,NU,X).
  !***LIBRARY   SLATEC
  !***CATEGORY  C3A2, C9
  !***TYPE      SINGLE PRECISION (XQMU-S, DXQMU-D)
  !***KEYWORDS  LEGENDRE FUNCTIONS
  !***AUTHOR  Smith, John M., (NBS and George Mason University)
  !***ROUTINES CALLED  XADD, XADJ, XPQNU
  !***REVISION HISTORY  (YYMMDD)
  !   820728  DATE WRITTEN
  !   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
  !   901019  Revisions to prologue.  (DWL and WRB)
  !   901106  Corrected order of sections in prologue and added TYPE
  !           section.  (WRB)
  !   920127  Revised PURPOSE section of prologue.  (DWL)
  !***END PROLOGUE  XQMU
  INTEGER Id, Ierror, ipq, ipq1, ipq2, Ipqa, k, mu, Mu1, Mu2
  DIMENSION Pqa(*), Ipqa(*)
  REAL dmu, nu, Nu1, Nu2, pq, Pqa, pq1, pq2, Sx, X, x1, x2
  REAL Theta
  !***FIRST EXECUTABLE STATEMENT  XQMU
  Ierror = 0
  mu = 0
  !
  !        CALL XPQNU TO OBTAIN Q(0.,NU1,X)
  !
  CALL XPQNU(Nu1,Nu2,mu,Theta,Id,Pqa,Ipqa,Ierror)
  IF ( Ierror/=0 ) RETURN
  pq2 = Pqa(1)
  ipq2 = Ipqa(1)
  mu = 1
  !
  !        CALL XPQNU TO OBTAIN Q(1.,NU1,X)
  !
  CALL XPQNU(Nu1,Nu2,mu,Theta,Id,Pqa,Ipqa,Ierror)
  IF ( Ierror/=0 ) RETURN
  nu = Nu1
  k = 0
  mu = 1
  dmu = 1.
  pq1 = Pqa(1)
  ipq1 = Ipqa(1)
  IF ( Mu1<=0 ) THEN
    k = k + 1
    Pqa(k) = pq2
    Ipqa(k) = ipq2
    IF ( Mu2<1 ) RETURN
  ENDIF
  IF ( Mu1<=1 ) THEN
    k = k + 1
    Pqa(k) = pq1
    Ipqa(k) = ipq1
    IF ( Mu2<=1 ) RETURN
  ENDIF
  DO
    !
    !        FORWARD RECURRENCE IN MU TO OBTAIN
    !                  Q(MU1,NU,X),Q(MU1+1,NU,X),....,Q(MU2,NU,X) USING
    !             Q(MU+1,NU,X)=-2.*MU*X*SQRT(1./(1.-X**2))*Q(MU,NU,X)
    !                               -(NU+MU)*(NU-MU+1.)*Q(MU-1,NU,X)
    !
    x1 = -2.*dmu*X*Sx*pq1
    x2 = (nu+dmu)*(nu-dmu+1.)*pq2
    CALL XADD(x1,ipq1,-x2,ipq2,pq,ipq,Ierror)
    IF ( Ierror/=0 ) RETURN
    CALL XADJ(pq,ipq,Ierror)
    IF ( Ierror/=0 ) RETURN
    pq2 = pq1
    ipq2 = ipq1
    pq1 = pq
    ipq1 = ipq
    mu = mu + 1
    dmu = dmu + 1.
    IF ( mu>=Mu1 ) THEN
      k = k + 1
      Pqa(k) = pq
      Ipqa(k) = ipq
      IF ( Mu2<=mu ) EXIT
    ENDIF
  ENDDO
  RETURN
END SUBROUTINE XQMU