!** DXPMU
SUBROUTINE DXPMU(Nu1,Nu2,Mu1,Mu2,Theta,X,Sx,Id,Pqa,Ipqa,Ierror)
  !> To compute the values of Legendre functions for DXLEGF.
  !            Method: backward mu-wise recurrence for P(-MU,NU,X) for
  !            fixed nu to obtain P(-MU2,NU1,X), P(-(MU2-1),NU1,X), ...,
  !            P(-MU1,NU1,X) and store in ascending mu order.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  C3A2, C9
  !***
  ! **Type:**      DOUBLE PRECISION (XPMU-S, DXPMU-D)
  !***
  ! **Keywords:**  LEGENDRE FUNCTIONS
  !***
  ! **Author:**  Smith, John M., (NBS and George Mason University)
  !***
  ! **Routines called:**  DXADD, DXADJ, DXPQNU

  !* REVISION HISTORY  (YYMMDD)
  !   820728  DATE WRITTEN
  !   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
  !   901019  Revisions to prologue.  (DWL and WRB)
  !   901106  Changed all specific intrinsics to generic.  (WRB)
  !           Corrected order of sections in prologue and added TYPE
  !           section.  (WRB)
  !   920127  Revised PURPOSE section of prologue.  (DWL)

  INTEGER :: Id, Ierror, Mu1, Mu2, Ipqa(Mu2-Mu1+1)
  REAL(DP) :: Pqa(Mu2-Mu1+1), Nu1, Nu2, X, Sx, Theta
  INTEGER :: ip0, j, mu, n
  REAL(DP) :: p0, x1, x2
  !
  !        CALL DXPQNU TO OBTAIN P(-MU2,NU,X)
  !
  !* FIRST EXECUTABLE STATEMENT  DXPMU
  Ierror = 0
  CALL DXPQNU(Nu1,Nu2,Mu2,Theta,Id,Pqa,Ipqa,Ierror)
  IF( Ierror/=0 ) RETURN
  p0 = Pqa(1)
  ip0 = Ipqa(1)
  mu = Mu2 - 1
  !
  !        CALL DXPQNU TO OBTAIN P(-MU2-1,NU,X)
  !
  CALL DXPQNU(Nu1,Nu2,mu,Theta,Id,Pqa,Ipqa,Ierror)
  IF( Ierror/=0 ) RETURN
  n = Mu2 - Mu1 + 1
  Pqa(n) = p0
  Ipqa(n) = ip0
  IF( n/=1 ) THEN
    Pqa(n-1) = Pqa(1)
    Ipqa(n-1) = Ipqa(1)
    IF( n/=2 ) THEN
      j = n - 2
      DO
        !
        !        BACKWARD RECURRENCE IN MU TO OBTAIN
        !              P(-MU2,NU1,X),P(-(MU2-1),NU1,X),....P(-MU1,NU1,X)
        !              USING
        !              (NU-MU)*(NU+MU+1.)*P(-(MU+1),NU,X)=
        !                2.*MU*X*SQRT((1./(1.-X**2))*P(-MU,NU,X)-P(-(MU-1),NU,X)
        !
        x1 = 2._DP*mu*X*Sx*Pqa(j+1)
        x2 = -(Nu1-mu)*(Nu1+mu+1._DP)*Pqa(j+2)
        CALL DXADD(x1,Ipqa(j+1),x2,Ipqa(j+2),Pqa(j),Ipqa(j),Ierror)
        IF( Ierror/=0 ) RETURN
        CALL DXADJ(Pqa(j),Ipqa(j),Ierror)
        IF( Ierror/=0 ) RETURN
        IF( j==1 ) EXIT
        j = j - 1
        mu = mu - 1
      END DO
    END IF
  END IF
END SUBROUTINE DXPMU
