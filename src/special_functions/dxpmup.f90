!** DXPMUP
SUBROUTINE DXPMUP(Nu1,Nu2,Mu1,Mu2,Pqa,Ipqa,Ierror)
  !>
  !  To compute the values of Legendre functions for DXLEGF.
  !            This subroutine transforms an array of Legendre functions
  !            of the first kind of negative order stored in array PQA
  !            into Legendre functions of the first kind of positive
  !            order stored in array PQA. The original array is destroyed.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  C3A2, C9
  !***
  ! **Type:**      DOUBLE PRECISION (XPMUP-S, DXPMUP-D)
  !***
  ! **Keywords:**  LEGENDRE FUNCTIONS
  !***
  ! **Author:**  Smith, John M., (NBS and George Mason University)
  !***
  ! **Routines called:**  DXADJ

  !* REVISION HISTORY  (YYMMDD)
  !   820728  DATE WRITTEN
  !   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
  !   901019  Revisions to prologue.  (DWL and WRB)
  !   901106  Changed all specific intrinsics to generic.  (WRB)
  !           Corrected order of sections in prologue and added TYPE
  !           section.  (WRB)
  !   920127  Revised PURPOSE section of prologue.  (DWL)

  INTEGER i, Ierror, Ipqa(*), iprod, j, k, l, mu, Mu1, Mu2, n
  REAL(8) :: dmu, nu, Nu1, Nu2, Pqa(*), prod
  !* FIRST EXECUTABLE STATEMENT  DXPMUP
  Ierror = 0
  nu = Nu1
  mu = Mu1
  dmu = mu
  n = INT(Nu2-Nu1+.1D0) + (Mu2-Mu1) + 1
  j = 1
  IF ( MOD(REAL(nu),1.)==0. ) THEN
    DO WHILE ( dmu>=nu+1.D0 )
      Pqa(j) = 0.D0
      Ipqa(j) = 0
      j = j + 1
      IF ( j>n ) RETURN
      !        INCREMENT EITHER MU OR NU AS APPROPRIATE.
      IF ( Nu2-Nu1>.5D0 ) nu = nu + 1.D0
      IF ( Mu2>Mu1 ) mu = mu + 1
    END DO
  END IF
  !
  !        TRANSFORM P(-MU,NU,X) TO P(MU,NU,X) USING
  !        P(MU,NU,X)=(NU-MU+1)*(NU-MU+2)*...*(NU+MU)*P(-MU,NU,X)*(-1)**MU
  !
  prod = 1.D0
  iprod = 0
  k = 2*mu
  IF ( k/=0 ) THEN
    DO l = 1, k
      prod = prod*(dmu-nu-l)
      CALL DXADJ(prod,iprod,Ierror)
    END DO
    IF ( Ierror/=0 ) RETURN
  END IF
  DO i = j, n
    IF ( mu/=0 ) THEN
      Pqa(i) = Pqa(i)*prod*(-1)**mu
      Ipqa(i) = Ipqa(i) + iprod
      CALL DXADJ(Pqa(i),Ipqa(i),Ierror)
      IF ( Ierror/=0 ) RETURN
    END IF
    IF ( Nu2-Nu1>.5D0 ) THEN
      prod = prod*(-dmu-nu-1.D0)/(dmu-nu-1.D0)
      CALL DXADJ(prod,iprod,Ierror)
      IF ( Ierror/=0 ) RETURN
      nu = nu + 1.D0
    ELSE
      prod = (dmu-nu)*prod*(-dmu-nu-1.D0)
      CALL DXADJ(prod,iprod,Ierror)
      IF ( Ierror/=0 ) RETURN
      mu = mu + 1
      dmu = dmu + 1.D0
    END IF
  END DO
END SUBROUTINE DXPMUP
