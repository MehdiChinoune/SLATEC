!** DXPNRM
SUBROUTINE DXPNRM(Nu1,Nu2,Mu1,Mu2,Pqa,Ipqa,Ierror)
  !> To compute the values of Legendre functions for DXLEGF.
  !            This subroutine transforms an array of Legendre functions
  !            of the first kind of negative order stored in array PQA
  !            into normalized Legendre polynomials stored in array PQA.
  !            The original array is destroyed.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  C3A2, C9
  !***
  ! **Type:**      DOUBLE PRECISION (XPNRM-S, DXPNRM-D)
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

  INTEGER :: i, Ierror, Ipqa(*), iprod, j, k, l, mu, Mu1, Mu2
  REAL(DP) :: c1, dmu, nu, Nu1, Nu2, Pqa(*), prod
  !* FIRST EXECUTABLE STATEMENT  DXPNRM
  Ierror = 0
  l = INT( (Mu2-Mu1) + (Nu2-Nu1+1.5_DP) )
  mu = Mu1
  dmu = Mu1
  nu = Nu1
  !
  !         IF MU >NU, NORM P =0.
  !
  j = 1
  DO WHILE( dmu>nu )
    Pqa(j) = 0._DP
    Ipqa(j) = 0
    j = j + 1
    IF( j>l ) RETURN
    !
    !        INCREMENT EITHER MU OR NU AS APPROPRIATE.
    !
    IF( Mu2>Mu1 ) dmu = dmu + 1._DP
    IF( Nu2-Nu1>.5_DP ) nu = nu + 1._DP
  END DO
  !
  !         TRANSFORM P(-MU,NU,X) INTO NORMALIZED P(MU,NU,X) USING
  !              NORM P(MU,NU,X)=
  !                 SQRT((NU+.5)*FACTORIAL(NU+MU)/FACTORIAL(NU-MU))
  !                              *P(-MU,NU,X)
  !
  prod = 1._DP
  iprod = 0
  k = 2*mu
  IF( k>0 ) THEN
    DO i = 1, k
      prod = prod*SQRT(nu+dmu+1._DP-i)
      CALL DXADJ(prod,iprod,Ierror)
    END DO
    IF( Ierror/=0 ) RETURN
  END IF
  DO i = j, l
    c1 = prod*SQRT(nu+.5_DP)
    Pqa(i) = Pqa(i)*c1
    Ipqa(i) = Ipqa(i) + iprod
    CALL DXADJ(Pqa(i),Ipqa(i),Ierror)
    IF( Ierror/=0 ) RETURN
    IF( Nu2-Nu1>.5_DP ) THEN
      prod = SQRT(nu+dmu+1._DP)*prod
      IF( nu/=dmu-1._DP ) prod = prod/SQRT(nu-dmu+1._DP)
      CALL DXADJ(prod,iprod,Ierror)
      IF( Ierror/=0 ) RETURN
      nu = nu + 1._DP
    ELSEIF( dmu>=nu ) THEN
      prod = 0._DP
      iprod = 0
      mu = mu + 1
      dmu = dmu + 1._DP
    ELSE
      prod = SQRT(nu+dmu+1._DP)*prod
      IF( nu>dmu ) prod = prod*SQRT(nu-dmu)
      CALL DXADJ(prod,iprod,Ierror)
      IF( Ierror/=0 ) RETURN
      mu = mu + 1
      dmu = dmu + 1._DP
    END IF
  END DO
END SUBROUTINE DXPNRM
