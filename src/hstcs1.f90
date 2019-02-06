!*==HSTCS1.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK HSTCS1
      SUBROUTINE HSTCS1(Intl,A,B,M,Mbdcnd,Bda,Bdb,C,D,N,Nbdcnd,Bdc,Bdd,Elmbda,F,
     &                  Idimf,Pertrb,Ierr1,Am,Bm,Cm,An,Bn,Cn,Snth,Rsq,Wrk)
      IMPLICIT NONE
!*--HSTCS16
!*** Start of declarations inserted by SPAG
      REAL A , a1 , a2 , a3 , Am , An , B , Bda , Bdb , Bdc , Bdd , Bm , Bn , 
     &     C , Cm , Cn , D , dr , dth , dthsq
      REAL Elmbda , F , Pertrb , Rsq , Snth , Wrk , x , y
      INTEGER i , Idimf , Ierr1 , Intl , isw , j , M , Mbdcnd , N , nb , Nbdcnd
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  HSTCS1
!***SUBSIDIARY
!***PURPOSE  Subsidiary to HSTCSP
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (HSTCS1-S)
!***AUTHOR  (UNKNOWN)
!***SEE ALSO  HSTCSP
!***ROUTINES CALLED  BLKTRI
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  HSTCS1
      DIMENSION Bda(*) , Bdb(*) , Bdc(*) , Bdd(*) , F(Idimf,*) , Am(*) , Bm(*) , 
     &          Cm(*) , An(*) , Bn(*) , Cn(*) , Snth(*) , Rsq(*) , Wrk(*)
!***FIRST EXECUTABLE STATEMENT  HSTCS1
      dth = (B-A)/M
      dthsq = dth*dth
      DO i = 1 , M
        Snth(i) = SIN(A+(i-0.5)*dth)
      ENDDO
      dr = (D-C)/N
      DO j = 1 , N
        Rsq(j) = (C+(j-0.5)*dr)**2
      ENDDO
!
!     MULTIPLY RIGHT SIDE BY R(J)**2
!
      DO j = 1 , N
        x = Rsq(j)
        DO i = 1 , M
          F(i,j) = x*F(i,j)
        ENDDO
      ENDDO
!
!      DEFINE COEFFICIENTS AM,BM,CM
!
      x = 1./(2.*COS(dth/2.))
      DO i = 2 , M
        Am(i) = (Snth(i-1)+Snth(i))*x
        Cm(i-1) = Am(i)
      ENDDO
      Am(1) = SIN(A)
      Cm(M) = SIN(B)
      DO i = 1 , M
        x = 1./Snth(i)
        y = x/dthsq
        Am(i) = Am(i)*y
        Cm(i) = Cm(i)*y
        Bm(i) = Elmbda*x*x - Am(i) - Cm(i)
      ENDDO
!
!     DEFINE COEFFICIENTS AN,BN,CN
!
      x = C/dr
      DO j = 1 , N
        An(j) = (x+j-1)**2
        Cn(j) = (x+j)**2
        Bn(j) = -(An(j)+Cn(j))
      ENDDO
      isw = 1
      nb = Nbdcnd
      IF ( C==0..AND.nb==2 ) nb = 6
!
!     ENTER DATA ON THETA BOUNDARIES
!
      SELECT CASE (Mbdcnd)
      CASE (3,4,8)
        Bm(1) = Bm(1) + Am(1)
        x = dth*Am(1)
        DO j = 1 , N
          F(1,j) = F(1,j) + x*Bda(j)
        ENDDO
      CASE (5,6,9)
      CASE DEFAULT
        Bm(1) = Bm(1) - Am(1)
        x = 2.*Am(1)
        DO j = 1 , N
          F(1,j) = F(1,j) - x*Bda(j)
        ENDDO
      END SELECT
      SELECT CASE (Mbdcnd)
      CASE (2,3,6)
        Bm(M) = Bm(M) + Cm(M)
        x = dth*Cm(M)
        DO j = 1 , N
          F(M,j) = F(M,j) - x*Bdb(j)
        ENDDO
      CASE (7,8,9)
      CASE DEFAULT
        Bm(M) = Bm(M) - Cm(M)
        x = 2.*Cm(M)
        DO j = 1 , N
          F(M,j) = F(M,j) - x*Bdb(j)
        ENDDO
      END SELECT
!
!     ENTER DATA ON R BOUNDARIES
!
      SELECT CASE (nb)
      CASE (3,4)
        Bn(1) = Bn(1) + An(1)
        x = dr*An(1)
        DO i = 1 , M
          F(i,1) = F(i,1) + x*Bdc(i)
        ENDDO
      CASE (5,6)
      CASE DEFAULT
        Bn(1) = Bn(1) - An(1)
        x = 2.*An(1)
        DO i = 1 , M
          F(i,1) = F(i,1) - x*Bdc(i)
        ENDDO
      END SELECT
      SELECT CASE (nb)
      CASE (2,3,6)
        Bn(N) = Bn(N) + Cn(N)
        x = dr*Cn(N)
        DO i = 1 , M
          F(i,N) = F(i,N) - x*Bdd(i)
        ENDDO
      CASE DEFAULT
        Bn(N) = Bn(N) - Cn(N)
        x = 2.*Cn(N)
        DO i = 1 , M
          F(i,N) = F(i,N) - x*Bdd(i)
        ENDDO
      END SELECT
!
!     CHECK FOR SINGULAR PROBLEM.  IF SINGULAR, PERTURB F.
!
      Pertrb = 0.
      SELECT CASE (Mbdcnd)
      CASE (1,2,4,5,7)
      CASE DEFAULT
        SELECT CASE (nb)
        CASE (1,2,4,5)
        CASE DEFAULT
          IF ( Elmbda<0 ) THEN
          ELSEIF ( Elmbda==0 ) THEN
            isw = 2
            DO i = 1 , M
              x = 0.
              DO j = 1 , N
                x = x + F(i,j)
              ENDDO
              Pertrb = Pertrb + x*Snth(i)
            ENDDO
            x = 0.
            DO j = 1 , N
              x = x + Rsq(j)
            ENDDO
            Pertrb = 2.*(Pertrb*SIN(dth/2.))/(x*(COS(A)-COS(B)))
            DO j = 1 , N
              x = Rsq(j)*Pertrb
              DO i = 1 , M
                F(i,j) = F(i,j) - x
              ENDDO
            ENDDO
          ELSE
            Ierr1 = 10
          ENDIF
        END SELECT
      END SELECT
      a2 = 0.
      DO i = 1 , M
        a2 = a2 + F(i,1)
      ENDDO
      a2 = a2/Rsq(1)
!
!     INITIALIZE BLKTRI
!
      IF ( Intl==0 ) CALL BLKTRI(0,1,N,An,Bn,Cn,1,M,Am,Bm,Cm,Idimf,F,Ierr1,Wrk)
!
!     CALL BLKTRI TO SOLVE SYSTEM OF EQUATIONS.
!
      CALL BLKTRI(1,1,N,An,Bn,Cn,1,M,Am,Bm,Cm,Idimf,F,Ierr1,Wrk)
      IF ( isw==2.AND.C==0..AND.Nbdcnd==2 ) THEN
        a1 = 0.
        a3 = 0.
        DO i = 1 , M
          a1 = a1 + Snth(i)*F(i,1)
          a3 = a3 + Snth(i)
        ENDDO
        a1 = a1 + Rsq(1)*a2/2.
        IF ( Mbdcnd==3 ) a1 = a1 + (SIN(B)*Bdb(1)-SIN(A)*Bda(1))/(2.*(B-A))
        a1 = a1/a3
        a1 = Bdc(1) - a1
        DO i = 1 , M
          DO j = 1 , N
            F(i,j) = F(i,j) + a1
          ENDDO
        ENDDO
      ENDIF
      END SUBROUTINE HSTCS1
