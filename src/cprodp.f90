!*==CPRODP.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK CPRODP
      SUBROUTINE CPRODP(Nd,Bd,Nm1,Bm1,Nm2,Bm2,Na,Aa,X,Yy,M,A,B,C,D,U,Y)
      IMPLICIT NONE
!*--CPRODP5
!*** Start of declarations inserted by SPAG
      REAL A , Aa , B , Bm1 , Bm2 , C , rt , X , Yy
      INTEGER ia , id , iflg , j , k , M , m1 , m2 , mm , mm2 , Na , Nd , Nm1 , 
     &        Nm2
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  CPRODP
!***SUBSIDIARY
!***PURPOSE  Subsidiary to BLKTRI
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (CPRODP-S, CPROCP-C)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
! PRODP applies a sequence of matrix operations to the vector X and
! stores the result in YY. (Periodic boundary conditions and COMPLEX
! case)
!
! BD,BM1,BM2     are arrays containing roots of certain B polynomials.
! ND,NM1,NM2     are the lengths of the arrays BD,BM1,BM2 respectively.
! AA             Array containing scalar multipliers of the vector X.
! NA             is the length of the array AA.
! X,YY      The matrix operations are applied to X and the result is YY.
! A,B,C          are arrays which contain the tridiagonal matrix.
! M              is the order of the matrix.
! D,U,Y          are working arrays.
! ISGN           determines whether or not a change in sign is made.
!
!***SEE ALSO  BLKTRI
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  CPRODP
!
      COMPLEX Y , D , U , v , den , bh , ym , am , y1 , y2 , yh , Bd , crt
      DIMENSION A(*) , B(*) , C(*) , X(*) , Y(*) , D(*) , U(*) , Bd(*) , Bm1(*)
     &          , Bm2(*) , Aa(*) , Yy(*)
!***FIRST EXECUTABLE STATEMENT  CPRODP
      DO j = 1 , M
        Y(j) = CMPLX(X(j),0.)
      ENDDO
      mm = M - 1
      mm2 = M - 2
      id = Nd
      m1 = Nm1
      m2 = Nm2
      ia = Na
 100  iflg = 0
      IF ( id>0 ) THEN
        crt = Bd(id)
        id = id - 1
        iflg = 1
!
! BEGIN SOLUTION TO SYSTEM
!
        bh = B(M) - crt
        ym = Y(M)
        den = B(1) - crt
        D(1) = C(1)/den
        U(1) = A(1)/den
        Y(1) = Y(1)/den
        v = CMPLX(C(M),0.)
        IF ( mm2>=2 ) THEN
          DO j = 2 , mm2
            den = B(j) - crt - A(j)*D(j-1)
            D(j) = C(j)/den
            U(j) = -A(j)*U(j-1)/den
            Y(j) = (Y(j)-A(j)*Y(j-1))/den
            bh = bh - v*U(j-1)
            ym = ym - v*Y(j-1)
            v = -v*D(j-1)
          ENDDO
        ENDIF
        den = B(M-1) - crt - A(M-1)*D(M-2)
        D(M-1) = (C(M-1)-A(M-1)*U(M-2))/den
        Y(M-1) = (Y(M-1)-A(M-1)*Y(M-2))/den
        am = A(M) - v*D(M-2)
        bh = bh - v*U(M-2)
        ym = ym - v*Y(M-2)
        den = bh - am*D(M-1)
        IF ( ABS(den)/=0 ) THEN
          Y(M) = (ym-am*Y(M-1))/den
        ELSE
          Y(M) = (1.,0.)
        ENDIF
        Y(M-1) = Y(M-1) - D(M-1)*Y(M)
        DO j = 2 , mm
          k = M - j
          Y(k) = Y(k) - D(k)*Y(k+1) - U(k)*Y(M)
        ENDDO
      ENDIF
      IF ( m1<=0 ) THEN
        IF ( m2<=0 ) THEN
          IF ( ia>0 ) THEN
            rt = Aa(ia)
            ia = ia - 1
            iflg = 1
!
! SCALAR MULTIPLICATION
!
            DO j = 1 , M
              Y(j) = rt*Y(j)
            ENDDO
          ENDIF
          IF ( iflg>0 ) GOTO 100
          DO j = 1 , M
            Yy(j) = REAL(Y(j))
          ENDDO
          GOTO 99999
        ELSE
          rt = Bm2(m2)
          m2 = m2 - 1
        ENDIF
      ELSEIF ( m2<=0 ) THEN
        rt = Bm1(m1)
        m1 = m1 - 1
      ELSEIF ( ABS(Bm1(m1))<=ABS(Bm2(m2)) ) THEN
        rt = Bm2(m2)
        m2 = m2 - 1
      ELSE
        rt = Bm1(m1)
        m1 = m1 - 1
      ENDIF
!
! MATRIX MULTIPLICATION
!
      yh = Y(1)
      y1 = (B(1)-rt)*Y(1) + C(1)*Y(2) + A(1)*Y(M)
      IF ( mm>=2 ) THEN
        DO j = 2 , mm
          y2 = A(j)*Y(j-1) + (B(j)-rt)*Y(j) + C(j)*Y(j+1)
          Y(j-1) = y1
          y1 = y2
        ENDDO
      ENDIF
      Y(M) = A(M)*Y(M-1) + (B(M)-rt)*Y(M) + C(M)*yh
      Y(M-1) = y1
      iflg = 1
      GOTO 100
99999 END SUBROUTINE CPRODP
