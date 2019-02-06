!*==CUNK1.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK CUNK1
      SUBROUTINE CUNK1(Z,Fnu,Kode,Mr,N,Y,Nz,Tol,Elim,Alim)
      IMPLICIT NONE
!*--CUNK15
!***BEGIN PROLOGUE  CUNK1
!***SUBSIDIARY
!***PURPOSE  Subsidiary to CBESK
!***LIBRARY   SLATEC
!***TYPE      ALL (CUNK1-A, ZUNK1-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CUNK1 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
!     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
!     UNIFORM ASYMPTOTIC EXPANSION.
!     MR INDICATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
!     NZ=-1 MEANS AN OVERFLOW WILL OCCUR
!
!***SEE ALSO  CBESK
!***ROUTINES CALLED  CS1S2, CUCHK, CUNIK, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CUNK1
      COMPLEX cfn , ck , cone , crsc , cs , cscl , csgn , cspn , csr , css , 
     &        cwrk , cy , czero , c1 , c2 , phi , rz , sum , s1 , s2 , Y , Z , 
     &        zeta1 , zeta2 , zr , phid , zeta1d , zeta2d , sumd
      REAL Alim , ang , aphi , asc , ascle , bry , cpn , c2i , c2m , c2r , 
     &     Elim , fmr , fn , fnf , Fnu , pi , rs1 , sgn , spn , Tol , x , R1MACH
      INTEGER i , ib , iflag , ifn , il , init , inu , iuf , k , kdflg , kflag , 
     &        kk , Kode , Mr , N , nw , Nz , j , ipard , initd , ic , m
      DIMENSION bry(3) , init(2) , Y(N) , sum(2) , phi(2) , zeta1(2) , zeta2(2)
     &          , cy(2) , cwrk(16,3) , css(3) , csr(3)
      DATA czero , cone/(0.0E0,0.0E0) , (1.0E0,0.0E0)/
      DATA pi/3.14159265358979324E0/
!***FIRST EXECUTABLE STATEMENT  CUNK1
      kdflg = 1
      Nz = 0
!-----------------------------------------------------------------------
!     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
!     THE UNDERFLOW LIMIT
!-----------------------------------------------------------------------
      cscl = CMPLX(1.0E0/Tol,0.0E0)
      crsc = CMPLX(Tol,0.0E0)
      css(1) = cscl
      css(2) = cone
      css(3) = crsc
      csr(1) = crsc
      csr(2) = cone
      csr(3) = cscl
      bry(1) = 1.0E+3*R1MACH(1)/Tol
      bry(2) = 1.0E0/bry(1)
      bry(3) = R1MACH(2)
      x = REAL(Z)
      zr = Z
      IF ( x<0.0E0 ) zr = -Z
      j = 2
      DO i = 1 , N
!-----------------------------------------------------------------------
!     J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J
!-----------------------------------------------------------------------
        j = 3 - j
        fn = Fnu + (i-1)
        init(j) = 0
        CALL CUNIK(zr,fn,2,0,Tol,init(j),phi(j),zeta1(j),zeta2(j),sum(j),
     &             cwrk(1,j))
        IF ( Kode==1 ) THEN
          s1 = zeta1(j) - zeta2(j)
        ELSE
          cfn = CMPLX(fn,0.0E0)
          s1 = zeta1(j) - cfn*(cfn/(zr+zeta2(j)))
        ENDIF
!-----------------------------------------------------------------------
!     TEST FOR UNDERFLOW AND OVERFLOW
!-----------------------------------------------------------------------
        rs1 = REAL(s1)
        IF ( ABS(rs1)<=Elim ) THEN
          IF ( kdflg==1 ) kflag = 2
          IF ( ABS(rs1)>=Alim ) THEN
!-----------------------------------------------------------------------
!     REFINE  TEST AND SCALE
!-----------------------------------------------------------------------
            aphi = ABS(phi(j))
            rs1 = rs1 + ALOG(aphi)
            IF ( ABS(rs1)>Elim ) GOTO 50
            IF ( kdflg==1 ) kflag = 1
            IF ( rs1>=0.0E0 ) THEN
              IF ( kdflg==1 ) kflag = 3
            ENDIF
          ENDIF
!-----------------------------------------------------------------------
!     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
!     EXPONENT EXTREMES
!-----------------------------------------------------------------------
          s2 = phi(j)*sum(j)
          c2r = REAL(s1)
          c2i = AIMAG(s1)
          c2m = EXP(c2r)*REAL(css(kflag))
          s1 = CMPLX(c2m,0.0E0)*CMPLX(COS(c2i),SIN(c2i))
          s2 = s2*s1
          IF ( kflag==1 ) THEN
            CALL CUCHK(s2,nw,bry(1),Tol)
            IF ( nw/=0 ) GOTO 50
          ENDIF
          cy(kdflg) = s2
          Y(i) = s2*csr(kflag)
          IF ( kdflg==2 ) GOTO 100
          kdflg = 2
          CYCLE
        ENDIF
 50     IF ( rs1>0.0E0 ) GOTO 600
!-----------------------------------------------------------------------
!     FOR X.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
!-----------------------------------------------------------------------
        IF ( x<0.0E0 ) GOTO 600
        kdflg = 1
        Y(i) = czero
        Nz = Nz + 1
        IF ( i/=1 ) THEN
          IF ( Y(i-1)/=czero ) THEN
            Y(i-1) = czero
            Nz = Nz + 1
          ENDIF
        ENDIF
      ENDDO
      i = N
 100  rz = CMPLX(2.0E0,0.0E0)/zr
      ck = CMPLX(fn,0.0E0)*rz
      ib = i + 1
      IF ( N<ib ) GOTO 300
!-----------------------------------------------------------------------
!     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW, SET SEQUENCE TO ZERO
!     ON UNDERFLOW
!-----------------------------------------------------------------------
      fn = Fnu + (N-1)
      ipard = 1
      IF ( Mr/=0 ) ipard = 0
      initd = 0
      CALL CUNIK(zr,fn,2,ipard,Tol,initd,phid,zeta1d,zeta2d,sumd,cwrk(1,3))
      IF ( Kode==1 ) THEN
        s1 = zeta1d - zeta2d
      ELSE
        cfn = CMPLX(fn,0.0E0)
        s1 = zeta1d - cfn*(cfn/(zr+zeta2d))
      ENDIF
      rs1 = REAL(s1)
      IF ( ABS(rs1)<=Elim ) THEN
        IF ( ABS(rs1)<Alim ) GOTO 200
!-----------------------------------------------------------------------
!     REFINE ESTIMATE AND TEST
!-----------------------------------------------------------------------
        aphi = ABS(phid)
        rs1 = rs1 + ALOG(aphi)
        IF ( ABS(rs1)<Elim ) GOTO 200
      ENDIF
      IF ( rs1>0.0E0 ) GOTO 600
!-----------------------------------------------------------------------
!     FOR X.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
!-----------------------------------------------------------------------
      IF ( x<0.0E0 ) GOTO 600
      Nz = N
      DO i = 1 , N
        Y(i) = czero
      ENDDO
      RETURN
!-----------------------------------------------------------------------
!     RECUR FORWARD FOR REMAINDER OF THE SEQUENCE
!-----------------------------------------------------------------------
 200  s1 = cy(1)
      s2 = cy(2)
      c1 = csr(kflag)
      ascle = bry(kflag)
      DO i = ib , N
        c2 = s2
        s2 = ck*s2 + s1
        s1 = c2
        ck = ck + rz
        c2 = s2*c1
        Y(i) = c2
        IF ( kflag<3 ) THEN
          c2r = REAL(c2)
          c2i = AIMAG(c2)
          c2r = ABS(c2r)
          c2i = ABS(c2i)
          c2m = MAX(c2r,c2i)
          IF ( c2m>ascle ) THEN
            kflag = kflag + 1
            ascle = bry(kflag)
            s1 = s1*c1
            s2 = c2
            s1 = s1*css(kflag)
            s2 = s2*css(kflag)
            c1 = csr(kflag)
          ENDIF
        ENDIF
      ENDDO
 300  IF ( Mr==0 ) RETURN
!-----------------------------------------------------------------------
!     ANALYTIC CONTINUATION FOR RE(Z).LT.0.0E0
!-----------------------------------------------------------------------
      Nz = 0
      fmr = Mr
      sgn = -SIGN(pi,fmr)
!-----------------------------------------------------------------------
!     CSPN AND CSGN ARE COEFF OF K AND I FUNCTIONS RESP.
!-----------------------------------------------------------------------
      csgn = CMPLX(0.0E0,sgn)
      inu = Fnu
      fnf = Fnu - inu
      ifn = inu + N - 1
      ang = fnf*sgn
      cpn = COS(ang)
      spn = SIN(ang)
      cspn = CMPLX(cpn,spn)
      IF ( MOD(ifn,2)==1 ) cspn = -cspn
      asc = bry(1)
      kk = N
      iuf = 0
      kdflg = 1
      ib = ib - 1
      ic = ib - 1
      DO k = 1 , N
        fn = Fnu + (kk-1)
!-----------------------------------------------------------------------
!     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
!     FUNCTION ABOVE
!-----------------------------------------------------------------------
        m = 3
        IF ( N>2 ) THEN
          IF ( (kk==N).AND.(ib<N) ) GOTO 350
          IF ( (kk/=ib).AND.(kk/=ic) ) THEN
            initd = 0
            GOTO 350
          ENDIF
        ENDIF
        initd = init(j)
        phid = phi(j)
        zeta1d = zeta1(j)
        zeta2d = zeta2(j)
        sumd = sum(j)
        m = j
        j = 3 - j
 350    CALL CUNIK(zr,fn,1,0,Tol,initd,phid,zeta1d,zeta2d,sumd,cwrk(1,m))
        IF ( Kode==1 ) THEN
          s1 = -zeta1d + zeta2d
        ELSE
          cfn = CMPLX(fn,0.0E0)
          s1 = -zeta1d + cfn*(cfn/(zr+zeta2d))
        ENDIF
!-----------------------------------------------------------------------
!     TEST FOR UNDERFLOW AND OVERFLOW
!-----------------------------------------------------------------------
        rs1 = REAL(s1)
        IF ( ABS(rs1)>Elim ) GOTO 450
        IF ( kdflg==1 ) iflag = 2
        IF ( ABS(rs1)>=Alim ) THEN
!-----------------------------------------------------------------------
!     REFINE  TEST AND SCALE
!-----------------------------------------------------------------------
          aphi = ABS(phid)
          rs1 = rs1 + ALOG(aphi)
          IF ( ABS(rs1)>Elim ) GOTO 450
          IF ( kdflg==1 ) iflag = 1
          IF ( rs1>=0.0E0 ) THEN
            IF ( kdflg==1 ) iflag = 3
          ENDIF
        ENDIF
        s2 = csgn*phid*sumd
        c2r = REAL(s1)
        c2i = AIMAG(s1)
        c2m = EXP(c2r)*REAL(css(iflag))
        s1 = CMPLX(c2m,0.0E0)*CMPLX(COS(c2i),SIN(c2i))
        s2 = s2*s1
        IF ( iflag==1 ) THEN
          CALL CUCHK(s2,nw,bry(1),Tol)
          IF ( nw/=0 ) s2 = CMPLX(0.0E0,0.0E0)
        ENDIF
 400    cy(kdflg) = s2
        c2 = s2
        s2 = s2*csr(iflag)
!-----------------------------------------------------------------------
!     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
!-----------------------------------------------------------------------
        s1 = Y(kk)
        IF ( Kode/=1 ) THEN
          CALL CS1S2(zr,s1,s2,nw,asc,Alim,iuf)
          Nz = Nz + nw
        ENDIF
        Y(kk) = s1*cspn + s2
        kk = kk - 1
        cspn = -cspn
        IF ( c2/=czero ) THEN
          IF ( kdflg==2 ) GOTO 500
          kdflg = 2
          CYCLE
        ELSE
          kdflg = 1
          CYCLE
        ENDIF
 450    IF ( rs1>0.0E0 ) GOTO 600
        s2 = czero
        GOTO 400
      ENDDO
      k = N
 500  il = N - k
      IF ( il==0 ) RETURN
!-----------------------------------------------------------------------
!     RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
!     K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP
!     INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES.
!-----------------------------------------------------------------------
      s1 = cy(1)
      s2 = cy(2)
      cs = csr(iflag)
      ascle = bry(iflag)
      fn = (inu+il)
      DO i = 1 , il
        c2 = s2
        s2 = s1 + CMPLX(fn+fnf,0.0E0)*rz*s2
        s1 = c2
        fn = fn - 1.0E0
        c2 = s2*cs
        ck = c2
        c1 = Y(kk)
        IF ( Kode/=1 ) THEN
          CALL CS1S2(zr,c1,c2,nw,asc,Alim,iuf)
          Nz = Nz + nw
        ENDIF
        Y(kk) = c1*cspn + c2
        kk = kk - 1
        cspn = -cspn
        IF ( iflag<3 ) THEN
          c2r = REAL(ck)
          c2i = AIMAG(ck)
          c2r = ABS(c2r)
          c2i = ABS(c2i)
          c2m = MAX(c2r,c2i)
          IF ( c2m>ascle ) THEN
            iflag = iflag + 1
            ascle = bry(iflag)
            s1 = s1*cs
            s2 = ck
            s1 = s1*css(iflag)
            s2 = s2*css(iflag)
            cs = csr(iflag)
          ENDIF
        ENDIF
      ENDDO
      RETURN
 600  Nz = -1
      END SUBROUTINE CUNK1
