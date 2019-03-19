!** DLSI
SUBROUTINE DLSI(W,Mdw,Ma,Mg,N,Prgopt,X,Rnorm,Mode,Ws,Ip)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to DLSEI
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (LSI-S, DLSI-D)
  !***
  ! **Author:**  Hanson, R. J., (SNLA)
  !***
  ! **Description:**
  !
  !     This is a companion subprogram to DLSEI.  The documentation for
  !     DLSEI has complete usage instructions.
  !
  !     Solve..
  !              AX = B,  A  MA by N  (least squares equations)
  !     subject to..
  !
  !              GX.GE.H, G  MG by N  (inequality constraints)
  !
  !     Input..
  !
  !      W(*,*) contains  (A B) in rows 1,...,MA+MG, cols 1,...,N+1.
  !                       (G H)
  !
  !     MDW,MA,MG,N
  !              contain (resp) var. dimension of W(*,*),
  !              and matrix dimensions.
  !
  !     PRGOPT(*),
  !              Program option vector.
  !
  !     OUTPUT..
  !
  !      X(*),RNORM
  !
  !              Solution vector(unless MODE=2), length of AX-B.
  !
  !      MODE
  !              =0   Inequality constraints are compatible.
  !              =2   Inequality constraints contradictory.
  !
  !      WS(*),
  !              Working storage of dimension K+N+(MG+2)*(N+7),
  !              where K=MAX(MA+MG,N).
  !      IP(MG+2*N+1)
  !              Integer working storage
  !
  !***
  ! **Routines called:**  D1MACH, DASUM, DAXPY, DCOPY, DDOT, DH12, DHFTI,
  !                    DLPDP, DSCAL, DSWAP

  !* REVISION HISTORY  (YYMMDD)
  !   790701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890618  Completely restructured and extensively revised (WRB & RWC)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   900604  DP version created from SP version.  (RWC)
  !   920422  Changed CALL to DHFTI to include variable MA.  (WRB)
  
  INTEGER Ip(*), Ma, Mdw, Mg, Mode, N
  REAL(8) :: Prgopt(*), Rnorm, W(Mdw,*), Ws(*), X(*)
  !
  EXTERNAL D1MACH, DASUM, DAXPY, DCOPY, DDOT, DH12, DHFTI, DLPDP, &
    DSCAL, DSWAP
  REAL(8) :: D1MACH, DASUM, DDOT
  !
  REAL(8) :: anorm, drelpr, fac, gam, rb, tau, tol, xnorm, temp(1)
  INTEGER i, j, k, key, krank, krm1, krp1, l, last, link, m, &
    map1, mdlpdp, minman, n1, n2, n3, next, np1
  LOGICAL cov, first, sclcov
  !
  SAVE drelpr, first
  DATA first/.TRUE./
  !
  !* FIRST EXECUTABLE STATEMENT  DLSI
  !
  !     Set the nominal tolerance used in the code.
  !
  IF ( first ) drelpr = D1MACH(4)
  first = .FALSE.
  tol = SQRT(drelpr)
  !
  Mode = 0
  Rnorm = 0.D0
  m = Ma + Mg
  np1 = N + 1
  krank = 0
  IF ( N>0.AND.m>0 ) THEN
    !
    !     To process option vector.
    !
    cov = .FALSE.
    sclcov = .TRUE.
    last = 1
    link = INT( Prgopt(1) )
    DO
      !
      IF ( link>1 ) THEN
        key = INT( Prgopt(last+1) )
        IF ( key==1 ) cov = Prgopt(last+2)/=0.D0
        IF ( key==10 ) sclcov = Prgopt(last+2)==0.D0
        IF ( key==5 ) tol = MAX(drelpr,Prgopt(last+2))
        next = INT( Prgopt(link) )
        last = link
        link = next
        CYCLE
      ENDIF
      !
      !     Compute matrix norm of least squares equations.
      !
      anorm = 0.D0
      DO j = 1, N
        anorm = MAX(anorm,DASUM(Ma,W(1,j),1))
      ENDDO
      !
      !     Set tolerance for DHFTI( ) rank test.
      !
      tau = tol*anorm
      !
      !     Compute Householder orthogonal decomposition of matrix.
      !
      CALL DCOPY(N,0.D0,0,Ws,1)
      CALL DCOPY(Ma,W(1,np1),1,Ws,1)
      k = MAX(m,N)
      minman = MIN(Ma,N)
      n1 = k + 1
      n2 = n1 + N
      CALL DHFTI(W,Mdw,Ma,N,Ws,Ma,1,tau,krank,temp,Ws(n2),Ws(n1),Ip)
      Rnorm = temp(1)
      fac = 1.D0
      gam = Ma - krank
      IF ( krank<Ma.AND.sclcov ) fac = Rnorm**2/gam
      !
      !     Reduce to DLPDP and solve.
      !
      map1 = Ma + 1
      !
      !     Compute inequality rt-hand side for DLPDP.
      !
      IF ( Ma<m ) THEN
        IF ( minman>0 ) THEN
          DO i = map1, m
            W(i,np1) = W(i,np1) - DDOT(N,W(i,1),Mdw,Ws,1)
          ENDDO
          !
          !           Apply permutations to col. of inequality constraint matrix.
          !
          DO i = 1, minman
            CALL DSWAP(Mg,W(map1,i),1,W(map1,Ip(i)),1)
          ENDDO
          !
          !           Apply Householder transformations to constraint matrix.
          !
          IF ( krank>0.AND.krank<N ) THEN
            DO i = krank, 1, -1
              CALL DH12(2,i,krank+1,N,W(i,1),Mdw,Ws(n1+i-1),W(map1,1),Mdw,1,&
                Mg)
            ENDDO
          ENDIF
          !
          !           Compute permuted inequality constraint matrix times r-inv.
          !
          DO i = map1, m
            DO j = 1, krank
              W(i,j) = (W(i,j)-DDOT(j-1,W(1,j),1,W(i,1),Mdw))/W(j,j)
            ENDDO
          ENDDO
        ENDIF
        !
        !        Solve the reduced problem with DLPDP algorithm,
        !        the least projected distance problem.
        !
        CALL DLPDP(W(map1,1),Mdw,Mg,krank,N-krank,Prgopt,X,xnorm,mdlpdp,&
          Ws(n2),Ip(N+1))
        !
        !        Compute solution in original coordinates.
        !
        IF ( mdlpdp==1 ) THEN
          DO i = krank, 1, -1
            X(i) = (X(i)-DDOT(krank-i,W(i,i+1),Mdw,X(i+1),1))/W(i,i)
          ENDDO
          !
          !           Apply Householder transformation to solution vector.
          !
          IF ( krank<N ) THEN
            DO i = 1, krank
              CALL DH12(2,i,krank+1,N,W(i,1),Mdw,Ws(n1+i-1),X,1,1,1)
            ENDDO
          ENDIF
          !
          !           Repermute variables to their input order.
          !
          IF ( minman>0 ) THEN
            DO i = minman, 1, -1
              CALL DSWAP(1,X(i),1,X(Ip(i)),1)
            ENDDO
            !
            !              Variables are now in original coordinates.
            !              Add solution of unconstrained problem.
            !
            DO i = 1, N
              X(i) = X(i) + Ws(i)
            ENDDO
            !
            !              Compute the residual vector norm.
            !
            Rnorm = SQRT(Rnorm**2+xnorm**2)
          ENDIF
        ELSE
          Mode = 2
        ENDIF
      ELSE
        CALL DCOPY(N,Ws,1,X,1)
      ENDIF
      !
      !     Compute covariance matrix based on the orthogonal decomposition
      !     from DHFTI( ).
      !
      IF ( .NOT.(.NOT.cov.OR.krank<=0) ) THEN
        krm1 = krank - 1
        krp1 = krank + 1
        !
        !     Copy diagonal terms to working array.
        !
        CALL DCOPY(krank,W,Mdw+1,Ws(n2),1)
        !
        !     Reciprocate diagonal terms.
        !
        DO j = 1, krank
          W(j,j) = 1.D0/W(j,j)
        ENDDO
        !
        !     Invert the upper triangular QR factor on itself.
        !
        IF ( krank>1 ) THEN
          DO i = 1, krm1
            DO j = i + 1, krank
              W(i,j) = -DDOT(j-i,W(i,i),Mdw,W(i,j),1)*W(j,j)
            ENDDO
          ENDDO
        ENDIF
        !
        !     Compute the inverted factor times its transpose.
        !
        DO i = 1, krank
          DO j = i, krank
            W(i,j) = DDOT(krank+1-j,W(i,j),Mdw,W(j,j),Mdw)
          ENDDO
        ENDDO
        !
        !     Zero out lower trapezoidal part.
        !     Copy upper triangular to lower triangular part.
        !
        IF ( krank<N ) THEN
          DO j = 1, krank
            CALL DCOPY(j,W(1,j),1,W(j,1),Mdw)
          ENDDO
          !
          DO i = krp1, N
            CALL DCOPY(i,0.D0,0,W(i,1),Mdw)
          ENDDO
          !
          !        Apply right side transformations to lower triangle.
          !
          n3 = n2 + krp1
          DO i = 1, krank
            l = n1 + i
            k = n2 + i
            rb = Ws(l-1)*Ws(k-1)
            !
            !           If RB.GE.0.D0, transformation can be regarded as zero.
            !
            IF ( rb<0.D0 ) THEN
              rb = 1.D0/rb
              !
              !              Store unscaled rank one Householder update in work array.
              !
              CALL DCOPY(N,0.D0,0,Ws(n3),1)
              l = n1 + i
              k = n3 + i
              Ws(k-1) = Ws(l-1)
              !
              DO j = krp1, N
                Ws(n3+j-1) = W(i,j)
              ENDDO
              !
              DO j = 1, N
                Ws(j) = rb*(DDOT(j-i,W(j,i),Mdw,Ws(n3+i-1),1)+DDOT(N-j+1,W(j&
                  ,j),1,Ws(n3+j-1),1))
              ENDDO
              !
              l = n3 + i
              gam = 0.5D0*rb*DDOT(N-i+1,Ws(l-1),1,Ws(i),1)
              CALL DAXPY(N-i+1,gam,Ws(l-1),1,Ws(i),1)
              DO j = i, N
                DO l = 1, i - 1
                  W(j,l) = W(j,l) + Ws(n3+j-1)*Ws(l)
                ENDDO
                !
                DO l = i, j
                  W(j,l) = W(j,l) + Ws(j)*Ws(n3+l-1) + Ws(l)*Ws(n3+j-1)
                ENDDO
              ENDDO
            ENDIF
          ENDDO
          !
          !        Copy lower triangle to upper triangle to symmetrize the
          !        covariance matrix.
          !
          DO i = 1, N
            CALL DCOPY(i,W(i,1),Mdw,W(1,i),1)
          ENDDO
        ENDIF
        !
        !     Repermute rows and columns.
        !
        DO i = minman, 1, -1
          k = Ip(i)
          IF ( i/=k ) THEN
            CALL DSWAP(1,W(i,i),1,W(k,k),1)
            CALL DSWAP(i-1,W(1,i),1,W(1,k),1)
            CALL DSWAP(k-i-1,W(i,i+1),Mdw,W(i+1,k),1)
            CALL DSWAP(N-k,W(i,k+1),Mdw,W(k,k+1),Mdw)
          ENDIF
        ENDDO
        !
        !     Put in normalized residual sum of squares scale factor
        !     and symmetrize the resulting covariance matrix.
        !
        DO j = 1, N
          CALL DSCAL(j,fac,W(1,j),1)
          CALL DCOPY(j,W(1,j),1,W(j,1),Mdw)
        ENDDO
      ENDIF
      EXIT
    ENDDO
  ENDIF
  !
  Ip(1) = krank
  Ip(2) = N + MAX(m,N) + (Mg+2)*(N+7)
END SUBROUTINE DLSI
