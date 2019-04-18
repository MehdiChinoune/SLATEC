!** HQR2
SUBROUTINE HQR2(Nm,N,Low,Igh,H,Wr,Wi,Z,Ierr)
  !>
  !  Compute the eigenvalues and eigenvectors of a real upper
  !            Hessenberg matrix using QR method.
  !***
  ! **Library:**   SLATEC (EISPACK)
  !***
  ! **Category:**  D4C2B
  !***
  ! **Type:**      SINGLE PRECISION (HQR2-S, COMQR2-C)
  !***
  ! **Keywords:**  EIGENVALUES, EIGENVECTORS, EISPACK
  !***
  ! **Author:**  Smith, B. T., et al.
  !***
  ! **Description:**
  !
  !     This subroutine is a translation of the ALGOL procedure HQR2,
  !     NUM. MATH. 16, 181-204(1970) by Peters and Wilkinson.
  !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
  !
  !     This subroutine finds the eigenvalues and eigenvectors
  !     of a REAL UPPER Hessenberg matrix by the QR method.  The
  !     eigenvectors of a REAL GENERAL matrix can also be found
  !     if  ELMHES  and  ELTRAN  or  ORTHES  and  ORTRAN  have
  !     been used to reduce this general matrix to Hessenberg form
  !     and to accumulate the similarity transformations.
  !
  !     On INPUT
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameters, H and Z, as declared in the calling
  !          program dimension statement.  NM is an INTEGER variable.
  !
  !        N is the order of the matrix H.  N is an INTEGER variable.
  !          N must be less than or equal to NM.
  !
  !        LOW and IGH are two INTEGER variables determined by the
  !          balancing subroutine  BALANC.  If  BALANC  has not been
  !          used, set LOW=1 and IGH equal to the order of the matrix, N.
  !
  !        H contains the upper Hessenberg matrix.  H is a two-dimensional
  !          REAL array, dimensioned H(NM,N).
  !
  !        Z contains the transformation matrix produced by  ELTRAN
  !          after the reduction by  ELMHES, or by  ORTRAN  after the
  !          reduction by  ORTHES, if performed.  If the eigenvectors
  !          of the Hessenberg matrix are desired, Z must contain the
  !          identity matrix.  Z is a two-dimensional REAL array,
  !          dimensioned Z(NM,M).
  !
  !     On OUTPUT
  !
  !        H has been destroyed.
  !
  !        WR and WI contain the real and imaginary parts, respectively,
  !          of the eigenvalues.  The eigenvalues are unordered except
  !          that complex conjugate pairs of values appear consecutively
  !          with the eigenvalue having the positive imaginary part first.
  !          If an error exit is made, the eigenvalues should be correct
  !          for indices IERR+1, IERR+2, ..., N.  WR and WI are one-
  !          dimensional REAL arrays, dimensioned WR(N) and WI(N).
  !
  !        Z contains the real and imaginary parts of the eigenvectors.
  !          If the J-th eigenvalue is real, the J-th column of Z
  !          contains its eigenvector.  If the J-th eigenvalue is complex
  !          with positive imaginary part, the J-th and (J+1)-th
  !          columns of Z contain the real and imaginary parts of its
  !          eigenvector.  The eigenvectors are unnormalized.  If an
  !          error exit is made, none of the eigenvectors has been found.
  !
  !        IERR is an INTEGER flag set to
  !          Zero       for normal return,
  !          J          if the J-th eigenvalue has not been
  !                     determined after a total of 30*N iterations.
  !                     The eigenvalues should be correct for indices
  !                     IERR+1, IERR+2, ..., N, but no eigenvectors are
  !                     computed.
  !
  !     Calls CDIV for complex division.
  !
  !     Questions and comments should be directed to B. S. Garbow,
  !     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
  !     ------------------------------------------------------------------
  !
  !***
  ! **References:**  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
  !                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
  !                 system Routines - EISPACK Guide, Springer-Verlag,
  !                 1976.
  !***
  ! **Routines called:**  CDIV

  !* REVISION HISTORY  (YYMMDD)
  !   760101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  !
  INTEGER i, j, k, l, m, N, en, ii, jj, ll, mm, na, Nm, nn
  INTEGER Igh, itn, its, Low, mp2, enm2, Ierr
  REAL H(Nm,*), Wr(*), Wi(*), Z(Nm,*)
  REAL p, q, r, s, t, w, x, y, ra, sa, vi, vr, zz, norm, s1, s2
  LOGICAL notlas
  !
  !* FIRST EXECUTABLE STATEMENT  HQR2
  Ierr = 0
  norm = 0.0E0
  k = 1
  !     .......... STORE ROOTS ISOLATED BY BALANC
  !                AND COMPUTE MATRIX NORM ..........
  DO i = 1, N
    !
    DO j = k, N
      norm = norm + ABS(H(i,j))
    END DO
    !
    k = i
    IF ( i<Low.OR.i>Igh ) THEN
      Wr(i) = H(i,i)
      Wi(i) = 0.0E0
    END IF
  END DO
  !
  en = Igh
  t = 0.0E0
  itn = 30*N
  !     .......... SEARCH FOR NEXT EIGENVALUES ..........
  100 CONTINUE
  IF ( en<Low ) THEN
    !     .......... ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND
    !                VECTORS OF UPPER TRIANGULAR FORM ..........
    IF ( norm/=0.0E0 ) THEN
      !     .......... FOR EN=N STEP -1 UNTIL 1 DO -- ..........
      DO nn = 1, N
        en = N + 1 - nn
        p = Wr(en)
        q = Wi(en)
        na = en - 1
        IF ( q<0 ) THEN
          !     .......... COMPLEX VECTOR ..........
          m = na
          !     .......... LAST VECTOR COMPONENT CHOSEN IMAGINARY SO THAT
          !                EIGENVECTOR MATRIX IS TRIANGULAR ..........
          IF ( ABS(H(en,na))<=ABS(H(na,en)) ) THEN
            CALL CDIV(0.0E0,-H(na,en),H(na,na)-p,q,H(na,na),H(na,en))
          ELSE
            H(na,na) = q/H(en,na)
            H(na,en) = -(H(en,en)-p)/H(en,na)
          END IF
          H(en,na) = 0.0E0
          H(en,en) = 1.0E0
          enm2 = na - 1
          IF ( enm2/=0 ) THEN
            !     .......... FOR I=EN-2 STEP -1 UNTIL 1 DO -- ..........
            DO ii = 1, enm2
              i = na - ii
              w = H(i,i) - p
              ra = 0.0E0
              sa = H(i,en)
              !
              DO j = m, na
                ra = ra + H(i,j)*H(j,na)
                sa = sa + H(i,j)*H(j,en)
              END DO
              !
              IF ( Wi(i)>=0.0E0 ) THEN
                m = i
                IF ( Wi(i)/=0.0E0 ) THEN
                  !     .......... SOLVE COMPLEX EQUATIONS ..........
                  x = H(i,i+1)
                  y = H(i+1,i)
                  vr = (Wr(i)-p)*(Wr(i)-p) + Wi(i)*Wi(i) - q*q
                  vi = (Wr(i)-p)*2.0E0*q
                  IF ( vr==0.0E0.AND.vi==0.0E0 ) THEN
                    s1 = norm*(ABS(w)+ABS(q)+ABS(x)+ABS(y)+ABS(zz))
                    vr = s1
                    DO
                      vr = 0.5E0*vr
                      IF ( s1+vr<=s1 ) THEN
                        vr = 2.0E0*vr
                        EXIT
                      END IF
                    END DO
                  END IF
                  CALL CDIV(x*r-zz*ra+q*sa,x*s-zz*sa-q*ra,vr,vi,H(i,na),H(i,en))
                  IF ( ABS(x)<=ABS(zz)+ABS(q) ) THEN
                    CALL CDIV(-r-y*H(i,na),-s-y*H(i,en),zz,q,H(i+1,na),H(i+1,en))
                  ELSE
                    H(i+1,na) = (-ra-w*H(i,na)+q*H(i,en))/x
                    H(i+1,en) = (-sa-w*H(i,en)-q*H(i,na))/x
                  END IF
                ELSE
                  CALL CDIV(-ra,-sa,w,q,H(i,na),H(i,en))
                END IF
              ELSE
                zz = w
                r = ra
                s = sa
              END IF
            END DO
          END IF
        ELSEIF ( q==0 ) THEN
          !     .......... REAL VECTOR ..........
          m = en
          H(en,en) = 1.0E0
          IF ( na/=0 ) THEN
            !     .......... FOR I=EN-1 STEP -1 UNTIL 1 DO -- ..........
            DO ii = 1, na
              i = en - ii
              w = H(i,i) - p
              r = H(i,en)
              IF ( m<=na ) THEN
                !
                DO j = m, na
                  r = r + H(i,j)*H(j,en)
                END DO
              END IF
              !
              IF ( Wi(i)>=0.0E0 ) THEN
                m = i
                IF ( Wi(i)/=0.0E0 ) THEN
                  !     .......... SOLVE REAL EQUATIONS ..........
                  x = H(i,i+1)
                  y = H(i+1,i)
                  q = (Wr(i)-p)*(Wr(i)-p) + Wi(i)*Wi(i)
                  t = (x*s-zz*r)/q
                  H(i,en) = t
                  IF ( ABS(x)<=ABS(zz) ) THEN
                    H(i+1,en) = (-s-y*t)/zz
                  ELSE
                    H(i+1,en) = (-r-w*t)/x
                  END IF
                ELSE
                  t = w
                  IF ( t==0.0E0 ) THEN
                    t = norm
                    DO
                      t = 0.5E0*t
                      IF ( norm+t<=norm ) THEN
                        t = 2.0E0*t
                        EXIT
                      END IF
                    END DO
                  END IF
                  H(i,en) = -r/t
                END IF
              ELSE
                zz = w
                s = r
              END IF
              !     .......... END REAL VECTOR ..........
            END DO
          END IF
        END IF
        !     .......... END COMPLEX VECTOR ..........
      END DO
      !     .......... END BACK SUBSTITUTION.
      !                VECTORS OF ISOLATED ROOTS ..........
      DO i = 1, N
        IF ( i<Low.OR.i>Igh ) THEN
          !
          DO j = i, N
            Z(i,j) = H(i,j)
          END DO
        END IF
        !
      END DO
      !     .......... MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
      !                VECTORS OF ORIGINAL FULL MATRIX.
      !                FOR J=N STEP -1 UNTIL LOW DO -- ..........
      DO jj = Low, N
        j = N + Low - jj
        m = MIN(j,Igh)
        !
        DO i = Low, Igh
          zz = 0.0E0
          !
          DO k = Low, m
            zz = zz + Z(i,k)*H(k,j)
          END DO
          !
          Z(i,j) = zz
        END DO
        !
      END DO
    END IF
    RETURN
  ELSE
    its = 0
    na = en - 1
    enm2 = na - 1
  END IF
  !     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
  !                FOR L=EN STEP -1 UNTIL LOW DO -- ..........
  200 CONTINUE
  DO ll = Low, en
    l = en + Low - ll
    IF ( l==Low ) EXIT
    s = ABS(H(l-1,l-1)) + ABS(H(l,l))
    IF ( s==0.0E0 ) s = norm
    s2 = s + ABS(H(l,l-1))
    IF ( s2==s ) EXIT
  END DO
  !     .......... FORM SHIFT ..........
  x = H(en,en)
  IF ( l==en ) THEN
    !     .......... ONE ROOT FOUND ..........
    H(en,en) = x + t
    Wr(en) = H(en,en)
    Wi(en) = 0.0E0
    en = na
  ELSE
    y = H(na,na)
    w = H(en,na)*H(na,en)
    IF ( l==na ) THEN
      !     .......... TWO ROOTS FOUND ..........
      p = (y-x)/2.0E0
      q = p*p + w
      zz = SQRT(ABS(q))
      H(en,en) = x + t
      x = H(en,en)
      H(na,na) = y + t
      IF ( q<0.0E0 ) THEN
        !     .......... COMPLEX PAIR ..........
        Wr(na) = x + p
        Wr(en) = x + p
        Wi(na) = zz
        Wi(en) = -zz
      ELSE
        !     .......... REAL PAIR ..........
        zz = p + SIGN(zz,p)
        Wr(na) = x + zz
        Wr(en) = Wr(na)
        IF ( zz/=0.0E0 ) Wr(en) = x - w/zz
        Wi(na) = 0.0E0
        Wi(en) = 0.0E0
        x = H(en,na)
        s = ABS(x) + ABS(zz)
        p = x/s
        q = zz/s
        r = SQRT(p*p+q*q)
        p = p/r
        q = q/r
        !     .......... ROW MODIFICATION ..........
        DO j = na, N
          zz = H(na,j)
          H(na,j) = q*zz + p*H(en,j)
          H(en,j) = q*H(en,j) - p*zz
        END DO
        !     .......... COLUMN MODIFICATION ..........
        DO i = 1, en
          zz = H(i,na)
          H(i,na) = q*zz + p*H(i,en)
          H(i,en) = q*H(i,en) - p*zz
        END DO
        !     .......... ACCUMULATE TRANSFORMATIONS ..........
        DO i = Low, Igh
          zz = Z(i,na)
          Z(i,na) = q*zz + p*Z(i,en)
          Z(i,en) = q*Z(i,en) - p*zz
          !
        END DO
      END IF
      en = enm2
    ELSEIF ( itn==0 ) THEN
      !     .......... SET ERROR -- NO CONVERGENCE TO AN
      !                EIGENVALUE AFTER 30*N ITERATIONS ..........
      Ierr = en
      RETURN
    ELSE
      IF ( its==10.OR.its==20 ) THEN
        !     .......... FORM EXCEPTIONAL SHIFT ..........
        t = t + x
        !
        DO i = Low, en
          H(i,i) = H(i,i) - x
        END DO
        !
        s = ABS(H(en,na)) + ABS(H(na,enm2))
        x = 0.75E0*s
        y = x
        w = -0.4375E0*s*s
      END IF
      its = its + 1
      itn = itn - 1
      !     .......... LOOK FOR TWO CONSECUTIVE SMALL
      !                SUB-DIAGONAL ELEMENTS.
      !                FOR M=EN-2 STEP -1 UNTIL L DO -- ..........
      DO mm = l, enm2
        m = enm2 + l - mm
        zz = H(m,m)
        r = x - zz
        s = y - zz
        p = (r*s-w)/H(m+1,m) + H(m,m+1)
        q = H(m+1,m+1) - zz - r - s
        r = H(m+2,m+1)
        s = ABS(p) + ABS(q) + ABS(r)
        p = p/s
        q = q/s
        r = r/s
        IF ( m==l ) EXIT
        s1 = ABS(p)*(ABS(H(m-1,m-1))+ABS(zz)+ABS(H(m+1,m+1)))
        s2 = s1 + ABS(H(m,m-1))*(ABS(q)+ABS(r))
        IF ( s2==s1 ) EXIT
      END DO
      !
      mp2 = m + 2
      !
      DO i = mp2, en
        H(i,i-2) = 0.0E0
        IF ( i/=mp2 ) H(i,i-3) = 0.0E0
      END DO
      !     .......... DOUBLE QR STEP INVOLVING ROWS L TO EN AND
      !                COLUMNS M TO EN ..........
      DO k = m, na
        notlas = k/=na
        IF ( k/=m ) THEN
          p = H(k,k-1)
          q = H(k+1,k-1)
          r = 0.0E0
          IF ( notlas ) r = H(k+2,k-1)
          x = ABS(p) + ABS(q) + ABS(r)
          IF ( x==0.0E0 ) CYCLE
          p = p/x
          q = q/x
          r = r/x
        END IF
        s = SIGN(SQRT(p*p+q*q+r*r),p)
        IF ( k==m ) THEN
          IF ( l/=m ) H(k,k-1) = -H(k,k-1)
        ELSE
          H(k,k-1) = -s*x
        END IF
        p = p + s
        x = p/s
        y = q/s
        zz = r/s
        q = q/p
        r = r/p
        !     .......... ROW MODIFICATION ..........
        DO j = k, N
          p = H(k,j) + q*H(k+1,j)
          IF ( notlas ) THEN
            p = p + r*H(k+2,j)
            H(k+2,j) = H(k+2,j) - p*zz
          END IF
          H(k+1,j) = H(k+1,j) - p*y
          H(k,j) = H(k,j) - p*x
        END DO
        !
        j = MIN(en,k+3)
        !     .......... COLUMN MODIFICATION ..........
        DO i = 1, j
          p = x*H(i,k) + y*H(i,k+1)
          IF ( notlas ) THEN
            p = p + zz*H(i,k+2)
            H(i,k+2) = H(i,k+2) - p*r
          END IF
          H(i,k+1) = H(i,k+1) - p*q
          H(i,k) = H(i,k) - p
        END DO
        !     .......... ACCUMULATE TRANSFORMATIONS ..........
        DO i = Low, Igh
          p = x*Z(i,k) + y*Z(i,k+1)
          IF ( notlas ) THEN
            p = p + zz*Z(i,k+2)
            Z(i,k+2) = Z(i,k+2) - p*r
          END IF
          Z(i,k+1) = Z(i,k+1) - p*q
          Z(i,k) = Z(i,k) - p
        END DO
        !
      END DO
      !
      GOTO 200
    END IF
  END IF
  GOTO 100
  RETURN
END SUBROUTINE HQR2
