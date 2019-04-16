!** HQR
SUBROUTINE HQR(Nm,N,Low,Igh,H,Wr,Wi,Ierr)
  !>
  !***
  !  Compute the eigenvalues of a real upper Hessenberg matrix
  !            using the QR method.
  !***
  ! **Library:**   SLATEC (EISPACK)
  !***
  ! **Category:**  D4C2B
  !***
  ! **Type:**      SINGLE PRECISION (HQR-S, COMQR-C)
  !***
  ! **Keywords:**  EIGENVALUES, EIGENVECTORS, EISPACK
  !***
  ! **Author:**  Smith, B. T., et al.
  !***
  ! **Description:**
  !
  !     This subroutine is a translation of the ALGOL procedure HQR,
  !     NUM. MATH. 14, 219-231(1970) by Martin, Peters, and Wilkinson.
  !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 359-371(1971).
  !
  !     This subroutine finds the eigenvalues of a REAL
  !     UPPER Hessenberg matrix by the QR method.
  !
  !     On INPUT
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameter, H, as declared in the calling program
  !          dimension statement.  NM is an INTEGER variable.
  !
  !        N is the order of the matrix H.  N is an INTEGER variable.
  !          N must be less than or equal to NM.
  !
  !        LOW and IGH are two INTEGER variables determined by the
  !          balancing subroutine  BALANC.  If  BALANC  has not been
  !          used, set LOW=1 and IGH equal to the order of the matrix, N.
  !
  !        H contains the upper Hessenberg matrix.  Information about
  !          the transformations used in the reduction to Hessenberg
  !          form by  ELMHES  or  ORTHES, if performed, is stored
  !          in the remaining triangle under the Hessenberg matrix.
  !          H is a two-dimensional REAL array, dimensioned H(NM,N).
  !
  !     On OUTPUT
  !
  !        H has been destroyed.  Therefore, it must be saved before
  !          calling  HQR  if subsequent calculation and back
  !          transformation of eigenvectors is to be performed.
  !
  !        WR and WI contain the real and imaginary parts, respectively,
  !          of the eigenvalues.  The eigenvalues are unordered except
  !          that complex conjugate pairs of values appear consecutively
  !          with the eigenvalue having the positive imaginary part first.
  !          If an error exit is made, the eigenvalues should be correct
  !          for indices IERR+1, IERR+2, ..., N.  WR and WI are one-
  !          dimensional REAL arrays, dimensioned WR(N) and WI(N).
  !
  !        IERR is an INTEGER flag set to
  !          Zero       for normal return,
  !          J          if the J-th eigenvalue has not been
  !                     determined after a total of 30*N iterations.
  !                     The eigenvalues should be correct for indices
  !                     IERR+1, IERR+2, ..., N.
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
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   760101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  !
  INTEGER i, j, k, l, m, N, en, ll, mm, na, Nm, Igh, itn, its, &
    Low, mp2, enm2, Ierr
  REAL H(Nm,*), Wr(*), Wi(*)
  REAL p, q, r, s, t, w, x, y, zz, norm, s1, s2
  LOGICAL notlas
  !
  !* FIRST EXECUTABLE STATEMENT  HQR
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
  IF ( en<Low ) RETURN
  its = 0
  na = en - 1
  enm2 = na - 1
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
    Wr(en) = x + t
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
      x = x + t
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
        DO j = k, en
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
        DO i = l, j
          p = x*H(i,k) + y*H(i,k+1)
          IF ( notlas ) THEN
            p = p + zz*H(i,k+2)
            H(i,k+2) = H(i,k+2) - p*r
          END IF
          H(i,k+1) = H(i,k+1) - p*q
          H(i,k) = H(i,k) - p
        END DO
        !
      END DO
      !
      GOTO 200
    END IF
  END IF
  GOTO 100
  RETURN
END SUBROUTINE HQR
