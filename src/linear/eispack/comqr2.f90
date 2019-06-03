!** COMQR2
SUBROUTINE COMQR2(Nm,N,Low,Igh,Ortr,Orti,Hr,Hi,Wr,Wi,Zr,Zi,Ierr)
  !>
  !  Compute the eigenvalues and eigenvectors of a complex upper
  !            Hessenberg matrix.
  !***
  ! **Library:**   SLATEC (EISPACK)
  !***
  ! **Category:**  D4C2B
  !***
  ! **Type:**      COMPLEX (HQR2-S, COMQR2-C)
  !***
  ! **Keywords:**  EIGENVALUES, EIGENVECTORS, EISPACK
  !***
  ! **Author:**  Smith, B. T., et al.
  !***
  ! **Description:**
  !
  !     This subroutine is a translation of a unitary analogue of the
  !     ALGOL procedure  COMLR2, NUM. MATH. 16, 181-204(1970) by Peters
  !     and Wilkinson.
  !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
  !     The unitary analogue substitutes the QR algorithm of Francis
  !     (COMP. JOUR. 4, 332-345(1962)) for the LR algorithm.
  !
  !     This subroutine finds the eigenvalues and eigenvectors
  !     of a COMPLEX UPPER Hessenberg matrix by the QR
  !     method.  The eigenvectors of a COMPLEX GENERAL matrix
  !     can also be found if  CORTH  has been used to reduce
  !     this general matrix to Hessenberg form.
  !
  !     On INPUT
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameters, HR, HI, ZR, and ZI, as declared in the
  !          calling program dimension statement.  NM is an INTEGER
  !          variable.
  !
  !        N is the order of the matrix H=(HR,HI).  N is an INTEGER
  !          variable.  N must be less than or equal to NM.
  !
  !        LOW and IGH are two INTEGER variables determined by the
  !          balancing subroutine  CBAL.  If  CBAL  has not been used,
  !          set LOW=1 and IGH equal to the order of the matrix, N.
  !
  !        ORTR and ORTI contain information about the unitary trans-
  !          formations used in the reduction by  CORTH, if performed.
  !          Only elements LOW through IGH are used.  If the eigenvectors
  !          of the Hessenberg matrix are desired, set ORTR(J) and
  !          ORTI(J) to 0.0E0 for these elements.  ORTR and ORTI are
  !          one-dimensional REAL arrays, dimensioned ORTR(IGH) and
  !          ORTI(IGH).
  !
  !        HR and HI contain the real and imaginary parts, respectively,
  !          of the complex upper Hessenberg matrix.  Their lower
  !          triangles below the subdiagonal contain information about
  !          the unitary transformations used in the reduction by  CORTH,
  !          if performed.  If the eigenvectors of the Hessenberg matrix
  !          are desired, these elements may be arbitrary.  HR and HI
  !          are two-dimensional REAL arrays, dimensioned HR(NM,N) and
  !          HI(NM,N).
  !
  !     On OUTPUT
  !
  !        ORTR, ORTI, and the upper Hessenberg portions of HR and HI
  !          have been destroyed.
  !
  !        WR and WI contain the real and imaginary parts, respectively,
  !          of the eigenvalues of the upper Hessenberg matrix.  If an
  !          error exit is made, the eigenvalues should be correct for
  !          indices IERR+1, IERR+2, ..., N.  WR and WI are one-
  !          dimensional REAL arrays, dimensioned WR(N) and WI(N).
  !
  !        ZR and ZI contain the real and imaginary parts, respectively,
  !          of the eigenvectors.  The eigenvectors are unnormalized.
  !          If an error exit is made, none of the eigenvectors has been
  !          found.  ZR and ZI are two-dimensional REAL arrays,
  !          dimensioned ZR(NM,N) and ZI(NM,N).
  !
  !        IERR is an INTEGER flag set to
  !          Zero       for normal return,
  !          J          if the J-th eigenvalue has not been
  !                     determined after a total of 30*N iterations.
  !                     The eigenvalues should be correct for indices
  !                     IERR+1, IERR+2, ..., N, but no eigenvectors are
  !                     computed.
  !
  !     Calls CSROOT for complex square root.
  !     Calls PYTHAG(A,B) for sqrt(A**2 + B**2).
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
  ! **Routines called:**  CDIV, CSROOT, PYTHAG

  !* REVISION HISTORY  (YYMMDD)
  !   760101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  !
  INTEGER i, j, k, l, m, N, en, ii, jj, ll, Nm, nn, Igh, ip1
  INTEGER itn, its, Low, lp1, enm1, iend, Ierr
  REAL(SP) Hr(Nm,*), Hi(Nm,*), Wr(*), Wi(*), Zr(Nm,*), Zi(Nm,*)
  REAL(SP) Ortr(*), Orti(*)
  REAL(SP) si, sr, ti, tr, xi, xr, yi, yr, zzi, zzr, norm, s1, s2
  !
  !* FIRST EXECUTABLE STATEMENT  COMQR2
  Ierr = 0
  !     .......... INITIALIZE EIGENVECTOR MATRIX ..........
  DO i = 1, N
    !
    DO j = 1, N
      Zr(i,j) = 0.0E0
      Zi(i,j) = 0.0E0
      IF ( i==j ) Zr(i,j) = 1.0E0
    END DO
  END DO
  !     .......... FORM THE MATRIX OF ACCUMULATED TRANSFORMATIONS
  !                FROM THE INFORMATION LEFT BY CORTH ..........
  iend = Igh - Low - 1
  IF ( iend<0 ) GOTO 100
  IF ( iend/=0 ) THEN
    !     .......... FOR I=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
    DO ii = 1, iend
      i = Igh - ii
      IF ( Ortr(i)/=0.0E0.OR.Orti(i)/=0.0E0 ) THEN
        IF ( Hr(i,i-1)/=0.0E0.OR.Hi(i,i-1)/=0.0E0 ) THEN
          !     .......... NORM BELOW IS NEGATIVE OF H FORMED IN CORTH ..........
          norm = Hr(i,i-1)*Ortr(i) + Hi(i,i-1)*Orti(i)
          ip1 = i + 1
          !
          DO k = ip1, Igh
            Ortr(k) = Hr(k,i-1)
            Orti(k) = Hi(k,i-1)
          END DO
          !
          DO j = i, Igh
            sr = 0.0E0
            si = 0.0E0
            !
            DO k = i, Igh
              sr = sr + Ortr(k)*Zr(k,j) + Orti(k)*Zi(k,j)
              si = si + Ortr(k)*Zi(k,j) - Orti(k)*Zr(k,j)
            END DO
            !
            sr = sr/norm
            si = si/norm
            !
            DO k = i, Igh
              Zr(k,j) = Zr(k,j) + sr*Ortr(k) - si*Orti(k)
              Zi(k,j) = Zi(k,j) + sr*Orti(k) + si*Ortr(k)
            END DO
            !
          END DO
        END IF
      END IF
      !
    END DO
  END IF
  !     .......... CREATE REAL SUBDIAGONAL ELEMENTS ..........
  l = Low + 1
  !
  DO i = l, Igh
    ll = MIN(i+1,Igh)
    IF ( Hi(i,i-1)/=0.0E0 ) THEN
      norm = PYTHAG(Hr(i,i-1),Hi(i,i-1))
      yr = Hr(i,i-1)/norm
      yi = Hi(i,i-1)/norm
      Hr(i,i-1) = norm
      Hi(i,i-1) = 0.0E0
      !
      DO j = i, N
        si = yr*Hi(i,j) - yi*Hr(i,j)
        Hr(i,j) = yr*Hr(i,j) + yi*Hi(i,j)
        Hi(i,j) = si
      END DO
      !
      DO j = 1, ll
        si = yr*Hi(j,i) + yi*Hr(j,i)
        Hr(j,i) = yr*Hr(j,i) - yi*Hi(j,i)
        Hi(j,i) = si
      END DO
      !
      DO j = Low, Igh
        si = yr*Zi(j,i) + yi*Zr(j,i)
        Zr(j,i) = yr*Zr(j,i) - yi*Zi(j,i)
        Zi(j,i) = si
      END DO
    END IF
    !
  END DO
  !     .......... STORE ROOTS ISOLATED BY CBAL ..........
  100 CONTINUE
  DO i = 1, N
    IF ( i<Low.OR.i>Igh ) THEN
      Wr(i) = Hr(i,i)
      Wi(i) = Hi(i,i)
    END IF
  END DO
  !
  en = Igh
  tr = 0.0E0
  ti = 0.0E0
  itn = 30*N
  !     .......... SEARCH FOR NEXT EIGENVALUE ..........
  200 CONTINUE
  IF ( en<Low ) THEN
    !     .......... ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND
    !                VECTORS OF UPPER TRIANGULAR FORM ..........
    norm = 0.0E0
    !
    DO i = 1, N
      !
      DO j = i, N
        norm = norm + ABS(Hr(i,j)) + ABS(Hi(i,j))
      END DO
    END DO
    !
    IF ( N/=1.AND.norm/=0.0E0 ) THEN
      !     .......... FOR EN=N STEP -1 UNTIL 2 DO -- ..........
      DO nn = 2, N
        en = N + 2 - nn
        xr = Wr(en)
        xi = Wi(en)
        enm1 = en - 1
        !     .......... FOR I=EN-1 STEP -1 UNTIL 1 DO -- ..........
        DO ii = 1, enm1
          i = en - ii
          zzr = Hr(i,en)
          zzi = Hi(i,en)
          IF ( i/=enm1 ) THEN
            ip1 = i + 1
            !
            DO j = ip1, enm1
              zzr = zzr + Hr(i,j)*Hr(j,en) - Hi(i,j)*Hi(j,en)
              zzi = zzi + Hr(i,j)*Hi(j,en) + Hi(i,j)*Hr(j,en)
            END DO
          END IF
          !
          yr = xr - Wr(i)
          yi = xi - Wi(i)
          IF ( yr==0.0E0.AND.yi==0.0E0 ) THEN
            yr = norm
            DO
              yr = 0.5E0*yr
              IF ( norm+yr<=norm ) THEN
                yr = 2.0E0*yr
                EXIT
              END IF
            END DO
          END IF
          CALL CDIV(zzr,zzi,yr,yi,Hr(i,en),Hi(i,en))
        END DO
        !
      END DO
      !     .......... END BACKSUBSTITUTION ..........
      enm1 = N - 1
      !     .......... VECTORS OF ISOLATED ROOTS ..........
      DO i = 1, enm1
        IF ( i<Low.OR.i>Igh ) THEN
          ip1 = i + 1
          !
          DO j = ip1, N
            Zr(i,j) = Hr(i,j)
            Zi(i,j) = Hi(i,j)
          END DO
        END IF
        !
      END DO
      !     .......... MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
      !                VECTORS OF ORIGINAL FULL MATRIX.
      !                FOR J=N STEP -1 UNTIL LOW+1 DO -- ..........
      DO jj = Low, enm1
        j = N + Low - jj
        m = MIN(j-1,Igh)
        !
        DO i = Low, Igh
          zzr = Zr(i,j)
          zzi = Zi(i,j)
          !
          DO k = Low, m
            zzr = zzr + Zr(i,k)*Hr(k,j) - Zi(i,k)*Hi(k,j)
            zzi = zzi + Zr(i,k)*Hi(k,j) + Zi(i,k)*Hr(k,j)
          END DO
          !
          Zr(i,j) = zzr
          Zi(i,j) = zzi
        END DO
        !
      END DO
    END IF
    RETURN
  ELSE
    its = 0
    enm1 = en - 1
  END IF
  !     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
  !                FOR L=EN STEP -1 UNTIL LOW DO -- ..........
  300 CONTINUE
  DO ll = Low, en
    l = en + Low - ll
    IF ( l==Low ) EXIT
    s1 = ABS(Hr(l-1,l-1)) + ABS(Hi(l-1,l-1)) + ABS(Hr(l,l)) + ABS(Hi(l,l))
    s2 = s1 + ABS(Hr(l,l-1))
    IF ( s2==s1 ) EXIT
  END DO
  !     .......... FORM SHIFT ..........
  IF ( l==en ) THEN
    !     .......... A ROOT FOUND ..........
    Hr(en,en) = Hr(en,en) + tr
    Wr(en) = Hr(en,en)
    Hi(en,en) = Hi(en,en) + ti
    Wi(en) = Hi(en,en)
    en = enm1
    GOTO 200
  ELSEIF ( itn==0 ) THEN
    !     .......... SET ERROR -- NO CONVERGENCE TO AN
    !                EIGENVALUE AFTER 30*N ITERATIONS ..........
    Ierr = en
  ELSE
    IF ( its==10.OR.its==20 ) THEN
      !     .......... FORM EXCEPTIONAL SHIFT ..........
      sr = ABS(Hr(en,enm1)) + ABS(Hr(enm1,en-2))
      si = 0.0E0
    ELSE
      sr = Hr(en,en)
      si = Hi(en,en)
      xr = Hr(enm1,en)*Hr(en,enm1)
      xi = Hi(enm1,en)*Hr(en,enm1)
      IF ( xr/=0.0E0.OR.xi/=0.0E0 ) THEN
        yr = (Hr(enm1,enm1)-sr)/2.0E0
        yi = (Hi(enm1,enm1)-si)/2.0E0
        CALL CSROOT(yr**2-yi**2+xr,2.0E0*yr*yi+xi,zzr,zzi)
        IF ( yr*zzr+yi*zzi<0.0E0 ) THEN
          zzr = -zzr
          zzi = -zzi
        END IF
        CALL CDIV(xr,xi,yr+zzr,yi+zzi,xr,xi)
        sr = sr - xr
        si = si - xi
      END IF
    END IF
    !
    DO i = Low, en
      Hr(i,i) = Hr(i,i) - sr
      Hi(i,i) = Hi(i,i) - si
    END DO
    !
    tr = tr + sr
    ti = ti + si
    its = its + 1
    itn = itn - 1
    !     .......... REDUCE TO TRIANGLE (ROWS) ..........
    lp1 = l + 1
    !
    DO i = lp1, en
      sr = Hr(i,i-1)
      Hr(i,i-1) = 0.0E0
      norm = PYTHAG(PYTHAG(Hr(i-1,i-1),Hi(i-1,i-1)),sr)
      xr = Hr(i-1,i-1)/norm
      Wr(i-1) = xr
      xi = Hi(i-1,i-1)/norm
      Wi(i-1) = xi
      Hr(i-1,i-1) = norm
      Hi(i-1,i-1) = 0.0E0
      Hi(i,i-1) = sr/norm
      !
      DO j = i, N
        yr = Hr(i-1,j)
        yi = Hi(i-1,j)
        zzr = Hr(i,j)
        zzi = Hi(i,j)
        Hr(i-1,j) = xr*yr + xi*yi + Hi(i,i-1)*zzr
        Hi(i-1,j) = xr*yi - xi*yr + Hi(i,i-1)*zzi
        Hr(i,j) = xr*zzr - xi*zzi - Hi(i,i-1)*yr
        Hi(i,j) = xr*zzi + xi*zzr - Hi(i,i-1)*yi
      END DO
      !
    END DO
    !
    si = Hi(en,en)
    IF ( si/=0.0E0 ) THEN
      norm = PYTHAG(Hr(en,en),si)
      sr = Hr(en,en)/norm
      si = si/norm
      Hr(en,en) = norm
      Hi(en,en) = 0.0E0
      IF ( en/=N ) THEN
        ip1 = en + 1
        !
        DO j = ip1, N
          yr = Hr(en,j)
          yi = Hi(en,j)
          Hr(en,j) = sr*yr + si*yi
          Hi(en,j) = sr*yi - si*yr
        END DO
      END IF
    END IF
    !     .......... INVERSE OPERATION (COLUMNS) ..........
    DO j = lp1, en
      xr = Wr(j-1)
      xi = Wi(j-1)
      !
      DO i = 1, j
        yr = Hr(i,j-1)
        yi = 0.0E0
        zzr = Hr(i,j)
        zzi = Hi(i,j)
        IF ( i/=j ) THEN
          yi = Hi(i,j-1)
          Hi(i,j-1) = xr*yi + xi*yr + Hi(j,j-1)*zzi
        END IF
        Hr(i,j-1) = xr*yr - xi*yi + Hi(j,j-1)*zzr
        Hr(i,j) = xr*zzr + xi*zzi - Hi(j,j-1)*yr
        Hi(i,j) = xr*zzi - xi*zzr - Hi(j,j-1)*yi
      END DO
      !
      DO i = Low, Igh
        yr = Zr(i,j-1)
        yi = Zi(i,j-1)
        zzr = Zr(i,j)
        zzi = Zi(i,j)
        Zr(i,j-1) = xr*yr - xi*yi + Hi(j,j-1)*zzr
        Zi(i,j-1) = xr*yi + xi*yr + Hi(j,j-1)*zzi
        Zr(i,j) = xr*zzr + xi*zzi - Hi(j,j-1)*yr
        Zi(i,j) = xr*zzi - xi*zzr - Hi(j,j-1)*yi
      END DO
      !
    END DO
    !
    IF ( si/=0.0E0 ) THEN
      !
      DO i = 1, en
        yr = Hr(i,en)
        yi = Hi(i,en)
        Hr(i,en) = sr*yr - si*yi
        Hi(i,en) = sr*yi + si*yr
      END DO
      !
      DO i = Low, Igh
        yr = Zr(i,en)
        yi = Zi(i,en)
        Zr(i,en) = sr*yr - si*yi
        Zi(i,en) = sr*yi + si*yr
        !
      END DO
    END IF
    GOTO 300
  END IF
  RETURN
END SUBROUTINE COMQR2
