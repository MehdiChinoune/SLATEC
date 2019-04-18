!** COMQR
SUBROUTINE COMQR(Nm,N,Low,Igh,Hr,Hi,Wr,Wi,Ierr)
  !>
  !  Compute the eigenvalues of complex upper Hessenberg matrix
  !            using the QR method.
  !***
  ! **Library:**   SLATEC (EISPACK)
  !***
  ! **Category:**  D4C2B
  !***
  ! **Type:**      COMPLEX (HQR-S, COMQR-C)
  !***
  ! **Keywords:**  EIGENVALUES, EIGENVECTORS, EISPACK
  !***
  ! **Author:**  Smith, B. T., et al.
  !***
  ! **Description:**
  !
  !     This subroutine is a translation of a unitary analogue of the
  !     ALGOL procedure  COMLR, NUM. MATH. 12, 369-376(1968) by Martin
  !     and Wilkinson.
  !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 396-403(1971).
  !     The unitary analogue substitutes the QR algorithm of Francis
  !     (COMP. JOUR. 4, 332-345(1962)) for the LR algorithm.
  !
  !     This subroutine finds the eigenvalues of a COMPLEX
  !     upper Hessenberg matrix by the QR method.
  !
  !     On INPUT
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameters, HR and HI, as declared in the calling
  !          program dimension statement.  NM is an INTEGER variable.
  !
  !        N is the order of the matrix H=(HR,HI).  N is an INTEGER
  !          variable.  N must be less than or equal to NM.
  !
  !        LOW and IGH are two INTEGER variables determined by the
  !          balancing subroutine  CBAL.  If  CBAL  has not been used,
  !          set LOW=1 and IGH equal to the order of the matrix, N.
  !
  !        HR and HI contain the real and imaginary parts, respectively,
  !          of the complex upper Hessenberg matrix.  Their lower
  !          triangles below the subdiagonal contain information about
  !          the unitary transformations used in the reduction by  CORTH,
  !          if performed.  HR and HI are two-dimensional REAL arrays,
  !          dimensioned HR(NM,N) and HI(NM,N).
  !
  !     On OUTPUT
  !
  !        The upper Hessenberg portions of HR and HI have been
  !          destroyed.  Therefore, they must be saved before calling
  !          COMQR  if subsequent calculation of eigenvectors is to
  !          be performed.
  !
  !        WR and WI contain the real and imaginary parts, respectively,
  !          of the eigenvalues of the upper Hessenberg matrix.  If an
  !          error exit is made, the eigenvalues should be correct for
  !          indices IERR+1, IERR+2, ..., N.  WR and WI are one-
  !          dimensional REAL arrays, dimensioned WR(N) and WI(N).
  !
  !        IERR is an INTEGER flag set to
  !          Zero       for normal return,
  !          J          if the J-th eigenvalue has not been
  !                     determined after a total of 30*N iterations.
  !                     The eigenvalues should be correct for indices
  !                     IERR+1, IERR+2, ..., N.
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
  INTEGER i, j, l, N, en, ll, Nm, Igh, itn, its, Low, lp1, enm1, Ierr
  REAL Hr(Nm,*), Hi(Nm,*), Wr(*), Wi(*)
  REAL si, sr, ti, tr, xi, xr, yi, yr, zzi, zzr, norm, s1, s2
  !
  !* FIRST EXECUTABLE STATEMENT  COMQR
  Ierr = 0
  IF ( Low/=Igh ) THEN
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
        DO j = i, Igh
          si = yr*Hi(i,j) - yi*Hr(i,j)
          Hr(i,j) = yr*Hr(i,j) + yi*Hi(i,j)
          Hi(i,j) = si
        END DO
        !
        DO j = Low, ll
          si = yr*Hi(j,i) + yi*Hr(j,i)
          Hr(j,i) = yr*Hr(j,i) - yi*Hi(j,i)
          Hi(j,i) = si
        END DO
      END IF
      !
    END DO
  END IF
  !     .......... STORE ROOTS ISOLATED BY CBAL ..........
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
  100 CONTINUE
  IF ( en<Low ) RETURN
  its = 0
  enm1 = en - 1
  !     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
  !                FOR L=EN STEP -1 UNTIL LOW E0 -- ..........
  200 CONTINUE
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
    Wr(en) = Hr(en,en) + tr
    Wi(en) = Hi(en,en) + ti
    en = enm1
    GOTO 100
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
      DO j = i, en
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
    END IF
    !     .......... INVERSE OPERATION (COLUMNS) ..........
    DO j = lp1, en
      xr = Wr(j-1)
      xi = Wi(j-1)
      !
      DO i = l, j
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
    END DO
    !
    IF ( si/=0.0E0 ) THEN
      !
      DO i = l, en
        yr = Hr(i,en)
        yi = Hi(i,en)
        Hr(i,en) = sr*yr - si*yi
        Hi(i,en) = sr*yi + si*yr
        !
      END DO
    END IF
    GOTO 200
  END IF
  RETURN
END SUBROUTINE COMQR
