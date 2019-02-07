!*==COMLR2.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK COMLR2
SUBROUTINE COMLR2(Nm,N,Low,Igh,Int,Hr,Hi,Wr,Wi,Zr,Zi,Ierr)
  IMPLICIT NONE
  !*--COMLR25
  !***BEGIN PROLOGUE  COMLR2
  !***PURPOSE  Compute the eigenvalues and eigenvectors of a complex upper
  !            Hessenberg matrix using the modified LR method.
  !***LIBRARY   SLATEC (EISPACK)
  !***CATEGORY  D4C2B
  !***TYPE      COMPLEX (COMLR2-C)
  !***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK, LR METHOD
  !***AUTHOR  Smith, B. T., et al.
  !***DESCRIPTION
  !
  !     This subroutine is a translation of the ALGOL procedure COMLR2,
  !     NUM. MATH. 16, 181-204(1970) by Peters and Wilkinson.
  !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
  !
  !     This subroutine finds the eigenvalues and eigenvectors
  !     of a COMPLEX UPPER Hessenberg matrix by the modified LR
  !     method.  The eigenvectors of a COMPLEX GENERAL matrix
  !     can also be found if  COMHES  has been used to reduce
  !     this general matrix to Hessenberg form.
  !
  !     On INPUT
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameters, HR, HI, ZR and ZI, as declared in the
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
  !        INT contains information on the rows and columns
  !          interchanged in the reduction by  COMHES, if performed.
  !          Only elements LOW through IGH are used.  If you want the
  !          eigenvectors of a complex general matrix, leave INT as it
  !          came from  COMHES.  If the eigenvectors of the Hessenberg
  !          matrix are desired, set INT(J)=J for these elements.  INT
  !          is a one-dimensional INTEGER array, dimensioned INT(IGH).
  !
  !        HR and HI contain the real and imaginary parts, respectively,
  !          of the complex upper Hessenberg matrix.  Their lower
  !          triangles below the subdiagonal contain the multipliers
  !          which were used in the reduction by  COMHES, if performed.
  !          If the eigenvectors of a complex general matrix are
  !          desired, leave these multipliers in the lower triangles.
  !          If the eigenvectors of the Hessenberg matrix are desired,
  !          these elements must be set to zero.  HR and HI are
  !          two-dimensional REAL arrays, dimensioned HR(NM,N) and
  !          HI(NM,N).
  !
  !     On OUTPUT
  !
  !        The upper Hessenberg portions of HR and HI have been
  !          destroyed, but the location HR(1,1) contains the norm
  !          of the triangularized matrix.
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
  !     Calls CDIV for complex division.
  !
  !     Questions and comments should be directed to B. S. Garbow,
  !     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
  !     ------------------------------------------------------------------
  !
  !***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
  !                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
  !                 system Routines - EISPACK Guide, Springer-Verlag,
  !                 1976.
  !***ROUTINES CALLED  CDIV, CSROOT
  !***REVISION HISTORY  (YYMMDD)
  !   760101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  COMLR2
  !
  INTEGER i, j, k, l, m, N, en, ii, jj, ll, mm, Nm, nn, Igh, &
    im1, ip1
  INTEGER itn, its, Low, mp1, enm1, iend, Ierr
  REAL Hr(Nm,*), Hi(Nm,*), Wr(*), Wi(*), Zr(Nm,*), Zi(Nm,*)
  REAL si, sr, ti, tr, xi, xr, yi, yr, zzi, zzr, norm, s1, s2
  INTEGER Int(*)
  !
  !***FIRST EXECUTABLE STATEMENT  COMLR2
  Ierr = 0
  !     .......... INITIALIZE EIGENVECTOR MATRIX ..........
  DO i = 1, N
    !
    DO j = 1, N
      Zr(i,j) = 0.0E0
      Zi(i,j) = 0.0E0
      IF ( i==j ) Zr(i,j) = 1.0E0
    ENDDO
  ENDDO
  !     .......... FORM THE MATRIX OF ACCUMULATED TRANSFORMATIONS
  !                FROM THE INFORMATION LEFT BY COMHES ..........
  iend = Igh - Low - 1
  IF ( iend>0 ) THEN
    !     .......... FOR I=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
    DO ii = 1, iend
      i = Igh - ii
      ip1 = i + 1
      !
      DO k = ip1, Igh
        Zr(k,i) = Hr(k,i-1)
        Zi(k,i) = Hi(k,i-1)
      ENDDO
      !
      j = Int(i)
      IF ( i/=j ) THEN
        !
        DO k = i, Igh
          Zr(i,k) = Zr(j,k)
          Zi(i,k) = Zi(j,k)
          Zr(j,k) = 0.0E0
          Zi(j,k) = 0.0E0
        ENDDO
        !
        Zr(j,i) = 1.0E0
      ENDIF
    ENDDO
  ENDIF
  !     .......... STORE ROOTS ISOLATED BY CBAL ..........
  DO i = 1, N
    IF ( i<Low.OR.i>Igh ) THEN
      Wr(i) = Hr(i,i)
      Wi(i) = Hi(i,i)
    ENDIF
  ENDDO
  !
  en = Igh
  tr = 0.0E0
  ti = 0.0E0
  itn = 30*N
  !     .......... SEARCH FOR NEXT EIGENVALUE ..........
  100 CONTINUE
  IF ( en<Low ) THEN
    !     .......... ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND
    !                VECTORS OF UPPER TRIANGULAR FORM ..........
    norm = 0.0E0
    !
    DO i = 1, N
      !
      DO j = i, N
        norm = norm + ABS(Hr(i,j)) + ABS(Hi(i,j))
      ENDDO
    ENDDO
    !
    Hr(1,1) = norm
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
            ENDDO
          ENDIF
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
              ENDIF
            ENDDO
          ENDIF
          CALL CDIV(zzr,zzi,yr,yi,Hr(i,en),Hi(i,en))
        ENDDO
        !
      ENDDO
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
          ENDDO
        ENDIF
        !
      ENDDO
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
          ENDDO
          !
          Zr(i,j) = zzr
          Zi(i,j) = zzi
        ENDDO
        !
      ENDDO
    ENDIF
    GOTO 99999
  ELSE
    its = 0
    enm1 = en - 1
  ENDIF
  !     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
  !                FOR L=EN STEP -1 UNTIL LOW DO -- ..........
  200 CONTINUE
  DO ll = Low, en
    l = en + Low - ll
    IF ( l==Low ) EXIT
    s1 = ABS(Hr(l-1,l-1)) + ABS(Hi(l-1,l-1)) + ABS(Hr(l,l)) + ABS(Hi(l,l))
    s2 = s1 + ABS(Hr(l,l-1)) + ABS(Hi(l,l-1))
    IF ( s2==s1 ) EXIT
  ENDDO
  !     .......... FORM SHIFT ..........
  IF ( l==en ) THEN
    !     .......... A ROOT FOUND ..........
    Hr(en,en) = Hr(en,en) + tr
    Wr(en) = Hr(en,en)
    Hi(en,en) = Hi(en,en) + ti
    Wi(en) = Hi(en,en)
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
      si = ABS(Hi(en,enm1)) + ABS(Hi(enm1,en-2))
    ELSE
      sr = Hr(en,en)
      si = Hi(en,en)
      xr = Hr(enm1,en)*Hr(en,enm1) - Hi(enm1,en)*Hi(en,enm1)
      xi = Hr(enm1,en)*Hi(en,enm1) + Hi(enm1,en)*Hr(en,enm1)
      IF ( xr/=0.0E0.OR.xi/=0.0E0 ) THEN
        yr = (Hr(enm1,enm1)-sr)/2.0E0
        yi = (Hi(enm1,enm1)-si)/2.0E0
        CALL CSROOT(yr**2-yi**2+xr,2.0E0*yr*yi+xi,zzr,zzi)
        IF ( yr*zzr+yi*zzi<0.0E0 ) THEN
          zzr = -zzr
          zzi = -zzi
        ENDIF
        CALL CDIV(xr,xi,yr+zzr,yi+zzi,xr,xi)
        sr = sr - xr
        si = si - xi
      ENDIF
    ENDIF
    !
    DO i = Low, en
      Hr(i,i) = Hr(i,i) - sr
      Hi(i,i) = Hi(i,i) - si
    ENDDO
    !
    tr = tr + sr
    ti = ti + si
    its = its + 1
    itn = itn - 1
    !     .......... LOOK FOR TWO CONSECUTIVE SMALL
    !                SUB-DIAGONAL ELEMENTS ..........
    xr = ABS(Hr(enm1,enm1)) + ABS(Hi(enm1,enm1))
    yr = ABS(Hr(en,enm1)) + ABS(Hi(en,enm1))
    zzr = ABS(Hr(en,en)) + ABS(Hi(en,en))
    !     .......... FOR M=EN-1 STEP -1 UNTIL L DO -- ..........
    DO mm = l, enm1
      m = enm1 + l - mm
      IF ( m==l ) EXIT
      yi = yr
      yr = ABS(Hr(m,m-1)) + ABS(Hi(m,m-1))
      xi = zzr
      zzr = xr
      xr = ABS(Hr(m-1,m-1)) + ABS(Hi(m-1,m-1))
      s1 = zzr/yi*(zzr+xr+xi)
      s2 = s1 + yr
      IF ( s2==s1 ) EXIT
    ENDDO
    !     .......... TRIANGULAR DECOMPOSITION H=L*R ..........
    mp1 = m + 1
    !
    DO i = mp1, en
      im1 = i - 1
      xr = Hr(im1,im1)
      xi = Hi(im1,im1)
      yr = Hr(i,im1)
      yi = Hi(i,im1)
      IF ( ABS(xr)+ABS(xi)>=ABS(yr)+ABS(yi) ) THEN
        CALL CDIV(yr,yi,xr,xi,zzr,zzi)
        Wr(i) = -1.0E0
      ELSE
        !     .......... INTERCHANGE ROWS OF HR AND HI ..........
        DO j = im1, N
          zzr = Hr(im1,j)
          Hr(im1,j) = Hr(i,j)
          Hr(i,j) = zzr
          zzi = Hi(im1,j)
          Hi(im1,j) = Hi(i,j)
          Hi(i,j) = zzi
        ENDDO
        !
        CALL CDIV(xr,xi,yr,yi,zzr,zzi)
        Wr(i) = 1.0E0
      ENDIF
      Hr(i,im1) = zzr
      Hi(i,im1) = zzi
      !
      DO j = i, N
        Hr(i,j) = Hr(i,j) - zzr*Hr(im1,j) + zzi*Hi(im1,j)
        Hi(i,j) = Hi(i,j) - zzr*Hi(im1,j) - zzi*Hr(im1,j)
      ENDDO
      !
    ENDDO
    !     .......... COMPOSITION R*L=H ..........
    DO j = mp1, en
      xr = Hr(j,j-1)
      xi = Hi(j,j-1)
      Hr(j,j-1) = 0.0E0
      Hi(j,j-1) = 0.0E0
      !     .......... INTERCHANGE COLUMNS OF HR, HI, ZR, AND ZI,
      !                IF NECESSARY ..........
      IF ( Wr(j)>0.0E0 ) THEN
        !
        DO i = 1, j
          zzr = Hr(i,j-1)
          Hr(i,j-1) = Hr(i,j)
          Hr(i,j) = zzr
          zzi = Hi(i,j-1)
          Hi(i,j-1) = Hi(i,j)
          Hi(i,j) = zzi
        ENDDO
        !
        DO i = Low, Igh
          zzr = Zr(i,j-1)
          Zr(i,j-1) = Zr(i,j)
          Zr(i,j) = zzr
          zzi = Zi(i,j-1)
          Zi(i,j-1) = Zi(i,j)
          Zi(i,j) = zzi
        ENDDO
      ENDIF
      !
      DO i = 1, j
        Hr(i,j-1) = Hr(i,j-1) + xr*Hr(i,j) - xi*Hi(i,j)
        Hi(i,j-1) = Hi(i,j-1) + xr*Hi(i,j) + xi*Hr(i,j)
      ENDDO
      !     .......... ACCUMULATE TRANSFORMATIONS ..........
      DO i = Low, Igh
        Zr(i,j-1) = Zr(i,j-1) + xr*Zr(i,j) - xi*Zi(i,j)
        Zi(i,j-1) = Zi(i,j-1) + xr*Zi(i,j) + xi*Zr(i,j)
      ENDDO
      !
    ENDDO
    !
    GOTO 200
  ENDIF
  99999 CONTINUE
  END SUBROUTINE COMLR2
