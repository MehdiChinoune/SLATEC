!*==COMLR.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK COMLR
SUBROUTINE COMLR(Nm,N,Low,Igh,Hr,Hi,Wr,Wi,Ierr)
  IMPLICIT NONE
  !*--COMLR5
  !***BEGIN PROLOGUE  COMLR
  !***PURPOSE  Compute the eigenvalues of a complex upper Hessenberg
  !            matrix using the modified LR method.
  !***LIBRARY   SLATEC (EISPACK)
  !***CATEGORY  D4C2B
  !***TYPE      COMPLEX (COMLR-C)
  !***KEYWORDS  EIGENVALUES, EISPACK, LR METHOD
  !***AUTHOR  Smith, B. T., et al.
  !***DESCRIPTION
  !
  !     This subroutine is a translation of the ALGOL procedure COMLR,
  !     NUM. MATH. 12, 369-376(1968) by Martin and Wilkinson.
  !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 396-403(1971).
  !
  !     This subroutine finds the eigenvalues of a COMPLEX
  !     UPPER Hessenberg matrix by the modified LR method.
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
  !          triangles below the subdiagonal contain the multipliers
  !          which were used in the reduction by  COMHES, if performed.
  !          HR and HI are two-dimensional REAL arrays, dimensioned
  !          HR(NM,N) and HI(NM,N).
  !
  !     On OUTPUT
  !
  !        The upper Hessenberg portions of HR and HI have been
  !          destroyed.  Therefore, they must be saved before calling
  !          COMLR  if subsequent calculation of eigenvectors is to
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
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  COMLR
  !
  INTEGER i , j , l , m , N , en , ll , mm , Nm , Igh , im1 , itn , its , &
    Low , mp1 , enm1 , Ierr
  REAL Hr(Nm,*) , Hi(Nm,*) , Wr(*) , Wi(*)
  REAL si , sr , ti , tr , xi , xr , yi , yr , zzi , zzr , s1 , s2
  !
  !***FIRST EXECUTABLE STATEMENT  COMLR
  Ierr = 0
  !     .......... STORE ROOTS ISOLATED BY CBAL ..........
  DO i = 1 , N
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
  100  IF ( en<Low ) GOTO 99999
  its = 0
  enm1 = en - 1
  !     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
  !                FOR L=EN STEP -1 UNTIL LOW E0 -- ..........
  200  DO ll = Low , en
  l = en + Low - ll
  IF ( l==Low ) EXIT
  s1 = ABS(Hr(l-1,l-1)) + ABS(Hi(l-1,l-1)) + ABS(Hr(l,l)) + ABS(Hi(l,l))
  s2 = s1 + ABS(Hr(l,l-1)) + ABS(Hi(l,l-1))
  IF ( s2==s1 ) EXIT
ENDDO
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
  DO i = Low , en
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
  DO mm = l , enm1
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
  DO i = mp1 , en
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
      DO j = im1 , en
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
    DO j = i , en
      Hr(i,j) = Hr(i,j) - zzr*Hr(im1,j) + zzi*Hi(im1,j)
      Hi(i,j) = Hi(i,j) - zzr*Hi(im1,j) - zzi*Hr(im1,j)
    ENDDO
    !
  ENDDO
  !     .......... COMPOSITION R*L=H ..........
  DO j = mp1 , en
    xr = Hr(j,j-1)
    xi = Hi(j,j-1)
    Hr(j,j-1) = 0.0E0
    Hi(j,j-1) = 0.0E0
    !     .......... INTERCHANGE COLUMNS OF HR AND HI,
    !                IF NECESSARY ..........
    IF ( Wr(j)>0.0E0 ) THEN
      !
      DO i = l , j
        zzr = Hr(i,j-1)
        Hr(i,j-1) = Hr(i,j)
        Hr(i,j) = zzr
        zzi = Hi(i,j-1)
        Hi(i,j-1) = Hi(i,j)
        Hi(i,j) = zzi
      ENDDO
    ENDIF
    !
    DO i = l , j
      Hr(i,j-1) = Hr(i,j-1) + xr*Hr(i,j) - xi*Hi(i,j)
      Hi(i,j-1) = Hi(i,j-1) + xr*Hi(i,j) + xi*Hr(i,j)
    ENDDO
    !
  ENDDO
  !
  GOTO 200
ENDIF
99999 END SUBROUTINE COMLR
