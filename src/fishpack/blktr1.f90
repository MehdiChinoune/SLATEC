!DECK BLKTR1
SUBROUTINE BLKTR1(N,An,Bn,Cn,M,Am,Bm,Cm,Idimy,Y,B,W1,W2,W3,Wd,Ww,Wu,PRDCT,CPRDCT)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  BLKTR1
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to BLKTRI
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (BLKTR1-S, CBLKT1-C)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  ! BLKTR1 solves the linear system set up by BLKTRI.
  !
  ! B  contains the roots of all the B polynomials.
  ! W1,W2,W3,WD,WW,WU  are all working arrays.
  ! PRDCT  is either PRODP or PROD depending on whether the boundary
  ! conditions in the M direction are periodic or not.
  ! CPRDCT is either CPRODP or CPROD which are the complex versions
  ! of PRODP and PROD. These are called in the event that some
  ! of the roots of the B sub P polynomial are complex.
  !
  !***SEE ALSO  BLKTRI
  !***ROUTINES CALLED  INDXA, INDXB, INDXC
  !***COMMON BLOCKS    CBLKT
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  BLKTR1
  REAL Am, An, B, Bm, Bn, Cm, Cn, CNV, dum, EPS, W1, W2, W3, Wd, Wu, Ww, Y
  INTEGER i, i1, i2, i3, i4, Idimy, idxa, idxc, if, ifd, IK, &
    im1, im2, im3, imi1, imi2, ip, ip1, ip2, ip3
  INTEGER ipi1, ipi2, ipi3, ir, irm1, iz, izr, j, K, kdo, l, ll, &
    M, N, na, nc, NCMplx, NM, nm1, nm2
  INTEGER nm3, np, np1, np2, np3, NPP, nz
  DIMENSION An(*), Bn(*), Cn(*), Am(*), Bm(*), Cm(*), B(*), W1(*), &
    W2(*), W3(*), Wd(*), Ww(*), Wu(*), Y(Idimy,*)
  COMMON /CBLKT / NPP, K, EPS, CNV, NM, NCMplx, IK
  !***FIRST EXECUTABLE STATEMENT  BLKTR1
  kdo = K - 1
  DO l = 1, kdo
    ir = l - 1
    i2 = 2**ir
    i1 = i2/2
    i3 = i2 + i1
    i4 = i2 + i2
    irm1 = ir - 1
    CALL INDXB(i2,ir,im2,nm2)
    CALL INDXB(i1,irm1,im3,nm3)
    CALL INDXB(i3,irm1,im1,nm1)
    CALL PRDCT(nm2,B(im2),nm3,B(im3),nm1,B(im1),0,dum,Y(1,i2),W3,M,Am,Bm,Cm,&
      Wd,Ww,Wu)
    if = 2**K
    DO i = i4, if, i4
      IF ( i<=NM ) THEN
        ipi1 = i + i1
        ipi2 = i + i2
        ipi3 = i + i3
        CALL INDXC(i,ir,idxc,nc)
        IF ( i<if ) THEN
          CALL INDXA(i,ir,idxa,na)
          CALL INDXB(i-i1,irm1,im1,nm1)
          CALL INDXB(ipi2,ir,ip2,np2)
          CALL INDXB(ipi1,irm1,ip1,np1)
          CALL INDXB(ipi3,irm1,ip3,np3)
          CALL PRDCT(nm1,B(im1),0,dum,0,dum,na,An(idxa),W3,W1,M,Am,Bm,Cm,Wd,&
            Ww,Wu)
          IF ( ipi2<=NM ) THEN
            CALL PRDCT(np2,B(ip2),np1,B(ip1),np3,B(ip3),0,dum,Y(1,ipi2),W3,&
              M,Am,Bm,Cm,Wd,Ww,Wu)
            CALL PRDCT(np1,B(ip1),0,dum,0,dum,nc,Cn(idxc),W3,W2,M,Am,Bm,Cm,&
              Wd,Ww,Wu)
          ELSE
            DO j = 1, M
              W3(j) = 0.
              W2(j) = 0.
            ENDDO
          ENDIF
          DO j = 1, M
            Y(j,i) = W1(j) + W2(j) + Y(j,i)
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  IF ( NPP==0 ) THEN
    !
    !     THE PERIODIC CASE IS TREATED USING THE CAPACITANCE MATRIX METHOD
    !
    if = 2**K
    i = if/2
    i1 = i/2
    CALL INDXB(i-i1,K-2,im1,nm1)
    CALL INDXB(i+i1,K-2,ip1,np1)
    CALL INDXB(i,K-1,iz,nz)
    CALL PRDCT(nz,B(iz),nm1,B(im1),np1,B(ip1),0,dum,Y(1,i),W1,M,Am,Bm,Cm,Wd,&
      Ww,Wu)
    izr = i
    DO j = 1, M
      W2(j) = W1(j)
    ENDDO
    DO ll = 2, K
      l = K - ll + 1
      ir = l - 1
      i2 = 2**ir
      i1 = i2/2
      i = i2
      CALL INDXC(i,ir,idxc,nc)
      CALL INDXB(i,ir,iz,nz)
      CALL INDXB(i-i1,ir-1,im1,nm1)
      CALL INDXB(i+i1,ir-1,ip1,np1)
      CALL PRDCT(np1,B(ip1),0,dum,0,dum,nc,Cn(idxc),W1,W1,M,Am,Bm,Cm,Wd,Ww,&
        Wu)
      DO j = 1, M
        W1(j) = Y(j,i) + W1(j)
      ENDDO
      CALL PRDCT(nz,B(iz),nm1,B(im1),np1,B(ip1),0,dum,W1,W1,M,Am,Bm,Cm,Wd,&
        Ww,Wu)
    ENDDO
    DO ll = 2, K
      l = K - ll + 1
      ir = l - 1
      i2 = 2**ir
      i1 = i2/2
      i4 = i2 + i2
      ifd = if - i2
      DO i = i2, ifd, i4
        IF ( i-i2==izr ) THEN
          IF ( i>NM ) EXIT
          CALL INDXA(i,ir,idxa,na)
          CALL INDXB(i,ir,iz,nz)
          CALL INDXB(i-i1,ir-1,im1,nm1)
          CALL INDXB(i+i1,ir-1,ip1,np1)
          CALL PRDCT(nm1,B(im1),0,dum,0,dum,na,An(idxa),W2,W2,M,Am,Bm,Cm,Wd,&
            Ww,Wu)
          DO j = 1, M
            W2(j) = Y(j,i) + W2(j)
          ENDDO
          CALL PRDCT(nz,B(iz),nm1,B(im1),np1,B(ip1),0,dum,W2,W2,M,Am,Bm,Cm,&
            Wd,Ww,Wu)
          izr = i
          IF ( i==NM ) GOTO 50
        ENDIF
      ENDDO
    ENDDO
    50 CONTINUE
    DO j = 1, M
      Y(j,NM+1) = Y(j,NM+1) - Cn(NM+1)*W1(j) - An(NM+1)*W2(j)
    ENDDO
    CALL INDXB(if/2,K-1,im1,nm1)
    CALL INDXB(if,K-1,ip,np)
    IF ( NCMplx/=0 ) THEN
      CALL CPRDCT(NM+1,B(ip),nm1,B(im1),0,dum,0,dum,Y(1,NM+1),Y(1,NM+1),M,&
        Am,Bm,Cm,W1,W3,Ww)
    ELSE
      CALL PRDCT(NM+1,B(ip),nm1,B(im1),0,dum,0,dum,Y(1,NM+1),Y(1,NM+1),M,Am,&
        Bm,Cm,Wd,Ww,Wu)
    ENDIF
    DO j = 1, M
      W1(j) = An(1)*Y(j,NM+1)
      W2(j) = Cn(NM)*Y(j,NM+1)
      Y(j,1) = Y(j,1) - W1(j)
      Y(j,NM) = Y(j,NM) - W2(j)
    ENDDO
    DO l = 1, kdo
      ir = l - 1
      i2 = 2**ir
      i4 = i2 + i2
      i1 = i2/2
      i = i4
      CALL INDXA(i,ir,idxa,na)
      CALL INDXB(i-i2,ir,im2,nm2)
      CALL INDXB(i-i2-i1,ir-1,im3,nm3)
      CALL INDXB(i-i1,ir-1,im1,nm1)
      CALL PRDCT(nm2,B(im2),nm3,B(im3),nm1,B(im1),0,dum,W1,W1,M,Am,Bm,Cm,Wd,&
        Ww,Wu)
      CALL PRDCT(nm1,B(im1),0,dum,0,dum,na,An(idxa),W1,W1,M,Am,Bm,Cm,Wd,Ww,&
        Wu)
      DO j = 1, M
        Y(j,i) = Y(j,i) - W1(j)
      ENDDO
    ENDDO
    !
    izr = NM
    DO l = 1, kdo
      ir = l - 1
      i2 = 2**ir
      i1 = i2/2
      i3 = i2 + i1
      i4 = i2 + i2
      irm1 = ir - 1
      DO i = i4, if, i4
        ipi1 = i + i1
        ipi2 = i + i2
        ipi3 = i + i3
        IF ( ipi2==izr ) THEN
          CALL INDXC(i,ir,idxc,nc)
          CALL INDXB(ipi2,ir,ip2,np2)
          CALL INDXB(ipi1,irm1,ip1,np1)
          CALL INDXB(ipi3,irm1,ip3,np3)
          CALL PRDCT(np2,B(ip2),np1,B(ip1),np3,B(ip3),0,dum,W2,W2,M,Am,Bm,&
            Cm,Wd,Ww,Wu)
          CALL PRDCT(np1,B(ip1),0,dum,0,dum,nc,Cn(idxc),W2,W2,M,Am,Bm,Cm,Wd,&
            Ww,Wu)
          DO j = 1, M
            Y(j,i) = Y(j,i) - W2(j)
          ENDDO
          izr = i
          EXIT
        ELSEIF ( i==izr ) THEN
          EXIT
        ENDIF
      ENDDO
    ENDDO
  ENDIF
  !
  ! BEGIN BACK SUBSTITUTION PHASE
  !
  DO ll = 1, K
    l = K - ll + 1
    ir = l - 1
    irm1 = ir - 1
    i2 = 2**ir
    i1 = i2/2
    i4 = i2 + i2
    ifd = if - i2
    DO i = i2, ifd, i4
      IF ( i<=NM ) THEN
        imi1 = i - i1
        imi2 = i - i2
        ipi1 = i + i1
        ipi2 = i + i2
        CALL INDXA(i,ir,idxa,na)
        CALL INDXC(i,ir,idxc,nc)
        CALL INDXB(i,ir,iz,nz)
        CALL INDXB(imi1,irm1,im1,nm1)
        CALL INDXB(ipi1,irm1,ip1,np1)
        IF ( i<=i2 ) THEN
          DO j = 1, M
            W1(j) = 0.
          ENDDO
        ELSE
          CALL PRDCT(nm1,B(im1),0,dum,0,dum,na,An(idxa),Y(1,imi2),W1,M,Am,&
            Bm,Cm,Wd,Ww,Wu)
        ENDIF
        IF ( ipi2<=NM ) THEN
          CALL PRDCT(np1,B(ip1),0,dum,0,dum,nc,Cn(idxc),Y(1,ipi2),W2,M,Am,&
            Bm,Cm,Wd,Ww,Wu)
        ELSE
          DO j = 1, M
            W2(j) = 0.
          ENDDO
        ENDIF
        DO j = 1, M
          W1(j) = Y(j,i) + W1(j) + W2(j)
        ENDDO
        CALL PRDCT(nz,B(iz),nm1,B(im1),np1,B(ip1),0,dum,W1,Y(1,i),M,Am,Bm,&
          Cm,Wd,Ww,Wu)
      ENDIF
    ENDDO
  ENDDO
END SUBROUTINE BLKTR1