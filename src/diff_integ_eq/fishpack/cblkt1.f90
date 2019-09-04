!** CBLKT1
PURE SUBROUTINE CBLKT1(An,Cn,M,Am,Bm,Cm,Idimy,Y,B,W1,W2,W3,Wd,Ww,Wu,PRDCT,CPRDCT)
  !> Subsidiary to CBLKTR
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      COMPLEX (BLKTR1-S, CBLKT1-C)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  ! CBLKT1 solves the linear system of routine CBLKTR.
  !
  ! B  contains the roots of all the B polynomials.
  ! W1,W2,W3,WD,WW,WU  are all working arrays.
  ! PRDCT is either PROCP or PROC depending on whether the boundary
  ! conditions in the M direction are periodic or not.
  ! CPRDCT is either CPROCP or CPROC which are called if some of the zeros
  ! of the B polynomials are complex.
  !
  !***
  ! **See also:**  CBLKTR
  !***
  ! **Routines called:**  INXCA, INXCB, INXCC
  !***
  ! COMMON BLOCKS    CCBLK

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  USE CCBLK, ONLY : k_com, ncmplx_com, nm_com, npp_com
  !
  INTERFACE
    PURE SUBROUTINE PRDCT(Nd,Bd,Nm1,Bm1,Nm2,Bm2,Na,Aa,X,Y,M,A,B,C,D,U,W)
      IMPORT SP
      INTEGER, INTENT(IN) :: M, Na, Nd, Nm1, Nm2
      REAL(SP), INTENT(IN) :: Aa(Na), Bd(Nd), Bm1(Nm1), Bm2(Nm2)
      COMPLEX(SP), INTENT(IN) :: X(M), A(M), B(M), C(M)
      COMPLEX(SP), INTENT(OUT) :: Y(M), D(M), W(M), U(M)
    END SUBROUTINE PRDCT
    PURE SUBROUTINE CPRDCT(Nd,Bd,Nm1,Bm1,Nm2,Bm2,Na,Aa,X,Y,M,A,B,C,D,W,Yy)
      IMPORT SP
      INTEGER, INTENT(IN) :: M, Na, Nd, Nm1, Nm2
      REAL(SP), INTENT(IN) :: Aa(Na), Bm1(Nm1), Bm2(Nm2)
      COMPLEX(SP), INTENT(IN) :: Bd(Nd), X(M), A(M), B(M), C(M)
      COMPLEX(SP), INTENT(OUT) :: Y(M), D(M), W(M), Yy(M)
    END SUBROUTINE CPRDCT
  END INTERFACE
  INTEGER, INTENT(IN) :: Idimy, M
  REAL(SP), INTENT(IN) :: An(nm_com), Cn(nm_com), B(:)
  COMPLEX(SP), INTENT(IN) :: Am(M), Bm(M), Cm(M)
  COMPLEX(SP), INTENT(INOUT) :: Y(Idimy,nm_com)
  COMPLEX(SP), INTENT(OUT) :: W1(M), W2(M), W3(M), Wd(M), Ww(M), Wu(M)
  !
  INTEGER :: i, i1, i2, i3, i4, idxa, idxc, if, ifd, im1, im2, im3, imi1, &
    imi2, ip, ip1, ip2, ip3, ipi1, ipi2, ipi3, ir, irm1, iz, izr, j, kdo, l, ll, &
    na, nc, nm1, nm2, nm3, np, np1, np2, np3, nz
  REAL(SP) :: dum(0)
  COMPLEX(SP) :: bc(nm_com/2), duma(M)
  !* FIRST EXECUTABLE STATEMENT  CBLKT1
  kdo = k_com - 1
  DO l = 1, kdo
    ir = l - 1
    i2 = 2**ir
    i1 = i2/2
    i3 = i2 + i1
    i4 = i2 + i2
    irm1 = ir - 1
    CALL INXCB(i2,ir,im2,nm2)
    CALL INXCB(i1,irm1,im3,nm3)
    CALL INXCB(i3,irm1,im1,nm1)
    CALL PRDCT(nm2,B(im2:im2+nm2-1),nm3,B(im3:im3+nm3-1),nm1,B(im1:im1+nm1-1),0,&
      dum,Y(1,i2),W3,M,Am,Bm,Cm,Wd,Ww,Wu)
    if = 2**k_com
    DO i = i4, if, i4
      IF( i<=nm_com ) THEN
        ipi1 = i + i1
        ipi2 = i + i2
        ipi3 = i + i3
        CALL INXCC(i,ir,idxc,nc)
        IF( i<if ) THEN
          CALL INXCA(i,ir,idxa,na)
          CALL INXCB(i-i1,irm1,im1,nm1)
          CALL INXCB(ipi2,ir,ip2,np2)
          CALL INXCB(ipi1,irm1,ip1,np1)
          CALL INXCB(ipi3,irm1,ip3,np3)
          CALL PRDCT(nm1,B(im1:im1+nm1-1),0,dum,0,dum,na,An(idxa),W3,W1,M,Am,Bm,&
            Cm,Wd,Ww,Wu)
          IF( ipi2<=nm_com ) THEN
            CALL PRDCT(np2,B(ip2:ip2+np2-1),np1,B(ip1:ip1+np1-1),np3,B(ip3:ip3+np3-1),0,&
              dum,Y(1,ipi2),W3,M,Am,Bm,Cm,Wd,Ww,Wu)
            CALL PRDCT(np1,B(ip1:ip1+np1-1),0,dum,0,dum,nc,Cn(idxc),W3,W2,M,Am,Bm,Cm,&
              Wd,Ww,Wu)
          ELSE
            DO j = 1, M
              W3(j) = (0._SP,0._SP)
              W2(j) = (0._SP,0._SP)
            END DO
          END IF
          DO j = 1, M
            Y(j,i) = W1(j) + W2(j) + Y(j,i)
          END DO
        END IF
      END IF
    END DO
  END DO
  IF( npp_com==0 ) THEN
    !
    !     THE PERIODIC CASE IS TREATED USING THE CAPACITANCE MATRIX METHOD
    !
    if = 2**k_com
    i = if/2
    i1 = i/2
    CALL INXCB(i-i1,k_com-2,im1,nm1)
    CALL INXCB(i+i1,k_com-2,ip1,np1)
    CALL INXCB(i,k_com-1,iz,nz)
    CALL PRDCT(nz,B(iz:iz+nz-1),nm1,B(im1:im1+nm1-1),np1,B(ip1:ip1+np1-1),0,dum,&
      Y(1,i),W1,M,Am,Bm,Cm,Wd,Ww,Wu)
    izr = i
    DO j = 1, M
      W2(j) = W1(j)
    END DO
    DO ll = 2, k_com
      l = k_com - ll + 1
      ir = l - 1
      i2 = 2**ir
      i1 = i2/2
      i = i2
      CALL INXCC(i,ir,idxc,nc)
      CALL INXCB(i,ir,iz,nz)
      CALL INXCB(i-i1,ir-1,im1,nm1)
      CALL INXCB(i+i1,ir-1,ip1,np1)
      duma = W1
      CALL PRDCT(np1,B(ip1:ip1+np1-1),0,dum,0,dum,nc,Cn(idxc),duma,W1,M,Am,Bm,Cm,Wd,Ww,Wu)
      DO j = 1, M
        W1(j) = Y(j,i) + W1(j)
      END DO
      duma = W1
      CALL PRDCT(nz,B(iz:iz+nz-1),nm1,B(im1:im1+nm1-1),np1,B(ip1:ip1+np1-1),0,&
        dum,duma,W1,M,Am,Bm,Cm,Wd,Ww,Wu)
    END DO
    DO ll = 2, k_com
      l = k_com - ll + 1
      ir = l - 1
      i2 = 2**ir
      i1 = i2/2
      i4 = i2 + i2
      ifd = if - i2
      DO i = i2, ifd, i4
        IF( i-i2==izr ) THEN
          IF( i>nm_com ) EXIT
          CALL INXCA(i,ir,idxa,na)
          CALL INXCB(i,ir,iz,nz)
          CALL INXCB(i-i1,ir-1,im1,nm1)
          CALL INXCB(i+i1,ir-1,ip1,np1)
          duma = W2
          CALL PRDCT(nm1,B(im1:im1+nm1-1),0,dum,0,dum,na,An(idxa),duma,W2,M,Am,Bm,&
            Cm,Wd,Ww,Wu)
          DO j = 1, M
            W2(j) = Y(j,i) + W2(j)
          END DO
          duma = W2
          CALL PRDCT(nz,B(iz:iz+nz-1),nm1,B(im1:im1+nm1-1),np1,B(ip1:ip1+np1-1),0,&
            dum,duma,W2,M,Am,Bm,Cm,Wd,Ww,Wu)
          izr = i
          IF( i==nm_com ) GOTO 50
        END IF
      END DO
    END DO
    50 CONTINUE
    DO j = 1, M
      Y(j,nm_com+1) = Y(j,nm_com+1) - Cn(nm_com+1)*W1(j) - An(nm_com+1)*W2(j)
    END DO
    CALL INXCB(if/2,k_com-1,im1,nm1)
    CALL INXCB(if,k_com-1,ip,np)
    IF( ncmplx_com/=0 ) THEN
      bc = [ ( CMPLX( B(ip+i), B(ip+i+1), SP ), i = 1, nm_com, 2 ) ]
      CALL CPRDCT(nm_com+1,bc,nm1,B(im1:im1+nm1-1),0,dum,0,dum,Y(1,nm_com+1),Y(1,nm_com+1),M,&
        Am,Bm,Cm,W1,W3,Ww)
    ELSE
      CALL PRDCT(nm_com+1,B(ip:ip+nm_com),nm1,B(im1:im1+nm1-1),0,dum,0,dum,Y(1,nm_com+1),&
        Y(1,nm_com+1),M,Am,Bm,Cm,Wd,Ww,Wu)
    END IF
    DO j = 1, M
      W1(j) = An(1)*Y(j,nm_com+1)
      W2(j) = Cn(nm_com)*Y(j,nm_com+1)
      Y(j,1) = Y(j,1) - W1(j)
      Y(j,nm_com) = Y(j,nm_com) - W2(j)
    END DO
    DO l = 1, kdo
      ir = l - 1
      i2 = 2**ir
      i4 = i2 + i2
      i1 = i2/2
      i = i4
      CALL INXCA(i,ir,idxa,na)
      CALL INXCB(i-i2,ir,im2,nm2)
      CALL INXCB(i-i2-i1,ir-1,im3,nm3)
      CALL INXCB(i-i1,ir-1,im1,nm1)
      duma = W1
      CALL PRDCT(nm2,B(im2:im2+nm2-1),nm3,B(im3:im3+nm3-1),nm1,B(im1:im1-nm1-1),0,&
        dum,duma,W1,M,Am,Bm,Cm,Wd,Ww,Wu)
      duma = W1
      CALL PRDCT(nm1,B(im1:im1+nm1-1),0,dum,0,dum,na,An(idxa),duma,W1,M,Am,Bm,Cm,Wd,Ww,Wu)
      DO j = 1, M
        Y(j,i) = Y(j,i) - W1(j)
      END DO
    END DO
    !
    izr = nm_com
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
        IF( ipi2==izr ) THEN
          CALL INXCC(i,ir,idxc,nc)
          CALL INXCB(ipi2,ir,ip2,np2)
          CALL INXCB(ipi1,irm1,ip1,np1)
          CALL INXCB(ipi3,irm1,ip3,np3)
          duma = W2
          CALL PRDCT(np2,B(ip2:ip2+np2-1),np1,B(ip1:ip1+np1-1),np3,B(ip3:ip3+np3-1),0,&
            dum,duma,W2,M,Am,Bm,Cm,Wd,Ww,Wu)
          duma = W2
          CALL PRDCT(np1,B(ip1:ip1+np1-1),0,dum,0,dum,nc,Cn(idxc),duma,W2,M,Am,Bm,Cm,Wd,&
            Ww,Wu)
          DO j = 1, M
            Y(j,i) = Y(j,i) - W2(j)
          END DO
          izr = i
          EXIT
        ELSEIF( i==izr ) THEN
          EXIT
        END IF
      END DO
    END DO
  END IF
  !
  ! BEGIN BACK SUBSTITUTION PHASE
  !
  DO ll = 1, k_com
    l = k_com - ll + 1
    ir = l - 1
    irm1 = ir - 1
    i2 = 2**ir
    i1 = i2/2
    i4 = i2 + i2
    ifd = if - i2
    DO i = i2, ifd, i4
      IF( i<=nm_com ) THEN
        imi1 = i - i1
        imi2 = i - i2
        ipi1 = i + i1
        ipi2 = i + i2
        CALL INXCA(i,ir,idxa,na)
        CALL INXCC(i,ir,idxc,nc)
        CALL INXCB(i,ir,iz,nz)
        CALL INXCB(imi1,irm1,im1,nm1)
        CALL INXCB(ipi1,irm1,ip1,np1)
        IF( i<=i2 ) THEN
          DO j = 1, M
            W1(j) = (0._SP,0._SP)
          END DO
        ELSE
          CALL PRDCT(nm1,B(im1:im1+nm1-1),0,dum,0,dum,na,An(idxa),Y(1,imi2),W1,M,Am,&
            Bm,Cm,Wd,Ww,Wu)
        END IF
        IF( ipi2<=nm_com ) THEN
          CALL PRDCT(np1,B(ip1:ip1+np1-1),0,dum,0,dum,nc,Cn(idxc),Y(1,ipi2),W2,M,Am,&
            Bm,Cm,Wd,Ww,Wu)
        ELSE
          DO j = 1, M
            W2(j) = (0._SP,0._SP)
          END DO
        END IF
        DO j = 1, M
          W1(j) = Y(j,i) + W1(j) + W2(j)
        END DO
        CALL PRDCT(nz,B(iz:iz+nz-1),nm1,B(im1:im1+nm1-1),np1,B(ip1:ip1+np1-1),0,&
          dum,W1,Y(1,i),M,Am,Bm,Cm,Wd,Ww,Wu)
      END IF
    END DO
  END DO
  !
END SUBROUTINE CBLKT1