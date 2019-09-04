!** RADFG
PURE SUBROUTINE RADFG(Ido,Ip,L1,Idl1,Cc,C1,C2,Ch,Ch2,Wa)
  !> Calculate the fast Fourier transform of subvectors of arbitrary length.
  !***
  ! **Library:**   SLATEC (FFTPACK)
  !***
  ! **Type:**      SINGLE PRECISION (RADFG-S)
  !***
  ! **Author:**  Swarztrauber, P. N., (NCAR)
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   830401  Modified to use SLATEC library source file format.
  !   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
  !           (a) changing dummy array size declarations (1) to (*),
  !           (b) changing references to intrinsic function FLOAT
  !               to REAL(SP), and
  !           (c) changing definition of variable TPI by using
  !               FORTRAN intrinsic function ATAN instead of a DATA
  !               statement.
  !   881128  Modified by Dick Valent to meet prologue standards.
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)

  INTEGER , INTENT(IN) :: Idl1, Ido, Ip, L1
  REAL(SP), INTENT(IN) :: Wa(:)
  REAL(SP), INTENT(INOUT) :: Ch2(Idl1,Ip)
  REAL(SP), INTENT(OUT) :: Cc(Ido,Ip,L1), C1(Ido,L1,Ip), C2(Idl1,Ip), Ch(Ido,L1,Ip)
  INTEGER :: i, ic, idij, idp2, ik, ipp2, ipph, is, j, j2, jc, k, l, lc, nbd
  REAL(SP) :: ai1, ai2, ar1, ar1h, ar2, ar2h, arg, dc2, dcp, ds2, dsp, tpi
  !* FIRST EXECUTABLE STATEMENT  RADFG
  tpi = 8._SP*ATAN(1._SP)
  arg = tpi/Ip
  dcp = COS(arg)
  dsp = SIN(arg)
  ipph = (Ip+1)/2
  ipp2 = Ip + 2
  idp2 = Ido + 2
  nbd = (Ido-1)/2
  IF( Ido==1 ) THEN
    DO ik = 1, Idl1
      C2(ik,1) = Ch2(ik,1)
    END DO
  ELSE
    DO ik = 1, Idl1
      Ch2(ik,1) = C2(ik,1)
    END DO
    DO j = 2, Ip
      DO k = 1, L1
        Ch(1,k,j) = C1(1,k,j)
      END DO
    END DO
    IF( nbd>L1 ) THEN
      is = -Ido
      DO j = 2, Ip
        is = is + Ido
        DO k = 1, L1
          idij = is
          DO i = 3, Ido, 2
            idij = idij + 2
            Ch(i-1,k,j) = Wa(idij-1)*C1(i-1,k,j) + Wa(idij)*C1(i,k,j)
            Ch(i,k,j) = Wa(idij-1)*C1(i,k,j) - Wa(idij)*C1(i-1,k,j)
          END DO
        END DO
      END DO
    ELSE
      is = -Ido
      DO j = 2, Ip
        is = is + Ido
        idij = is
        DO i = 3, Ido, 2
          idij = idij + 2
          DO k = 1, L1
            Ch(i-1,k,j) = Wa(idij-1)*C1(i-1,k,j) + Wa(idij)*C1(i,k,j)
            Ch(i,k,j) = Wa(idij-1)*C1(i,k,j) - Wa(idij)*C1(i-1,k,j)
          END DO
        END DO
      END DO
    END IF
    IF( nbd<L1 ) THEN
      DO j = 2, ipph
        jc = ipp2 - j
        DO i = 3, Ido, 2
          DO k = 1, L1
            C1(i-1,k,j) = Ch(i-1,k,j) + Ch(i-1,k,jc)
            C1(i-1,k,jc) = Ch(i,k,j) - Ch(i,k,jc)
            C1(i,k,j) = Ch(i,k,j) + Ch(i,k,jc)
            C1(i,k,jc) = Ch(i-1,k,jc) - Ch(i-1,k,j)
          END DO
        END DO
      END DO
    ELSE
      DO j = 2, ipph
        jc = ipp2 - j
        DO k = 1, L1
          DO i = 3, Ido, 2
            C1(i-1,k,j) = Ch(i-1,k,j) + Ch(i-1,k,jc)
            C1(i-1,k,jc) = Ch(i,k,j) - Ch(i,k,jc)
            C1(i,k,j) = Ch(i,k,j) + Ch(i,k,jc)
            C1(i,k,jc) = Ch(i-1,k,jc) - Ch(i-1,k,j)
          END DO
        END DO
      END DO
    END IF
  END IF
  DO j = 2, ipph
    jc = ipp2 - j
    DO k = 1, L1
      C1(1,k,j) = Ch(1,k,j) + Ch(1,k,jc)
      C1(1,k,jc) = Ch(1,k,jc) - Ch(1,k,j)
    END DO
  END DO
  !
  ar1 = 1._SP
  ai1 = 0._SP
  DO l = 2, ipph
    lc = ipp2 - l
    ar1h = dcp*ar1 - dsp*ai1
    ai1 = dcp*ai1 + dsp*ar1
    ar1 = ar1h
    DO ik = 1, Idl1
      Ch2(ik,l) = C2(ik,1) + ar1*C2(ik,2)
      Ch2(ik,lc) = ai1*C2(ik,Ip)
    END DO
    dc2 = ar1
    ds2 = ai1
    ar2 = ar1
    ai2 = ai1
    DO j = 3, ipph
      jc = ipp2 - j
      ar2h = dc2*ar2 - ds2*ai2
      ai2 = dc2*ai2 + ds2*ar2
      ar2 = ar2h
      DO ik = 1, Idl1
        Ch2(ik,l) = Ch2(ik,l) + ar2*C2(ik,j)
        Ch2(ik,lc) = Ch2(ik,lc) + ai2*C2(ik,jc)
      END DO
    END DO
  END DO
  DO j = 2, ipph
    DO ik = 1, Idl1
      Ch2(ik,1) = Ch2(ik,1) + C2(ik,j)
    END DO
  END DO
  !
  IF( Ido<L1 ) THEN
    DO i = 1, Ido
      DO k = 1, L1
        Cc(i,1,k) = Ch(i,k,1)
      END DO
    END DO
  ELSE
    DO k = 1, L1
      DO i = 1, Ido
        Cc(i,1,k) = Ch(i,k,1)
      END DO
    END DO
  END IF
  DO j = 2, ipph
    jc = ipp2 - j
    j2 = j + j
    DO k = 1, L1
      Cc(Ido,j2-2,k) = Ch(1,k,j)
      Cc(1,j2-1,k) = Ch(1,k,jc)
    END DO
  END DO
  IF( Ido==1 ) RETURN
  IF( nbd<L1 ) THEN
    DO j = 2, ipph
      jc = ipp2 - j
      j2 = j + j
      DO i = 3, Ido, 2
        ic = idp2 - i
        DO k = 1, L1
          Cc(i-1,j2-1,k) = Ch(i-1,k,j) + Ch(i-1,k,jc)
          Cc(ic-1,j2-2,k) = Ch(i-1,k,j) - Ch(i-1,k,jc)
          Cc(i,j2-1,k) = Ch(i,k,j) + Ch(i,k,jc)
          Cc(ic,j2-2,k) = Ch(i,k,jc) - Ch(i,k,j)
        END DO
      END DO
    END DO
    RETURN
  END IF
  DO j = 2, ipph
    jc = ipp2 - j
    j2 = j + j
    DO k = 1, L1
      DO i = 3, Ido, 2
        ic = idp2 - i
        Cc(i-1,j2-1,k) = Ch(i-1,k,j) + Ch(i-1,k,jc)
        Cc(ic-1,j2-2,k) = Ch(i-1,k,j) - Ch(i-1,k,jc)
        Cc(i,j2-1,k) = Ch(i,k,j) + Ch(i,k,jc)
        Cc(ic,j2-2,k) = Ch(i,k,jc) - Ch(i,k,j)
      END DO
    END DO
  END DO
  !
  RETURN
END SUBROUTINE RADFG