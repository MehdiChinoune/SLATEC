!** PASSF
SUBROUTINE PASSF(Nac,Ido,Ip,L1,Idl1,Cc,C1,C2,Ch,Ch2,Wa)
  !>
  !  Calculate the fast Fourier transform of subvectors of
  !            arbitrary length.
  !***
  ! **Library:**   SLATEC (FFTPACK)
  !***
  ! **Type:**      SINGLE PRECISION (PASSF-S)
  !***
  ! **Author:**  Swarztrauber, P. N., (NCAR)
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   830401  Modified to use SLATEC library source file format.
  !   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
  !           changing dummy array size declarations (1) to (*).
  !   881128  Modified by Dick Valent to meet prologue standards.
  !   890831  Modified array declarations.  (WRB)
  !   891009  Removed unreferenced variable.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)

  INTEGER :: Nac, Idl1, Ido, Ip, L1
  REAL(SP) :: C1(Ido,L1,Ip), C2(Idl1,Ip), Cc(Ido,Ip,L1), Ch(Ido,L1,Ip), Ch2(Idl1,Ip), Wa(:)
  INTEGER :: i, idij, idj, idl, idlj, idot, idp, ik, inc, ipp2, ipph, j, jc, k, l, lc
  REAL(SP) :: wai, war
  !* FIRST EXECUTABLE STATEMENT  PASSF
  idot = Ido/2
  ipp2 = Ip + 2
  ipph = (Ip+1)/2
  idp = Ip*Ido
  !
  IF ( Ido<L1 ) THEN
    DO j = 2, ipph
      jc = ipp2 - j
      DO i = 1, Ido
        DO k = 1, L1
          Ch(i,k,j) = Cc(i,j,k) + Cc(i,jc,k)
          Ch(i,k,jc) = Cc(i,j,k) - Cc(i,jc,k)
        END DO
      END DO
    END DO
    DO i = 1, Ido
      DO k = 1, L1
        Ch(i,k,1) = Cc(i,1,k)
      END DO
    END DO
  ELSE
    DO j = 2, ipph
      jc = ipp2 - j
      DO k = 1, L1
        DO i = 1, Ido
          Ch(i,k,j) = Cc(i,j,k) + Cc(i,jc,k)
          Ch(i,k,jc) = Cc(i,j,k) - Cc(i,jc,k)
        END DO
      END DO
    END DO
    DO k = 1, L1
      DO i = 1, Ido
        Ch(i,k,1) = Cc(i,1,k)
      END DO
    END DO
  END IF
  idl = 2 - Ido
  inc = 0
  DO l = 2, ipph
    lc = ipp2 - l
    idl = idl + Ido
    DO ik = 1, Idl1
      C2(ik,l) = Ch2(ik,1) + Wa(idl-1)*Ch2(ik,2)
      C2(ik,lc) = -Wa(idl)*Ch2(ik,Ip)
    END DO
    idlj = idl
    inc = inc + Ido
    DO j = 3, ipph
      jc = ipp2 - j
      idlj = idlj + inc
      IF ( idlj>idp ) idlj = idlj - idp
      war = Wa(idlj-1)
      wai = Wa(idlj)
      DO ik = 1, Idl1
        C2(ik,l) = C2(ik,l) + war*Ch2(ik,j)
        C2(ik,lc) = C2(ik,lc) - wai*Ch2(ik,jc)
      END DO
    END DO
  END DO
  DO j = 2, ipph
    DO ik = 1, Idl1
      Ch2(ik,1) = Ch2(ik,1) + Ch2(ik,j)
    END DO
  END DO
  DO j = 2, ipph
    jc = ipp2 - j
    DO ik = 2, Idl1, 2
      Ch2(ik-1,j) = C2(ik-1,j) - C2(ik,jc)
      Ch2(ik-1,jc) = C2(ik-1,j) + C2(ik,jc)
      Ch2(ik,j) = C2(ik,j) + C2(ik-1,jc)
      Ch2(ik,jc) = C2(ik,j) - C2(ik-1,jc)
    END DO
  END DO
  Nac = 1
  IF ( Ido==2 ) RETURN
  Nac = 0
  DO ik = 1, Idl1
    C2(ik,1) = Ch2(ik,1)
  END DO
  DO j = 2, Ip
    DO k = 1, L1
      C1(1,k,j) = Ch(1,k,j)
      C1(2,k,j) = Ch(2,k,j)
    END DO
  END DO
  IF ( idot>L1 ) THEN
    idj = 2 - Ido
    DO j = 2, Ip
      idj = idj + Ido
      DO k = 1, L1
        idij = idj
        DO i = 4, Ido, 2
          idij = idij + 2
          C1(i-1,k,j) = Wa(idij-1)*Ch(i-1,k,j) + Wa(idij)*Ch(i,k,j)
          C1(i,k,j) = Wa(idij-1)*Ch(i,k,j) - Wa(idij)*Ch(i-1,k,j)
        END DO
      END DO
    END DO
    RETURN
  END IF
  idij = 0
  DO j = 2, Ip
    idij = idij + 2
    DO i = 4, Ido, 2
      idij = idij + 2
      DO k = 1, L1
        C1(i-1,k,j) = Wa(idij-1)*Ch(i-1,k,j) + Wa(idij)*Ch(i,k,j)
        C1(i,k,j) = Wa(idij-1)*Ch(i,k,j) - Wa(idij)*Ch(i-1,k,j)
      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE PASSF
