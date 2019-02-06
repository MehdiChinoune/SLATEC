!*==PASSB.f90  processed by SPAG 6.72Dc at 10:55 on  6 Feb 2019
!DECK PASSB
      SUBROUTINE PASSB(Nac,Ido,Ip,L1,Idl1,Cc,C1,C2,Ch,Ch2,Wa)
      IMPLICIT NONE
!*--PASSB5
!*** Start of declarations inserted by SPAG
      REAL C1 , C2 , Cc , Ch , Ch2 , Wa , wai , war
      INTEGER i , idij , idj , idl , Idl1 , idlj , Ido , idot , idp , ik , inc , 
     &        Ip , ipp2 , ipph , j , jc , k , l , L1 , lc
      INTEGER Nac
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  PASSB
!***SUBSIDIARY
!***PURPOSE  Calculate the fast Fourier transform of subvectors of
!            arbitrary length.
!***LIBRARY   SLATEC (FFTPACK)
!***TYPE      SINGLE PRECISION (PASSB-S)
!***AUTHOR  Swarztrauber, P. N., (NCAR)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   830401  Modified to use SLATEC library source file format.
!   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
!           changing dummy array size declarations (1) to (*).
!   881128  Modified by Dick Valent to meet prologue standards.
!   890831  Modified array declarations.  (WRB)
!   891009  Removed unreferenced variable.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  PASSB
      DIMENSION Ch(Ido,L1,*) , Cc(Ido,Ip,*) , C1(Ido,L1,*) , Wa(*) , C2(Idl1,*)
     &          , Ch2(Idl1,*)
!***FIRST EXECUTABLE STATEMENT  PASSB
      idot = Ido/2
      ipp2 = Ip + 2
      ipph = (Ip+1)/2
      idp = Ip*Ido
!
      IF ( Ido<L1 ) THEN
        DO j = 2 , ipph
          jc = ipp2 - j
          DO i = 1 , Ido
!DIR$ IVDEP
            DO k = 1 , L1
              Ch(i,k,j) = Cc(i,j,k) + Cc(i,jc,k)
              Ch(i,k,jc) = Cc(i,j,k) - Cc(i,jc,k)
            ENDDO
          ENDDO
        ENDDO
        DO i = 1 , Ido
!DIR$ IVDEP
          DO k = 1 , L1
            Ch(i,k,1) = Cc(i,1,k)
          ENDDO
        ENDDO
      ELSE
        DO j = 2 , ipph
          jc = ipp2 - j
          DO k = 1 , L1
!DIR$ IVDEP
            DO i = 1 , Ido
              Ch(i,k,j) = Cc(i,j,k) + Cc(i,jc,k)
              Ch(i,k,jc) = Cc(i,j,k) - Cc(i,jc,k)
            ENDDO
          ENDDO
        ENDDO
        DO k = 1 , L1
!DIR$ IVDEP
          DO i = 1 , Ido
            Ch(i,k,1) = Cc(i,1,k)
          ENDDO
        ENDDO
      ENDIF
      idl = 2 - Ido
      inc = 0
      DO l = 2 , ipph
        lc = ipp2 - l
        idl = idl + Ido
!DIR$ IVDEP
        DO ik = 1 , Idl1
          C2(ik,l) = Ch2(ik,1) + Wa(idl-1)*Ch2(ik,2)
          C2(ik,lc) = Wa(idl)*Ch2(ik,Ip)
        ENDDO
        idlj = idl
        inc = inc + Ido
        DO j = 3 , ipph
          jc = ipp2 - j
          idlj = idlj + inc
          IF ( idlj>idp ) idlj = idlj - idp
          war = Wa(idlj-1)
          wai = Wa(idlj)
!DIR$ IVDEP
          DO ik = 1 , Idl1
            C2(ik,l) = C2(ik,l) + war*Ch2(ik,j)
            C2(ik,lc) = C2(ik,lc) + wai*Ch2(ik,jc)
          ENDDO
        ENDDO
      ENDDO
      DO j = 2 , ipph
!DIR$ IVDEP
        DO ik = 1 , Idl1
          Ch2(ik,1) = Ch2(ik,1) + Ch2(ik,j)
        ENDDO
      ENDDO
      DO j = 2 , ipph
        jc = ipp2 - j
!DIR$ IVDEP
        DO ik = 2 , Idl1 , 2
          Ch2(ik-1,j) = C2(ik-1,j) - C2(ik,jc)
          Ch2(ik-1,jc) = C2(ik-1,j) + C2(ik,jc)
          Ch2(ik,j) = C2(ik,j) + C2(ik-1,jc)
          Ch2(ik,jc) = C2(ik,j) - C2(ik-1,jc)
        ENDDO
      ENDDO
      Nac = 1
      IF ( Ido==2 ) RETURN
      Nac = 0
      DO ik = 1 , Idl1
        C2(ik,1) = Ch2(ik,1)
      ENDDO
      DO j = 2 , Ip
!DIR$ IVDEP
        DO k = 1 , L1
          C1(1,k,j) = Ch(1,k,j)
          C1(2,k,j) = Ch(2,k,j)
        ENDDO
      ENDDO
      IF ( idot>L1 ) THEN
        idj = 2 - Ido
        DO j = 2 , Ip
          idj = idj + Ido
          DO k = 1 , L1
            idij = idj
!DIR$ IVDEP
            DO i = 4 , Ido , 2
              idij = idij + 2
              C1(i-1,k,j) = Wa(idij-1)*Ch(i-1,k,j) - Wa(idij)*Ch(i,k,j)
              C1(i,k,j) = Wa(idij-1)*Ch(i,k,j) + Wa(idij)*Ch(i-1,k,j)
            ENDDO
          ENDDO
        ENDDO
        GOTO 99999
      ENDIF
      idij = 0
      DO j = 2 , Ip
        idij = idij + 2
        DO i = 4 , Ido , 2
          idij = idij + 2
!DIR$ IVDEP
          DO k = 1 , L1
            C1(i-1,k,j) = Wa(idij-1)*Ch(i-1,k,j) - Wa(idij)*Ch(i,k,j)
            C1(i,k,j) = Wa(idij-1)*Ch(i,k,j) + Wa(idij)*Ch(i-1,k,j)
          ENDDO
        ENDDO
      ENDDO
      RETURN
99999 END SUBROUTINE PASSB
