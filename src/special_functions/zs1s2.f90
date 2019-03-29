!** ZS1S2
SUBROUTINE ZS1S2(Zrr,Zri,S1r,S1i,S2r,S2i,Nz,Ascle,Alim,Iuf)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to ZAIRY and ZBESK
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CS1S2-A, ZS1S2-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     ZS1S2 TESTS FOR A POSSIBLE UNDERFLOW RESULTING FROM THE
  !     ADDITION OF THE I AND K FUNCTIONS IN THE ANALYTIC CON-
  !     TINUATION FORMULA WHERE S1=K FUNCTION AND S2=I FUNCTION.
  !     ON KODE=1 THE I AND K FUNCTIONS ARE DIFFERENT ORDERS OF
  !     MAGNITUDE, BUT FOR KODE=2 THEY CAN BE OF THE SAME ORDER
  !     OF MAGNITUDE AND THE MAXIMUM MUST BE AT LEAST ONE
  !     PRECISION ABOVE THE UNDERFLOW LIMIT.
  !
  !***
  ! **See also:**  ZAIRY, ZBESK
  !***
  ! **Routines called:**  ZABS, ZEXP, ZLOG

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !   930122  Added ZEXP and ZLOG to EXTERNAL statement.  (RWC)

  !     COMPLEX CZERO,C1,S1,S1D,S2,ZR
  REAL(8) :: aa, Alim, aln, Ascle, as1, as2, c1i, c1r, s1di, s1dr, S1i, S1r, &
    S2i, S2r, Zri, Zrr
  INTEGER Iuf, idum, Nz
  REAL(8), EXTERNAL :: ZABS
  EXTERNAL :: ZEXP, ZLOG
  REAL(8), PARAMETER :: zeror = 0.0D0, zeroi = 0.0D0
  !* FIRST EXECUTABLE STATEMENT  ZS1S2
  Nz = 0
  as1 = ZABS(S1r,S1i)
  as2 = ZABS(S2r,S2i)
  IF ( S1r/=0.0D0.OR.S1i/=0.0D0 ) THEN
    IF ( as1/=0.0D0 ) THEN
      aln = -Zrr - Zrr + LOG(as1)
      s1dr = S1r
      s1di = S1i
      S1r = zeror
      S1i = zeroi
      as1 = zeror
      IF ( aln>=(-Alim) ) THEN
        CALL ZLOG(s1dr,s1di,c1r,c1i,idum)
        c1r = c1r - Zrr - Zrr
        c1i = c1i - Zri - Zri
        CALL ZEXP(c1r,c1i,S1r,S1i)
        as1 = ZABS(S1r,S1i)
        Iuf = Iuf + 1
      ENDIF
    ENDIF
  ENDIF
  aa = MAX(as1,as2)
  IF ( aa>Ascle ) RETURN
  S1r = zeror
  S1i = zeroi
  S2r = zeror
  S2i = zeroi
  Nz = 1
  Iuf = 0
END SUBROUTINE ZS1S2
