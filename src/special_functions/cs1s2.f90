!** CS1S2
SUBROUTINE CS1S2(Zr,S1,S2,Nz,Ascle,Alim,Iuf)
  !>
  !  Subsidiary to CAIRY and CBESK
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CS1S2-A, ZS1S2-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     CS1S2 TESTS FOR A POSSIBLE UNDERFLOW RESULTING FROM THE
  !     ADDITION OF THE I AND K FUNCTIONS IN THE ANALYTIC CON-
  !     TINUATION FORMULA WHERE S1=K FUNCTION AND S2=I FUNCTION.
  !     ON KODE=1 THE I AND K FUNCTIONS ARE DIFFERENT ORDERS OF
  !     MAGNITUDE, BUT FOR KODE=2 THEY CAN BE OF THE SAME ORDER
  !     OF MAGNITUDE AND THE MAXIMUM MUST BE AT LEAST ONE
  !     PRECISION ABOVE THE UNDERFLOW LIMIT.
  !
  !***
  ! **See also:**  CAIRY, CBESK
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)

  COMPLEX c1, S1, s1d, S2, Zr
  REAL aa, Alim, aln, Ascle, as1, as2, xx
  INTEGER Iuf, Nz
  COMPLEX, PARAMETER :: czero  = (0.0E0,0.0E0)
  !* FIRST EXECUTABLE STATEMENT  CS1S2
  Nz = 0
  as1 = ABS(S1)
  as2 = ABS(S2)
  aa = REAL(S1)
  aln = AIMAG(S1)
  IF ( aa/=0.0E0.OR.aln/=0.0E0 ) THEN
    IF ( as1/=0.0E0 ) THEN
      xx = REAL(Zr)
      aln = -xx - xx + ALOG(as1)
      s1d = S1
      S1 = czero
      as1 = 0.0E0
      IF ( aln>=(-Alim) ) THEN
        c1 = CLOG(s1d) - Zr - Zr
        S1 = CEXP(c1)
        as1 = ABS(S1)
        Iuf = Iuf + 1
      END IF
    END IF
  END IF
  aa = MAX(as1,as2)
  IF ( aa>Ascle ) RETURN
  S1 = czero
  S2 = czero
  Nz = 1
  Iuf = 0
END SUBROUTINE CS1S2
