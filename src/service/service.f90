MODULE service
  use ISO_FORTRAN_ENV, ONLY : INPUT_UNIT, OUTPUT_UNIT, ERROR_UNIT, &
    NUMERIC_STORAGE_SIZE, CHARACTER_STORAGE_SIZE
  IMPLICIT NONE
  INTEGER, PARAMETER :: SP = SELECTED_REAL_KIND(6,37)
  INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(15,307)
  INTEGER, PARAMETER :: QP = SELECTED_REAL_KIND(33,4931)
  !   I/O unit numbers:
  !     I1MACH( 1) = the standard input unit.
  !     I1MACH( 2) = the standard output unit.
  !     I1MACH( 3) = the standard punch unit.
  !     I1MACH( 4) = the standard error message unit.
  !
  !   Words:
  !     I1MACH( 5) = the number of bits per integer storage unit.
  !     I1MACH( 6) = the number of characters per integer storage unit.
  !
  !   Integers:
  !     assume integers are represented in the S-digit, base-A form
  !
  !                sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
  !
  !                where 0 <= X(I) < A for I=0,...,S-1.
  !     I1MACH( 7) = A, the base.
  !     I1MACH( 8) = S, the number of base-A digits.
  !     I1MACH( 9) = A**S - 1, the largest magnitude.
  !
  !   Floating-Point Numbers:
  !     Assume floating-point numbers are represented in the T-digit,
  !     base-B form
  !                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
  !
  !                where 0 <= X(I) < B for I=1,...,T,
  !                0 < X(1), and EMIN <= E <= EMAX.
  !     I1MACH(10) = B, the base.
  !
  !   Single-Precision:
  !     I1MACH(11) = T, the number of base-B digits.
  !     I1MACH(12) = EMIN, the smallest exponent E.
  !     I1MACH(13) = EMAX, the largest exponent E.
  !
  !   Double-Precision:
  !     I1MACH(14) = T, the number of base-B digits.
  !     I1MACH(15) = EMIN, the smallest exponent E.
  !     I1MACH(16) = EMAX, the largest exponent E.
  INTEGER, PARAMETER :: I1MACH(16) = [ INPUT_UNIT, OUTPUT_UNIT, OUTPUT_UNIT, &
    ERROR_UNIT, NUMERIC_STORAGE_SIZE, CHARACTER_STORAGE_SIZE, &
    RADIX(1), DIGITS(1), HUGE(1), RADIX(1._SP), &
    DIGITS(1._SP), MINEXPONENT(1._SP), MAXEXPONENT(1._SP), &
    DIGITS(1._DP), MINEXPONENT(1._DP), MAXEXPONENT(1._DP) ]
  !   R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
  !   R1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
  !   R1MACH(3) = B**(-T), the smallest relative spacing.
  !   R1MACH(4) = B**(1-T), the largest relative spacing.
  !   R1MACH(5) = LOG10(B)
  !
  !   Assume single precision numbers are represented in the T-digit,
  !   base-B form
  !
  !              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
  !
  !   where 0 <= X(I) < B for I=1,...,T, 0 < X(1), and
  !   EMIN <= E <= EMAX.
  REAL(SP), PARAMETER :: R1MACH(5) = [ TINY(1._SP), HUGE(1._SP), SPACING(0.5_SP), &
    SPACING(1._SP), LOG10( REAL( RADIX(1._SP), SP ) ) ]
  !   D1MACH( 1) = B**(EMIN-1), the smallest positive magnitude.
  !   D1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
  !   D1MACH( 3) = B**(-T), the smallest relative spacing.
  !   D1MACH( 4) = B**(1-T), the largest relative spacing.
  !   D1MACH( 5) = LOG10(B)
  !
  !   Assume double precision numbers are represented in the T-digit,
  !   base-B form
  !
  !              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
  !
  !   where 0 <= X(I) < B for I=1,...,T, 0 < X(1), and
  !   EMIN <= E <= EMAX.
  REAL(DP), PARAMETER :: D1MACH(5) = [ TINY(1._DP), HUGE(1._DP), SPACING(0.5_DP), &
    SPACING(1._DP), LOG10( REAL( RADIX(1._DP), DP ) ) ]

CONTAINS
  include"fdump.f90"
  include"j4save.f90"
  include"numxer.f90"
  include"xerbla.f90"
  include"xerclr.f90"
  include"xerdmp.f90"
  include"xermax.f90"
  include"xermsg.f90"
  include"xerprn.f90"
  include"xersve.f90"
  include"xgetf.f90"
  include"xgetua.f90"
  include"xgetun.f90"
  include"xsetf.f90"
  include"xsetua.f90"
  include"xsetun.f90"
END MODULE service