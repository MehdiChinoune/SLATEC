
      REAL FUNCTION CABS1(zdum)
        COMPLEX, INTENT(IN) :: zdum

        CABS1 = ABS(REAL(zdum)) + ABS(AIMAG(zdum))

      END FUNCTION CABS1
