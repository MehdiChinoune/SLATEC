
REAL FUNCTION DELX(xa,xb)
  COMPLEX, INTENT(IN) :: xa, xb

  DELX = ABS(REAL(xa-xb)) + ABS(AIMAG(xa-xb))

END FUNCTION DELX
