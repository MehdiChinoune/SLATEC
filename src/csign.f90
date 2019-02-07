
COMPLEX FUNCTION CSIGN(zdum1,zdum2)
  COMPLEX, INTENT(IN) :: zdum1 , zdum2

  CSIGN = ABS(zdum1)*(zdum2/ABS(zdum2))

END FUNCTION CSIGN
