
COMPLEX(SP) ELEMENTAL FUNCTION CSIGN1(zdum1,zdum2)
  USE blas, ONLY : SCABS1
  COMPLEX(SP), INTENT(IN) :: zdum1, zdum2

  CSIGN1 = SCABS1(zdum1)*(zdum2/SCABS1(zdum2))

END FUNCTION CSIGN1