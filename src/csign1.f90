
      COMPLEX FUNCTION CSIGN1(zdum1,zdum2)
        COMPLEX, INTENT(IN) :: zdum1 , zdum2
        REAL, EXTERNAL :: CABS1

        CSIGN1 = CABS1(zdum1)*(zdum2/CABS1(zdum2))

      END FUNCTION CSIGN1
