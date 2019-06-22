MODULE common_mod
  USE service, ONLY: SP
  IMPLICIT NONE

CONTAINS
  include"get_argument.f90"
  include"ismpl.f90"
  include"pass.f90"
END MODULE common_mod