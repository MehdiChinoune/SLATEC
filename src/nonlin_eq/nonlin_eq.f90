MODULE nonlin_eq
  use service
  use linear
  IMPLICIT NONE

CONTAINS
  include"cpevl.f90"
  include"cpevlr.f90"
  include"cpqr79.f90"
  include"cpzero.f90"
  include"dfzero.f90"
  include"dsos.f90"
  include"dsoseq.f90"
  include"dsossl.f90"
  include"fzero.f90"
  include"rpqr79.f90"
  include"rpzero.f90"
  include"sos.f90"
  include"soseqs.f90"
  include"sossol.f90"
END MODULE nonlin_eq