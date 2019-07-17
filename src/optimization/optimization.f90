include"splp/la05ds.f90"
include"splp/la05dd.f90"

MODULE optimization
  USE service, ONLY : SP, DP
  IMPLICIT NONE

CONTAINS
  include"splp/dfulmt.f90"
  include"splp/dpchng.f90"
  include"splp/dpincw.f90"
  include"splp/dpinit.f90"
  include"splp/dpintm.f90"
  include"splp/dplpce.f90"
  include"splp/dplpdm.f90"
  include"splp/dplpfe.f90"
  include"splp/dplpfl.f90"
  include"splp/dplpmn.f90"
  include"splp/dplpmu.f90"
  include"splp/dplpup.f90"
  include"splp/dpnnzr.f90"
  include"splp/dpopt.f90"
  include"splp/dsplp.f90"
  include"splp/dusrmt.f90"
  include"splp/fulmat.f90"
  include"splp/idloc.f90"
  include"splp/iploc.f90"
  include"splp/la05ad.f90"
  include"splp/la05as.f90"
  include"splp/la05bd.f90"
  include"splp/la05bs.f90"
  include"splp/la05cd.f90"
  include"splp/la05cs.f90"
  include"splp/la05ed.f90"
  include"splp/la05es.f90"
  include"splp/mc20ad.f90"
  include"splp/mc20as.f90"
  include"splp/pchngs.f90"
  include"splp/pinitm.f90"
  include"splp/pnnzrs.f90"
  include"splp/spincw.f90"
  include"splp/spinit.f90"
  include"splp/splp.f90"
  include"splp/splpce.f90"
  include"splp/splpdm.f90"
  include"splp/splpfe.f90"
  include"splp/splpfl.f90"
  include"splp/splpmn.f90"
  include"splp/splpmu.f90"
  include"splp/splpup.f90"
  include"splp/spopt.f90"
  include"splp/usrmat.f90"
END MODULE optimization