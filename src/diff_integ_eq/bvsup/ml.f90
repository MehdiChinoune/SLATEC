MODULE ML
  USE service, ONLY : SP, eps_sp, huge_sp
  IMPLICIT NONE
  !
  REAL(SP), PARAMETER :: uro_com = eps_sp, sru_com = SQRT(uro_com), &
    twou_com = 2._SP*uro_com, fouru_com = 4._SP*uro_com, dd = -LOG10(uro_com), &
    sqovfl_com = SQRT(huge_sp)
  INTEGER, PARAMETER :: lpar_com = INT( 0.5_SP*dd ), ke = INT( 0.5_SP + 0.75_SP*dd )
  REAL(SP), PARAMETER :: eps_com = 10._SP**(-2*ke)
  !
  REAL(SP) :: c_com, xsav_com, px_com, pwcnd_com, tnd_com, x_com, xbeg_com, xend_com, &
    xot_com, xop_com, ae_com, re_com, tol_com
  INTEGER :: igofx_com, inhomo_com, ivp_com, ncomp_com, nfc_com, info_com(15), &
    istkop_com, knswot_com, kop_com, lotjp_com, mnswot_com, nswot_com, nxpts_com, &
    nic_com, nopg_com, mxnon_com, ndisk_com, ntape_com, neq_com, indpvt_com, &
    integ_com, nps_com, ntp_com, neqivp_com, numort_com, nfcc_com, icoco_com, &
    kkkzpw_com, needw_com, neediw_com, k1_com, k2_com, k3_com, k4_com, k5_com, &
    k6_com, k7_com, k8_com, k9_com, k10_com, k11_com, l1_com, l2_com, kkkint_com, &
    lllint_com, nofst_com
END MODULE ML