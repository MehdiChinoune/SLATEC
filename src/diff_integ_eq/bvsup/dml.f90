MODULE DML
  USE service, ONLY : DP, eps_dp, huge_dp
  IMPLICIT NONE
  !
  REAL(DP), PARAMETER :: uro_com = eps_dp, sru_com = SQRT(uro_com), &
    twou_com = 2._DP*uro_com, fouru_com = 4._DP*uro_com, dd = -LOG10(uro_com), &
    sqovfl_com = SQRT(huge_dp)
  INTEGER, PARAMETER :: lpar_com = INT( 0.5_DP*dd ), ke = INT( 0.5_DP + 0.75_DP*dd )
  REAL(DP), PARAMETER :: eps_com = 10._DP**(-2*ke)
  !
  REAL(DP) :: c_com, xsav_com, px_com, pwcnd_com, tnd_com, x_com, xbeg_com, xend_com, &
    xot_com, xop_com, ae_com, re_com, tol_com
  INTEGER :: igofx_com, inhomo_com, ivp_com, ncomp_com, nfc_com, info_com(15), &
    istkop_com, knswot_com, kop_com, lotjp_com, mnswot_com, nswot_com, nxpts_com, &
    nic_com, nopg_com, mxnon_com, ndisk_com, ntape_com, neq_com, indpvt_com, &
    integ_com, nps_com, ntp_com, neqivp_com, numort_com, nfcc_com, icoco_com, &
    kkkzpw_com, needw_com, neediw_com, k1_com, k2_com, k3_com, k4_com, k5_com, &
    k6_com, k7_com, k8_com, k9_com, k10_com, k11_com, l1_com, l2_com, kkkint_com, &
    lllint_com, nofst_com
END MODULE DML