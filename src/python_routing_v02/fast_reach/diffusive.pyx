import cython

from fortran_wrappers cimport c_diffnw

cdef struct QE:
    double q_ev_g,
    double elv_ev_g

@cython.boundscheck(False)
cdef void diffnw(double dtini_g,
             double t0_g,
             double tfin_g,
             double saveinterval_g,
             double saveinterval_ev_g,
             double dt_ql_g,
             double dt_ub_g,
             double dt_db_g,
             int nts_ql_g,
             int nts_ub_g,
             int nts_db_g,
             int mxncomp_g,
             int nrch_g,
             double z_ar_g,
             double bo_ar_g,
             double traps_ar_g,
             double tw_ar_g,
             double twcc_ar_g,
             double mann_ar_g,
             double manncc_ar_g,
             double so_ar_g,
             double dx_ar_g,
             int nhincr_m_g,
             int nhincr_f_g,
             double ufhlt_m_g,
             double ufqlt_m_g,
             double ufhlt_f_g,
             double ufqlt_f_g,
             int frnw_col,
             int frnw_g,
             double qlat_g,
             double ubcd_g,
             double dbcd_g,
             double cfl_g,
             double theta_g,
             int tzeq_flag_g,
             int y_opt_g,
             double so_llm_g,
             int ntss_ev_g,
             QE *rv) nogil:
    
    cdef:
        double q_ev_g
        double elv_ev_g
        
    c_diffnw(
            &dtini_g,
            &t0_g,
            &tfin_g,
            &saveinterval_g,
            &saveinterval_ev_g,
            &dt_ql_g,
            &dt_ub_g,
            &dt_db_g,
            &nts_ql_g,
            &nts_ub_g,
            &nts_db_g,
            &mxncomp_g,
            &nrch_g,
            &z_ar_g,
            &bo_ar_g,
            &traps_ar_g,
            &tw_ar_g,
            &twcc_ar_g,
            &mann_ar_g,
            &manncc_ar_g,
            &so_ar_g,
            &dx_ar_g,
            &nhincr_m_g,
            &nhincr_f_g,
            &ufhlt_m_g,
            &ufqlt_m_g,
            &ufhlt_f_g,
            &ufqlt_f_g,
            &frnw_col,
            &frnw_g,
            &qlat_g,
            &ubcd_g,
            &dbcd_g,
            &cfl_g,
            &theta_g,
            &tzeq_flag_g,
            &y_opt_g,
            &so_llm_g,
            &ntss_ev_g,
            &q_ev_g,
            &elv_ev_g)
    
    rv.q_ev_g = q_ev_g
    rv.elv_ev_g = elv_ev_g
    
    