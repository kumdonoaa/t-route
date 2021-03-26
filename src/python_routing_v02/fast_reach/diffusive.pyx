import cython
import numpy as np
cimport numpy as np

from fortran_wrappers cimport c_diffnw

# TO DO load some example inputs to test the module

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
             double[::1,:] z_ar_g,
             double[::1,:] bo_ar_g,
             double[::1,:] traps_ar_g,
             double[::1,:] tw_ar_g,
             double[::1,:] twcc_ar_g,
             double[::1,:] mann_ar_g,
             double[::1,:] manncc_ar_g,
             double[::1,:] so_ar_g,
             double[::1,:] dx_ar_g,
             int nhincr_m_g,
             int nhincr_f_g,
             double[::1,:,:] ufhlt_m_g,
             double[::1,:,:] ufqlt_m_g,
             double[::1,:,:] ufhlt_f_g,
             double[::1,:,:] ufqlt_f_g,
             int frnw_col,
             double[::1,:] frnw_g,
             double[::1,:,:] qlat_g,
             double[::1,:] ubcd_g,
             double[::1] dbcd_g,
             double cfl_g,
             double theta_g,
             int tzeq_flag_g,
             int y_opt_g,
             double so_llm_g,
             int ntss_ev_g,
             double[:,:,:] out_q,
             double[:,:,:] out_elv):
        
    cdef:
        double[::1,:,:] q_ev_g = np.empty([ntss_ev_g,mxncomp_g,nrch_g], dtype = np.double, order = 'F') 
        double[::1,:,:] elv_ev_g = np.empty([ntss_ev_g,mxncomp_g,nrch_g], dtype = np.double, order = 'F') 
    
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
            &z_ar_g[0,0],
            &bo_ar_g[0,0],
            &traps_ar_g[0,0],
            &tw_ar_g[0,0],
            &twcc_ar_g[0,0],
            &mann_ar_g[0,0],
            &manncc_ar_g[0,0],
            &so_ar_g[0,0],
            &dx_ar_g[0,0],
            &nhincr_m_g,
            &nhincr_f_g,
            &ufhlt_m_g[0,0,0],
            &ufqlt_m_g[0,0,0],
            &ufhlt_f_g[0,0,0],
            &ufqlt_f_g[0,0,0],
            &frnw_col,
            &frnw_g[0,0],
            &qlat_g[0,0,0],
            &ubcd_g[0,0],
            &dbcd_g[0],
            &cfl_g,
            &theta_g,
            &tzeq_flag_g,
            &y_opt_g,
            &so_llm_g,
            &ntss_ev_g,
            &q_ev_g[0,0,0],
            &elv_ev_g[0,0,0])
    
    # copy data from Fortran to Python memory view
    out_q[:,:,:] = q_ev_g[::1,:,:]
    out_elv[:,:,:] = elv_ev_g[::1,:,:]
    
cpdef object compute_diffusive(double dtini_g,
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
                 double[::1,:] z_ar_g,
                 double[::1,:] bo_ar_g,
                 double[::1,:] traps_ar_g,
                 double[::1,:] tw_ar_g,
                 double[::1,:] twcc_ar_g,
                 double[::1,:] mann_ar_g,
                 double[::1,:] manncc_ar_g,
                 double[::1,:] so_ar_g,
                 double[::1,:] dx_ar_g,
                 int nhincr_m_g,
                 int nhincr_f_g,
                 double[::1,:,:] ufhlt_m_g,
                 double[::1,:,:] ufqlt_m_g,
                 double[::1,:,:] ufhlt_f_g,
                 double[::1,:,:] ufqlt_f_g,
                 int frnw_col,
                 double[::1,:] frnw_g,
                 double[::1,:,:] qlat_g,
                 double[::1,:] ubcd_g,
                 double[::1] dbcd_g,
                 double cfl_g,
                 double theta_g,
                 int tzeq_flag_g,
                 int y_opt_g,
                 double so_llm_g,
                 int ntss_ev_g):
    
    cdef:
        double[:,:,:] out_q = np.empty([ntss_ev_g,mxncomp_g,nrch_g], dtype = np.double)
        double[:,:,:] out_elv = np.empty([ntss_ev_g,mxncomp_g,nrch_g], dtype = np.double)
    
    diffnw(dtini_g,
         t0_g,
         tfin_g,
         saveinterval_g,
         saveinterval_ev_g,
         dt_ql_g,
         dt_ub_g,
         dt_db_g,
         nts_ql_g,
         nts_ub_g,
         nts_db_g,
         mxncomp_g,
         nrch_g,
         z_ar_g,
         bo_ar_g,
         traps_ar_g,
         tw_ar_g,
         twcc_ar_g,
         mann_ar_g,
         manncc_ar_g,
         so_ar_g,
         dx_ar_g,
         nhincr_m_g,
         nhincr_f_g,
         ufhlt_m_g,
         ufqlt_m_g,
         ufhlt_f_g,
         ufqlt_f_g,
         frnw_col,
         frnw_g,
         qlat_g,
         ubcd_g,
         dbcd_g,
         cfl_g,
         theta_g,
         tzeq_flag_g,
         y_opt_g,
         so_llm_g,
         ntss_ev_g,
         out_q,
         out_elv)
    
    return np.asarray(out_q), np.asarray(out_elv)