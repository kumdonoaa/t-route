import cython
import numpy as np
import pandas as pd
cimport numpy as np

from .fortran_wrappers cimport c_diffnw_cnt
from .. import diffusive_utils as diff_utils
import troute.nhd_network_utilities_v02 as nnu
import troute.nhd_network as nhd_network

# TO DO load some example inputs to test the module

@cython.boundscheck(False)
cdef void diffnw_cnt(double dtini_g,
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
             double[::1,:] iniq,    
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

    c_diffnw_cnt(
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
            &iniq[0,0],
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

cpdef object compute_diffusive_tst(
    int nsteps,
    int qts_subdivisions,
    list reaches_wTypes, # a list of tuples
    dict rconn,
    const long[:] data_idx,
    object[:] data_cols,
    const float[:,:] data_values,
    const float[:,:] initial_conditions,
    const float[:,:] qlat_values,
    list lake_numbers_col,
    const double[:,:] wbody_cols,
    dict waterbody_parameters,
    const int[:,:] reservoir_types,
    bint reservoir_type_specified,
    str model_start_time,
    const float[:,:] usgs_values,
    const int[:] usgs_positions_list,
    const float[:,:] lastobs_values,
    dict upstream_results={},
    bint assume_short_ts=False,
    bint return_courant=False,
    dict diffusive_parameters=False
    ):

    # segment connections dictionary
    connections = nhd_network.reverse_network(rconn)

    # network tailwater
    tw = list(nhd_network.headwaters(rconn))[0]

    # network reaches
    reach_list = []
    for i in reaches_wTypes:
        reach_list.append(i[0])

    # generate diffusive inputs
    diff_inputs = diff_utils.diffusive_input_data_v02(
        tw,
        connections,
        rconn,
        reach_list,
        diffusive_parameters,
        np.asarray(data_cols),
        np.asarray(data_idx),
        np.asarray(data_values),
        np.asarray(qlat_values),
        np.asarray(initial_conditions),
        upstream_results,
        qts_subdivisions,
        nsteps
        )

    # unpack/declare diffusive input variables
    cdef:
        double dtini_g = diff_inputs["dtini_g"]
        double t0_g = diff_inputs["t0_g"]
        double tfin_g  = diff_inputs["tfin_g"]
        double saveinterval_g = diff_inputs["saveinterval_g"]
        double saveinterval_ev_g = diff_inputs["saveinterval_ev_g"]
        double dt_ql_g = diff_inputs["dt_ql_g"]
        double dt_ub_g = diff_inputs["dt_ub_g"]
        double dt_db_g = diff_inputs["dt_db_g"]
        int nts_ql_g = diff_inputs["nts_ql_g"]
        int nts_ub_g = diff_inputs["nts_ub_g"]
        int nts_db_g = diff_inputs["nts_db_g"]
        int mxncomp_g = diff_inputs["mxncomp_g"]
        int nrch_g = diff_inputs["nrch_g"]
        double[::1,:] z_ar_g = np.asfortranarray(diff_inputs["z_ar_g"])
        double[::1,:] bo_ar_g = np.asfortranarray(diff_inputs["bo_ar_g"])
        double[::1,:] traps_ar_g = np.asfortranarray(diff_inputs["traps_ar_g"])
        double[::1,:] tw_ar_g = np.asfortranarray(diff_inputs["tw_ar_g"])
        double[::1,:] twcc_ar_g = np.asfortranarray(diff_inputs["twcc_ar_g"])
        double[::1,:] mann_ar_g = np.asfortranarray(diff_inputs["mann_ar_g"])
        double[::1,:] manncc_ar_g = np.asfortranarray(diff_inputs["manncc_ar_g"])
        double[::1,:] so_ar_g = np.asfortranarray(diff_inputs["so_ar_g"])
        double[::1,:] dx_ar_g = np.asfortranarray(diff_inputs["dx_ar_g"])
        int nhincr_m_g = diff_inputs["nhincr_m_g"]
        int nhincr_f_g = diff_inputs["nhincr_f_g"]
        double[::1,:,:] ufhlt_m_g = np.asfortranarray(diff_inputs["ufhlt_m_g"])
        double[::1,:,:] ufqlt_m_g = np.asfortranarray(diff_inputs["ufqlt_m_g"])
        double[::1,:,:] ufhlt_f_g = np.asfortranarray(diff_inputs["ufhlt_f_g"])
        double[::1,:,:] ufqlt_f_g = np.asfortranarray(diff_inputs["ufqlt_f_g"])
        int frnw_col = diff_inputs["frnw_col"]
        double[::1,:] frnw_g = np.asfortranarray(diff_inputs["frnw_g"], dtype = np.double)
        double[::1,:,:] qlat_g = np.asfortranarray(diff_inputs["qlat_g"])
        double[::1,:] ubcd_g = np.asfortranarray(diff_inputs["ubcd_g"])
        double[::1] dbcd_g = np.asfortranarray(diff_inputs["dbcd_g"])
        double cfl_g = diff_inputs["cfl_g"]
        double theta_g = diff_inputs["theta_g"]
        int tzeq_flag_g = diff_inputs["tzeq_flag_g"]
        int y_opt_g = diff_inputs["y_opt_g"]
        double so_llm_g = diff_inputs["so_llm_g"]
        int ntss_ev_g = diff_inputs["ntss_ev_g"]
        double[::1,:] iniq = np.asfortranarray(diff_inputs["iniq"])
        double[:,:,:] out_q = np.empty([ntss_ev_g,mxncomp_g,nrch_g], dtype = np.double)
        double[:,:,:] out_elv = np.empty([ntss_ev_g,mxncomp_g,nrch_g], dtype = np.double)

    # call diffusive compute kernel
    diffnw_cnt(dtini_g,
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
     iniq,
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

    # re-format outputs
    index_array, flowvelelv = diff_utils.unpack_output(
                                diff_inputs["pynw"],
                                diff_inputs["ordered_reaches"],
                                out_q,
                                out_elv
                                )

    return index_array, flowvelelv
