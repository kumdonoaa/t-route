import numpy as np
import pathlib
import sys

ENV_IS_CL = False
if ENV_IS_CL:
    root = pathlib.Path("/", "content", "t-route")
elif not ENV_IS_CL:
    root = pathlib.Path("../..").resolve()
    sys.path.append("fast_reach")

## network and reach utilities
import troute.nhd_network_utilities_v02 as nnu
import diffusive

# variables needed to run diffusive model
dtini_g = 1.0  #inital time interval [secs]
t0_g = 0.0 #simulation starting time [hrs]
tfin_g = 1.0 #simulation ending time [hrs]
saveinterval_g = 1.0 #output publishing time interval [secs]
saveinterval_ev_g = saveinterval_g

dt_ql_g = 1.0 #time interval of wrf-hydro qlat input data [sec]
dt_ub_g = 1.0 #time interval of input data for upstream boundary condition [sec]
dt_db_g = 1.0 #time interval of input data for downstream boundary condition [sec]
            
            
nts_ql_g = 1 #number of simulation timesteps per qlateral time interval [-]
nts_ub_g = 1 #number of simulation timesteps per upstream boundary time interval [-]
nts_db_g = 1 #number of simulation timesteps per downstream boundary time interval [-]

mxncomp_g = 4 #maximum number of nodes in a reach
nrch_g = 6 #number network reaches

z_ar_g = np.empty([mxncomp_g,nrch_g], dtype = np.double) #altitude of each segment of each reach
bo_ar_g = np.empty([mxncomp_g,nrch_g], dtype = np.double) #channel bottom width
traps_ar_g = np.empty([mxncomp_g,nrch_g], dtype = np.double) #side slope of main channel bank
tw_ar_g = np.empty([mxncomp_g,nrch_g], dtype = np.double) #top width of main channel
twcc_ar_g = np.empty([mxncomp_g,nrch_g], dtype = np.double) #top width of the floodplain
mann_ar_g = np.empty([mxncomp_g,nrch_g], dtype = np.double) #manning's n for main channel
manncc_ar_g = np.empty([mxncomp_g,nrch_g], dtype = np.double) #manning's n for floodplain
so_ar_g = np.empty([mxncomp_g,nrch_g], dtype = np.double) #channel bottom slope
dx_ar_g = np.empty([mxncomp_g,nrch_g], dtype = np.double) #segment length

nhincr_m_g = 4 # number of water table elevation increments used for making uniform flow lookup table for MC
nhincr_f_g = 4 # number of water table elevation increments used for making uniform flow lookup table for floodplain
ufhlt_m_g = np.empty([mxncomp_g, nrch_g, nhincr_m_g], dtype = np.double) #water elevation array for uniform flow lookup table
ufqlt_m_g = np.empty([mxncomp_g, nrch_g, nhincr_m_g], dtype = np.double) #discharge array for uniform flow lookup table
ufhlt_f_g = np.empty([mxncomp_g, nrch_g, nhincr_f_g], dtype = np.double) #water elevation array for uniform flow lookup table
ufqlt_f_g = np.empty([mxncomp_g, nrch_g, nhincr_f_g], dtype = np.double) #water elevation array for uniform flow lookup table

frnw_col = 4 #numer of columns in fortran network map
frnw_g = np.empty([nrch_g, frnw_col], dtype = 'int32') # Fortran network map


qlat_g = np.empty([nts_ql_g, mxncomp_g, nrch_g], dtype = np.double) 
ubcd_g = np.empty([nts_ub_g, nrch_g], dtype = np.double) 
dbcd_g = np.empty([nts_ub_g], dtype = np.double) 

cfl_g = 0.98
theta_g = 0.25
tzeq_flag_g = 1
y_opt_g = 0
so_llm_g = 0.00001
ntss_ev_g = 10

out_q, out_elv = diffusive.compute_diffusive(dtini_g,
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
                                np.asfortranarray(z_ar_g),
                                np.asfortranarray(bo_ar_g),
                                np.asfortranarray(traps_ar_g),
                                np.asfortranarray(tw_ar_g),
                                np.asfortranarray(twcc_ar_g),
                                np.asfortranarray(mann_ar_g),
                                np.asfortranarray(manncc_ar_g),
                                np.asfortranarray(so_ar_g),
                                np.asfortranarray(dx_ar_g),
                                nhincr_m_g,
                                nhincr_f_g,
                                np.asfortranarray(ufhlt_m_g),
                                np.asfortranarray(ufqlt_m_g),
                                np.asfortranarray(ufhlt_f_g),
                                np.asfortranarray(ufqlt_f_g),
                                frnw_col,
                                np.asfortranarray(frnw_g),
                                np.asfortranarray(qlat_g),
                                np.asfortranarray(ubcd_g),
                                np.asfortranarray(dbcd_g),
                                cfl_g,
                                theta_g,
                                tzeq_flag_g,
                                y_opt_g,
                                so_llm_g,
                                ntss_ev_g,
                               )

