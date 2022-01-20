import numpy as np
from functools import partial, reduce
import troute.nhd_network as nhd_network
import os
import netCDF4 as nc
import pandas as pd

def adj_alt1(
    mx_jorder, ordered_reaches, geo_cols, geo_index, geo_data, dbfksegID, z_all
):
    """
    Adjust reach altitude data so that altitude of last node in reach is equal to that of the head segment of
    the neighboring downstream reach. 
    
    Parameters
    ----------
    mx_jorder -- (int) maximum network reach order
    geo_cols -- (ndarray of strs) column headers for geomorphic parameters data array (geo_data)
    geo_index -- (ndarray of int64s) row indices for geomorphic parameters data array (geo_data)
    geo_data --(ndarray of float32s) geomorphic parameters data arra
    dbfksegID -- (int) segment ID of fake node (=bottom node) of TW reach that hosts downstream boundary 
                        condition for a network that is being routed. 
    z_all -- (dict) adjusted altitude dictionary with placeholder values to be replaced
    
    Returns
    ----------
    z_all -- (dict) adjusted altitude dictionary
    """

    for x in range(mx_jorder, -1, -1):
        for head_segment, reach in ordered_reaches[x]:
            seg_list = reach["segments_list"]
            ncomp = reach["number_segments"]
            for seg in range(0, ncomp):
                segID = seg_list[seg]
                if seg == ncomp - 1 and seg_list.count(dbfksegID) == 0:
                    # At junction, the altitude of fake segment (=bottom node)  of an upstream reach
                    # is equal to that of the first segment (=top node) of the downstream reach

                    # head segment id of downstream reach from a junction
                    dsrchID = reach["downstream_head_segment"]

                    idx_dsrchID = np.where(geo_index == dsrchID)
                    idx_alt = np.where(geo_cols == "alt")
                    z_all[segID]["adj.alt"][0] = geo_data[idx_dsrchID, idx_alt]

                elif seg == ncomp - 1 and seg_list.count(dbfksegID) > 0:
                    # Terminal downstream fakesegment
                    ## AD HOC: need to be corrected later
                    segID2 = seg_list[seg - 1]
                    idx_segID2 = np.where(geo_index == segID2)
                    idx_so = np.where(geo_cols == "s0")
                    idx_dx = np.where(geo_cols == "dx")

                    So = geo_data[idx_segID2, idx_so]
                    dx = geo_data[idx_segID2, idx_dx]
                    z_all[segID]["adj.alt"][0] = z_all[segID2]["adj.alt"][0] - So * dx
                else:
                    idx_segID = np.where(geo_index == segID)
                    idx_alt = np.where(geo_cols == "alt")
                    z_all[segID]["adj.alt"][0] = geo_data[idx_segID, idx_alt]

    return z_all


def fp_mainstem_network_map(
     mainstem_headseg_list,mx_jorder, ordered_reaches, rchbottom_reaches, nrch_g, frnw_col, dbfksegID, pynw
):
    """
    Channel network mapping between Python and Fortran
    
    Parameters
    ----------
    mainstem_headseg_list -- (int) a list of link IDs of headsegs of related mainstem reaches 
    mx_jorder -- (int) maximum network reach order
    ordered_reaches -- (dict) reaches and reach metadata by junction order
    rchbottom_reaches -- (dict) reaches and reach metadata keyed by reach tail segment
    nrch_g -- (int) number of reaches in the network
    frnw_col -- (int) number of columns in the fortran network map
    dbfskegID -- (str) segment ID of fake (ghost) node at network downstream boundary 
    pynw -- (dict) ordered reach head segments
    
    Returns
    -------
    frnw_g -- (nparray of int) Fortran-Python network mapping array
    """

    #  Store headwater reach and upstream reaches above a junction
    #  as well as downstream reach after a junction
    #  into python-extension-fortran variables.
    frnw_g = np.zeros((nrch_g, frnw_col), dtype=int)
    frj = -1
    for x in range(mx_jorder, -1, -1):
        for head_segment, reach in ordered_reaches[x]:
            seg_list = reach["segments_list"]
            ncomp = reach["number_segments"]
            frj = frj + 1
            frnw_g[frj, 0] = ncomp

            if not reach["upstream_bottom_segments"]:
                # headwater reach
                frnw_g[frj, 2] = 0  # the number of upstream reaches
            else:
                # reaches before a junction
                nusrch = len(reach["upstream_bottom_segments"])
                frnw_g[frj, 2] = nusrch  # the number of upstream reaches
                usrch_bseg_list = list(reach["upstream_bottom_segments"])
                usrch_hseg_mainstem_j=-100 # ini.value of frj* of reach on mainstem that is just upstream of the current reach of frj
                i = 0
                for usrch in range(0, nusrch):
                    usrch_bseg_id = usrch_bseg_list[
                        usrch
                    ]  # upstream reach's bottom segment
                    usrch_hseg_id = rchbottom_reaches[usrch_bseg_id]["segments_list"][0]
                    # find Fortran js corresponding to individual usrchid
                    for j, sid in pynw.items():
                        if sid == usrch_hseg_id:
                            i = i + 1
                            frnw_g[frj, 2 + i] = j
                            # find if the selected headseg ID of an upstream reach belong to mainstem segment
                            if mainstem_headseg_list.count(usrch_hseg_id) > 0:
                                usrch_hseg_mainstem_j=j
                # store frj of mainstem headseg's reach in the just upstream of the current reach of frj
                frnw_g[frj, 2 + nusrch + 1] = usrch_hseg_mainstem_j
                            

            if seg_list.count(dbfksegID) > 0:
                # a reach where downstream boundary condition is set.
                frnw_g[
                    frj, 1
                ] = -100  # head_segment ID that is in terminal downstream reach.
                # That is, -100 indicates the reach of the head segment is
                # terminal downstream reach where ds.bcond. happens.
            else:
                # reach after a junction
                dsrch_hseg_id = reach["downstream_head_segment"]
                # fortran j index equivalent to dsrchID.
                frnw_g[frj, 1] = [
                    j for j, sid in pynw.items() if sid == dsrch_hseg_id[0]
                ][0]

    # Adust frnw_g element values according to Fortran-Python index relationship, that is Python i = Fortran i+1
    for frj in range(0, nrch_g):
        frnw_g[frj, 1] = frnw_g[frj, 1] + 1  # downstream reach index for frj reach
        if frnw_g[frj, 2] > 0:
            nusrch = frnw_g[frj, 2]
            for i in range(0, nusrch+1):  #<- +1 is for additional column of frj of mainstem headseg's reach in the upstream of the current reach of frj
                frnw_g[frj, 3 + i] = (
                    frnw_g[frj, 3 + i] + 1
                )  # upstream reach indicds for frj reach

    return frnw_g


def fp_chgeo_map(
    mx_jorder, ordered_reaches, geo_cols, geo_index, geo_data, z_all, mxncomp_g, nrch_g
):
    """
    Channel geometry data mapping between Python and Fortran
    
    Parameters
    ----------
    mx_jorder -- (int) maximum network reach order
    ordered_reaches -- (dict) reaches and reach metadata by junction order
    geo_cols -- (ndarray of strs) column headers for geomorphic parameters data array (geo_data)
    geo_index -- (ndarray of int64s) row indices for geomorphic parameters data array (geo_data)
    geo_data --(ndarray of float32s) geomorphic parameters data array
    z_all -- (dict) adjusted altitude dictionary
    mxncomp_g -- (int) maximum number of nodes in a reach
    nrch_g -- (int) number of reaches in the network
    
    Returns
    -------
    z_ar_g -- (numpy of float64s) altitude (meters)
    bo_ar_g -- (numpy of float64s) bottom width (meters)
    traps_ar_g -- (numpy of float64s) sideslope (m/m)
    tw_ar_g -- (numpy of float64s) top width (meters)
    twcc_ar_g -- (numpy of float64s) top width of compound channel (meters)
    mann_ar_g -- (numpy of float64s) manning's roughness (seconds/meters^1/3)
    mancc_ar_g -- (numpy of float64s) manning's roughness compound channel (seconds/meters^(1/3))
    so_ar_g -- (numpy of float64s) bottom slope (m/m)
    dx_ar_g -- (numpy of float64s) segment length (meters)
    """

    z_ar_g = np.zeros((mxncomp_g, nrch_g))
    bo_ar_g = np.zeros((mxncomp_g, nrch_g))
    traps_ar_g = np.zeros((mxncomp_g, nrch_g))
    tw_ar_g = np.zeros((mxncomp_g, nrch_g))
    twcc_ar_g = np.zeros((mxncomp_g, nrch_g))
    mann_ar_g = np.zeros((mxncomp_g, nrch_g))
    manncc_ar_g = np.zeros((mxncomp_g, nrch_g))
    so_ar_g = np.zeros((mxncomp_g, nrch_g))
    dx_ar_g = np.zeros((mxncomp_g, nrch_g))
    frj = -1
    for x in range(mx_jorder, -1, -1):
        for head_segment, reach in ordered_reaches[x]:
            seg_list = reach["segments_list"]
            ncomp = reach["number_segments"]
            frj = frj + 1
            for seg in range(0, ncomp):
                if seg == ncomp - 1:
                    segID = seg_list[seg - 1]
                else:
                    segID = seg_list[seg]

                idx_segID = np.where(geo_index == segID)

                idx_par = np.where(geo_cols == "bw")
                bo_ar_g[seg, frj] = geo_data[idx_segID, idx_par]

                idx_par = np.where(geo_cols == "cs")
                traps_ar_g[seg, frj] = geo_data[idx_segID, idx_par]

                idx_par = np.where(geo_cols == "tw")
                tw_ar_g[seg, frj] = geo_data[idx_segID, idx_par]

                idx_par = np.where(geo_cols == "twcc")
                twcc_ar_g[seg, frj] = geo_data[idx_segID, idx_par]

                idx_par = np.where(geo_cols == "n")
                mann_ar_g[seg, frj] = geo_data[idx_segID, idx_par]

                idx_par = np.where(geo_cols == "ncc")
                manncc_ar_g[seg, frj] = geo_data[idx_segID, idx_par]

                idx_par = np.where(geo_cols == "s0")
                so_ar_g[seg, frj] = geo_data[idx_segID, idx_par]

                idx_par = np.where(geo_cols == "dx")
                dx_ar_g[seg, frj] = geo_data[idx_segID, idx_par]

                segID1 = seg_list[seg]
                z_ar_g[seg, frj] = z_all[segID1]["adj.alt"][0]

    return (
        z_ar_g,
        bo_ar_g,
        traps_ar_g,
        tw_ar_g,
        twcc_ar_g,
        mann_ar_g,
        manncc_ar_g,
        so_ar_g,
        dx_ar_g,
    )


def fp_qlat_map(
    mx_jorder,
    ordered_reaches,
    nts_ql_g,
    geo_cols,
    geo_index,
    geo_data,
    qlat_data,
    qlat_g,
):
    """
    lateral inflow mapping between Python and Fortran
    
    Parameters
    ----------
    mx_jorder -- (int) maximum network reach order
    ordered_reaches -- (dict) reaches and reach metadata by junction order
    nts_ql_g -- (int) numer of qlateral timesteps
    geo_cols -- (ndarray of strs) column headers for geomorphic parameters data array (geo_data)
    geo_index -- (ndarray of int64) row indices for geomorphic parameters data array (geo_data)
    geo_data -- (ndarray of float32) geomorphic parameters data array
    qlat_data -- (ndarray of float32) qlateral data (m3/sec)
    qlat_g -- (ndarray of float32) empty qlateral array to be filled
    
    Returns
    -------
    qlat_g -- (ndarray of float32) qlateral array (m3/sec/m)
    
    Notes
    -----
    data in qlat_g are normalized by segment length with units of m2/sec = m3/sec/m
    """

    frj = -1
    for x in range(mx_jorder, -1, -1):
        for head_segment, reach in ordered_reaches[x]:
            seg_list = reach["segments_list"]
            ncomp = reach["number_segments"]
            frj = frj + 1
            for seg in range(0, ncomp):
                segID = seg_list[seg]
                for tsi in range(0, nts_ql_g):
                    if seg < ncomp - 1:

                        idx_segID = np.where(geo_index == segID)
                        idx_par = np.where(geo_cols == "dx")

                        tlf = qlat_data[idx_segID, tsi]  # [m^3/sec]
                        dx = geo_data[idx_segID, idx_par]  # [meter]
                        qlat_g[tsi, seg, frj] = tlf / dx  # [m^2/sec]

                    else:
                        qlat_g[
                            tsi, seg, frj
                        ] = 0.0  # seg=ncomp is actually for bottom node in Fotran code.
                        # And, lateral flow enters between adjacent nodes.

    return qlat_g


def fp_ubcd_map(frnw_g, pynw, nts_ub_g, nrch_g, ds_seg, upstream_inflows):
    """
    Upstream boundary condition mapping between Python and Fortran
    
    Parameters
    ----------
    frnw_g -- (nparray of int) Fortran-Python network mapping array
    pynw -- (dict) ordered reach head segments
    nrch_g -- (int) number of reaches in the network
    ds_seg -- (ndarray of int64) row indices for downstream segments recieving flows in upstream_flows
    upstream_inflows (ndarray of float32) upstream_inflows (m3/sec) to segments in ds_seg
    
    Returns
    -------
    ubcd_g -- (ndarray of float32) upstream boundary data (m3/sec)
    """

    ubcd_g = np.zeros((nts_ub_g, nrch_g))
    frj = -1

    #loop over every segment in network
    for frj in range(nrch_g):

        #if this is a head segment
        if frnw_g[frj, 2] == 0:  # the number of upstream reaches is zero.
            
            head_segment = pynw[frj]
            
            if head_segment in set(ds_seg):
                
                idx_dssegID = np.where(np.asarray(list(set(ds_seg))) == head_segment)
                for tsi in range(0, nts_ub_g):
                    ubcd_g[tsi, frj] = upstream_inflows[idx_dssegID, tsi]  # [m^3/s]

    return ubcd_g


def fp_dbcd_map(usgsID2tw=None, usgssDT=None, usgseDT=None, usgspCd=None):
    """
    Downststream boundary condition mapping between Python and Fortran using USGS stage observations
    
    Parameters
    ----------
    usgsID2tw -- (str) usgs site ID
    usgssDT -- (str) start data request date (yyyy-mm-dd) 
    usgseDT -- (str) end data request date (yyyy-mm-dd) 
    usgspCd --  (list of str) usgs parameter code 
    
    Returns
    -------
    nts_db_g -- (int) number of timesteps in downstream boundary data
    dbcd_g -- (ndarray of float64) downstream boundary data (meters)
    """

    # ** 1) downstream stage (here, lake elevation) boundary condition
    # from nwis_client.iv import IVDataService
    # install via: pip install hydrotools.nwis_client
    if usgsID2tw:
        try:
            from hydrotools.nwis_client.iv import IVDataService
        except ImportError as err:
            print(err, end="... ")
            print(
                "Please install hydrotools.nwis_client via: `pip install hydrotools.nwis_client`"
            )
            raise  # ensures program exit by re-raising the error.

        # from evaluation_tools.nwis_client.iv import IVDataService
        # Retrieve streamflow and stage data from two sites
        # Note: 1. Retrieved data all are based on UTC time zone (UTC is 4 hours ahead of Eastern Time during
        #          daylight saving time and 5 hours ahead during standard time)
        #       2. Retrieved data are always 1 hour ahead of stated startDT and 1 hour ahead of endDT,
        #          where starDT or endDT equal to yyyy-mm-dd 00:00.
        #       3. Also, retrieved data in 15 min so there are always four more data before startDT
        #          and four less data before endDT.
        #          For example, startDT='2018-08-01' and endDT='2018-09-01', then the retrieved data starts by
        #          2018-07-31-23:00:00
        #          2018-07-31-23:15:00
        #          2018-07-31-23:30:00
        #          2018-07-31-23:45:00
        #          2018-08-01-00:00:00
        #               .......
        #          2018-08-31-22:00:00
        #          2018-08-31-22:15:00
        #          2018-08-31-22:30:00
        #          2018-08-31-22:45:00
        #          2018-08-31-23:00:00
        #       4. '00060' for discharge [ft^3/s]
        #          '00065' for stage [ft]
        #          '62614' for Elevation, lake/res,NGVD29 [ft]
        ivds = IVDataService()
        observations_data = ivds.get(
            sites=usgsID2tw,  # sites='01646500,0208758850',
            startDT=usgssDT,  #'2018-08-01',
            endDT=usgseDT,  #'2020-09-01',
            # parameterCd='62614'
            parameterCd=usgspCd,
        )
        nts_db_g = len(observations_data) - 4
        # ** 4 is added to make data used here has its date time as from startDT 00:00 to (endDT-1day) 23:00, UTC
        # ** usgs data at this site uses NGVD1929 feet datum while 'alt' of RouteLink uses NAD88 meter datum.
        #  -> Has to convert accordingly !!!!
        #
        # source: https://pubs.usgs.gov/sir/2010/5040/section.html
        # Over most USGS study area it is used that NGVD = NAVD88 - 3.6 feet
        dbcd_g = np.zeros(nts_db_g)
        for tsi in range(0, nts_db_g):
            i = tsi + 4
            dmy = observations_data.iloc[i, 4]
            dbcd_g[tsi] = 0.3048 * (
                dmy + 3.6
            )  # accuracy with +-0.5feet for 95 percent of USGS study area.
            # 0.3048 to covert ft to meter. [meter]
    else:
        nts_db_g = 1
        dbcd_g = -np.ones(nts_db_g)
        
    return nts_db_g, dbcd_g

def fp_naturalxsec_map(
                ordered_reaches,
                mainstem_headseg_list, 
                inland_bathyNC, 
                comid_bathy, 
                geo_index, 
                geo_cols, 
                geo_data, 
                mx_jorder, 
                mxncomp_g, 
                nrch_g,
                dbfksegID):
    """
    natural cross section mapping between Python and Fortran using eHydro_ned_cross_sections data
    
    Parameters
    ----------
    ordered_reaches -- (dict) reaches and reach metadata by junction order
    mainstem_headseg_list -- (int) a list of link IDs of headsegs of related mainstem reaches 
    inland_bathyNC --  bathymetry NC data of channel cross section
    geo_index -- (ndarray of int64s) row indices for geomorphic parameters data array (geo_data)
    mx_jorder -- (int) max junction order
    mxncomp_g -- (int) maximum number of nodes in a reach
    nrch_g -- (int) number of reaches in the network
    dbfksegID -- (int) segment ID of fake node (=bottom node) of TW reach that hosts downstream boundary 
                        condition for a network that is being routed. 
    
    Returns
    -------
    x_bathy_g -- (numpy of float64s) lateral distance of bathy data points
    z_bathy_g -- (numpy of float64s) elevation of bathy data points
    size_bathy_g -- (integer) the nubmer of bathy data points of each cross section
    mxnbathy_g -- (integer) maximum size of bathy data points
    """    
    # take only related data from bathy NC file
    inland_bathyDF = pd.DataFrame()
    inland_bathyDF['comid'] =  pd.DataFrame(np.ma.filled(inland_bathyNC['comid'][:]))
    # this is lateral distance from x=0 where the NWM line crosses the cross section later line. 
    # For left side of x=0 (when looking upstream-to-downstream), the lateral distance takes negative values 
    # for left side of main channel or left floodplains whiel positive for right side of main ch. or right floodplains.
    inland_bathyDF['x_bathy'] =  pd.DataFrame(np.ma.filled(inland_bathyNC['xid_d'][:]))
    inland_bathyDF['z_bathy'] =   pd.DataFrame(np.ma.filled(inland_bathyNC['z'][:]))
    inland_bathyDF['n'] =   pd.DataFrame(np.ma.filled(inland_bathyNC['n'][:]))
    
    # find maximum size of x or z bathy data
    sizelist=[]
    ncomid= len(comid_bathy)
    for i in range(0,ncomid): # for bathy data of arbitarily selected three segments
        inland_bathyDF_subset = inland_bathyDF[inland_bathyDF['comid'] == comid_bathy[i]]    
        size= inland_bathyDF_subset['x_bathy'].size
        sizelist.append(size)
    mxnbathy_g=max(sizelist)

    x_bathy_g = np.zeros((mxnbathy_g, mxncomp_g, nrch_g))
    z_bathy_g = np.zeros((mxnbathy_g, mxncomp_g, nrch_g))
    mann_bathy_g = np.zeros((mxnbathy_g, mxncomp_g, nrch_g))
    size_bathy_g= np.zeros((mxncomp_g, nrch_g), dtype=int)
    
    #import pdb; pdb.set_trace()
    
    frj = -1
    for x in range(mx_jorder, -1, -1):
        for head_segment, reach in ordered_reaches[x]:
            seg_list = reach["segments_list"]
            ncomp = reach["number_segments"]
            frj = frj + 1
            headsegID= seg_list[0]                  
    ## **** for now, bathy data here doesn't cover Cape fear river we're experimenting with. So, we picked
    ## **** arbitrarily three comid from the bathy NC data and get the bathy data to pass to three segments
    ## **** we're considering now. Basically, three mainstem reaches with each having only one segment.                
            if mainstem_headseg_list.count(headsegID) > 0:  
                # when selected segID is mainstem:
                #import pdb; pdb.set_trace()                   
                for seg in range(0,ncomp):
                    if seg==ncomp-1 and x > 0:
                        #In node-configuration, bottom node takes bathy of the first segment of the downtream reach
                        # except TW reach where the bottom node bathy is interpolated by bathy of the last segment 
                        # with so*0.5*dx of 
                        segID= reach["downstream_head_segment"]
                    else:
                        # In node configuration, for example, for 2 segments of a reach, 
                        # 1st node(=top node) takes bathy of the first segment and
                        # 2nd node takes bathy of the second segment
                        # 3rd node takes bathy of, as explained right above, the firt segment 
                        segID= seg_list[seg]
                    
                     ## **** all hypothetical coding
                    if x==2 and seg==0:
                        idx_comid=0
                    if x==2 and seg==1:
                        idx_comid=1                        
                    if x==1 and seg==0:
                        idx_comid=1
                    if x==1 and seg==1:
                        idx_comid=2
                    if x==0 and seg==0:
                        idx_comid=2
                    if x==0 and seg==1:
                        #when tail water node, z_bathy will be adjusted by so*dx of related segment
                        idx_comid=2                                      
                    
                    inland_bathyDF_subset = inland_bathyDF[inland_bathyDF['comid'] == comid_bathy[idx_comid]]
                    size_bathy_g[seg, frj]= inland_bathyDF_subset['x_bathy'].size
                    
                    for idp in range(0, size_bathy_g[seg, frj]):
                        x_bathy_g[idp, seg, frj]  = inland_bathyDF_subset [['x_bathy']].iloc[idp]
                        z_bathy_g[idp, seg, frj] =  inland_bathyDF_subset [['z_bathy']].iloc[idp]
                        mann_bathy_g[idp, seg, frj] = inland_bathyDF_subset [['n']].iloc[idp]
                        
                if seg == ncomp - 1 and seg_list.count(dbfksegID) > 0:
                    # At tailwater bottom node (as # of segment = # node +1), z_bathy need to be adjusted by so*dx of the last segment
                    segID2 = seg_list[seg - 1]
                    idx_segID2 = np.where(geo_index == segID2)
                    idx_so = np.where(geo_cols == "s0")
                    idx_dx = np.where(geo_cols == "dx")

                    So = geo_data[idx_segID2, idx_so]
                    dx = geo_data[idx_segID2, idx_dx]
                    
                    for idp in range(0, size_bathy_g[seg, frj]):
                        z_bathy_g[idp, seg, frj] = z_bathy_g[idp, seg, frj] - So * dx   
                    #print(So)
                    #print(dx)
                    #print(So*dx)

   
    return x_bathy_g, z_bathy_g, mann_bathy_g, size_bathy_g, mxnbathy_g


def diffusive_input_data_v01(
    tw,
    connections,
    rconn,
    reach_list,
    mainstem_headseg_list,
    q_wrf_data,    
    diffusive_parameters,
    geo_cols,
    geo_index,
    geo_data,
    qlat_data,
    initial_conditions,
    upstream_results,
    qts_subdivisions,
    nsteps,
    dt,
    lake_segs,
    wbody_params,
    comid_bathy,
    inland_bathyNC,
):
    """
    Build input data objects for diffusive wave model
    
    Parameters
    ----------
    tw -- (int) Tailwater segment ID
    connections -- (dict) donwstream connections for each segment in the network
    rconn -- (dict) upstream connections for each segment in the network
    reach_list -- (list of lists) lists of segments comprising different reaches in the network
    diffusive_parametters -- (dict) Diffusive wave model parameters
    geo_cols -- (ndarray of strs) column headers for geomorphic parameters data array (geo_data)
    geo_index -- (ndarray of int64s) row indices for geomorphic parameters data array (geo_data)
    geo_data --(ndarray of float32s) geomorphic parameters data array
    qlat_data -- (ndarray of float32) qlateral data (m3/sec)
    initial_conditions -- (ndarray of float32) initial flow (m3/sec) and depth (m above ch bottom) states for network nodes
    upstream_results -- (dict) with values of 1d arrays upstream flow, velocity, and depth   
    qts_subdivisions -- (int) number of qlateral timestep subdivisions 
 
    Returns
    -------
    diff_ins -- (dict) formatted inputs for diffusive wave model
    """
    
    dt_ql_g = 3600 #dt * qts_subdivisions    # lateral inflow timestep (sec)    
    dt_ub_g = 3600 #dt    # upstream boundary condition timestep (sec)    
    dt_db_g = 3600 #dt * qts_subdivisions    # downstream boundary condition timestep (sec)    
    saveinterval_tu = 60.0    # time interval at which flow and depth simulations are written out by Tulane diffusive model    
    saveinterval_cnt = 12*dt  #* qts_subdivisions    # time interval at which depth is written out by cnt model
    dtini_g = dt    # simulation time step (for cn-mod, it changes but fixed for cnt and cnx)
    dt_qtrib_g = 60.0  #3600.0    #tributary (including mainstem upstream boundary) data time step [sec]
    t0_g = 0.0     # simulation start hr **set to zero for Fortran computation
    tfin_g = 72.0  # (dt * nsteps)/60/60 # [hr]
    
    # USGS data related info.
    #usgsID = diffusive_parameters.get("usgsID", None)
    #seg2usgsID = diffusive_parameters.get("link2usgsID", None)
    #usgssDT = diffusive_parameters.get("usgs_start_date", None)
    #usgseDT = diffusive_parameters.get("usgs_end_date", None)
    #usgspCd = diffusive_parameters.get("usgs_parameterCd", None)
    
    # CNX: store the time variable values in the above to an array
    '''
    timestep_ar_g = np.zeros(8)
    timestep_ar_g[0]= dtini_g
    timestep_ar_g[1]= t0_g
    timestep_ar_g[2]= tfin_g
    timestep_ar_g[3]= 0.0 # not used  
    timestep_ar_g[4]= dt_ql_g
    timestep_ar_g[5]= dt_ub_g
    timestep_ar_g[6]= dt_db_g
    timestep_ar_g[7]= dt_qtrib_g   
    '''
    '''
    # CN-mod: store the time variable values in the above to an array
    timestep_ar_g = np.zeros(8)
    timestep_ar_g[0]= dtini_g
    timestep_ar_g[1]= t0_g
    timestep_ar_g[2]= tfin_g
    timestep_ar_g[3]= saveinterval_tu
    timestep_ar_g[4]= dt_ql_g
    timestep_ar_g[5]= dt_ub_g
    timestep_ar_g[6]= dt_db_g
    timestep_ar_g[7]= dt_qtrib_g   
    '''
    # CNT: store the time variable values in the above to an array
    timestep_ar_g = np.zeros(8)
    timestep_ar_g[0]= dtini_g
    timestep_ar_g[1]= t0_g
    timestep_ar_g[2]= tfin_g
    timestep_ar_g[3]= saveinterval_cnt
    timestep_ar_g[4]= dt_ql_g
    timestep_ar_g[5]= dt_ub_g
    timestep_ar_g[6]= dt_db_g
    timestep_ar_g[7]= dt_qtrib_g   
    #'''
    
    
    
    # diffusive parameters
    #cfl_g = diffusive_parameters.get("courant_number_upper_limit", None)
    #theta_g = diffusive_parameters.get("theta_parameter", None)
    #tzeq_flag_g = diffusive_parameters.get("chgeo_computation_flag", None)
    #y_opt_g = diffusive_parameters.get("water_elevation_computation_flag", None)
    #so_llm_g = diffusive_parameters.get("bed_slope_lower_limit", None) 
    
    # CNX Parameters
    '''   
    paradim=6 
    para_ar_g= np.zeros(paradim)
    para_ar_g[0]= 0.95 #lower limit of Courant number.
    para_ar_g[1]= 0.05  # lower limit of celerity.
    para_ar_g[2]= 0.002    # lower limit of discharge.
    para_ar_g[3]= 0.00001    # lower limit of channel bed slope.
    para_ar_g[4]= 12 # upper limit of the number of uniformly distributed sub-nodes, including ghost nodes, on each stream segment.
    para_ar_g[5]= 2  # number of sub-nodes added to existing uniformly distributed sub-nodes, used to provide adequate discharge downstream boundary condition. 
    '''
    '''
    # CN-mod parameters
    paradim=10
    para_ar_g= np.zeros(paradim)
    para_ar_g[0]= 0.95    # Courant number (default: 0.95)
    para_ar_g[1]= 0.5     # lower limit of celerity (default: 0.5)
    para_ar_g[2]= 50.0    # lower limit of diffusivity (default: 50)
    para_ar_g[3]= 5000.0  # upper limit of diffusivity (default: 1000)
    para_ar_g[4]= -15.0   # lower limit of dimensionless diffusivity, used to determine b/t normal depth and diffusive depth
    para_ar_g[5]= -10.0   #upper limit of dimensionless diffusivity, used to determine b/t normal depth and diffusive depth
    para_ar_g[6]= 1.0     # 0:run Bisection to compute water level; 1: Newton Raphson (default: 1.0)
    para_ar_g[7]= 0.02831   # lower limit of discharge (default: 0.02831 cms)
    para_ar_g[8]= 0.0001    # lower limit of channel bed slope (default: 0.0001)
    para_ar_g[9]= 1.0     # weight in numerically computing 2nd derivative: 0: explicit, 1: implicit (default: 1.0)
    '''    
    # CNT parameters
    paradim=9
    para_ar_g= np.zeros(paradim)
    para_ar_g[0]= 0.5    # lower limit of celerity
    para_ar_g[1] = 50.0   # lower limit of diffusivity
    para_ar_g[2] = 400.0  # upper limit of diffusivity
    para_ar_g[3] = 0.1    # lower limit of dimensionless space step, used to reduce ocillation of hydrograph
    para_ar_g[4] = 0.3    # lower limit of dimensionless time step, used to reduce ociallation of hydrograph
    para_ar_g[5] = 2  # the number of ghost time node, used to get adequate dischage boundary condition in time.
    para_ar_g[6] = 0.02831 # lower limit of discharge
    para_ar_g[7] = 0.0001 # lower limit of channel bed slope
    para_ar_g[8] =  3  # 1: bisection, 2: simple iterative with contribution from dU/dX, 3: Newton Raphson for water depth

    # number of reaches in network
    nrch_g = len(reach_list)

    # maximum number of nodes in a reach
    mxncomp_g = 0
    for r in reach_list:
        nnodes = len(r) + 1
        if nnodes > mxncomp_g:
            mxncomp_g = nnodes
    
    ds_seg = []
    offnet_segs = []
    upstream_flow_array = np.zeros((len(ds_seg), nsteps+1))
    if upstream_results:
        
        # create a list of segments downstream of offnetwork upstreams [ds_seg]
        # and a list of offnetwork upstream segments [offnet_segs]
        inv_map = nhd_network.reverse_network(rconn)
        for seg in upstream_results:
            ds_seg.append(inv_map[seg][0])
            offnet_segs.append(seg)
        
        # populate an array of upstream flows (boundary condtions)
        upstream_flow_array = np.zeros((len(set(ds_seg)), nsteps+1))
        for j, seg in enumerate(set(ds_seg)):
            
            # offnetwork-upstream connections
            us_segs = rconn[seg]
            
            # sum upstream flows and initial conditions
            usq = np.zeros((len(us_segs), nsteps))
            us_iniq = 0
            for k, s in enumerate(us_segs):
                usq[k] = upstream_results[s]['results'][::3]
                
                if s in lake_segs:
                    # initial conditions from wbody_param array
                    idx_segID = np.where(np.asarray(lake_segs) == s)
                    us_iniq += wbody_params[idx_segID,9]
                else:
                    # initial conditions from initial_conditions array
                    idx_segID = np.where(geo_index == s)
                    us_iniq += initial_conditions[idx_segID,0]
            
            # write upstream flows to upstream_flow_array
            upstream_flow_array[j,1:] = np.sum(usq, axis = 0)
            upstream_flow_array[j,0] = us_iniq
    
    # Order reaches by junction depth
    path_func = partial(nhd_network.split_at_waterbodies_and_junctions, set(offnet_segs), rconn)
    tr = nhd_network.dfs_decomposition_depth_tuple(rconn, path_func)    
    
    jorder_reaches = sorted(tr, key=lambda x: x[0])
    mx_jorder = max(jorder_reaches)[0]  # maximum junction order of subnetwork of TW
    
    ordered_reaches = {}
    rchhead_reaches = {}
    rchbottom_reaches = {}
    z_all = {}
    for o, rch in jorder_reaches:

        # add one more segment(fake) to the end of a list of segments to account for node configuration.
        fksegID = int(str(rch[-1]) + str(2))
        rch.append(fksegID)

        # additional segment(fake) to upstream bottom segments
        if any(j in rconn[rch[0]] for j in offnet_segs):
            fk_usbseg = []
        else:
            fk_usbseg = [int(str(x) + str(2)) for x in rconn[rch[0]]]            

        if o not in ordered_reaches:
            ordered_reaches.update({o: []})
        ordered_reaches[o].append(
            [
                rch[0],
                {
                    "number_segments": len(rch),
                    "segments_list": rch,
                    "upstream_bottom_segments": fk_usbseg,
                    "downstream_head_segment": connections[rch[-2]],
                },
            ]
        )

        if rch[0] not in rchhead_reaches:
            # a list of segments for a given head segment
            rchhead_reaches.update(
                {rch[0]: {"number_segments": len(rch), "segments_list": rch}}
            )
            # a list of segments for a given bottom segment
            rchbottom_reaches.update(
                {rch[-1]: {"number_segments": len(rch), "segments_list": rch}}
            )
        # for channel altitude adjustment
        z_all.update({seg: {"adj.alt": np.zeros(1)} for seg in rch})

        # cahnnel geometry data
        a = np.where(geo_cols == "cs")
        geo_data[:, a] = 1.0 / geo_data[:, a]
    

    # --------------------------------------------------------------------------------------
    #                                 Step 0-3
    #    Adjust altitude so that altitude of the last sement of a reach is equal to that
    #    of the first segment of its downstream reach right after their common junction.
    # --------------------------------------------------------------------------------------
    dbfksegID = int(str(tw) + str(2))

    adj_alt1(
        mx_jorder, ordered_reaches, geo_cols, geo_index, geo_data, dbfksegID, z_all
    )
    
    # --------------------------------------------------------------------------------------
    #                                 Step 0-4
    #     Make Fortran-Python channel network mapping variables.
    # --------------------------------------------------------------------------------------

    # build a list of head segments in descending reach order [headwater -> tailwater]
    pynw = {}
    frj = -1
    for x in range(mx_jorder, -1, -1):
        for head_segment, reach in ordered_reaches[x]:
            frj = frj + 1
            pynw[frj] = head_segment
    

    frnw_col = diffusive_parameters.get("fortran_nework_map_col_number", None)
    frnw_g = fp_mainstem_network_map(
        mainstem_headseg_list, mx_jorder, ordered_reaches, rchbottom_reaches, nrch_g, frnw_col, dbfksegID, pynw
    )

    # covert data type from integer to float for frnw
    dfrnw_g = np.zeros((nrch_g, frnw_col), dtype=float)
    for j in range(0, nrch_g):
        for col in range(0, frnw_col):
            dfrnw_g[j, col] = float(frnw_g[j, col])
    
    np.savetxt("./temp/frnw_ar.txt", frnw_g, fmt='%10i')
    #np.savetxt("./temp/pynw_ar.txt", pynw, fmt='%20i')
    # ---------------------------------------------------------------------------------
    #                              Step 0-5
    #                  Prepare channel geometry data
    # ---------------------------------------------------------------------------------
    (
        z_ar_g,
        bo_ar_g,
        traps_ar_g,
        tw_ar_g,
        twcc_ar_g,
        mann_ar_g,
        manncc_ar_g,
        so_ar_g,
        dx_ar_g,
    ) = fp_chgeo_map(
        mx_jorder,
        ordered_reaches,
        geo_cols,
        geo_index,
        geo_data,
        z_all,
        mxncomp_g,
        nrch_g,
    )
    # saved txt file has row for segment and colunum for frj
    np.savetxt("./temp/z_ar.txt", z_ar_g, fmt='%15.5f')
    np.savetxt("./temp/bo_ar.txt", bo_ar_g, fmt='%15.5f')
    np.savetxt("./temp/traps_ar.txt", traps_ar_g, fmt='%15.5f')
    np.savetxt("./temp/tw_ar.txt", tw_ar_g, fmt='%15.5f')
    np.savetxt("./temp/twcc_ar.txt", twcc_ar_g, fmt='%15.5f')
    np.savetxt("./temp/mann_ar.txt", mann_ar_g, fmt='%15.5f')
    np.savetxt("./temp/manncc_ar.txt", manncc_ar_g, fmt='%15.5f')
    np.savetxt("./temp/so_ar.txt", so_ar_g, fmt='%15.6f') 
    np.savetxt("./temp/dx_ar.txt", dx_ar_g, fmt='%15.5f')   
    
    # ---------------------------------------------------------------------------------
    #                              Step 0-6
    #                  Prepare initial conditions data
    # ---------------------------------------------------------------------------------
    output_path='./temp'
    
    iniq = np.zeros((mxncomp_g, nrch_g))
    frj = -1
    for x in range(mx_jorder, -1, -1):
        for head_segment, reach in ordered_reaches[x]:
            seg_list = reach["segments_list"]
            ncomp = reach["number_segments"]
            frj = frj + 1
            for seg in range(0, ncomp):
                if seg == ncomp - 1:
                    segID = seg_list[seg - 1]
                else:
                    segID = seg_list[seg]
                    
                idx_segID = np.where(geo_index == segID)
                iniq[seg, frj] = initial_conditions[idx_segID, 0]
                if iniq[seg, frj]<0.0001:
                    iniq[seg, frj]=0.0001
    
    frj= -1
    with open("./temp/iniq_ar.txt",'a') as pyqini:    
        for x in range(mx_jorder,-1,-1):   
            for head_segment, reach in ordered_reaches[x]:                  
                seg_list= reach["segments_list"]
                ncomp= reach["number_segments"]             
                frj= frj+1     
                for seg in range(0,ncomp):                               
                    segID= seg_list[seg]                  
                    dmy1= iniq[seg, frj]
                    pyqini.write("%s %s %s\n" %\
                                (str(frj+1), str(seg+1), str(dmy1)) )                   
                                # <- +1 is added to frj for accommodating Fortran indexing. 

    # ---------------------------------------------------------------------------------
    #                              Step 0-7

    #                  Prepare lateral inflow data
    # ---------------------------------------------------------------------------------
    nts_ql_g = (
        int((tfin_g - t0_g) * 3600.0 / dt_ql_g)
    )  # the number of the entire time steps of lateral flow data

    qlat_g = np.zeros((nts_ql_g, mxncomp_g, nrch_g))

    fp_qlat_map(
        mx_jorder,
        ordered_reaches,
        nts_ql_g,
        geo_cols,
        geo_index,
        geo_data,
        qlat_data,
        qlat_g,
    )
       
    frj= -1
    with open("./temp/qlat_ar.txt",'a') as pyql:
    #if 1==1:
        for x in range(mx_jorder,-1,-1):   
            for head_segment, reach in ordered_reaches[x]:                  
                seg_list= reach["segments_list"]
                ncomp= reach["number_segments"]             
                frj= frj+1     
                for seg in range(0,ncomp):                               
                    segID= seg_list[seg]
                    for tsi in range (0,nts_ql_g):
                        t_min= dt_ql_g/60.0*float(tsi)                        
                        dmy1= qlat_g[tsi, seg, frj]
                        pyql.write("%s %s %s %s\n" %\
                                (str(frj+1), str(seg+1), str(t_min), str(dmy1)) )                   
                                # <- +1 is added to frj for accommodating Fortran indexing. 
    # ---------------------------------------------------------------------------------
    #                              Step 0-8

    #       Prepare upstream boundary (top segments of head basin reaches) data
    # ---------------------------------------------------------------------------------
    #nts_ub_g = upstream_flow_array.shape[1]
    nts_ub_g = int((tfin_g - t0_g) * 3600.0 / dt_ub_g)
    ubcd_g = fp_ubcd_map(frnw_g, pynw, nts_ub_g, nrch_g, ds_seg, upstream_flow_array)


        
    # ---------------------------------------------------------------------------------
    #                              Step 0-9

    #       Prepare downstrea boundary (bottom segments of TW reaches) data
    # ---------------------------------------------------------------------------------
    '''
    if seg2usgsID:
        if tw in seg2usgsID:
            ipos = seg2usgsID.index(tw)
            usgsID2tw = usgsID[ipos]
        else:
            usgsID2tw=None
    else:
        usgsID2tw=None
    
    nts_db_g, dbcd_g = fp_dbcd_map(usgsID2tw, usgssDT, usgseDT, usgspCd)
    '''
    nts_db_g = 1
    dbcd_g = -np.ones(nts_db_g)   

    # ---------------------------------------------------------------------------------------------
    #                              Step 0-9-2

    #       Prepare tributary q time series data generated by MC that flow into a juction boundary  
    # ---------------------------------------------------------------------------------------------
    #dt_qtrib_g=3600.0 # MC chrout file time step in [sec]  ** TO DO: move to yaml or else.
    nts_qtrib_g = int((tfin_g - t0_g) * 3600.0 / dt_qtrib_g)
    
    qtrib_g = np.zeros((nts_qtrib_g, nrch_g))
    frj = -1
    for x in range(mx_jorder, -1, -1):
        for head_segment, reach in ordered_reaches[x]:
            frj = frj + 1
            if mainstem_headseg_list.count(head_segment) == 0:                
                idx_trib = np.where(geo_index == head_segment)               
                #qtrib_g[:,frj]= q_wrf_data[idx_trib, 0:nts_qtrib_g] 
                qtrib_g[:,frj]= 0.0    #test  
                #qtrib=[value for count, value in enumerate(q_wrf_data[idx_trib]) if count<nts_qtrib_g]
                #qtrib2=np.asarray(qtrib)
                #qtrib_g[:,frj]= qtrib2
                #for i in range(0,nts_qtrib_g):
                #    qtrib_g[i,frj]=q_wrf_data[idx_trib,i]
    # synthetic inflow at mainstem upper boundary
    import math
    qsi= 100.0
    wsi= 80.0
    zsi= 0.6
    for n in range(0,nts_qtrib_g):
        tmin= dt_qtrib_g/60.0*n
        if n==0:
            tmin= 0.1
        qtrib_g[n,1]= qsi*math.exp(zsi*(2.0 - tmin/wsi - wsi/tmin))/((tmin/wsi)**(3.0/2.0))   
    
    frj= -1
    with open(os.path.join(output_path,"qtrib_ar.txt"),'a') as pyqtrib:
        for x in range(mx_jorder,-1,-1):   
            for head_segment, reach in ordered_reaches[x]:                  
                seg_list= reach["segments_list"]
                ncomp= reach["number_segments"]             
                frj= frj+1     
                for tsi in range (0, nts_qtrib_g):
                    t_min=  dt_qtrib_g/60.0*float(tsi)                        
                    dmy1= qtrib_g[tsi, frj]
                    pyqtrib.write("%s %s %s\n" %\
                                 (str(frj+1), str(t_min), str(dmy1)) )  
                                 # <- +1 is added to frj for accommodating Fortran indexing.      
    # ---------------------------------------------------------------------------------
    #                              Step 0-10

    #                 Prepare cross section bathymetry data
    # ---------------------------------------------------------------------------------
    
    x_bathy_g, z_bathy_g, mann_bathy_g, size_bathy_g, mxnbathy_g = fp_naturalxsec_map(        
                                                                   ordered_reaches,                                             
                                                                   mainstem_headseg_list, 
                                                                   inland_bathyNC, 
                                                                   comid_bathy, 
                                                                   geo_index, 
                                                                   geo_cols, 
                                                                   geo_data, 
                                                                   mx_jorder,
                                                                   mxncomp_g, 
                                                                   nrch_g,
                                                                   dbfksegID)
    
    #import pdb;  pdb.set_trace()
    frj= -1
    with open(os.path.join(output_path,"bathy_data.txt"),'a') as pybathy:
        for x in range(mx_jorder,-1,-1):   
            for head_segment, reach in ordered_reaches[x]:                  
                seg_list= reach["segments_list"]
                ncomp= reach["number_segments"]             
                frj= frj+1
                headsegID= seg_list[0]                 
                if mainstem_headseg_list.count(headsegID) > 0:
                    for seg in range(0,ncomp):  
                        #import pdb; pdb.set_trace()
                        for idp in range(0, size_bathy_g[seg, frj]):
                            dmy1 = x_bathy_g[idp, seg, frj]
                            dmy2 = z_bathy_g[idp, seg, frj]
                            dmy3 = mann_bathy_g[idp, seg, frj]                            
                            pybathy.write("%s %s %s %s %s %s\n" %\
                                 (str(idp+1), str(seg+1), str(frj+1), str(dmy1), str(dmy2), str(dmy3))) 
    # this doesn't work:  
    #np.savetxt("./temp/x_bathy_g.txt", x_bathy_g, fmt='%15.5f')  
    #np.savetxt("./temp/z_bathy_g.txt", z_bathy_g, fmt='%15.5f') 
    #np.savetxt("./temp/mann_bathy_g.txt", mann_bathy_g, fmt='%15.5f') 
    #np.savetxt("./temp/size_bathy_g.txt", size_bathy_g, fmt='%15.5f')  
 
    
    # TODO: Call uniform flow lookup table creation kernel    
    # ---------------------------------------------------------------------------------
    #                              Step 0-11

    #                       Build input dictionary
    # ---------------------------------------------------------------------------------
    ntss_ev_g = int((tfin_g - t0_g) * 3600.0 / dt) + 1

    # build a dictionary of diffusive model inputs and helper variables
    diff_ins = {}

    # model input parameters
        # model input parameters
    diff_ins["timestep_ar_g"]= timestep_ar_g  
    diff_ins["nts_ql_g"] = nts_ql_g
    diff_ins["nts_ub_g"] = nts_ub_g
    diff_ins["nts_db_g"] = nts_db_g
    diff_ins["ntss_ev_g"] = ntss_ev_g
    
    diff_ins["nts_qtrib_g"] = nts_qtrib_g
    
    diff_ins["mxncomp_g"] = mxncomp_g
    diff_ins["nrch_g"] = nrch_g
    diff_ins["z_ar_g"] = z_ar_g
    diff_ins["bo_ar_g"] = bo_ar_g
    diff_ins["traps_ar_g"] = traps_ar_g
    diff_ins["tw_ar_g"] = tw_ar_g
    diff_ins["twcc_ar_g"] = twcc_ar_g
    diff_ins["mann_ar_g"] = mann_ar_g
    diff_ins["manncc_ar_g"] = manncc_ar_g
    diff_ins["so_ar_g"] = so_ar_g
    diff_ins["dx_ar_g"] = dx_ar_g
    diff_ins["frnw_col"] = frnw_col
    diff_ins["frnw_g"] = dfrnw_g
    diff_ins["qlat_g"] = qlat_g
    diff_ins["ubcd_g"] = ubcd_g
    diff_ins["dbcd_g"] = dbcd_g
    
    diff_ins["qtrib_g"] = qtrib_g
    
    diff_ins["paradim"] = paradim 
    diff_ins["para_ar_g"] = para_ar_g   
    diff_ins["iniq"] = iniq
    # python-fortran crosswalk data
    diff_ins["pynw"] = pynw
    diff_ins["ordered_reaches"] = ordered_reaches
    # bathymetry data
    diff_ins["x_bathy_g"]= x_bathy_g
    diff_ins["z_bathy_g"]= z_bathy_g
    diff_ins["mann_bathy_g"]= mann_bathy_g
    diff_ins["size_bathy_g"]= size_bathy_g
    diff_ins["mxnbathy_g"]= mxnbathy_g

    
    
    '''    
    diff_ins["dtini_g"] = dtini_g
    diff_ins["t0_g"] = t0_g
    diff_ins["tfin_g"] = tfin_g
    diff_ins["saveinterval_tu"] = saveinterval_tu
    diff_ins["saveinterval_cnt"] = saveinterval_cnt
    diff_ins["dt_ql_g"] = dt_ql_g
    diff_ins["dt_ub_g"] = dt_ub_g
    diff_ins["dt_db_g"] = dt_db_g
    diff_ins["nts_ql_g"] = nts_ql_g
    diff_ins["nts_ub_g"] = nts_ub_g
    diff_ins["nts_db_g"] = nts_db_g
    diff_ins["mxncomp_g"] = mxncomp_g
    diff_ins["nrch_g"] = nrch_g
    diff_ins["z_ar_g"] = z_ar_g
    diff_ins["bo_ar_g"] = bo_ar_g
    diff_ins["traps_ar_g"] = traps_ar_g
    diff_ins["tw_ar_g"] = tw_ar_g
    diff_ins["twcc_ar_g"] = twcc_ar_g
    diff_ins["mann_ar_g"] = mann_ar_g
    diff_ins["manncc_ar_g"] = manncc_ar_g
    diff_ins["so_ar_g"] = so_ar_g
    diff_ins["dx_ar_g"] = dx_ar_g
    #diff_ins["nhincr_m_g"] = nhincr_m_g
    #diff_ins["nhincr_f_g"] = nhincr_f_g
    #diff_ins["ufhlt_m_g"] = ufhlt_m_g
    #diff_ins["ufqlt_m_g"] = ufqlt_m_g
    #diff_ins["ufhlt_f_g"] = ufhlt_f_g
    #diff_ins["ufqlt_f_g"] = ufqlt_f_g
    diff_ins["frnw_col"] = frnw_col
    diff_ins["frnw_g"] = dfrnw_g
    diff_ins["qlat_g"] = qlat_g
    diff_ins["ubcd_g"] = ubcd_g
    diff_ins["dbcd_g"] = dbcd_g
    diff_ins["cfl_g"] = cfl_g
    diff_ins["theta_g"] = theta_g
    diff_ins["tzeq_flag_g"] = tzeq_flag_g
    diff_ins["y_opt_g"] = y_opt_g
    diff_ins["so_llm_g"] = so_llm_g
    diff_ins["ntss_ev_g"] = ntss_ev_g
    diff_ins["iniq"] = iniq

    # python-fortran crosswalk data
    diff_ins["pynw"] = pynw
    diff_ins["ordered_reaches"] = ordered_reaches    
    '''
    #import pdb; pdb.set_trace()
    
    return diff_ins


def unpack_output(pynw, ordered_reaches, out_q, out_elv):
    """
    Unpack diffusive wave output arrays
    
    Parameters
    ----------
    pynw -- (dict) ordered reach head segments
    ordered_reaches -- 
    out_q -- (diffusive._memoryviewslice of float64) diffusive wave model flow output (m3/sec)
    out_elv -- (diffusive._memoryviewslice of float64) diffusive wave model water surface elevation output (meters)
    
    Returns
    -------
    np.asarray(rch_list, dtype=np.intp) - segment indices 
    np.asarray(dat_all, dtype = 'float32') - flow, velocity, elevation array
    """

    reach_heads = list(pynw.values())
    nts = len(out_q[:, 0, 0])

    i = 1
    rch_list = []
    for o in ordered_reaches.keys():
        for rch in ordered_reaches[o]:

            rch_segs = rch[1]["segments_list"]
            rch_list.extend(rch_segs[:-1])

            j = reach_heads.index(rch[0])

            if i == 1:
                dat_all = np.empty((len(rch_segs)-1, nts * 3))
                dat_all[:] = np.nan
                # flow result
                dat_all[:, ::3] = np.transpose(np.array(out_q[:, 1 : len(rch_segs), j]))
                # elevation result
                dat_all[:, 2::3] = np.transpose(
                    np.array(out_elv[:, 1 : len(rch_segs), j])
                )

            else:
                dat_all_c = np.empty((len(rch_segs)-1, nts * 3))
                dat_all_c[:] = np.nan
                # flow result
                dat_all_c[:, ::3] = np.transpose(
                    np.array(out_q[:, 1 : len(rch_segs), j])
                )
                # elevation result
                dat_all_c[:, 2::3] = np.transpose(
                    np.array(out_elv[:, 1 : len(rch_segs), j])
                )
                # concatenate
                dat_all = np.concatenate((dat_all, dat_all_c))

            i += 1

    return np.asarray(rch_list, dtype=np.intp), np.asarray(dat_all, dtype="float32")