#from libcpp.map cimport map
#from libcpp.vector cimport vector

cdef (float, float, float, float) simple_da(
    const float timestep,
    const float routing_period,
    const float decay_coeff,
    const float gage_maxtimestep,
    const float target_val,
    const float model_val,
    float lastobs_time,
    float lastobs_val,
    bint da_check_gage=*,
) nogil


cdef float simple_da_with_decay(
    const float last_valid_obs,
    const float model_val,
    const float minutes_since_last_valid,
    const float decay_coeff,
) nogil


cdef float obs_persist_shift(
    const float last_valid_obs,
    const float model_val,
    const float minutes_since_last_valid,
    const float decay_coeff,
) nogil


cdef (float, float, int) simple_scaling(
    float timestep,
    float routing_period,
    float gage_maxtimestep,
    float target_val,
    float model_val,
    float lastobs_time,
    float lastobs_val,
    int segment_idx,
    dict upstream_connections,
    dict segment_id_to_idx,
    #map[int, vector[int]]& upstream_connections,
    #map[int, int]& segment_id_to_idx,
    const float[:] totaldasqkm,
    float[:] q_values_timestep,
    float[:] simple_scaled_Q,
    int[:] simple_scaled_segment_idx,      # buffer of segment idx
    int insert_pos                         # number of currently filled entries
) 