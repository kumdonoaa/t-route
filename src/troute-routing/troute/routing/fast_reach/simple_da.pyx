from libc.math cimport exp, isnan, NAN, fabs, pow
from libc.stdio cimport printf
#from libcpp.vector cimport vector
#from libcpp.map cimport map
#from libcpp.set cimport set
#from libcpp.utility cimport pair

cpdef float simple_da_with_decay_py(
    const float last_valid_obs,
    const float model_val,
    const float minutes_since_last_valid,
    const float decay_coeff,
):
    """
    pass-through for using pytest with `simple_da_with_decay`
    """
    return simple_da_with_decay(
        last_valid_obs,
        model_val,
        minutes_since_last_valid,
        decay_coeff,
    )


cdef (float, float, float, float) simple_da(
    const float timestep,
    const float routing_period,
    const float decay_coeff,
    const float gage_maxtimestep,
    const float target_val,
    const float model_val,
    float lastobs_time,
    float lastobs_val,
    bint da_check_gage = 0,
) nogil:
    """
    wrapper function to compute all DA elements
    """
    cdef float replacement_val, nudge_val, da_weighted_shift, da_decay_minutes,
    # cdef float lastobs_timestep, lastobs_value,

    # TODO: It is possible to remove the following branching logic if
    # we just loop over the timesteps during DA and post-DA, if that
    # is a major performance optimization. On the flip side, it would
    # probably introduce unwanted code complexity.
    if isnan(target_val):
        if da_check_gage:
            printf("THIS IS A NAN\t")
    # If we are still within the DA timeseries and the value is not a NaN,
    # then build the DA from the incoming value and update the lastobs arrays.
    if ((timestep <= gage_maxtimestep) and not (isnan(target_val))):
        if da_check_gage:
            printf("replace\t")
        replacement_val = target_val
        nudge_val = target_val - model_val
        # add/update lastobs_timestep
        lastobs_time = (timestep) * routing_period
        lastobs_val = target_val
    # In the unusual case that the observation is missing and the lastobs is also
    # missing, pass along the modeled value and flag the lastobs as still NaN.
    elif ((isnan(target_val)) and (isnan(lastobs_val))):
        replacement_val = model_val
        nudge_val = 0.0
        lastobs_val = NAN
        lastobs_time = NAN
    # When we are outside of the DA timeseries and/or the gage value is NaN
    # AND the lastobs is not null, use the decay calculation to estimate the
    # replacement value and leave the lastobs unmodified.
    else:
        if da_check_gage:
            printf("So we are fixing that...\t")
        if da_check_gage:
            printf("decay; target: %g\t", target_val)
        da_decay_minutes = ((timestep) * routing_period - lastobs_time) / 60 # seconds to minutes
        da_weighted_shift = obs_persist_shift(lastobs_val, model_val, da_decay_minutes, decay_coeff)
        nudge_val = da_weighted_shift
        # TODO: we need to export these values
        # replacement_val = simple_da_with_decay(lastobs_val, model_val, da_decay_minutes, decay_coeff)
        replacement_val = model_val + da_weighted_shift

        if da_check_gage:
            printf("a: %g\t", decay_coeff)
            printf("ts: %g\t", timestep)
            printf("dt: %g\t", routing_period)
            printf("min: %g\t", da_decay_minutes)
            printf("lo_t: %f\t", lastobs_time)
            printf("lov: %g\t", lastobs_val)
            printf("ndg: %g\t", da_weighted_shift)
            printf("orig: %g\t", model_val)
            printf("new: %g\t", replacement_val)
            printf("\n")

    return replacement_val, nudge_val, lastobs_time, lastobs_val,


cdef float simple_da_with_decay(
    const float last_valid_obs,
    const float model_val,
    const float minutes_since_last_valid,
    const float decay_coeff,
) nogil:
    """
    pass-through for computing final value instead of just 'nudge'.
    """
    return model_val + obs_persist_shift(
        last_valid_obs,
        model_val,
        minutes_since_last_valid,
        decay_coeff,
    )


cdef float obs_persist_shift(
    const float last_valid_obs,
    const float model_val,
    const float minutes_since_last_valid,
    const float decay_coeff,
) nogil:
    """
    Given a modeled value, last valid observation,
    time since that observation, and an exponential
    decay_coefficient, compute the 'nudge' value
    """

    cdef float da_weight, da_shift, da_weighted_shift
    da_weight = exp(fabs(minutes_since_last_valid)/-decay_coeff)  # TODO: This could be pre-calculated knowing when obs finish relative to simulation time
    # TODO: we need to be able to export these values to compute the 'Nudge'
    # One possibility would be to return only the nudge from this function...
    da_shift = last_valid_obs - model_val
    da_weighted_shift = da_shift * da_weight
    # printf("t: %g\ta: %g\t%g %g\tlo: %g\torig: %g --> new:%g\n", minutes_since_last_valid, decay_coeff, da_shift, da_weighted_shift, last_valid_obs, model_val, model_val + da_weighted_shift)
    return da_weighted_shift


cdef track_upstream_from_segment_idx(int segment_idx,
                                    dict upstream_connections,
                                    dict segment_id_to_idx):
    """
    Given a segment index (e.g., 6), find the corresponding segment ID (e.g., 115657),
    then track upstream segment IDs via upstream_connections,
    and return:
        - unordered list of upstream segment IDs (e.g., 115651, 115652, 115653, ...)
        - corresponding segment idx (e.g, 0, 2, 1, ...)

    Parameters
    ----------
    segment_idx : int
        The model index (e.g., 6)
    upstream_connections : dict[int, list[int]]
        Maps segment ID to upstream segment IDs
    segment_id_to_idx : dict[int, int]
        Maps segment ID to segment index

    Returns
    -------
    tuple of (list[int], list[int])
        Upstream segment IDs, and corresponding segment indices
    """
    cdef dict idx_to_segment_id = {}
    cdef int segment_id, idx

    # Build reverse mapping: segment_idx → segment_id
    for segment_id in segment_id_to_idx:
        idx = segment_id_to_idx[segment_id]
        idx_to_segment_id[idx] = segment_id

    if segment_idx not in idx_to_segment_id:
        raise ValueError(f"Segment index {segment_idx} not found in mapping.")

    cdef int start_segment_id = idx_to_segment_id[segment_idx]

    # DFS traversal
    cdef set visited = set()
    cdef list stack = [start_segment_id]
    cdef int current
    cdef list upstreams
    cdef int upstream

    cdef list upstream_segment_ids = []
    cdef list upstream_segment_idxs = []

    while stack:
        current = stack.pop()
        upstreams = upstream_connections.get(current, [])
        for upstream in upstreams:
            if upstream not in visited:
                visited.add(upstream)
                stack.append(upstream)

    for segment_id in visited:
        upstream_segment_ids.append(segment_id)
        if segment_id in segment_id_to_idx:
            upstream_segment_idxs.append(segment_id_to_idx[segment_id])

    return upstream_segment_ids, upstream_segment_idxs


cdef (float, float, int) simple_scaling(
    const float timestep,
    const float routing_period,
    const float gage_maxtimestep,
    const float target_val,
    const float model_val,
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
    int[:] simple_scaled_segment_idx,  
    int insert_pos
):
    """
    wrapper function to compute all DA elements
    """
    cdef float replacement_val, nudge_val
    cdef list upstream_segment_ids, upstream_segment_idxs
    
    cdef int segment_id, idx, i, j
    cdef float dQ0, dQ1, totalarea0, totalarea1
    cdef int ndim = q_values_timestep.shape[0]
    cdef bint already_scaled

    nudge_val = target_val - model_val
    lastobs_time = (timestep) * routing_period
    lastobs_val = target_val


    # Create a Python list initially with the same values as q_values_timestep.
    # The initial values will be updated by simple scaling for selected stream segments
    for i in range(ndim):
        simple_scaled_Q[i] = q_values_timestep[i]

    # Identify upstream segment IDs and their correspoinding indexes for a given segment
    # identified by its segment index. 
    upstream_segment_ids, upstream_segment_idxs = track_upstream_from_segment_idx(
        segment_idx,
        upstream_connections,
        segment_id_to_idx
    )
    print(f"given segment index: {segment_idx}")
    print("upstream_segment_ids:", " ".join([str(x) for x in upstream_segment_ids]))
    print("upstream_segment_idxs:", " ".join([str(x) for x in upstream_segment_idxs]))

    simple_scaled_Q[segment_idx] = target_val
    dQ0 = target_val - model_val
    totalarea0 = totaldasqkm[segment_idx]
    print(f"dQ0: {dQ0} totalarea0: {totalarea0}")

    # Insert initial segment_idx
    simple_scaled_segment_idx[insert_pos] = segment_idx
    insert_pos += 1

    for i in range(len(upstream_segment_idxs)):
        idx = upstream_segment_idxs[i]
        # Manual check for already scaled segment indices
        already_scaled = False
        for j in range(insert_pos):
            if simple_scaled_segment_idx[j] == idx:
                already_scaled = True
                break
        
        if not already_scaled:
            print(f"segment index in scaling: {idx}")
            totalarea1 = totaldasqkm[idx]
            simple_scaled_Q[idx] =  q_values_timestep[idx] + dQ0*pow(totalarea1/totalarea0, 0.78)
            simple_scaled_segment_idx[insert_pos] = idx
            insert_pos += 1
    
    print(f"insert_pos: {insert_pos}")
    print("simple_scaled_segment_idx:", 
      " ".join([str(simple_scaled_segment_idx[i]) for i in range(insert_pos)]))
    print("simple_scaled_Q:", " ".join([f"{simple_scaled_Q[i]:.3f}" for i in range(simple_scaled_Q.shape[0])]))

    return lastobs_time, lastobs_val, insert_pos


