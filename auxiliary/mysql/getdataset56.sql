SELECT 	param.time_step as step, 
	param.time, 
	prob.time_step as step, 
	prob.time_value as time, 
	prob.result_type as type,
	result.node_no, 
	result.node_idx, 
	result.result
FROM 	Dataset56 as ds56, 
	Dataset56_result as result, 
	Dataset56_trans_param as param, 
	Dataset56_trans_prob as prob

WHERE 	ds56.result_idx	= RESULTIDX

	AND	param.idx 	= ds56.dataset56_trans_param_idx
	AND	prob.idx	= ds56.dataset56_trans_prob_idx
	AND	ds56.idx	= result.dataset56_idx
