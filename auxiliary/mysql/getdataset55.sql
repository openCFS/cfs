SELECT 	param.time_step as step, 
	param.time, 
	prob.time_step as step, 
	prob.time_value as time, 
	prob.result_type as type,
	result.node_no, 
	result.node_idx as dof, 
	result.result
FROM 	Dataset55 as ds55, 
	Dataset55_result as result, 
	Dataset55_trans_param as param, 
	Dataset55_trans_prob as prob

WHERE 	ds55.result_idx	= RESULTIDX

	AND	param.idx 	= ds55.dataset55_trans_param_idx
	AND	prob.idx	= ds55.dataset55_trans_prob_idx
	AND	ds55.idx	= result.dataset55_idx
