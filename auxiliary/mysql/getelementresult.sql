SELECT 	trans.time_step as TimeStep, 
	trans.time_value as TimeValue, 
	trans.result_type as Type,
	value.elem_no as NodeNo, 
	value.dof as DOF, 
	value.result as Value
FROM 	Element_result as er, 
	Element_result_value as value, 
	Element_result_trans as trans

WHERE 	er.result_idx	= RESULTIDX
	AND	value.element_result_idx    = er.idx
	AND 	er.element_result_trans_idx = trans.idx
