SELECT  trans.time_step as Timestep, 
        trans.time_value as TimeValue, 
        trans.result_type as Type,
        result.node_no as NodeNo, 
        result.dof as DOF, 
        result.result as Value
FROM    Nodal_result as ns, 
	Nodal_result_trans as trans,
        Nodal_result_value as result 

WHERE   	ns.result_idx		 = RESULTIDX
        AND     ns.idx        		 = result.nodal_result_idx
	AND	ns.nodal_result_trans_idx= trans.idx

