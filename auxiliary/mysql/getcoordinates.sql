SELECT 	coord.node_label as node,
	coord.x_coord as x,
	coord.y_coord as y,
	coord.z_coord as z

FROM 	Node as node,
	Node_coordinates as coord

WHERE	node.result_idx		= RESULTIDX
	AND coord.node_idx	= node.idx
