SELECT 	coord.node_label as node,
	coord.x_coord as x,
	coord.y_coord as y,
	coord.z_coord as z
FROM 	Dataset781 as ds781,
	Dataset781_node as coord
WHERE	ds781.result_idx	= RESULTIDX
	AND coord.dataset781_idx=ds781.idx
