SELECT 	element.elem_label as Label,
	node.localidx as idx,
	element.elem_type_geo as geo,
	element.elem_grp_no as grp_no,
	node.node
FROM 	Element as element,
	Element_nodes as node
WHERE
	element.result_idx   = RESULTIDX
	AND node.element_idx = element.idx
