SELECT 	ds780.elem_label as elem_label,
	element.localidx as idx,
	ds780.elem_type_geo as geo,
	ds780.elem_grp_no as grp_no,
	element.node
FROM Dataset780 as ds780,
     Dataset780_node as element
WHERE
	ds780.result_idx = RESULTIDX
	AND element.dataset780_idx = ds780.idx
