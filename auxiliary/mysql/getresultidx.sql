SELECT	DATE_FORMAT(c.stamp,"%d.%m.%y   %H:%i:%S"),
	r.idx 
FROM	Calculation as c,
	Result as r
WHERE 	r.calculation_idx	= c.idx
