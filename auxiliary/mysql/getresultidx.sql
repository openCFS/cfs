SELECT	DATE_FORMAT(c.stamp,"%d.%m.%y   %H:%i:%S") as TimeStamp,
	i.filename,
	r.idx 
FROM	Calculation as c,
	Result as r,
	InputParam as i
WHERE 	r.calculation_idx	= c.idx
  AND	i.idx 			= c.inputparam_idx
