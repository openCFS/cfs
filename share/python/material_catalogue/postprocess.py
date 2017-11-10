#!/usr/bin/env python


# modifies Ngocs material catalog
# adds solid material and void material to catalog
mat = []
#vol = []
for line in open('detailed_stats_ortho_fcrc_bracket').readlines():
    mat.append(line.split())
#for line in open('detailed_stats_vol_ortho_fcrc_bracket').readlines():
#    vol.append(line.split())
mat[1][3:len(mat[0])] = ['1.346154e-06', '5.769231e-07', '5.769231e-07', '1.346154e-06', '5.769231e-07', '1.346154e-06', '3.846154e-07', '3.846154e-07', '3.846154e-07']
#vol[1][3]= '0.'
for i in range(len(mat)):
    if mat[i][0] == '10' or mat[i][1] == '10' or mat[i][2] == '10':
       #mat[i][3:len(mat[0])] = ['1.430976e+00', '6.734007e-01+00', '6.734007e-01', '1.430976e+00', '6.734007e-01','1.430976e+00', '3.787879e-01', '3.787879e-01', '3.787879e-01'] #Lufo Ti
       mat[i][3:len(mat[0])] = ['1.346154e+00', '5.769231e-01','5.769231e-01', '1.346154e+00', '5.769231e-01', '1.346154e+00', '3.846154e-01', '3.846154e-01', '3.846154e-01'] #99 lines

#for i in range(len(vol)):
#    if vol[i][0] == '10' or vol[i][1] == '10' or vol[i][2] == '10':
#        vol[i][3] = '1.'        
output = open('detailed_stats_ortho_fcrc_bracket_99lines.txt', 'w')
for i in range(len(mat)):
    for j in range(len(mat[0])):
        output.write("%s " % mat[i][j])
    output.write("\n")
output.close()
#output = open('detailed_stats_vol_ortho_fcrc_bracket.txt', 'w')
#for i in range(len(vol)):
#    for j in range(len(vol[0])):
  #      output.write("%s " % vol[i][j])
 #   output.write("\n")
#output.close()
