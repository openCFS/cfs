#!/bin/bash

# ANSREL should be 12.0, 12.1 or 13.0

ANSREL=13.0

# General Documentation
desc[1]="ANSYS, Inc. Installation and Licensing Guides"
pdfs[1]=ai_instl.pdf
desc[2]="ANSYS, Inc. Release Notes"
pdfs[2]=ai_rn.pdf
desc[3]="Legal Notice - ANSYS"
pdfs[3]=ai_ginfo.pdf

# Workbench Documentation
desc[4]="ANSYS Workbench User's Guide"
pdfs[4]=wb2_help.pdf
desc[5]="ANSYS Workbench Scripting Guide"
pdfs[5]=wb2_js.pdf
desc[6]="Design Exploration Application"
pdfs[6]=wb_dx.pdf
desc[7]="DesignModeler Application"
pdfs[7]=wb_dm.pdf
desc[8]="EKM Desktop User's Guide"
pdfs[8]=wb_ekmdt.pdf
desc[9]="Engineering Data Application"
pdfs[9]=wb_eda.pdf
desc[10]="External Connection User's Guide"
pdfs[10]=wb2_exc.pdf
desc[11]="FE Modeler Application"
pdfs[11]=wb_fem.pdf
desc[12]="Remote Solve Manager (RSM)"
pdfs[12]=wb_rsm.pdf

# Structural Mechanics Documentation
desc[13]="Mechanical APDL Verification Manual"
pdfs[13]=ans_vm.pdf
desc[14]="Mechanical (formerly Simulation)"
pdfs[14]=wb_sim.pdf
desc[15]="Command Reference"
pdfs[15]=ans_cmd.pdf
desc[16]="Element Reference"
pdfs[16]=ans_elem.pdf
desc[17]="Mechanical APDL Context Sensitive Help"
pdfs[17]=ans_wid.pdf
desc[18]="Theory Reference for the Mechanical APDL and Mechanical"
pdfs[18]=ans_thry.pdf
desc[19]="Mechanical APDL Programmer's Manual"
pdfs[19]=ans_prog.pdf
desc[20]="Advanced Analysis Techniques Guide"
pdfs[20]=ans_adv.pdf
desc[21]="Connection User's Guide"
pdfs[21]=ans_con.pdf
desc[22]="Basic Analysis Guide"
pdfs[22]=ans_bas.pdf
desc[23]="Coupled-Field Analysis Guide"
pdfs[23]=ans_cou.pdf
desc[24]="Contact Technology Guide"
pdfs[24]=ans_ctec.pdf
desc[25]="Distributed ANSYS Guide"
pdfs[25]=ans_dan.pdf
desc[26]="Fluids Analysis Guide"
pdfs[26]=ans_flu.pdf
desc[27]="Low-Frequency Electromagnetic Analysis Guide"
pdfs[27]=ans_lof.pdf
desc[28]="High-Frequency Electromagnetic Analysis Guide"
pdfs[28]=ans_hif.pdf
desc[29]="ANSYS LS-DYNA User's Guide"
pdfs[29]=ans_lsd.pdf
desc[30]="Modeling and Meshing Guide"
pdfs[30]=ans_mod.pdf
desc[31]="Operations Guide"
pdfs[31]=ans_ope.pdf
desc[32]="Performance Guide"
pdfs[32]=ans_per.pdf
desc[33]="Structural Analysis Guide"
pdfs[33]=ans_str.pdf
desc[34]="Rotordynamic Analysis Guide - ANSYS"
pdfs[34]=ans_rot.pdf
desc[35]="Technology Demonstration Guide"
pdfs[35]=ans_tec.pdf
desc[36]="Thermal Analysis Guide"
pdfs[36]=ans_the.pdf
desc[37]="Troubleshooting Guide"
pdfs[37]=ans_tbl.pdf
desc[38]="Multibody Analysis Guide"
pdfs[38]=ans_mul.pdf

# Meshing Documentation
desc[39]="CFX-Mesh"
pdfs[39]=wb_cm.pdf
desc[40]="Meshing User's Guide"
pdfs[40]=wb_msh.pdf
desc[41]="ANSYS ICEM CFD Help Manual"
pdfs[41]=icemcfd/icmhelp.pdf
desc[42]="ANSYS ICEM CFD User Manual"
pdfs[42]=icemcfd/icmuser.pdf
desc[43]="ANSYS ICEM CFD Tutorial Manual"
pdfs[43]=icemcfd/icmtut.pdf



# my.fit.edu/itresources/manuals/fluent13.0
desc[50]="FLUENT in Workbench Tutorials"
pdfs[50]=flu_wb.pdf
desc[51]="ANSYS FLUENT User's Guide"
pdfs[51]=flu_ug.pdf
desc[52]="ANSYS FLUENT UDF Manual"	
pdfs[52]=flu_udf.pdf
desc[53]="ANSYS FLUENT Theory Guide"
pdfs[53]=flu_th.pdf
desc[54]="ANSYS FLUENT Tutorial Guide"	
pdfs[54]=flu_tg.pdf
desc[55]="ANSYS FLUENT Text Command List"	
pdfs[55]=flu_tcl.pdf
desc[56]="Running ANSYS FLUENT Under SGE"
pdfs[56]=flu_sge.pdf
desc[57]="Running ANSYS FLUENT Under PBS Professional"
pdfs[57]=flu_pbs.pdf
desc[58]="ANSYS FLUENT Population Balance Module Manual"
pdfs[58]=flu_pbm.pdf
desc[59]="ANSYS FLUENT Magnetohydrodynamics (MHD) Module Manual"
pdfs[59]=flu_mhd.pdf
desc[60]="Running ANSYS FLUENT Under LSF"
pdfs[60]=flu_lsf.pdf
desc[61]="ANSYS FLUENT Getting Started Guide"
pdfs[61]=flu_gs.pdf	
desc[62]="ANSYS FLUENT Continuous Fiber Module Manual"
pdfs[62]=flu_fib.pdf
desc[63]="ANSYS FLUENT Fuel Cells Module Manual"
pdfs[63]=flu_fcm.pdf

desc[64]="CFX Tutorials"
pdfs[64]=cfx_tutr.pdf
desc[65]="CFX-Solver Theory Guide"
pdfs[65]=cfx_thry.pdf
desc[66]="CFX-Solver Manager User's Guide"
pdfs[66]=cfx_solv.pdf
desc[67]="CFX Reference Guide"
pdfs[67]=cfx_ref.pdf
desc[68]="CFX-Pre User's Guide"	
pdfs[68]=cfx_pre.pdf


#130
if [ "$ANSREL" = 13.0 ]; then
    desc[69]="CFD-Post User's Guide - ANSYS"
    pdfs[69]="cfd_post.pdf"

    desc[44]="ANSYS TurboGrid User's Guide"
    pdfs[44]=tg_user.pdf
    desc[45]="ANSYS TurboGrid Tutorials"
    pdfs[45]=tg_tutr.pdf
    desc[46]="TurboGrid Reference Guide"
    pdfs[46]=tg_ref.pdf
    desc[47]="ANSYS TurboGrid Introduction"
    pdfs[47]=tg_intr.pdf

#    desc[48]="ANSYS Icepack User's Guide"
#    pdfs[48]=ice_ug.pdf	
#    desc[49]="ANSYS Icepack Tutorials"
#    pdfs[49]=ice_tut.pdf

else
#12x
    desc[69]="ANSYS CFD-Post User's Guide"
    pdfs[69]=cfx/xpost.pdf

    desc[44]="ANSYS TurboGrid User's Guide"
    pdfs[44]=cfx/xtuser.pdf
    desc[45]="ANSYS TurboGrid Tutorials"
    pdfs[45]=cfx/xttutr.pdf
    desc[46]="TurboGrid Reference Guide"
    pdfs[46]=cfx/xtref.pdf
    desc[47]="ANSYS TurboGrid Introduction"
    pdfs[47]=cfx/xtintr.pdf

#    if [ "$ANSREL" = 12.1 ]; then
#	desc[48]="ANSYS Icepack User's Guide"
#	pdfs[48]=icepak/iceug.pdf
#	desc[49]="ANSYS Icepack Tutorials"
#	pdfs[49]=icepak/icetut.pdf
#    fi
fi

desc[70]="CFX-Solver Modeling Guide"
pdfs[70]=cfx_mod.pdf	
desc[71]="CFX Introduction"
pdfs[71]=cfx_intr.pdf
desc[72]="CFD-Post Tutorials"
pdfs[72]=cfd_posttutr.pdf

desc[73]="AUTODYN User Subroutines Tutorial"
pdfs[73]=adyn_usrsub.pdf
desc[74]="AUTODYN Parallel Processing Tutorial"
pdfs[74]=adyn_paratutorial.pdf
desc[75]="AUTODYN Composite Modelling"
pdfs[75]=adyn_composite.pdf

index=1
for i in "${desc[@]}"
do
    pdf=${pdfs[$index]}
    dir=$(dirname $pdf)
    base=$(basename $pdf)
    mkdir -p $dir
    echo "Downloading $i for version $ANSREL..."
    if [ ! -f $pdf ]; then
	wget http://www1.ansys.com/customer/content/documentation/$(echo $ANSREL | sed 's/\.//')/$pdf
	if [ ! "$?" = "0" ]; then
	    exit 1;
	fi
	mv $base $dir
    fi
    ((index++))
done