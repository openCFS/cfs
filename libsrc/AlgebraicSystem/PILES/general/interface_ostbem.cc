#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>
#include <fstream.h>
#include <math.h>

#include <general.hh>
#include <multigrid.hh>

namespace CoupledField
{

// OSTBEM stuff

void MakeElement(int numnode, int connect[], double coord[], double des[])
{
  disk_element_2d(&numnode, &connect[0], &coord[0], &des[0]);
}

void MakeSingleLaplace(Integer numelem, Integer numnode, Integer connect[],Double coord[], 
		       Integer ga, Integer nu, Double mat[], Double epar[])
{
  disk_single_laplace_2d(&numelem, &numnode, &connect[0], &coord[0], &ga, &nu, &mat[0], &epar[0]);
}

void MakeDoubleLaplace(Integer numelem1, Integer numnode1, Integer connect1[], Double coord1[],
		       Integer numelem2, Integer numnode2, Integer connect2[], Double coord2[],
		       Integer ga, Integer nu, Integer mu, Double mat[], Double epar[])
{
  disk_double_laplace_2d(&numelem1, &numnode1, &connect1[0], &coord1[0], 
			 &numelem2, &numnode2, &connect2[0], &coord2[0],
			 &ga, &nu, &mu, &mat[0], &epar[0]);
}

void MakeHyperLaplace(Integer numelem, Integer numnode, Integer type, Integer nu, 
		      Integer connect[], Double mat1[], Double mat2[], Double des[], 
		      Double epar[])
{
  make_hyper_laplace_2d(&numelem, &numnode, &type, &nu, &connect[0], &mat1[0], &mat2[0],
			&des[0], &epar[0]);
}

void MakeMasse(Integer numelem, Integer nu, Integer ga, Integer connect1[], Double coord1[],
	       Double des[], Integer connect2[], Double coord2[], Double mat[], Integer numnode)
{
  disk_masse_2d(&numelem, &nu, &ga, &connect1[0], &coord1[0], &des[0], &connect2[0], &coord2[0], 
		&mat[0], &numnode);
}

void CalcData(Integer numelem1, Integer numnode1, Integer connect1[], Double coord1[], 
	      Integer numelem2, Integer numnode2, Integer connect2[], Double coord2[],
	      Integer type, Integer nu, Integer mu, Double dir[], Double neu[], Double epar[], 
	      Double point[], Double val[])
{
//   darst2d(&numelem1, &numnode1, &connect1[0], &coord1[0], 
// 	  &numelem2, &numnode2, &connect2[0], &coord2[0], 
// 	  &type, &nu, &mu, &dir[0], &neu[0], &epar[0],
// 	  &point[0], &val[0]);
}

}
