#ifndef _BEM_H
#define _BEM_H
/*
 * Set *ONE* of these to 1 according to the format the FORTRAN compiler
 * uses for external names(e.g. gcc -DFORTRANDOUBLEUNDERSCORE)
 *
 * FORTRANNOUNDERSCORE		- no underscore appended
 * FORTRANCAPS			- all capitals
 * FORTRANUNDERSCORE		- one underscore appended
 * FORTRANDOUBLEUNDERSCORE	- two underscores appended if name contains
 *				  an underscore else one appended (g77 default)
*/

/// F77 Interface to matrix generation
/// ifdef  F77BEM: F77 libs of ostbem
/// ifndef F77BEM: obem.cc (no genration)
/// the library ost* is limited in size (2048!!!)
#ifdef F77BEM

#define  FORTRANUNDERSCORE
#ifdef FORTRANDOUBLEUNDERSCORE
#define FORTRANUNDERSCORE
#endif

#ifdef FORTRANUNDERSCORE


#ifdef FORTRANDOUBLEUNDERSCORE
#define disk_element		disk_element__
#define emt_parameter           emt_parameter__
#define disk_masse              disk_masse__
#define disk_single_laplace	disk_single_laplace__
#define disk_double_laplace	disk_double_laplace__
#define make_hyper_laplace	make_hyper_laplace__
#define darst3d                 darst3d__
#define disk_element_2d		disk_element_2d__
#define disk_masse_2d           disk_masse_2d__
#define disk_single_laplace_2d	disk_single_laplace_2d__
#define disk_double_laplace_2d	disk_double_laplace_2d__
#define make_hyper_laplace_2d	make_hyper_laplace_2d__
#define darst2d                 darst2d__
#else
#define disk_element		disk_element_
#define emt_parameter           emt_parameter_
#define disk_masse              disk_masse_
#define disk_single_laplace	disk_single_laplace_
#define disk_double_laplace	disk_double_laplace_
#define make_hyper_laplace	make_hyper_laplace_
#define darst3d                 darst3d_
#define disk_element_2d		disk_element_2d_
#define disk_masse_2d           disk_masse_2d_
#define disk_single_laplace_2d	disk_single_laplace_2d_
#define disk_double_laplace_2d	disk_double_laplace_2d_
#define make_hyper_laplace_2d	make_hyper_laplace_2d_
#define darst2d                 darst2d_
#endif						/* FORTRANDOUBLEUNDERSCORE */
#endif 						/* FORTRANUNDERSCORE */
#endif 						/* F77BEM */

extern "C"
{
  // N   number of elements 
  // IE  element to node list 
  // X   list of coordinates
  // KP  disk-element sizes
  extern int disk_element(int *, int *, double *, double *);

  // N   number of elements 
  // IE  element to node list 
  // X   list of coordinates
  // HE  element parameter
  extern int emt_parameter(int *, int *, double *, double *);


  // N   number of elements 
  // IE  element to node list 
  // IES element to node list 
  // X   list of coordinates
  // KP  disk-element sizes
  // XS  list of coordinates (dirichlet)
  // KM  mnass matrix
  // NP  number of nodes
  extern int disk_masse(int *, int *,int *,double *, double *,double *, double *, int *);

  // N   number of elements  (dirichlet)
  //-M   number of nodes (dirichlet) 
  // IE  element to node list (dirichlet)
  // X   list of point coordinates (dirichlet)
  //-D   2: Galerkin
  //+D   4: 7 point int. (or 1)
  // KE  single layer matrix
  //+HE  element paramter
  //-NU  0: pw costant Neumann data
  // extern int disk_single_laplace(int *, int *, int *, double *, int *, double *, int *);
  extern int disk_single_laplace(int *, int *, double *, int *, double *, double *);


  // NBN  number of elements (neumann)
  // MBN  number of nodes (neumann)
  // IEN  element to node list (neumann)
  // XN   list of point coordinates (neumann)
  //-NBD  number of elements  (dirichlet)
  //-MBD  number of nodes (dirichlet)
  //-IED  element to node list (dirichlet)
  //-XD   list of point coordinates (dirichlet)
  //+D    4: 7-point int. (or 1)
  // KD   double layer matrix
  //-D    2: Galerkin
  //-NU   0: pw costant Neumann data
  //+HE   element parameter
	 //extern int disk_double_laplace(int *, int *, int *, double *,
	 //				 int *, int *, int *, double *,
	 //				 double *, int *, int *);
  extern int disk_double_laplace(int *, int *, int *, double *, int *,
				 double *, double *);
  // N   number of elements 
  // M   number of nodes 
  // IE  element to node list 
  // X   list of coordinates
  // KE  single layer matrix
  // KH  hypersingular matrix
  extern int make_hyper_laplace(int *, int *, int *, double *, double *, double *);

  // N   number of elements
  // M   number of nodes
  // IE  element to node list
  // X   list of point coordinates
  // P   1: Laplace
  //-NU  0: pw costant Neumann data
  //-MU  1: pw linear dirichlet data
  // DD  vector of dirichlet data
  // ND  vector of neumann data
  // MAT material coeff (not used for laplace)
  // XC  point coordinates (x,y)
  // U   result
  //extern int darst3d(int *, int *, int *, double *, int *, int *, int *,
  //		     double *, double *, double *, double *, double *);
  extern int darst3d(int *, int *, int *, double *, int *, 
		     double *, double *, double *, double *, double *);
  // N   number of elements 
  // IE  element to node list 
  // X   list of coordinates
  // KP  disk-element sizes
  extern int disk_element_2d(int *, int *, double *, double *);


  // N   number of elements 
  // NU  0: pw costant Neumann data
  // D   2: Galerkin
  // IEM element to node list  (neumann)
  // XM  list of coordinates (neumann)
  // KP  disk-element sizes
  // IES element to node list 
  // X   list of coordinates
  // KM  mnass matrix
  // NP  number of nodes
  extern int disk_masse_2d(int *, int *,int *, int *,
			   double *, double *,int *, double *, double *, int *);  


  // N   number of elements  (dirichlet)
  // M   number of nodes (dirichlet)
  // IE  element to node list (dirichlet)
  // X   list of point coordinates (dirichlet)
  // D   2: Galerkin
  // NU  0: pw costant Neumann data
  // KE  single layer matrix
  // MAT material coeff (USED for laplace)
  extern int disk_single_laplace_2d(int *, int *, int *, double *, 
				    int *, int *,
				    double *, double *);

  // N   number of elements (neumann)
  // M   number of nodes (neumann)
  // IE  element to node list (neumann)
  // X   list of point coordinates (neumann)
  // NM  number of elements  (dirichlet)
  // MM  number of nodes (dirichlet)
  // IEM element to node list (dirichlet)
  // XM  list of point coordinates (dirichlet)
  // D   2: Galerkin
  // NU  0: pw costant Neumann data
  // MU  1: pw linear dirichlet data
  // KD  double layer matrix
  // MAT material coeff (USED for laplace)
  extern int disk_double_laplace_2d(int *, int *, int *, double *,
				    int *, int *, int *, double *, 
				    int *, int *, int *,
				    double *, double *);

  // N   number of elements 
  // M   number of nodes 
  // P   1: Laplace, 2: elasticity
  // NU  0: pw costant Neumann data
  // IE  element to node list 
  // KE  single layer matrix
  // KH  hypersingular matrix
  // KP  disk-element sizes
  // MAT material coeff (not used for laplace)
  extern int make_hyper_laplace_2d(int *, int *, int *, int *, int *,
				   double *, double *, double *, double *);


  // N   number of elements (neumann)
  // M   number of nodes (neumann)
  // IE  element to node list (neumann)
  // X   list of point coordinates (neumann)
  // NM   number of elements  (dirichlet)
  // MM   number of nodes (dirichlet)
  // IEM  element to node list (dirichlet)
  // XM   list of point coordinates (dirichlet)
  // P   1: Laplace, 2: elasticity
  // NU  0: pw costant Neumann data
  // MU  1: pw linear dirichlet data
  // DD  vector of dirichlet data
  // ND  vector of neumann data
  // MAT material coeff (not used for laplace)
  // XC  point coordinates (x,y)
  // U   result
  extern int darst2d(int *, int *, int *, double *, 
		     int *, int *, int *, double *, 
		     int *, int *, int *,
		     double *, double *, double *, double *, double *);
}  


#endif 						/* _BEM_H */


