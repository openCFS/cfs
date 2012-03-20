#ifndef SCPIP30_H_
#define SCPIP30_H_

/** Declare a C header for the SCPIP fortran implementation by Ch. Zillober */ 

extern "C" 
{
  /** SCPIP is a FORTRAN77 subroutine. It contains all necessary subroutines besides some linear algebra which
      will be mentioned below. To execute SCPIP, the corresponding file has to be compiled and linked with the
      object codes for function and gradient evaluations provided by the user. All calculations within the subrou-
      tines of this file are performed in double precision arithmetic.
   * @param n Number of variables, at least ¸ 1; not altered (NA) (IN) 
   * @param mie Number of inequality constraints (NA) (IN)
   * @param meq Number of equality constraints (NA) (IN)
   * @param iemax Dimension of inequality dependent arrays H ORG, Y IE, ACTIVE. Must be at least MIE and ¸ 1 (NA) (IN) 
   * @param eqmax Dimension of equality dependent arrays G ORG, Y EQ. Must be at least MEQ and ¸ 1 (NA) (IN)
   * @param x0 (N) X0 has to contain an initial guess of the solution. (IN)
   *        On return: X0 is replaced by the last computed iterate *
   * @param x_l (N) Lower bounds on the variables (NA) (IN)
   * @param x_u (N) Upper bounds on the variables (NA) (IN)
   * @param f_org Contains the objective function value of the last computed iterate
   * @param h_org (IEMAX) contains the values of the inequality constraints at the last computed iterate
   * @param g_org (EQMAX) contains the values of the equality constraints at the last computed iterate
   * @param df (N)contains the gradient of the objective function at X0
   * @param y_ie (IEMAX) contains the Lagrange multipliers for the inequality constraints at the last computed iterate (OUT)
   * @param y_eq (EQMAX) contains the Lagrange multipliers for the equality constraints at the last computed iterate (OUT)
   * @param y_l (N) contains the Lagrange multipliers for the lower bounds on the variables at the last computed iterate (OUT)
   * @param y_u (N) contains the Lagrange multipliers for the upper bounds on the variables at the last computed iterate (OUT)
   * @param icntl (13) Integer array of dimension 13 to be set by the user, but all components have reasonable defaults (NA). 
   *        A value 0 indicates that the defaults should be chosen (IN)
   * @param rcntl (6) double precision array of dimension 6 to be set by the user, 
   *        but all components have reasonable defaults (NA). 
   *        A value 0 indicates that the defaults should be chosen (IN)
   * @param info (23) Integer array of dimension 7[sic!] containing some problem information (OUT)
   * @param rinfo (5) Double precision array of dimension 5 containing some problem information (OUT)  
   * 
   * */
# if defined(_WIN32)
  void SCPIP40I(const int* n,            // ON INPUT (INTEGER = I): NUMBER OF VARIABLES, AT LEAST >=1. NOT ALTERED (NA)
                const int* mie,          // ON INPUT (I): NUMBER OF INEQUALITY CONSTRAINTS (NA)
                const int* meq,          // ON INPUT (I): NUMBER OF EQUALITY CONSTRAINTS (NA)
                const int* iemax,        // ON INPUT (I): DIMENSION OF INEQUALITY DEPENDENT ARRAYS H_ORG, Y_IE,ACTIVE. MUST BE AT LEAST MIE AND >=1. (NA)
                const int* eqmax,        // ON INPUT (I): DIMENSION OF EQUALITY DEPENDENT ARRAYS G_ORG, Y_EQ. MUST BE AT LEAST MEQ AND >=1. (NA) 
                double* x0,                   // ON INPUT (DOUBLE PRECISION = D): X0 HAS TO CONTAIN AN INITIAL GUESS OF THE SOLUTION 
                                               // ON RETURN: X0 IS REPLACED BY THE LAST COMPUTED ITERATE
                const double* x_l,      // ON INPUT (D): LOWER BOUNDS ON THE VARIABLES (NA)  
                const double* x_u,      // ON INPUT (D): UPPER BOUNDS ON THE VARIABLES (NA)
                double* f_org,                // ON RETURN (D), F_ORG CONTAINS THE OBJECTIVE FUNCTION VALUE OF THE LAST COMPUTED ITERATE
                double* h_org,                // ON RETURN (D), H_ORG CONTAINS THE VALUES OF THE INEQUALITY CONSTRAINTS AT THE LAST COMPUTED ITERATE
                double* g_org,                // ON RETURN (D), G_ORG CONTAINS THE VALUES OF THE EQUALITY CONSTRAINTS AT THE LAST COMPUTED ITERATE
                double* df,                   // ON RETURN (D), DF CONTAINS THE GRADIENT OF THE OBJECTIVE FUNCTION AT X0
                double* y_ie,                 // ON RETURN (D), Y_IE CONTAINS THE LAGRANGE MULTIPLIERS FOR THE INEQUALITY CONSTRAINTS AT THE LAST COMPUTED ITERATE  
                double* y_eq,                 // ON RETURN (D), Y_EQ CONTAINS THE LAGRANGE MULTIPLIERS FOR THE EQUALITY CONSTRAINTS AT THE LAST COMPUTED ITERATE
                double* y_l,                  // ON RETURN (D), Y_L CONTAINS THE LAGRANGE MULTIPLIERS FOR THE LOWER BOUNDS ON THE VARIABLES AT THE LAST COMPUTED ITERATE
                double* y_u,                  // ON RETURN (D), Y_U CONTAINS THE LAGRANGE MULTIPLIERS FOR THE UPPER BOUNDS ON THE VARIABLES AT THE LAST COMPUTED ITERATE
                const int* icntl,       // ON INPUT (I): INTEGER ARRAY OF DIMENSION .. TO BE SET BY THE USER, BUT ALL COMPONENTS HAVE REASONABLE DEFAULTS (NA).
                const double* rcntl,    // ON INPUT (D): DOUBLE PRECISION ARRAY OF DIMENSION 6 TO BE SET BY THE USER, BUT ALL COMPONENTS HAVE REASONABLE DEFAULTS (NA). 
                int* info,                    // ON RETURN (I): INTEGER ARRAY CONTAINING SOME PROBLEM INFORMATION.
                double* rinfo,                // ON RETURN (D): DOUBLE PRECISION ARRAY OF DIMENSION 5. CONTAINS INFORMATION ON SOME DATA
                const int* nout,         // ON INPUT (I), INTEGER, INDICATING OUTPUT UNIT NUMBER (NA)
                double* r_scp,                // (D) REAL WORKING ARRAY OF DIMENSION AT LEAST RDIM
                const int* rdim,         // ON INPUT (I): MUST BE AT LEAST 30*N+11*IEMAX+8+10*EQMAX (NA)
                double* r_sub,                // (D) REAL WORKING ARRAY OF DIMENSION AT LEAST RSUBDIM
                const int* rsubdim,      // ON INPUT (I): MUST BE AT LEAST 22*N+41*IEMAX+27*EQMAX+2*IELPAR+EQLPAR (NA)
                int* i_scp,                   // (I) INTEGER WORKING ARRAY OF DIMENSION AT LEAST IDIM
                const int* idim,         // ON INPUT (I): MUST BE AT LEAST 5*N+5*IEMAX+2*EQMAX+3 (NA)
                int* i_sub,                   // (I) INTEGER WORKING ARRAY OF DIMENSION AT LEAST ISUBDIM
                const int* isubdim,      // ON INPUT (I): MUST BE AT LEAST 2*N+3*IEMAX+2*EQMAX+IELPAR (NA)
                int* active,                  // ON RETURN (I): INTEGER ARRAY INDICATING THE INEQUALITY CONSTRAINTS ...
                                               // ON INPUT (I): A USER IS ALLOWED TO CHANGE ACTIVE BEFORE COMPUTING ...
                const int* mode,         // ON INPUT (I): (NA)
                                               // 1: USUAL SUBROUTINE CALL FOR FUNCTION VALUES AND GRADIENTS
                                               // 2: REVERSE COMMUNICATION
                int* ierr,                    // ON (VERY FIRST) INPUT (I): HAS INITIALLY TO BE 0. DO NOT ALTER IERR OUTSIDE THE ROUTINE!
                                               // ON RETURN:
                                               // <0 : REVERSE COMMUNICATION:
                                               // -1: FUNCTION VALUES ARE REQUESTED
                                               // -2: GRADIENTS ARE REQUESTED
                                               // -3: INFO(6) AND INFO(7) ARE COMPUTED. THE USER CAN NOW DIMENSION SPIW AND SPDW DYNAMICALLY. 
                                               //     IF THIS IS NOT DESIRED, JUST REJUMP INTO SCPIP AFTER IERR = -3.
                                               // 0: SUCCESSFUL COMPUTATION ...
                int* iern,                    // ON INPUT (I): IERN CONTAINS THE ROW NUMBERS, I.E. THE INDEX OF THE COMPONENT OF THE ENTRIES 
                                               // OF THE JACOBIAN OF THE INEQUALITY CONSTRAINTS (NA)
                int* iecn,                    // ON INPUT (I): IECN CONTAINS THE COLUMN NUMBERS, I.E. THE NUMBER OF THE
                                               // INEQUALITY OF THE ENTRIES OF THE JACOBIAN OF THE INEQUALITY CONSTRAINTS (NA)
                double* iederv,                  // ON INPUT (D): IEDERV CONTAINS THE VALUES OF THE ENTRIES OF THE JACOBIAN OF THE INEQUALITY CONSTRAINTS. ...
                const int* ielpar,       // ON INPUT (I): DIMENSION OF ARRAYS IERN, IECN AND IEDERV. MUST BE AT LEAST IELENG AND >=1. (NA)
                const int* ieleng,       // ON INPUT (I): ACTUAL NUMBER OF ENTRIES IN IEDERV (NA)
                int* eqrn,                    // ON INPUT (I): EQRN CONTAINS THE ROW NUMBERS, I.E. THE INDEX OF THE COMPONENT OF THE ENTRIES OF 
                                               // THE JACOBIAN OF THE EQUALITYCONSTRAINTS (NA)
                int* eqcn,                    // ON INPUT (I): EQCN CONTAINS THE COLUMN NUMBERS, I.E. THE NUMBER OF THE 
                                               // EQUALITY OF THE ENTRIES OF THE JACOBIAN OF THE EQUALITY CONSTRAINTS (NA)
                double* eqcoef,               // ON INPUT (D): EQCOEF CONTAINS THE VALUES OF THE ENTRIES OF THE JACOBIAN OF THE EQUALITY CONSTRAINTS.
                const int* eqlpar,       // ON INPUT (I): DIMENSION OF ARRAYS EQRN, EQCN AND EQCOEF. MUST BE AT LEAST EQLENG AND >=1. (NA)
                const int* eqleng,       // ON INPUT (I): ACTUAL NUMBER OF ENTRIES IN EQCOEF (NA)
                int* mactiv,                   // ON RETURN (I): NUMBER OF CONSTRAINTS (EQUALITIES AND INEQUALITIES)... 
                                               // THIS AFFECTS THE REQUIRED WORKING SPACE, SEE SPDWDIM BELOW, AND THUS THE SUBPROBLEM SOLUTION METHOD
                int* spiw,                    // (I) INTEGER WORKING ARRAY OF DIMENSION AT LEAST SPIWDIM. 
                                               // AFTER IERR = -3 ONE CAN DIMENSION SPIW AND SPDW DYNAMICALLY WITH THE VALUES OF INFO(6) AND INFO(7).
                const int* spiwdim,      // ON INPUT (I): HAS TO BE AT LEAST:
                                               // SPSTRAT=1: LINSYS=1,3:  1
                                               //            LINSYS=2  :  136*(IELENG+EQLENG)+32*(MIE+MEQ)+3*N+22
                                               //            IF YOU USE DYNAMIC ALLOCATION (ICNTL(13)=1) THEN YOU CAN USE LOWER INITIAL VALUES FOR SPIWDIM 
                                               //            AND LET YOUR CALLING ROUTINE DO THE REST.
                                               // SPSTRAT=2: LINSYS=1,3:  1
                                               //            LINSYS=2  : (IELENG)*5+2*MACTIV+1 (ESTIMATION)
                double* spdw,                 // (D) DOUBLE PRECISION WORKING ARRAY OF DIMENSION AT LEAST SPDWDIM. FOR
                const int* spdwdim,      // ON INPUT (I): HAS TO BE AT LEAST:
                                               // SPSTRAT=1, LINSYS=1:  (MACTIV+MEQ)**2
                                               //            LINSYS=2:  7*PRODLENG+2*(MIE+MEQ)   (ESTIMATION)
                                               //            LINSYS=3:  1
                                               //            IF YOU USE DYNAMIC ALLOCATION (ICNTL(13)=1) THEN YOU CAN USE ...
                                               // SPSTRAT=2, LINSYS=1:  N**2
                                               //            LINSYS=3:  1
                const int* spstrat,      // ON INPUT (I): SUBPROBLEM SOLUTION STATEGY:
                                               // 0: AUTOMATIC CHOICE BY THE PROGRAM
                                               // 1: APPROACH (S3) OF THE REFERENCE 4, ABOVE
                                               // 2: APPROACH (S4) OF THE REFERENCE 4,
                const int* linsys,       // ON INPUT (I): LINEAR SYSTEMS SOLVER THAT SHOULD BE APPLIED.
                                               // 1: DENSE CHOLESKY SOLVER
                                               // 2: SPARSE CHOLESKY SOLVER.
                                               // 3: CG SOLVER
               int* linear,             // ON INPUT (I): =0 , CONSTRAINT IS LIENAR, 1 NONLINEAR
			   const int* mf,		    // ON INPUT (I): NUMBER OF FEASIBILITY CONSTRAINTS
			   const int* lact,			// ON INPUT (I): NUMBER OF CONSTRAINTS THAT HAVE TO BE ACTIVE IN EACH ITERATE
			   int* setact);			// ON INPUT (I): ARRAY OF CONSTRAINTS WHICH HAVE TO BE ACTIVE IN EACH ITERATION


	 void CONCAVITY_D_2X2(
		 const double* x_init,
		 const int* m_elements,
		 double* m_h_org,
		 const double* eps,
		 const int* m_nie,
		 int* flag);

	 void CONCAVITY_GRAD_2X2(
		 double* m_iederv,
		 int* m_iern,
		 int* m_iecn,
		 int* m_ieleng,
		 const int* m_active,
		 const int* m_nie,
		 int* m_mactiv,
		 const double* m_xinit,
		 const int* m_elements,
		 const double* eps);


     void CONCAVITY_D_3X3(
		 const double* x_init,
		 const int* m_elements,
		 double* m_h_org,
		 const double* eps,
		 const int* m_nie,
		 int* flag);

	 void CONCAVITY_GRAD_3X3(
		 double* m_iederv,
		 int* m_iern,
		 int* m_iecn,
		 int* m_ieleng,
		 const int* m_active,
		 const int* m_nie,
		 int* m_mactiv,
		 const double* m_xinit,
		 const int* m_elements,
		 const double* eps);

     void CONCAVITY_D_6X6(	 
		 const double* x_init,
		 const int* m_elements,
		 double* m_h_org,
		 const double* eps,
		 const int* m_nie,
		 int* flag);

	 void CONCAVITY_GRAD_6X6(
		 double* m_iederv,
		 int* m_iern,
		 int* m_iecn,
		 int* m_ieleng,
		 const int* m_active,
		 const int* m_nie,
		 int* m_mactiv,
		 const double* m_xinit,
		 const int* m_elements,
		 const double* eps);

     void XEND(	 
		 const double* x_init,
		 const int* m_nvars);

# else
  void scpip40i_(const int* n,            // ON INPUT (INTEGER = I): NUMBER OF VARIABLES, AT LEAST >=1. NOT ALTERED (NA)
                const int* mie,          // ON INPUT (I): NUMBER OF INEQUALITY CONSTRAINTS (NA)
                const int* meq,          // ON INPUT (I): NUMBER OF EQUALITY CONSTRAINTS (NA)
                const int* iemax,        // ON INPUT (I): DIMENSION OF INEQUALITY DEPENDENT ARRAYS H_ORG, Y_IE,ACTIVE. MUST BE AT LEAST MIE AND >=1. (NA)
                const int* eqmax,        // ON INPUT (I): DIMENSION OF EQUALITY DEPENDENT ARRAYS G_ORG, Y_EQ. MUST BE AT LEAST MEQ AND >=1. (NA) 
                double* x0,                   // ON INPUT (DOUBLE PRECISION = D): X0 HAS TO CONTAIN AN INITIAL GUESS OF THE SOLUTION 
                                               // ON RETURN: X0 IS REPLACED BY THE LAST COMPUTED ITERATE
                const double* x_l,      // ON INPUT (D): LOWER BOUNDS ON THE VARIABLES (NA)  
                const double* x_u,      // ON INPUT (D): UPPER BOUNDS ON THE VARIABLES (NA)
                double* f_org,                // ON RETURN (D), F_ORG CONTAINS THE OBJECTIVE FUNCTION VALUE OF THE LAST COMPUTED ITERATE
                double* h_org,                // ON RETURN (D), H_ORG CONTAINS THE VALUES OF THE INEQUALITY CONSTRAINTS AT THE LAST COMPUTED ITERATE
                double* g_org,                // ON RETURN (D), G_ORG CONTAINS THE VALUES OF THE EQUALITY CONSTRAINTS AT THE LAST COMPUTED ITERATE
                double* df,                   // ON RETURN (D), DF CONTAINS THE GRADIENT OF THE OBJECTIVE FUNCTION AT X0
                double* y_ie,                 // ON RETURN (D), Y_IE CONTAINS THE LAGRANGE MULTIPLIERS FOR THE INEQUALITY CONSTRAINTS AT THE LAST COMPUTED ITERATE  
                double* y_eq,                 // ON RETURN (D), Y_EQ CONTAINS THE LAGRANGE MULTIPLIERS FOR THE EQUALITY CONSTRAINTS AT THE LAST COMPUTED ITERATE
                double* y_l,                  // ON RETURN (D), Y_L CONTAINS THE LAGRANGE MULTIPLIERS FOR THE LOWER BOUNDS ON THE VARIABLES AT THE LAST COMPUTED ITERATE
                double* y_u,                  // ON RETURN (D), Y_U CONTAINS THE LAGRANGE MULTIPLIERS FOR THE UPPER BOUNDS ON THE VARIABLES AT THE LAST COMPUTED ITERATE
                const int* icntl,       // ON INPUT (I): INTEGER ARRAY OF DIMENSION .. TO BE SET BY THE USER, BUT ALL COMPONENTS HAVE REASONABLE DEFAULTS (NA).
                const double* rcntl,    // ON INPUT (D): DOUBLE PRECISION ARRAY OF DIMENSION 6 TO BE SET BY THE USER, BUT ALL COMPONENTS HAVE REASONABLE DEFAULTS (NA). 
                int* info,                    // ON RETURN (I): INTEGER ARRAY CONTAINING SOME PROBLEM INFORMATION.
                double* rinfo,                // ON RETURN (D): DOUBLE PRECISION ARRAY OF DIMENSION 5. CONTAINS INFORMATION ON SOME DATA
                const int* nout,         // ON INPUT (I), INTEGER, INDICATING OUTPUT UNIT NUMBER (NA)
                double* r_scp,                // (D) REAL WORKING ARRAY OF DIMENSION AT LEAST RDIM
                const int* rdim,         // ON INPUT (I): MUST BE AT LEAST 30*N+11*IEMAX+8+10*EQMAX (NA)
                double* r_sub,                // (D) REAL WORKING ARRAY OF DIMENSION AT LEAST RSUBDIM
                const int* rsubdim,      // ON INPUT (I): MUST BE AT LEAST 22*N+41*IEMAX+27*EQMAX+2*IELPAR+EQLPAR (NA)
                int* i_scp,                   // (I) INTEGER WORKING ARRAY OF DIMENSION AT LEAST IDIM
                const int* idim,         // ON INPUT (I): MUST BE AT LEAST 5*N+5*IEMAX+2*EQMAX+3 (NA)
                int* i_sub,                   // (I) INTEGER WORKING ARRAY OF DIMENSION AT LEAST ISUBDIM
                const int* isubdim,      // ON INPUT (I): MUST BE AT LEAST 2*N+3*IEMAX+2*EQMAX+IELPAR (NA)
                int* active,                  // ON RETURN (I): INTEGER ARRAY INDICATING THE INEQUALITY CONSTRAINTS ...
                                               // ON INPUT (I): A USER IS ALLOWED TO CHANGE ACTIVE BEFORE COMPUTING ...
                const int* mode,         // ON INPUT (I): (NA)
                                               // 1: USUAL SUBROUTINE CALL FOR FUNCTION VALUES AND GRADIENTS
                                               // 2: REVERSE COMMUNICATION
                int* ierr,                    // ON (VERY FIRST) INPUT (I): HAS INITIALLY TO BE 0. DO NOT ALTER IERR OUTSIDE THE ROUTINE!
                                               // ON RETURN:
                                               // <0 : REVERSE COMMUNICATION:
                                               // -1: FUNCTION VALUES ARE REQUESTED
                                               // -2: GRADIENTS ARE REQUESTED
                                               // -3: INFO(6) AND INFO(7) ARE COMPUTED. THE USER CAN NOW DIMENSION SPIW AND SPDW DYNAMICALLY. 
                                               //     IF THIS IS NOT DESIRED, JUST REJUMP INTO SCPIP AFTER IERR = -3.
                                               // 0: SUCCESSFUL COMPUTATION ...
                int* iern,                    // ON INPUT (I): IERN CONTAINS THE ROW NUMBERS, I.E. THE INDEX OF THE COMPONENT OF THE ENTRIES 
                                               // OF THE JACOBIAN OF THE INEQUALITY CONSTRAINTS (NA)
                int* iecn,                    // ON INPUT (I): IECN CONTAINS THE COLUMN NUMBERS, I.E. THE NUMBER OF THE
                                               // INEQUALITY OF THE ENTRIES OF THE JACOBIAN OF THE INEQUALITY CONSTRAINTS (NA)
                double* iederv,                  // ON INPUT (D): IEDERV CONTAINS THE VALUES OF THE ENTRIES OF THE JACOBIAN OF THE INEQUALITY CONSTRAINTS. ...
                const int* ielpar,       // ON INPUT (I): DIMENSION OF ARRAYS IERN, IECN AND IEDERV. MUST BE AT LEAST IELENG AND >=1. (NA)
                const int* ieleng,       // ON INPUT (I): ACTUAL NUMBER OF ENTRIES IN IEDERV (NA)
                int* eqrn,                    // ON INPUT (I): EQRN CONTAINS THE ROW NUMBERS, I.E. THE INDEX OF THE COMPONENT OF THE ENTRIES OF 
                                               // THE JACOBIAN OF THE EQUALITYCONSTRAINTS (NA)
                int* eqcn,                    // ON INPUT (I): EQCN CONTAINS THE COLUMN NUMBERS, I.E. THE NUMBER OF THE 
                                               // EQUALITY OF THE ENTRIES OF THE JACOBIAN OF THE EQUALITY CONSTRAINTS (NA)
                double* eqcoef,               // ON INPUT (D): EQCOEF CONTAINS THE VALUES OF THE ENTRIES OF THE JACOBIAN OF THE EQUALITY CONSTRAINTS.
                const int* eqlpar,       // ON INPUT (I): DIMENSION OF ARRAYS EQRN, EQCN AND EQCOEF. MUST BE AT LEAST EQLENG AND >=1. (NA)
                const int* eqleng,       // ON INPUT (I): ACTUAL NUMBER OF ENTRIES IN EQCOEF (NA)
                int* mactiv,                   // ON RETURN (I): NUMBER OF CONSTRAINTS (EQUALITIES AND INEQUALITIES)... 
                                               // THIS AFFECTS THE REQUIRED WORKING SPACE, SEE SPDWDIM BELOW, AND THUS THE SUBPROBLEM SOLUTION METHOD
                int* spiw,                    // (I) INTEGER WORKING ARRAY OF DIMENSION AT LEAST SPIWDIM. 
                                               // AFTER IERR = -3 ONE CAN DIMENSION SPIW AND SPDW DYNAMICALLY WITH THE VALUES OF INFO(6) AND INFO(7).
                const int* spiwdim,      // ON INPUT (I): HAS TO BE AT LEAST:
                                               // SPSTRAT=1: LINSYS=1,3:  1
                                               //            LINSYS=2  :  136*(IELENG+EQLENG)+32*(MIE+MEQ)+3*N+22
                                               //            IF YOU USE DYNAMIC ALLOCATION (ICNTL(13)=1) THEN YOU CAN USE LOWER INITIAL VALUES FOR SPIWDIM 
                                               //            AND LET YOUR CALLING ROUTINE DO THE REST.
                                               // SPSTRAT=2: LINSYS=1,3:  1
                                               //            LINSYS=2  : (IELENG)*5+2*MACTIV+1 (ESTIMATION)
                double* spdw,                 // (D) DOUBLE PRECISION WORKING ARRAY OF DIMENSION AT LEAST SPDWDIM. FOR
                const int* spdwdim,      // ON INPUT (I): HAS TO BE AT LEAST:
                                               // SPSTRAT=1, LINSYS=1:  (MACTIV+MEQ)**2
                                               //            LINSYS=2:  7*PRODLENG+2*(MIE+MEQ)   (ESTIMATION)
                                               //            LINSYS=3:  1
                                               //            IF YOU USE DYNAMIC ALLOCATION (ICNTL(13)=1) THEN YOU CAN USE ...
                                               // SPSTRAT=2, LINSYS=1:  N**2
                                               //            LINSYS=3:  1
                const int* spstrat,      // ON INPUT (I): SUBPROBLEM SOLUTION STATEGY:
                                               // 0: AUTOMATIC CHOICE BY THE PROGRAM
                                               // 1: APPROACH (S3) OF THE REFERENCE 4, ABOVE
                                               // 2: APPROACH (S4) OF THE REFERENCE 4,
                const int* linsys,       // ON INPUT (I): LINEAR SYSTEMS SOLVER THAT SHOULD BE APPLIED.
                                               // 1: DENSE CHOLESKY SOLVER
                                               // 2: SPARSE CHOLESKY SOLVER.
                                               // 3: CG SOLVER
               int* linear,             // =0 , constraint is linear, 1 else
			   const int* mf,		    // ON INPUT (I): NUMBER OF FEASIBILITY CONSTRAINTS
			   const int* lact,			// ON INPUT (I): NUMBER OF CONSTRAINTS THAT HAVE TO BE ACTIVE IN EACH ITERATE
			   int* setact);			// ON INPUT (I): ARRAY OF CONSTRAINTS WHICH HAVE TO BE ACTIVE IN EACH ITERATION

     void concavity_d_2x2_(
		 const double* x_init,
		 const int* m_elements,
		 double* m_h_org,
		 const double* eps,
		 const int* m_nie,
		 int* flag);

	 void concavity_grad_2x2_(
		 double* m_iederv,
		 int* m_iern,
		 int* m_iecn,
		 int* m_ieleng,
		 const int* m_active,
		 const int* m_nie,
		 int* m_mactiv,
		 const double* m_xinit,
		 const int* m_elements,
		 const double* eps);

     void concavity_d_3x3_(
		 const double* x_init,
		 const int* m_elements,
		 double* m_h_org,
		 const double* eps,
		 const int* m_nie,
		 int* flag);

	 void concavity_grad_3x3_(
		 double* m_iederv,
		 int* m_iern,
		 int* m_iecn,
		 int* m_ieleng,
		 const int* m_active,
		 const int* m_nie,
		 int* m_mactiv,
		 const double* m_xinit,
		 const int* m_elements,
		 const double* eps);
 
	 void concavity_d_6x6_(
		 const double* x_init,
		 const int* m_elements,
		 double* m_h_org,
		 const double* eps,
		 const int* m_nie,
		 int* flag);

	 void concavity_grad_6x6_(
		 double* m_iederv,
		 int* m_iern,
		 int* m_iecn,
		 int* m_ieleng,
		 const int* m_active,
		 const int* m_nie,
		 int* m_mactiv,
		 const double* m_xinit,
		 const int* m_elements,
		 const double* eps);

# endif
} // end of extern_C




#endif /*SCPIP30_H_*/
