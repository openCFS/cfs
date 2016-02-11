/* This UDFs are for 2D computations onely. They are not fully tested so be
 * careful to use them.
 *
 * Computation of time derivative of pressure, laplacian of pressure, 
 * second spatial derivative of the lighthill tensor and convective term.
 * 
 * Be careful not to overwrite UDMs in different functions.
 * (Maybe specify them in on_loading!!!!!!!!!)
 * 
 *                 UDMs  UDSs
 * --------------------------
 * laplace_p          1     2
 * DivDivLighthill    1     2
 * convective_term    3     - 
 * dpdt               6     -
 * --------------------------
 * so in total       11     4
 *
 * Make sure to allocate all of these variables and hooking them before
 * running the simulation. Don't hook the functions later, this would cause
 * trouble with the number of time samples for mean velocity.
 *
 */

#include "udf.h"
#include "mem.h"

int init_ts;

DEFINE_EXECUTE_ON_LOADING(on_loading, libname)
{

    Message("\n::: Usage of sourceterms2D.c :::::::::::::::::::\n"
        "1. Make sure 11 user defined memories are activated at Define/User-Defined/...\n"
        "2. Hook all functions at Define/User-Defined/Funktion Hooks\n"
        "3. Run simulation\n\n"
        "Results are:\n"
        "UDM 0-4 previous pressure fileds\n"
        "UDM   5 dpdt\n"
        "UDM 6-7 velocity component sums\n"
        "UDM   8 convective term\n"
        "UDM   9 Laplacian of pressure\n"
        "UDM  10 Second derivative of LighthillTensor\n"
	"UDS 0-1 Gradient of Pressure\n"
	"UDS 2-3 Divergence of Lighthill tensor\n"
        "::::::::::::::::::::::::::::::::::::::::::::::\n");

    /* Check allocated variables */

    if (N_UDS < 4)
    {
        Message("There are just %d UDS specified. "
                "You must at least specify 4 UDS \n",N_UDS);
    }
    if (N_UDM < 11)
    {
        Message("::::::::::::::::::::::::::::::::::::::::::::::\n"
            "There are just %d UDMs specified.\n"
            "You must at least specify 11 UDMs\n"
            "::::::::::::::::::::::::::::::::::::::::::::::\n",N_UDM);
    }

/*
	#if RP_3D
    if (N_UDM < 12 || N_UDS < 7)
    {
        Message("::::::::::::::::::::::::::::::::::::::::::::::\n"
            "There are just %d UDMs specified.\n"
            "You must at least specify 11 UDMs\n"
            "::::::::::::::::::::::::::::::::::::::::::::::\n",N_UDM);
    }
	#else
	if (N_UDM < 11 || N_UDS < 5)
    {
        Message("::::::::::::::::::::::::::::::::::::::::::::::\n"
            "There are just %d UDMs specified.\n"
            "You must at least specify 11 UDMs\n"
            "::::::::::::::::::::::::::::::::::::::::::::::\n",N_UDM);
    }
	#endif
*/
    /* get initial time step number for the calculation of sampling time */
    init_ts = N_TIME;

	/* TIME DERIVATIVE */
	Set_User_Memory_Name(0,"P_N_-1");
	Set_User_Memory_Name(1,"P_N_-2");
	Set_User_Memory_Name(2,"P_N_-3");
	Set_User_Memory_Name(3,"P_N_-4");
	Set_User_Memory_Name(4,"P_N_-5");
	Set_User_Memory_Name(5,"dpdt");
	
	/*MeanFlow */
	Set_User_Memory_Name(6,"UmeanAccum");
	Set_User_Memory_Name(7,"VmeanAccum");
	
	/*Convection and Laplacian*/
	Set_User_Scalar_Name(0,"dP_dx");
	Set_User_Scalar_Name(1,"dP_dy");	
	Set_User_Memory_Name(8,"UmeanInNablaP");

	Set_User_Memory_Name(9,"DeltaP");
	
	
        /*Lighthill Sourceterm*/
	Set_User_Scalar_Name(2,"DivLHTensorX");
	Set_User_Scalar_Name(3,"DivLHTensorY");
	Set_User_Memory_Name(10,"DivDivLHTensor");
	
}

DEFINE_ADJUST(laplace_p,d)
{
	Thread *ct;
    cell_t c;
    
    thread_loop_c(ct,d)
    {
        begin_c_loop(c,ct)
        {
            C_UDSI(c,ct,0) = C_P_G(c,ct)[0];
            C_UDSI(c,ct,1) = C_P_G(c,ct)[1];
            C_UDMI(c,ct,9) = C_UDSI_G(c,ct,0)[0]+C_UDSI_G(c,ct,1)[1];
        }
        end_c_loop(c,ct)
    }
}

DEFINE_ADJUST(DivDivLighthill,d)
{
    Thread *ct;
    cell_t c;
    
    thread_loop_c(ct,d)
    {
        begin_c_loop(c,ct)
        {	
	   C_UDSI(c,ct,2) = 2*C_U(c,ct) * C_U_G(c,ct)[0] + C_U(c,ct) * C_V_G(c,ct)[1] + C_V(c,ct) * C_U_G(c,ct)[1];
							  
	   C_UDSI(c,ct,3) =  C_V(c,ct) * C_U_G(c,ct)[0] + C_U(c,ct) * C_V_G(c,ct)[0] + 2*C_V(c,ct) * C_V_G(c,ct)[1];
			
	   C_UDMI(c,ct,10) = C_UDSI_G(c,ct,2)[0] + C_UDSI_G(c,ct,3)[1];
        }
        end_c_loop(c,ct)
    }
}

/* convective_term computes (u_mean*nabla)p */
DEFINE_EXECUTE_AT_END(convective_term)
{
    Domain *d;
    Thread *ct;
    cell_t c;
    int num_samples;

    d = Get_Domain(1);

    /* calculate number of samples for mean velocity */
    num_samples = N_TIME - init_ts;
/*
    Message("Initial number of steps %i \n",init_ts);
    Message("Current number of steps %i \n",N_TIME);
    Message("Number of samples %i \n",num_samples);
*/

    thread_loop_c(ct,d)
    {

        begin_c_loop(c,ct)
        {
            /* calculating sum of velcities */
            C_UDMI(c,ct,6) =   C_UDMI(c,ct,6)  + C_U(c,ct); /* sum u */
            C_UDMI(c,ct,7) =   C_UDMI(c,ct,7)  + C_V(c,ct); /* sum v */

            /* computing convective term */
            C_UDMI(c,ct,8) = (C_UDMI(c,ct,6)/num_samples*C_P_G(c,ct)[0]
                            + C_UDMI(c,ct,7)/num_samples*C_P_G(c,ct)[1]);
        } 
        end_c_loop(c,ct)
    }
}


/* dpdt uses a one sided differtiatior to compute a pressure time derivative of order 5
 * http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/#deriv
 */
DEFINE_EXECUTE_AT_END(dpdt)
{
    Domain *d;
    Thread *ct;
    cell_t c;
    
    d = Get_Domain(1);

    thread_loop_c(ct,d)
    {
       begin_c_loop(c,ct)
        {
            /* Calculation of the pressure time derivative O5 */
            C_UDMI(c,ct,5) = 1.0/(16.0*CURRENT_TIMESTEP) * (C_P(c,ct) + 3.0*C_UDMI(c,ct,0) 
+ 2.0*C_UDMI(c,ct,1) - 2.0*C_UDMI(c,ct,2) - 3.0*C_UDMI(c,ct,3) - C_UDMI(c,ct,4));

            /* Updating of the old pressure fields */	
	    C_UDMI(c,ct,4) = C_UDMI(c,ct,3); /* n-5 */
            C_UDMI(c,ct,3) = C_UDMI(c,ct,2); /* n-4 */
            C_UDMI(c,ct,2) = C_UDMI(c,ct,1); /* n-3 */
            C_UDMI(c,ct,1) = C_UDMI(c,ct,0); /* n-2 */
            C_UDMI(c,ct,0) = C_P(c,ct);      /* n-1 */
        }
        end_c_loop(c,ct)
    }
}
