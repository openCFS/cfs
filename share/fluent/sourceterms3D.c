/* This UDFs are for 2D and 3D computation, but they are not fully.
 * Transformation is just implemented for rotation around z-axis!
 * 
 * Computation of time derivative of pressure, convective term and
 * transformation term of the derivative.
 * 
 * Be careful not to overwrite UDMs in different functions.
 * (Maybe specify them in on_loading)
 *
 *                UDMs
 * -------------------
 * convective_term   4
 * dpdt              6
 * transform         1
 * -------------------
 * so in total      11
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
    Message0("\n::: Usage of sourceterms.c :::::::::::::::::::\n"
        "1. Make sure 11 user defined memories are activated at Define/User-Defined/...\n"
        "2. Hook all functions at Define/User-Defined/Funktion Hooks\n"
        "3. Run simulation\n"
        "\n Transformation is just implemented for rotation around z-axis!\n\n"
        "Results are:\n"
        "UDM 0-2 velocity component sums\n"
        "UDM   3 convective term\n"
        "UDM 4-8 previous pressure fileds\n"
        "UDM   9 dpdt\n"
        "UDM  10 transformation term\n"
        "::::::::::::::::::::::::::::::::::::::::::::::\n");

    /* Check allocated variables */
/*
    if (N_UDS < 3)
    {
        Message("There are just %d UDS specified. "
                "You must at least specify 3 UDS \n",N_UDS);
    }
*/
    if (N_UDM < 11)
    {
        Message0("::::::::::::::::::::::::::::::::::::::::::::::\n"
            "There are just %d UDMs specified.\n"
            "You must at least specify 11 UDMs\n"
            "::::::::::::::::::::::::::::::::::::::::::::::\n",N_UDM);
    }

    /* get initial time step number for the calculation of sampling time */
    init_ts = N_TIME;

}

/* convective_term computes (u_mean*nabla)p */
DEFINE_EXECUTE_AT_END(convective_term)
{
    Domain *d;
    Thread *ct;
    cell_t c;
    
    d = Get_Domain(1);

    int num_samples;

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
            C_UDMI(c,ct,0) = C_UDMI(c,ct,0) + C_U(c,ct); /* sum u */
            C_UDMI(c,ct,1) = C_UDMI(c,ct,1) + C_V(c,ct); /* sum v */
            #if RP_3D
            C_UDMI(c,ct,2) = C_UDMI(c,ct,2) + C_W(c,ct); /* sum w */
            #endif

            /* computing convective term */
            #if RP_2D
            C_UDMI(c,ct,3) = (C_UDMI(c,ct,0)/num_samples*C_P_G(c,ct)[0]
                            + C_UDMI(c,ct,1)/num_samples*C_P_G(c,ct)[1]);
            #endif
            #if RP_3D
            C_UDMI(c,ct,3) = (C_UDMI(c,ct,0)/num_samples*C_P_G(c,ct)[0]
                            + C_UDMI(c,ct,1)/num_samples*C_P_G(c,ct)[1]
                            + C_UDMI(c,ct,2)/num_samples*C_P_G(c,ct)[2]);
            #endif
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
            C_UDMI(c,ct,9) = 1.0/(16.0*CURRENT_TIMESTEP) * (C_P(c,ct) + 3.0*C_UDMI(c,ct,4) 
+ 2.0*C_UDMI(c,ct,5) - 2.0*C_UDMI(c,ct,6) - 3.0*C_UDMI(c,ct,7) - C_UDMI(c,ct,8));

            /* Updating of the old pressure fields */
            C_UDMI(c,ct,8) = C_UDMI(c,ct,7); /* n-5 */
            C_UDMI(c,ct,7) = C_UDMI(c,ct,6); /* n-4 */
            C_UDMI(c,ct,6) = C_UDMI(c,ct,5); /* n-3 */
            C_UDMI(c,ct,5) = C_UDMI(c,ct,4); /* n-2 */
            C_UDMI(c,ct,4) = C_P(c,ct); /* n-1 */
        }
        end_c_loop(c,ct)
    }
}

/* transform computes (u_r*nabla)p */
DEFINE_EXECUTE_AT_END(transform)
{
    Domain *d;
    Thread *ct;
    cell_t c;
    
    d = Get_Domain(1);
    real omega,x[ND_ND];
    int i;
    /* real X[ND_ND]; */


    thread_loop_c(ct,d)
    {
/* to access fluid rotation data check out:
 * https://support.ansys.com/AnsysCustomerPortal/en_us/Knowledge%20Resources/Solutions/FLUENT/2006435
 */
   
      /* get rotation properties of current thread */
/*
      for(i=0;i<ND_ND;i++)
        {
            Message0("\n for thread %i ", ct);
            
            X[i]=THREAD_VAR(ct).fluid.origin[i];
            Message0("\n fluid_origin %i \t %g", i, X[i]);
            X[i]=THREAD_VAR(ct).fluid.axis[i];
            Message0("\n fluid_axis %i \t %g", i, X[i]);
            X[i]=THREAD_VAR(ct).fluid.velocity[i];
            Message0("\n fluid_velocity %i \t %g \n", i, X[i]);
        }
*/
        omega=THREAD_VAR(ct).fluid.omega;
        /* Message0("\n omega =  %g \n",omega); */

        begin_c_loop(c,ct)
        {
            
            /* computing transformation term */
            /* u_r = omega * r */

            C_CENTROID(x,c,ct) /* getting vector x to cell cetroid */

            /* correction if rotation is not around origin */
            for(i=0;i<ND_ND;i++)
            {
                x[i] = x[i] - THREAD_VAR(ct).fluid.origin[i];
            }

            /* compution (u_r * nabla) p */
            #if RP_2D
            C_UDMI(c,ct,10) = (omega*x[0]*C_P_G(c,ct)[1]
                             - omega*x[1]*C_P_G(c,ct)[0]);
            #endif
                                                         
            #if RP_3D
            /* So far just rotation around z axis implemented! */
            C_UDMI(c,ct,10) = (omega*x[0]*C_P_G(c,ct)[1]
                             - omega*x[1]*C_P_G(c,ct)[0]);
 
            /* correction for the rotation around a different axis */
            /* ... */
/*
            C_UDSI(c,ct,0) = (u_r*C_P_G(c,ct)[0]
                            + v_r*C_P_G(c,ct)[1]
                            + w_r*C_P_G(c,ct)[2]);
*/
           #endif

        }
        end_c_loop(c,ct)
    }
}
