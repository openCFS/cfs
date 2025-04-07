1

   Abaqus Student Edition 2018                                  Date 04-Jän-2020   Time 13:18:35
   

                                                                                 
                                                                                 
 STEP    1     INCREMENT     1     STEP TIME    0.00    


                        S T E P       1     C A L C U L A T I O N   O F   E I G E N V A L U E S

                                            F O R   B U C K L I N G   P R E D I C T I O N



                                                                                          

     THE SUBSPACE ITERATION METHOD IS USED FOR THIS ANALYSIS
     NUMBER OF EIGENVALUES                      5
     MAXIMUM NUMBER OF ITERATIONS                     30
     NUMBER OF VECTORS IN ITERATION                   10
                              ONLY INITIAL STRESS EFFECTS ARE INCLUDED IN THE
                              STIFFNESS MATRIX

     THIS IS A LINEAR PERTURBATION STEP.
     ALL LOADS ARE DEFINED AS CHANGE IN LOAD TO THE REFERENCE STATE

     EXTRAPOLATION WILL NOT BE USED

     CHARACTERISTIC ELEMENT LENGTH      2.74    

     DETAILS REGARDING ACTUAL SOLUTION WAVEFRONT REQUESTED

     DETAILED OUTPUT OF DIAGNOSTICS TO DATABASE REQUESTED

     PRINT OF INCREMENT NUMBER, TIME, ETC., TO THE MESSAGE FILE EVERY     1  INCREMENTS

  NUMBER OF EQUATIONS =      2976     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.21E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 
     ITERATION 1
      1       1.339188e+03     
      2       1.967448e+04     
      3       5.216100e+04     
      4       2.914920e+05     
      5       4.457789e+05     
      6       5.739771e+05     
      7       -1.263785e+06    
      8       -1.696845e+06    
      9       -3.826620e+06    
      10      6.664558e+06     
 
     NUMBER OF EIGENVALUES CONVERGED UPTO THIS ITERATION = 0
 
     ITERATION 2
      1       2.231476e+02     
      2       3.964552e+02     
      3       2.008446e+03     
      4       3.571516e+03     
      5       5.639380e+03     
      6       1.053760e+04     
      7       1.140517e+04     
      8       1.964881e+04     
      9       2.499309e+04     
      10      4.207240e+04     
 
     NUMBER OF EIGENVALUES CONVERGED UPTO THIS ITERATION = 0
 
     ITERATION 3
      1       2.231449e+02     
      2       3.964160e+02     
      3       2.005081e+03     
      4       3.552036e+03     
      5       5.553310e+03     
      6       9.821533e+03     
      7       1.087541e+04     
      8       1.800133e+04     
      9       2.047445e+04     
      10      2.831965e+04     
 
     NUMBER OF EIGENVALUES CONVERGED UPTO THIS ITERATION = 0
 
     ITERATION 4
      1       2.231449e+02     
      2       3.964160e+02     
      3       2.005072e+03     
      4       3.551884e+03     
      5       5.552188e+03     
      6       9.783465e+03     
      7       1.083863e+04     
      8       1.788941e+04     
      9       1.951022e+04     
      10      2.718541e+04     
 
     NUMBER OF EIGENVALUES CONVERGED UPTO THIS ITERATION = 2
 
     ITERATION 5
      1       2.231449e+02     
      2       3.964160e+02     
      3       2.005072e+03     
      4       3.551882e+03     
      5       5.552158e+03     
      6       9.780664e+03     
      7       1.083373e+04     
      8       1.784153e+04     
      9       1.914988e+04     
      10      2.688197e+04     
 
     NUMBER OF EIGENVALUES CONVERGED UPTO THIS ITERATION = 4
 
     ITERATION 6
      1       2.231449e+02     
      2       3.964160e+02     
      3       2.005072e+03     
      4       3.551882e+03     
      5       5.552157e+03     
      6       9.780443e+03     
      7       1.083309e+04     
      8       1.781932e+04     
      9       1.901252e+04     
      10      2.674267e+04     
 
     NUMBER OF EIGENVALUES CONVERGED UPTO THIS ITERATION = 5
1

   Abaqus Student Edition 2018                                  Date 04-Jän-2020   Time 13:18:37
   

                                                                                 
                                                                                 
 STEP    2     INCREMENT     1     STEP TIME    0.00    


                        S T E P       2     S T A T I C   A N A L Y S I S


                                                                                          

          FIXED TIME INCREMENTS
          TIME INCREMENT IS                                     1.00    
          TIME PERIOD IS                                        1.00    

     LINEAR EQUATION SOLVER TYPE         DIRECT SPARSE

 CONVERGENCE TOLERANCE PARAMETERS FOR FORCE    
     CRITERION FOR RESIDUAL FORCE     FOR A NONLINEAR PROBLEM          5.000E-03
     CRITERION FOR DISP.    CORRECTION IN A NONLINEAR PROBLEM          1.000E-02
     INIT. VALUE OF TIME AVG. FORCE     IS TIME AVG. FORCE     IN PREVIOUS STEP
     AVERAGE FORCE     IS TIME AVERAGE FORCE    
     ALTERNATE CRIT. FOR RESIDUAL FORCE     FOR A NONLINEAR PROBLEM    2.000E-02
     CRITERION FOR ZERO FORCE     RELATIVE TO TIME AVRG. FORCE         1.000E-05
     CRITERION FOR RESIDUAL FORCE     WHEN THERE IS ZERO FLUX          1.000E-05
     CRITERION FOR DISP.    CORRECTION WHEN THERE IS ZERO FLUX         1.000E-03
     CRITERION FOR RESIDUAL FORCE     FOR A LINEAR INCREMENT           1.000E-08
     FIELD CONVERSION RATIO                                             1.00    
     CRITERION FOR ZERO FORCE     REL. TO TIME AVRG. MAX. FORCE        1.000E-05
     CRITERION FOR ZERO DISP.    RELATIVE TO CHARACTERISTIC LENGTH     1.000E-08

     VOLUMETRIC STRAIN COMPATIBILITY TOLERANCE FOR HYBRID SOLIDS       1.000E-05
     AXIAL STRAIN COMPATIBILITY TOLERANCE FOR HYBRID BEAMS             1.000E-05
     TRANS. SHEAR STRAIN COMPATIBILITY TOLERANCE FOR HYBRID BEAMS      1.000E-05
     SOFT CONTACT CONSTRAINT COMPATIBILITY TOLERANCE FOR P>P0          5.000E-03
     SOFT CONTACT CONSTRAINT COMPATIBILITY TOLERANCE FOR P=0.0         0.100    
     CONTACT FORCE ERROR TOLERANCE FOR CONVERT SDI=YES                 1.00    
     DISPLACEMENT COMPATIBILITY TOLERANCE FOR DCOUP ELEMENTS           1.000E-05
     ROTATION COMPATIBILITY TOLERANCE FOR DCOUP ELEMENTS               1.000E-05

 EQUILIBRIUM WILL BE CHECKED FOR SEVERE DISCONTINUITY ITERATIONS

 ITERATION CONTROL PARAMETERS:
     FIRST EQUILIBRIUM ITERATION FOR CONSECUTIVE DIVERGENCE CHECK              4
     EQUILIBRIUM ITERATION AFTER WHICH ALTERNATE RESIDUAL IS USED              9
     MAXIMUM EQUILIBRIUM ITERATIONS ALLOWED                                   16
     MAXIMUM ITERATIONS FOR SEVERE DISCONTINUITIES                            12
     GLOBAL STABILIZATION CONTROL IS NOT USED

          PRINT OF INCREMENT NUMBER, TIME, ETC., EVERY    1  INCREMENTS

     THE MAXIMUM NUMBER OF INCREMENTS IN THIS STEP IS                     100

     LINEAR EXTRAPOLATION WILL BE USED

     CHARACTERISTIC ELEMENT LENGTH      2.74    

     DETAILS REGARDING ACTUAL SOLUTION WAVEFRONT REQUESTED

     DETAILED OUTPUT OF DIAGNOSTICS TO DATABASE REQUESTED

     PRINT OF INCREMENT NUMBER, TIME, ETC., TO THE MESSAGE FILE EVERY     1  INCREMENTS

     COLLECTING MODEL CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS

     COLLECTING STEP CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS


  INCREMENT     1 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT   1.00    

  NUMBER OF EQUATIONS =      2976     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.21E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =  0.10000    
       ELAPSED TOTAL CPU TIME (SEC) =  0.10000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.285E-03   TIME AVG. FORCE       4.285E-03
 LARGEST RESIDUAL FORCE            -5.656E-14   AT NODE        962   DOF  1
 LARGEST INCREMENT OF DISP.        -8.314E-05   AT NODE          2   DOF  1
 LARGEST CORRECTION TO DISP.       -8.314E-05   AT NODE          2   DOF  1
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED   1.00    ,  FRACTION OF STEP COMPLETED   1.00    
 STEP TIME COMPLETED        1.00    ,  TOTAL TIME COMPLETED         1.00    


          THE ANALYSIS HAS BEEN COMPLETED



     ANALYSIS SUMMARY:
     TOTAL OF          2  INCREMENTS
                       0  CUTBACKS IN AUTOMATIC INCREMENTATION
                       2  ITERATIONS INCLUDING CONTACT ITERATIONS IF PRESENT
                       2  PASSES THROUGH THE EQUATION SOLVER OF WHICH 
                       2  INVOLVE MATRIX DECOMPOSITION, INCLUDING
                       0  DECOMPOSITION(S) OF THE MASS MATRIX
                       1  REORDERING OF EQUATIONS TO MINIMIZE WAVEFRONT
                       0  ADDITIONAL RESIDUAL EVALUATIONS FOR LINE SEARCHES
                       0  ADDITIONAL OPERATOR EVALUATIONS FOR LINE SEARCHES
                       1  WARNING MESSAGES DURING USER INPUT PROCESSING
                       0  WARNING MESSAGES DURING ANALYSIS
                       0  ANALYSIS WARNINGS ARE NUMERICAL PROBLEM MESSAGES
                       0  ANALYSIS WARNINGS ARE NEGATIVE EIGENVALUE MESSAGES
                       0  ERROR MESSAGES



     JOB TIME SUMMARY
       USER TIME (SEC)      =   1.0000    
       SYSTEM TIME (SEC)    =  0.20000    
       TOTAL CPU TIME (SEC) =   1.2000    
       WALLCLOCK TIME (SEC) =          2
