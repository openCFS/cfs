1

   Abaqus Student Edition 2018                                  Date 04-Jän-2020   Time 12:56:17
   

                                                                                 
                                                                                 
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

     CHARACTERISTIC ELEMENT LENGTH      2.12    

     DETAILS REGARDING ACTUAL SOLUTION WAVEFRONT REQUESTED

     DETAILED OUTPUT OF DIAGNOSTICS TO DATABASE REQUESTED

     PRINT OF INCREMENT NUMBER, TIME, ETC., TO THE MESSAGE FILE EVERY     1  INCREMENTS

  NUMBER OF EQUATIONS =      1950     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 
     ITERATION 1
      1       1.023445e+05     
      2       -9.280161e+05    
      3       -9.297450e+05    
      4       -9.307332e+05    
      5       -9.324669e+05    
      6       -9.398199e+05    
      7       -9.627406e+05    
      8       -9.887388e+05    
      9       -1.039888e+06    
      10      -1.282108e+06    
 
     NUMBER OF EIGENVALUES CONVERGED UPTO THIS ITERATION = 0
 
     ITERATION 2
      1       2.222190e+02     
      2       2.014242e+03     
      3       5.849875e+03     
      4       1.388931e+04     
      5       3.768682e+04     
      6       2.046473e+05     
      7       5.683868e+05     
      8       -1.132024e+06    
      9       -1.354215e+06    
      10      -2.681391e+06    
 
     NUMBER OF EIGENVALUES CONVERGED UPTO THIS ITERATION = 0
 
     ITERATION 3
      1       2.222043e+02     
      2       1.994503e+03     
      3       5.511119e+03     
      4       1.072282e+04     
      5       1.763969e+04     
      6       2.663555e+04     
      7       3.648788e+04     
      8       5.926973e+04     
      9       6.991975e+04     
      10      1.253597e+05     
 
     NUMBER OF EIGENVALUES CONVERGED UPTO THIS ITERATION = 0
 
     ITERATION 4
      1       2.222043e+02     
      2       1.994502e+03     
      3       5.511024e+03     
      4       1.071742e+04     
      5       1.753937e+04     
      6       2.588770e+04     
      7       3.566720e+04     
      8       4.744113e+04     
      9       5.949008e+04     
      10      7.764319e+04     
 
     NUMBER OF EIGENVALUES CONVERGED UPTO THIS ITERATION = 2
 
     ITERATION 5
      1       2.222043e+02     
      2       1.994502e+03     
      3       5.511024e+03     
      4       1.071740e+04     
      5       1.753646e+04     
      6       2.587310e+04     
      7       3.562288e+04     
      8       4.687493e+04     
      9       5.910064e+04     
      10      7.508424e+04     
 
     NUMBER OF EIGENVALUES CONVERGED UPTO THIS ITERATION = 3
 
     ITERATION 6
      1       2.222043e+02     
      2       1.994502e+03     
      3       5.511024e+03     
      4       1.071740e+04     
      5       1.753634e+04     
      6       2.587212e+04     
      7       3.561653e+04     
      8       4.671766e+04     
      9       5.896399e+04     
      10      7.392010e+04     
 
     NUMBER OF EIGENVALUES CONVERGED UPTO THIS ITERATION = 4
 
     ITERATION 7
      1       2.222043e+02     
      2       1.994502e+03     
      3       5.511024e+03     
      4       1.071740e+04     
      5       1.753633e+04     
      6       2.587204e+04     
      7       3.561559e+04     
      8       4.666994e+04     
      9       5.890436e+04     
      10      7.327253e+04     
 
     NUMBER OF EIGENVALUES CONVERGED UPTO THIS ITERATION = 5
1

   Abaqus Student Edition 2018                                  Date 04-Jän-2020   Time 12:56:17
   

                                                                                 
                                                                                 
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

     CHARACTERISTIC ELEMENT LENGTH      2.12    

     DETAILS REGARDING ACTUAL SOLUTION WAVEFRONT REQUESTED

     DETAILED OUTPUT OF DIAGNOSTICS TO DATABASE REQUESTED

     PRINT OF INCREMENT NUMBER, TIME, ETC., TO THE MESSAGE FILE EVERY     1  INCREMENTS

     COLLECTING MODEL CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS

     COLLECTING STEP CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS


  INCREMENT     1 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT   1.00    

  NUMBER OF EQUATIONS =      1950     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      9.931E-03   TIME AVG. FORCE       9.931E-03
 LARGEST RESIDUAL FORCE            -4.037E-13   AT NODE        195   DOF  1
 LARGEST INCREMENT OF DISP.        -8.329E-05   AT NODE         39   DOF  1
 LARGEST CORRECTION TO DISP.       -8.329E-05   AT NODE         39   DOF  1
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
                       2  WARNING MESSAGES DURING USER INPUT PROCESSING
                       0  WARNING MESSAGES DURING ANALYSIS
                       0  ANALYSIS WARNINGS ARE NUMERICAL PROBLEM MESSAGES
                       0  ANALYSIS WARNINGS ARE NEGATIVE EIGENVALUE MESSAGES
                       0  ERROR MESSAGES



     JOB TIME SUMMARY
       USER TIME (SEC)      =  0.60000    
       SYSTEM TIME (SEC)    =  0.10000    
       TOTAL CPU TIME (SEC) =  0.70000    
       WALLCLOCK TIME (SEC) =          0
