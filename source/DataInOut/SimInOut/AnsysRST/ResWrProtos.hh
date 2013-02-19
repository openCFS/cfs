// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_SIM_RESWRPROTOS_HH
#define FILE_CFS_SIM_RESWRPROTOS_HH

namespace CoupledField
{

  // external Fortran functions

#if defined (__cplusplus)
  extern "C" {
 #endif

    // primary function:    Open result file
    //  input arguments:
    //     Nunit    (int,sc,in)       - Fortran Unit number for file (ANSYS uses 12)
    //     Lunit    (int,sc,in)       - Current print output unit (usually 6)
    //     Fname    (ch*(ncFname),sc,in) - The name (with extension) for the file
    //     ncFname  (int,sc,in)       - Number of characters in Fname
    //     Title    (ch*80,ar(2),in)  - Title and First subtitle
    //     JobName  (ch*8,sc,in)      - Jobname for file
    //     Units    (int,sc,in)       - 0, unknown  1, SI  2, CSG 
    //                                  3, U.S. Customary - foot
    //                                  4, U.S. Customary - inch
    //                                  6, MPA
    //     NumDOF   (int,sc,in)       - Number of DOF per node
    //     DOF      (int,ar(*),in)    - The DOFs per node
    //     UserCode (int,sc,in)       - Code for this application (user defined for
    //                                  now)
    //     MaxNode  (int,sc,in)       - Maximum node number used
    //     NumNode  (int,sc,in)       - Number of nodes used
    //     MaxElem  (int,sc,in)       - Maximum element number used
    //     NumElem  (int,sc,in)       - Number of elements used
    //     MaxResultSet (int,sc,in)   - Maximum number of result sets (usually 1000)

    //  output arguments:
    //     ResWrBegin (int,sc,out)    - 0, successful  other, error in file open

    extern Integer reswrbegin_(Integer* Nunit,
                               Integer* Lunit,
                               char* Fname,
                               Integer* ncFname,
                               char* Title,
                               char* JobName,
                               Integer* Units,
                               Integer* NumDOF,
                               Integer* DOF,
                               Integer* UserCode,
                               Integer* MaxNode,
                               Integer* NumNode,
                               Integer* MaxElem,
                               Integer* NumElem,
                               Integer* MaxResultSet,
                               Integer* kan,
                               Integer lenFname,
                               Integer lenTitle,
                               Integer lenJobName);

    // primary function:    Write Geometry Header Record

    //  input arguments:
    //     MaxType  (int,sc,in)       - Maximum element type
    //     NumType  (int,sc,in)       - Number of element types
    //     MaxReal  (int,sc,in)       - Maximum real constant set number
    //     NumReal  (int,sc,in)       - Number of real constant sets
    //     MaxCsys  (int,sc,in)       - Maximum coordinate system number
    //     NumCsys  (int,sc,in)       - Number of user coordinate systems

    //  output arguments:  none

    extern void reswrgeombegin_(Integer* MaxType,
                                Integer* NumType,
                                Integer* MaxReal,
                                Integer* NumReal,
                                Integer* MaxCsys,
                                Integer* NumCsys);
  
    // primary function:    Start element type records
    extern void reswrtypebegin_();
  
    // primary function:    Write an element type record

    //  input arguments:
    //     itype    (int,sc,in)         - Element type number
    //     iel//     (int,ar(IELCSZ),in) - Element characteristics
    extern void reswrtype_(Integer* itype,
                           Integer ielc[150]);

    extern void reswrtypeend_();

    // primary function:    Start real constant records
    extern void reswrrealbegin_();
  

    // primary function:    Write an real constant record
    //  input arguments:
    //     iReal    (int,sc,in)       - Real set number
    //     n        (int,sc,in)       - Number of real constants in set
    //     Rcon     (dp,ar(n),in)     - Real Constants

    //  output arguments:  none

    extern void reswrreal_(Integer* iReal,
                           Integer* n,
                           Double* Rcon);
  
    extern void reswrrealend_();
    
    // primary function:    Start coordinate system records
    extern void reswrcsysbegin_();

    // primary function:    Write an coordinate system record
    //  input arguments:
    //     iCsys    (int,sc,in)       - Coordinate system number
    //     Csys     (dp,ar(24),in)    - Coordinate system description

    extern void reswrcsys_(Integer* iCsys,
                           Double* Csys);

    // primary function:    Write out coordinate system pointers
    extern void reswrcsysend_();
    
    // primary function:    Start node records
    extern void reswrnodebegin_();
    
    // primary function:    Store a node
    //  input arguments:
    //     iNode    (int,sc,in)       - The node number
    //     xyzang   (dp,ar(6),in)     - x,y,z,thxy,thyz,thzx for node

    extern void reswrnode_(Integer* iNode,
                           Double xyzang[6]);

    // primary function:    Finish Node Storage
    extern void reswrnodeend_();

    // primary function:   Start element output
    extern void reswrelembegin_();

    // primary function:    Write out an element
    //  input arguments:
    //     iElem    (int,sc,in)       - The element number
    //     n        (int,sc,in)       - Number of nodes
    //     nodes    (int,ar(n),in)    - Element nodes
    //     ElemData (int,ar(10),in)   - Element information
    //                                    mat    - material reference number
    //                                    type   - element type number
    //                                    real   - real constant reference number
    //                                    elnum  - element number
    //                                    esys   - element coordinate system
    //                                    death  - death flag
    //                                             = 0 - alive
    //                                             = 1 - dead
    //                                    solidm - solid model reference
    //                                    shape  - coded shape key
    //                                    pexcl  - P-Method exclude key
    //                                    0      - Not used
    extern void reswrelem_(Integer* iElem,
                           Integer* n,
                           Integer* nodes,
                           Integer ElemData[10]);
    
    // primary function:   Finish element writing
    extern void reswrelemend_();

    // primary function:    Finish Geometry writing
    extern void reswrgeomend_();

    // primary function:    Write the solution header records
    //  input arguments:
    //     lstep    (int,sc,in)       - Load step number
    //     substep  (int,sc,in)       - Substep of this load step
    //     ncumit   (int,sc,in)       - Cumulative iteration number
    //     time     (dp,sc,in)       - Current solution time
    //     Title    (ch*80,ar(5),in)  - Title and 4 subtitles
    //     DofLab   (ch*4,ar(nDOF),in)- Labels for DOFs
    extern void reswrsolbegin_ (Integer* lstep,
                                Integer* substep,
                                Integer* ncumit,
                                Double* time,
                                Integer* kcmplx,
                                char Title[400],
                                char* DofLab,
                                //                                Integer lenTitle,
                                Integer lenDofLab);

    // primary function:     Begin output of displacement vector
    extern void reswrdispbegin_();

    // primary function:    Store nodal displacements
    //  input arguments:
    //     node     (int,sc,in)       - Node number
    //     Disp     (dp,ar(nDOF),in)  - Displacements
    extern void reswrdisp_(Integer* node,
                           Double* Disp);

    // primary function:    Put displacement vector on result file
    extern void reswrdispend_();

    // nData     (dp,ar(numNodes*nDOF),in)  - Displacements
    extern void reswrnodaldata_(Double* nData);

    // primary function:   Start reaction force output
    //  input arguments:
    //     nRforce  (int,sc,in)       - Number of reaction forces
    extern void reswrrforbegin_(Integer* nRforce);

    // primary function:    Store a reaction force
    //  input arguments:
    //     node     (int,sc,in)       - External node number
    //     idof     (int,sc,in)       - Internal dof number
    //     value    (dp,sc,in)        - Value of reaction force
    extern void reswrrfor_(Integer* node,
                           Integer* idof,
                           Double* value);


    // primary function:   Write out the reaction forces
    extern void reswrrforend_();
    
    // primary function:    Header for Boundary Conditions
    //  input arguments:
    //     nFixed   (int,sc,in)       - Number of Constraints
    //     nForce   (int,sc,in)       - Number of Forces
    extern void reswrbcbegin_(Integer* nFixed,
                              Integer* nForce);

    // primary function:    Start constraint output
    //  input arguments:
    //     nFixed   (int,sc,in)          - Number of constraints
    extern void reswrfixbegin_(Integer* nFixed);

    // primary function:    Store a constraint value
    //  input arguments:
    //     node     (int,sc,in)       - External node number
    //     idof     (int,sc,in)       - Internal dof number
    //     value    (dp,ar(4),in)     - Real,Imag, RealOld,ImagOld
    extern void reswrfix_(Integer* node,
                          Integer* idof,
                          Double* value);

    // primary function:    Write out the constraints
    extern void reswrfixend_();

    // primary function:    Start applied force output
    extern void reswrforcbegin_(Integer* nForce);

    //  input arguments:
    //     node     (int,sc,in)       - External node number
    //     idof     (int,sc,in)       - Internal dof number
    //     value    (dp,ar(4),in)     - Real,Imag, RealOld,ImagOld
    extern void reswrforc_(Integer* node,
                           Integer* idof,
                           Double* value);

    // primary function:    End of applied forces
    extern void reswrforcend_();

    // primary function:   End of BC output
    extern void reswrbcend_();

    // primary function:    Set up index array for element results
    extern void reswreresbegin_();

    // primary function:    Start of element results
    //  input arguments:
    //     iElem    (int,sc,in)        - Element Number
    extern void reswrestrbegin_(Integer* iElem);
                                
    // primary function:   Store element stress results
    //  input arguments:        
    //     iStr    (int,sc,in)       - Stress record number (1-25)
    //     nStr    (int,sc,in)       - Number of data words
    //     Str     (dp,ar(nStr),in)  - Stress values
    extern void reswrestr_(Integer* iStr,
                           Integer* nStr,
                           Double* Str);

    // primary function:     End of output for an element
    extern void reswrestrend_();

    // primary function:    Write out element index
    extern void reswreresend_();

    // primary function:    End of a Solution Result Set
    extern void reswrsolend_();

    // primary function:    Close result file
    extern void reswrend_();

    // Write just a .rst header for determining binlib version
    int reswrtest_(char* FName, int* ncFName, int lenFName);
  
#if defined (__cplusplus)
}
#endif
  
} // end of namespace

#endif
