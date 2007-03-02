#include <iostream>
#include <string>
#include <algorithm>

#include <def_use_pardiso.hh>
#include <def_use_metis.hh>

#include "utils/Exception.hh"
#include "external/ilupack/Ilupack.hh"
#include "matvec/matvec.hh"
#include "utils/utils.hh"
#include "utils/Enum.hh"
#include "solver/solver.hh"
#include "precond/precond.hh"


namespace OLAS 
{

template<typename T>
Ilupack<T>::Ilupack(OLAS_Params *myParams, OLAS_Report *myReport, MatrixEntryType type) 
{
    // Set pointers to communication objects
    myParams_  = myParams;
    myReport_  = myReport;
    
    if(type != COMPLEX && type != DOUBLE) throw Exception("unhandled type", type, __LINE__, __PRETTY_FUNCTION__);
    isComplex_ = type == COMPLEX;
    
    // set this to zero to signal the destructor if Setup() was called
    mat_.a = NULL;
    mat_.ia = NULL;
    mat_.ja = NULL;
    
    precondAge_ = 0;
    setupRuns_ = 0;
    
    // set the convenience pointers for the D* instances of ilupack structs.
    dmat_ptr_ = &mat_;
    zmat_ptr_ = reinterpret_cast<Zmat*>(&mat_);
    
    dparam_ptr_ = &param_;
    zparam_ptr_ = reinterpret_cast<ZILUPACKparam*>(&param_);

    dprecond_ptr_ = &precond_;
    zprecond_ptr_ = reinterpret_cast<ZAMGlevelmat*>(&precond_);
        
    // init the enums
    SetEnums();
    
    // set up the reference to the ilupack permutation functions */
    SetPermutations();
   
}



template<typename T>
Ilupack<T>::~Ilupack() 
{
    ReleaseIlupackMemory();
}

template<typename T>
void Ilupack<T>::ReleaseIlupackMemory() 
{
    // call the ilupack delete method only if Setup() which sets mat_ was called
    if(mat_.a == NULL) return;

    switch(matrix_)
    {
        case GNL: if(isComplex_) ZGNLAMGdelete(*zmat_ptr_, *zprecond_ptr_, nlev_, zparam_ptr_);
                            else DGNLAMGdelete(*dmat_ptr_, *dprecond_ptr_, nlev_, dparam_ptr_);
                  break;
        
        case SYM: if(isComplex_) ZSYMAMGdelete(*zmat_ptr_, *zprecond_ptr_, nlev_, zparam_ptr_);
                            else DSYMAMGdelete(*dmat_ptr_, *dprecond_ptr_, nlev_, dparam_ptr_);
                  break;

        case PD:  if(isComplex_) ZHPDAMGdelete(*zmat_ptr_, *zprecond_ptr_, nlev_, zparam_ptr_);
                            else DSPDAMGdelete(*dmat_ptr_, *dprecond_ptr_, nlev_, dparam_ptr_);
                  break;

        case HER: if(isComplex_) ZHERAMGdelete(*zmat_ptr_, *zprecond_ptr_, nlev_, zparam_ptr_);
                  break;
                  
        default: Exception("invalid matrix type", matrix_, __LINE__, __PRETTY_FUNCTION__);          
    }   
    
    mat_.a = NULL;           
}    

template<typename T>
void Ilupack<T>::SetMatrix(const BaseMatrix &base_mat)
{     
    const SparseOLASMatrix<T>& som = dynamic_cast<const SparseOLASMatrix<T>&>(base_mat);              

    // delete the old values before we allocate the new space - if could be faster! 
    if(mat_.a != NULL)
    {
        delete mat_.a;
        delete mat_.ia;
        delete mat_.ja;
    }

    UInt elements = 0;

    // pointers will point to the first actual element    
    const int* row_ptr = NULL;
    const int* col_ptr = NULL;
    const T*   val_ptr = NULL;

    if(som.GetStorageType() == SPARSE_SYM) 
    {
       const SCRS_Matrix<T>& scrs = dynamic_cast<const SCRS_Matrix<T>&>(som);        
       // Nnz is the total number = Diagonal + 2 * triagonals. But one is skipped!
       elements= scrs.GetNnz() - (UInt) (0.5 * (double) (scrs.GetNnz() - scrs.GetNrows()));
        
       row_ptr = scrs.GetRowPointer() + 1;
       col_ptr = scrs.GetColPointer() + 1;
       val_ptr = scrs.GetDataPointer() + 1;
    }   
    else 
    {
       const CRS_Matrix<T>& crs = dynamic_cast<const CRS_Matrix<T>&>(som);

       elements= crs.GetNnz();

       row_ptr = crs.GetRowPointer() + 1;
       col_ptr = crs.GetColPointer() + 1;
       val_ptr = crs.GetDataPointer() + 1;
    }   

    // allocate and copy the stuff
      
    // rows are simple but we have to handle the tailing element
    mat_.ia = new int[som.GetNrows() + 1]; // one plus due to CRS tail
    // ignore the olas 0th element garbage and add the one from above to the end
    std::copy(row_ptr, row_ptr + som.GetNrows() + 1, mat_.ia);  

    mat_.ja = new int[elements];
    std::copy(col_ptr, col_ptr + elements, mat_.ja);

    // values is brutal stuff, because binary the std::complex<double>* IS a double*!
    mat_.a = reinterpret_cast<double*>(new T[elements]);
    T* val_src = reinterpret_cast<T*>(mat_.a);
    std::copy(val_ptr, val_ptr + elements, val_src);
    
    mat_.nr = som.GetNrows();
    mat_.nc = som.GetNcols();
}


template<typename T>
void Ilupack<T>::Setup(BaseMatrix &sysMat) 
{
    // determine the matrix type. Symmetric/nonsymmetric, positive definite, ...
    // it is optional given in the xml file.
    
    // GNL, SYM, PD, HER, ....
    DetermineMatrixType(sysMat);

    // in case we already run release memory - it's save if first run
    ReleaseIlupackMemory(); 

    // set the ilupack matrix mat_ from the olas matrix - complex can be casted to double
    SetMatrix(sysMat);

    // gain and modify the default settings, does logging uses mat_ and sets param_ 
    InitParameters();

    // factorize
    nlev_ = 0; // initialized in the example, but why?
    
    // result
    int r;

    switch(matrix_)
    {
        case GNL: if(isComplex_) r = ZGNLAMGfactor(zmat_ptr_, zprecond_ptr_, &nlev_, zparam_ptr_, initial_.cplxFunc, regular_.cplxFunc, final_.cplxFunc);
                            else r = DGNLAMGfactor(dmat_ptr_, dprecond_ptr_, &nlev_, dparam_ptr_, initial_.realFunc, regular_.realFunc, final_.realFunc);
                  break;                             

        case SYM: if(isComplex_) r = ZSYMAMGfactor(zmat_ptr_, zprecond_ptr_, &nlev_, zparam_ptr_, initial_.cplxFunc, regular_.cplxFunc, final_.cplxFunc);
                            else r = DSYMAMGfactor(dmat_ptr_, dprecond_ptr_, &nlev_, dparam_ptr_, initial_.realFunc, regular_.realFunc, final_.realFunc); 
                  break;

        case PD : if(isComplex_) r = ZHPDAMGfactor(zmat_ptr_, zprecond_ptr_, &nlev_, zparam_ptr_, initial_.cplxFunc, regular_.cplxFunc, final_.cplxFunc);
                            else r = DSPDAMGfactor(dmat_ptr_, dprecond_ptr_, &nlev_, dparam_ptr_, initial_.realFunc, regular_.realFunc, final_.realFunc); 
                  break;

        case HER: if(isComplex_) r = ZHERAMGfactor(zmat_ptr_, zprecond_ptr_, &nlev_, zparam_ptr_, initial_.cplxFunc, regular_.cplxFunc, final_.cplxFunc);
                  break;
    }              
    
    if(r != 0){ 
        (*cla) << " error factorizing the ilupack preconditioner " << PrecondError(r, param_, nlev_) << std::endl;
        throw Exception("ilupack preconditioner failed", PrecondError(r, param_, nlev_), nlev_, __LINE__, __PRETTY_FUNCTION__);    
    }
    (*cla) << "factorizing got " << r << " " << (r != 0 ? PrecondError(r, param_, nlev_) : " ok ") 
           << " with " << nlev_ << " levels." << std::endl;
              
    // reset counter which holds the number of rhs's solved with this precond
    precondAge_ = 0;
    setupRuns_++;          
} 



template<typename T>
void Ilupack<T>::Solve(const BaseMatrix &base_mat, const BasePrecond &base_precond,
                             const BaseVector &base_rhs, BaseVector &base_sol, bool recursive_call )
{
    // the preconditioner sets the ilupack matrix
    if(mat_.a == NULL) throw Exception("Setup() not called before Solve()",__LINE__, __PRETTY_FUNCTION__);

    // we gain the solution pointer and compensate the first element of the OLAS one based stuff
    // Note, that we do this with the template form, such for complex, the 0 complex element is skipped.
    T* sol_ptr = dynamic_cast<Vector<T>&>(base_sol).GetPointer() + 1;
    
    // additionally remove the constness of the rhs, ilupack knows no const.
    T* rhs_ptr = const_cast<T*>(dynamic_cast<const Vector<T>&>(base_rhs).GetPointer() + 1);

    double *d_sol, *d_rhs;
    // is a ilupack type
    doublecomplex *z_sol, *z_rhs;
    
    d_sol = isComplex_ ? NULL : reinterpret_cast<double*>(sol_ptr);
    d_rhs = isComplex_ ? NULL : reinterpret_cast<double*>(rhs_ptr);
    z_sol = isComplex_ ? reinterpret_cast<doublecomplex*>(sol_ptr) : NULL;
    z_rhs = isComplex_ ? reinterpret_cast<doublecomplex*>(rhs_ptr) : NULL;

    // result
    int r;

    switch(matrix_)
    {
        case GNL: if(isComplex_) r = ZGNLAMGsolver(*zmat_ptr_, *zprecond_ptr_, nlev_, zparam_ptr_, z_sol, z_rhs);
                            else r = DGNLAMGsolver(*dmat_ptr_, *dprecond_ptr_, nlev_, dparam_ptr_, d_sol, d_rhs);
                  break;
        
        case SYM: if(isComplex_) r = ZSYMAMGsolver(*zmat_ptr_, *zprecond_ptr_, nlev_, zparam_ptr_, z_sol, z_rhs);
                            else r = DSYMAMGsolver(*dmat_ptr_, *dprecond_ptr_, nlev_, dparam_ptr_, d_sol, d_rhs);
                  break;

        case PD:  if(isComplex_) r = ZHPDAMGsolver(*zmat_ptr_, *zprecond_ptr_, nlev_, zparam_ptr_, z_sol, z_rhs);
                            else r = DSPDAMGsolver(*dmat_ptr_, *dprecond_ptr_, nlev_, dparam_ptr_, d_sol, d_rhs);
                  break;

        case HER: if(isComplex_) r = ZHERAMGsolver(*zmat_ptr_, *zprecond_ptr_, nlev_, zparam_ptr_, z_sol, z_rhs);
                  break;
    }              

    if(r != 0) 
    {
      (*cla) << " error from ilupack solver: " << SolverError(r) << std::endl;

      // try again if this is the first try but the setup() is some solvings ago
      if(recursive_call || precondAge_ == 0)
        throw Exception("error from ilupack solver", SolverError(r), __LINE__, __PRETTY_FUNCTION__, Exception::NO_SEGFAULT); 

      // try once more 
      std::cout << " -> call preconditioner and try once more" << std::endl;  
      (*cla) << " -> call preconditioner and try once more" << std::endl;
     
      Setup(const_cast<BaseMatrix&>(base_mat));

      Solve(base_mat, base_precond, base_rhs, base_sol, true);
    }    

    DumpSolverStatus(*cla);
    
    // solved another rhs with the current preconditioner
    precondAge_++;
}


template<typename T>
void Ilupack<T>::InitParameters()
{
    // output to cla or NULL
    std::ostringstream tmp;
    std::ostream& out = setupRuns_ == 0 ? *cla : tmp; 
  
    out << std::endl << "Ilupack: Available permutations for initial, regular and final reordering:" << std::endl;
    DumpPermutations(out);

    // set the permutations
    GetPermutation(INITIAL, initial_); 
    GetPermutation(REGULAR, regular_);
    GetPermutation(FINAL, final_);

    out << " Selected permutations:" << std::endl;
    out << " initial permutation: " << initial_.description << std::endl;
    out << " regular permutation: " << regular_.description << std::endl;
    out << " final   permutation: " << final_.description << std::endl;        

  
    // initializes the parameter block
    switch(matrix_)
    {
        case GNL: if(isComplex_) ZGNLAMGinit(*zmat_ptr_, zparam_ptr_);
                            else DGNLAMGinit(*dmat_ptr_, dparam_ptr_);
                  break;
        
        case SYM: if(isComplex_) ZSYMAMGinit(*zmat_ptr_, zparam_ptr_);
                            else DSYMAMGinit(*dmat_ptr_, dparam_ptr_);
                  break;

        case PD:  if(isComplex_) ZHPDAMGinit(*zmat_ptr_, zparam_ptr_);
                            else DSPDAMGinit(*dmat_ptr_, dparam_ptr_);
                  break;

        case HER: if(isComplex_) ZHERAMGinit(*zmat_ptr_, zparam_ptr_);
                            else throw Exception("ilupack matrix is set to hermitian but not complex", __LINE__, __PRETTY_FUNCTION__);
                  break;
        
        default: throw Exception("matrix type not implemented",  __LINE__, __PRETTY_FUNCTION__);          
    }              
        
    // read the default settings. I guess this simply extracts the values from magic indices in [i/f]parm
    int flag, elbow, max_it, new_int, nrestart = -1;
    double droptols[2], condest, restol, new_real;

    switch(matrix_)
    {
        case GNL: if(isComplex_) ZGNLAMGgetparams(*zparam_ptr_, &flag, &elbow, droptols, &condest, &restol, &max_it, &nrestart);
                            else DGNLAMGgetparams(*dparam_ptr_, &flag, &elbow, droptols, &condest, &restol, &max_it, &nrestart);
                  break;
        
        case SYM: if(isComplex_) ZSYMAMGgetparams(*zparam_ptr_, &flag, &elbow, droptols, &condest, &restol, &max_it);
                            else DSYMAMGgetparams(*dparam_ptr_, &flag, &elbow, droptols, &condest, &restol, &max_it);
                  break;

        case PD:  if(isComplex_) ZHPDAMGgetparams(*zparam_ptr_, &flag, &elbow, droptols, &condest, &restol, &max_it);
                            else DSPDAMGgetparams(*dparam_ptr_, &flag, &elbow, droptols, &condest, &restol, &max_it);
                  break;

        case HER: ZHERAMGgetparams(*zparam_ptr_, &flag, &elbow, droptols, &condest, &restol, &max_it);
                  break;
    }              

    out << std::endl << " Ilupack default and overwritten values: " << std::endl
                     << "======================================== " << std::endl;
    // it is magic stuff taken from the ilupack samples that we need index 1 ??:(
    out << " drop tolerance          : " << droptols[1];
    if(myParams_->HasKey("ilupack_dropTol")) {
        new_real = myParams_->GetDoubleValue("ilupack_dropTol");
        out << (new_real == droptols[1] ? " == " : " -> ") << new_real; 
        droptols[1] = new_real;
    }

    out << std::endl << " residual tolerance      : " << restol;
    if(myParams_->HasKey("ilupack_residualTol")) {
        new_real = myParams_->GetDoubleValue("ilupack_residualTol");
        out << (new_real == restol ? " == " : " -> ") << new_real; 
        restol = new_real;

    }

    out << std::endl << " elbow space             : " << elbow;
    if(myParams_->HasKey("ilupack_elbowSpace")) {
        new_int = myParams_->GetIntValue("ilupack_elbowSpace");
        out << " " << (new_int ==  elbow ? " == " : " -> ") << new_int;
        elbow = new_int;
    }

    out << std::endl << " condest (bound [L/U]^-1): " << condest;
    if(myParams_->HasKey("ilupack_condest")) {
        new_int = myParams_->GetIntValue("ilupack_condest");
        out << " " << (new_int == condest ? " ==" : "-> ") << new_int;
        condest = new_int;
    }

    out << std::endl << " max iterations          : " << max_it;
    if(myParams_->HasKey("ilupack_maxIter")) {
        new_int = myParams_->GetIntValue("ilupack_maxIter");
        out << " " << (new_int == max_it ? " == " : " -> ") << new_int;
        max_it = new_int;
    }
    
    // n restart only for GNL and up to now it cannot be set
    if(nrestart != -1) out << std::endl << " nrestart                : " << nrestart;

    // the solver is conditionally set!
    int org_solver = param_.ipar[SOLVER_TYPE];
    out << std::endl << " solver                  : " << solver.ToString(org_solver);
    if(myParams_->HasKey("ilupack_engine"))
    {
        new_int = solver.Parse(myParams_->GetStringValue("ilupack_engine"));
        out << " " << (org_solver == new_int ? " == " : " -> ") << solver.ToString(new_int);
        param_.ipar[SOLVER_TYPE] = new_int;
    }    
        
    // now the flag stuff 
    int org_flag = flag;
    // apply the settings from the xml file
    ApplyFlagSettings(flag);

    const int size = 26;

    std::map<int, std::string>::iterator iter;
    for(iter = flags.map.begin(); iter != flags.map.end(); iter++)
    {
        int f = iter->first;
        bool now = (flag & f) > 0;
        bool was = (org_flag & f) > 0;

        out << std::endl << " " << iter->second;
        out.width(size - iter->second.length());
        out << ": " << (was ? "on" : "off");
        
        // do we overwrite this setting?
        if(myParams_->HasContent(iter->second))
           out << (now == was ? " == " : " -> ") << (now ? "on" : "off");
    }
   
    out << std::endl << std::endl;

    // rewrite the updated parameters
    switch(matrix_)
    {
        case GNL: if(isComplex_) ZGNLAMGsetparams(*zmat_ptr_, zparam_ptr_, flag, elbow, droptols, condest, restol, max_it, nrestart);
                            else DGNLAMGsetparams(*dmat_ptr_, dparam_ptr_, flag, elbow, droptols, condest, restol, max_it, nrestart);
                  break;
        
        case SYM: if(isComplex_) ZSYMAMGsetparams(*zmat_ptr_, zparam_ptr_, flag, elbow, droptols, condest, restol, max_it);
                            else DSYMAMGsetparams(*dmat_ptr_, dparam_ptr_, flag, elbow, droptols, condest, restol, max_it); 
                  break;

        case PD:  if(isComplex_) ZHPDAMGsetparams(*zmat_ptr_, zparam_ptr_, flag, elbow, droptols, condest, restol, max_it);
                            else DSPDAMGsetparams(*dmat_ptr_, dparam_ptr_, flag, elbow, droptols, condest, restol, max_it);
                  break;

        case HER: ZHERAMGsetparams(*zmat_ptr_, zparam_ptr_, flag, elbow, droptols, condest, restol, max_it);
                  break;
    }        
 
}

template<typename T>
std::string Ilupack<T>::PrecondError(int result, const DILUPACKparam& param, int level)
{ 
    std::ostringstream os;

    // the text is more or less a copy & pase from the ilupack samples
    switch (result)
    {
        case  0: os << "No error";
                 break;
  
        case -1: os << "Input matrix may be wrong. ";
                 os << "(The elimination process has generated a row in L or U whose length is .gt. n";
                 break;                  

        case -2: os << "The matrix L overflows the array alu. Out of memory?";
                 break;

        case -3: os << "The matrix U overflows the array alu. Out of memory?"; 
                 break;

        case -4: os << "Illegal value for lfil"; 
                 break;

        case -5: os << "Zero row encountered"; 
                 break;
                 
        case -6: os << "Zero column encountered."; 
                 break;

        case -7: os << "Buffers too small. Check AMGSetup for nibuff and ndbuff.";
                 os << " Increase buffers to " << param.nibuff << " (integer)";
                 os << " and " << param.ndbuff << " (double)";
                 break;
                 
        default: os << "Zero pivot encountered at step number " << result << ".";
    }
    os << " (level=" << level <<")";
    return os.str();   
}

template<typename T>
std::string Ilupack<T>::SolverError(int result)
{ 
    std::ostringstream os;

    // the text is more or less a copy & pase from the ilupack samples
    switch (result)
    {
        case  0: os << "No error";
                 break;
  
        case -1: os << "Number of iteration steps exceeds its limit";
                 break;
                 
        case -2: os << "Not enough work space provided";
                 break;
                 
        case -3: os << "Algorithm breaks down";
                 break;

        default: os << "Solver exited with unspecified error code " << result;
    }
 
    return os.str();
}

template<typename T>
void Ilupack<T>::DetermineMatrixType(BaseMatrix &sysMat)
{ 
    // first determine it manually for some checking
    if(sysMat.GetStructureType() != STDMATRIX) throw Exception("Sorry, excpect StdMatrix", sysMat.GetStructureType(), __LINE__, __PRETTY_FUNCTION__); 

    const StdMatrix& stdMat = dynamic_cast<const StdMatrix&>(sysMat);
    MatrixStorageType mst = stdMat.GetStorageType(); 
    if(mst != SPARSE_SYM && mst != SPARSE_NONSYM) throw Exception("Sorry, expect sparse matrix", mst, __LINE__, __PRETTY_FUNCTION__);   
    
    if(myParams_->HasKey("ilupack_matrix")) {
        matrix_ = static_cast<Matrix>(matrix.Parse(myParams_->GetStringValue("ilupack_matrix")));
        // plausibility check -- killme: what is with hermitian?
        if(mst != SPARSE_SYM && matrix_ != GNL) 
            throw Exception("Matrix storrage is unsymmetric, so given ilupack_matrix is invalid", matrix.ToString(matrix_));
    } else {
        // ignore PD and HER
        matrix_ = mst == SPARSE_SYM ? SYM : GNL;
    }
    
    if(setupRuns_ == 0) (*cla) << " Use Ilupack matrix type " << matrix.ToString(matrix_) << std::endl;
}    
    


template<typename T>
void Ilupack<T>::SetEnums()
{
   solver.SetName("Ilupack:Solver");
   solver.Add(PCG,    "pcg");
   solver.Add(SBCG,   "sbcg");
   solver.Add(BCG,    "bcg");
   solver.Add(SQMR,   "sqmr");
   solver.Add(BCGSTAB,"bcgstab");
   solver.Add(TFQMR,  "tfqmr");
   solver.Add(FOM,    "fom");
   solver.Add(GMRES,  "gmres");
   solver.Add(FGMRES, "fgmres");
   solver.Add(DQGMRES,"dqgmres");

   matrix.SetName("Ilupack::Matrix");
   matrix.Add(GNL, "gnl");
   matrix.Add(SYM, "sym");
   matrix.Add(PD,  "pd");
   matrix.Add(HER, "her");   
   
   ordering.SetName("Ilupack::Ordering");
   ordering.Add(NONE, "null");
   ordering.Add(ND,"nd");
   ordering.Add(RCM,"rcm");
   ordering.Add(AMF,"amf");
   ordering.Add(AMD,"amd");
   ordering.Add(MMD,"mmd");
   ordering.Add(METIS_E,"metisE");
   ordering.Add(METIS_N,"metisN");
   ordering.Add(INDSET,"indset");
   ordering.Add(PP,"pp");
   ordering.Add(PQ,"pq");
   ordering.Add(FC,"fc");

   matching.SetName("Ilupack::Matching");
   matching.Add(PURE, "none");
   matching.Add(MC64, "mc64");
   matching.Add(MWM, "mwm"); 

   flags.SetName("Ilupack::Flags");
   flags.Add(FL_DROP_INVERSE, "DropInverse");
   flags.Add(FL_NO_SHIFT, "NoShift");
   flags.Add(FL_TISMENETSKY_SC,"TismenetskySC");
   flags.Add(FL_REPEAT_FACT,"RepeatFact");      
   flags.Add(FL_IMPROVED_ESTIMATE,"ImprovedEstimate");
   flags.Add(FL_DIAGONAL_COMPENSATION,"DiagonalCompensation");
   flags.Add(FL_COARSE_REDUCE,"CoarseReduce");      
   flags.Add(FL_FINAL_PIVOTING,"FinalPivoting");
   flags.Add(FL_ENSURE_SPD,"EnsureSPD");
   flags.Add(FL_SIMPLE_SC,"SimpleSC");      
   flags.Add(FL_PREPROCESS_INITIAL_SYSTEM,"PreprocessInitialSystem");
   flags.Add(FL_PREPROCESS_SUBSYSTEMS,"PreprocessSubsystem");
   flags.Add(FL_MULTI_PILUC,"MultiPiluc");      
   flags.Add(FL_RE_FACTOR,"ReFactor");
   flags.Add(FL_AGGRESSIVE_DROPPING,"AggressiveDropping");
   flags.Add(FL_DISCARD_MATRIX,"DiscardMatrix");
   flags.Add(FL_SYMMETRIC_STRUCTURE, "SymmetricStructure");
   
   permutationRole.SetName("Ilupack::PermuationRole");
   permutationRole.Add(INITIAL, "initial");
   permutationRole.Add(REGULAR, "regular");
   permutationRole.Add(FINAL,   "final");          
}


template<typename T>
void Ilupack<T>::AddPermutation(Matrix matrix, Ordering ordering, 
  Matching matching,RealPermFunc realFunc, ComplexPermFunc cplxFunc, 
  const std::string& ilupack_name, const std::string& code, const std::string& description)
{
    Permutation perm;
    perm.matrix = matrix;
    perm.ordering = ordering;
    perm.matching = matching;
    perm.realFunc = realFunc;
    perm.cplxFunc = cplxFunc;
    perm.ilupack_name = ilupack_name;
    perm.code    = code;
    perm.description = description;
    permutations.push_back(perm); 
}  
  
template<typename T>
void Ilupack<T>::SetPermutations()
{
   Permutation perm;

   // ------------------------------------- 
   // GNL plain permutaions - no extra libs
   AddPermutation(GNL, NONE, PURE, DGNLperm_null, ZGNLperm_null, "GNLperm_null", "", "no permutation");
   AddPermutation(GNL, ND, PURE, DGNLperm_nd, ZGNLperm_nd, "GNLperm_nd", "nd", "scaling + nested dissection ordering");  
   AddPermutation(GNL, RCM, PURE, DGNLperm_rcm, ZGNLperm_rcm, "GNLperm_rcm", "rcm", "scaling + reverse Cuthill-McKee ordering");    
   AddPermutation(GNL, MMD, PURE, DGNLperm_mmd, ZGNLperm_mmd, "GNLperm_mmd", "mmd", "scaling + minimum degree ordering");    
   AddPermutation(GNL, AMF, PURE, DGNLperm_amf, ZGNLperm_amf, "GNLperm_amf", "amf", "scaling + AMF ordering");
   AddPermutation(GNL, AMD, PURE, DGNLperm_amd, ZGNLperm_amd, "GNLperm_amd", "amd", "scaling + approximate minimum degree ordering");
   AddPermutation(GNL, PQ, PURE, DGNLperm_pq, ZGNLperm_pq, "GNLperm_pq", "PQs", "scaling + ddPQ");    
   AddPermutation(GNL, FC, PURE, DGNLperm_fc, ZGNLperm_fc, "GNLperm_fc", "FCs", "scaling + fine grid/coarse grid partitioning");
   AddPermutation(GNL, INDSET, PURE, DGNLperm_indset, ZGNLperm_indset, "GNLperm_indset", "inds", "scaling + independent set ordering");
   // GNL with METIS package 
#ifdef USE_METIS
   AddPermutation(GNL, METIS_E, PURE, DGNLperm_metis_e, ZGNLperm_metis_e, "GNLperm_metis_e", "mes", "scaling + MeTiS Edge nested dissection ordering");    
   AddPermutation(GNL, METIS_N, PURE, DGNLperm_metis_n, ZGNLperm_metis_n, "GNLperm_metis_n", "mns", "scaling + MeTiS Node nested dissection ordering");
#endif
   // GNL permuations based on PARDISO
#ifdef USE_PARDISO
   AddPermutation(GNL, RCM, MWM, DGNLperm_mwm_rcm, ZGNLperm_mwm_rcm, "GNLperm_mwm_rcm", "mwrc", "scaling + PARDISO max. weight matching + reverse Cuthill-McKee ordering");    
   AddPermutation(GNL, MMD, MWM, DGNLperm_mwm_mmd, ZGNLperm_mwm_mmd, "GNLperm_mwm_mmd", "mwad", "scaling + PARDISO max. weight matching + minimum degree ordering");
   AddPermutation(GNL, AMF, MWM, DGNLperm_mwm_amf, ZGNLperm_mwm_amf, "GNLperm_mwm_amf", "mwaf", "scaling + PARDISO max. weight matching + AMF ordering");
   AddPermutation(GNL, AMD, MWM, DGNLperm_mwm_amd, ZGNLperm_mwm_amd, "GNLperm_mwm_amd", "mwad", "scaling + PARDISO max. weight matching + approximate minimum degree ordering");    
   AddPermutation(GNL, METIS_E, MWM, DGNLperm_mwm_metis_e, ZGNLperm_mwm_metis_e, "GNLperm_mwm_metis_e", "mwme", "scaling + PARDISO max. weight matching + MeTiS Edge ND ordering");
   AddPermutation(GNL, METIS_N, MWM, DGNLperm_mwm_metis_n, ZGNLperm_mwm_metis_n, "GNLperm_mwm_metis_n", "mwmn", "scaling + PARDISO max. weight matching + MeTiS Node ND ordering");
#endif // USE_PARDISO
   // GNL permuations based on MC64
#ifdef USE_MC64
   AddPermutation(GNL, RCM, MC64, DGNLperm_mc64_rcm, ZGNLperm_mc64_rcm, "GNLperm_mc64_rcm", "mwrc", "scaling + MC64 max. weight matching + reverse Cuthill-McKee ordering");    
   AddPermutation(GNL, MMD, MC64, DGNLperm_mc64_mmd, ZGNLperm_mc64_mmd, "GNLperm_mc64_mmd", "mwad", "scaling + MC64 max. weight matching + minimum degree ordering");
   AddPermutation(GNL, AMF, MC64, DGNLperm_mc64_amf, ZGNLperm_mc64_amf, "GNLperm_mc64_amf", "mwaf", "scaling + MC64 max. weight matching + AMF ordering");
   AddPermutation(GNL, AMD, MC64, DGNLperm_mc64_amd, ZGNLperm_mc64_amd, "GNLperm_mc64_amd", "mwad", "scaling + MC64 max. weight matching + approximate minimum degree ordering");    
   AddPermutation(GNL, METIS_E, MC64, DGNLperm_mc64_metis_e, ZGNLperm_mc64_metis_e, "GNLperm_mc64_metis_e", "mwme", "scaling + MC64 max. weight matching + MeTiS Edge ND ordering");
   AddPermutation(GNL, METIS_N, MC64, DGNLperm_mc64_metis_n, ZGNLperm_mc64_metis_n, "GNLperm_mc64_metis_n", "mwmn", "scaling + MC64 max. weight matching + MeTiS Node ND ordering");
#endif // USE_MC64

   // ------------------------------------- 
   // SYM plain permutaions - no extra libs
   AddPermutation(SYM, FC, PURE, DSYMperm_fc, ZSYMperm_fc, "SYMperm_fc", "FCs", "scaling + fine grid/coarse grid partitioning");
   // SYM permuations based on PARDISO
#ifdef USE_PARDISO
   AddPermutation(SYM, RCM, MWM, DSYMperm_mwm_rcm, ZSYMperm_mwm_rcm, "SYMperm_mwm_rcm", "mwrc", "scaling + PARDISO sym. max. weight matching + reverse Cuthill-McKee ordering");    
   AddPermutation(SYM, MMD, MWM, DSYMperm_mwm_mmd, ZSYMperm_mwm_mmd, "SYMperm_mwm_mmd", "mwad", "scaling + PARDISO sym. max. weight matching + minimum degree ordering");
   AddPermutation(SYM, AMF, MWM, DSYMperm_mwm_amf, ZSYMperm_mwm_amf, "SYMperm_mwm_amf", "mwaf", "scaling + PARDISO sym. max. weight matching + AMF ordering");
   AddPermutation(SYM, AMD, MWM, DSYMperm_mwm_amd, ZSYMperm_mwm_amd, "SYMperm_mwm_amd", "mwad", "scaling + PARDISO sym. max. weight matching + approximate minimum degree ordering");    
   AddPermutation(SYM, METIS_E, MWM, DSYMperm_mwm_metis_e, ZSYMperm_mwm_metis_e, "SYMperm_mwm_metis_e", "mwme", "scaling + PARDISO sym. max. weight matching + MeTiS Edge ND ordering");
   AddPermutation(SYM, METIS_N, MWM, DSYMperm_mwm_metis_n, ZSYMperm_mwm_metis_n, "SYMperm_mwm_metis_n", "mwmn", "scaling + PARDISO sym. max. weight matching + MeTiS Node ND ordering");
#endif // USE_PARDISO
   // SYM permuations based on MC64
#ifdef USE_MC64
   AddPermutation(SYM, RCM, MC64, DSYMperm_mc64_rcm, ZSYMperm_mc64_rcm, "SYMperm_mc64_rcm", "mwrc", "scaling + MC64 sym. max. weight matching + reverse Cuthill-McKee ordering");    
   AddPermutation(SYM, MMD, MC64, DSYMperm_mc64_mmd, ZSYMperm_mc64_mmd, "SYMperm_mc64_mmd", "mwad", "scaling + MC64 sym. max. weight matching + minimum degree ordering");
   AddPermutation(SYM, AMF, MC64, DSYMperm_mc64_amf, ZSYMperm_mc64_amf, "SYMperm_mc64_amf", "mwaf", "scaling + MC64 sym. max. weight matching + AMF ordering");
   AddPermutation(SYM, AMD, MC64, DSYMperm_mc64_amd, ZSYMperm_mc64_amd, "SYMperm_mc64_amd", "mwad", "scaling + MC64 sym. max. weight matching + approximate minimum degree ordering");    
   AddPermutation(SYM, METIS_E, MC64, DSYMperm_mc64_metis_e, ZSYMperm_mc64_metis_e, "SYMperm_mc64_metis_e", "mwme", "scaling + MC64 sym. max. weight matching + MeTiS Edge ND ordering");
   AddPermutation(SYM, METIS_N, MC64, DSYMperm_mc64_metis_n, ZSYMperm_mc64_metis_n, "SYMperm_mc64_metis_n", "mwmn", "scaling + MC64 sym. max. weight matching + MeTiS Node ND ordering");
#endif // USE_MC64

   // ------------------------------------- 
   // Positive Definite plain permutaions - no extra libs
   AddPermutation(PD, PP, PURE, DSPDperm_pp, ZHPDperm_pp, "[S/H]PDperm_pp", "PPs", "scaling + symmetrized version of PQ");

   // ------------------------------------- 
   // Hermitian permutaions based on PARDISO
#ifdef USE_PARDISO
   AddPermutation(HER, RCM, MWM, NULL, ZHERperm_mwm_rcm, "HERperm_mwm_rcm", "mwrc", "scaling + PARDISO sym. max. weight matching + reverse Cuthill-McKee ordering");    
   AddPermutation(HER, MMD, MWM, NULL, ZHERperm_mwm_mmd, "HERperm_mwm_mmd", "mwad", "scaling + PARDISO sym. max. weight matching + minimum degree ordering");
   AddPermutation(HER, AMF, MWM, NULL, ZHERperm_mwm_amf, "HERperm_mwm_amf", "mwaf", "scaling + PARDISO sym. max. weight matching + AMF ordering");
   AddPermutation(HER, AMD, MWM, NULL, ZHERperm_mwm_amd, "HERperm_mwm_amd", "mwad", "scaling + PARDISO sym. max. weight matching + approximate minimum degree ordering");    
   AddPermutation(HER, METIS_E, MWM, NULL, ZHERperm_mwm_metis_e, "HERperm_mwm_metis_e", "mwme", "scaling + PARDISO sym. max. weight matching + MeTiS Edge ND ordering");
   AddPermutation(HER, METIS_N, MWM, NULL, ZHERperm_mwm_metis_n, "HERperm_mwm_metis_n", "mwmn", "scaling + PARDISO sym. max. weight matching + MeTiS Node ND ordering");
#endif // USE_PARDISO
   // Hermitian permuations based on MC64
#ifdef USE_MC64
   AddPermutation(HER, RCM, MC64, NULL, ZHERperm_mc64_rcm, "HERperm_mc64_rcm", "mwrc", "scaling + MC64 sym. max. weight matching + reverse Cuthill-McKee ordering");    
   AddPermutation(HER, MMD, MC64, NULL, ZHERperm_mc64_mmd, "HERperm_mc64_mmd", "mwad", "scaling + MC64 sym. max. weight matching + minimum degree ordering");
   AddPermutation(HER, AMF, MC64, NULL, ZHERperm_mc64_amf, "HERperm_mc64_amf", "mwaf", "scaling + MC64 sym. max. weight matching + AMF ordering");
   AddPermutation(HER, AMD, MC64, NULL, ZHERperm_mc64_amd, "HERperm_mc64_amd", "mwad", "scaling + MC64 sym. max. weight matching + approximate minimum degree ordering");    
   AddPermutation(HER, METIS_E, MC64, NULL, ZHERperm_mc64_metis_e, "HERperm_mc64_metis_e", "mwme", "scaling + MC64 sym. max. weight matching + MeTiS Edge ND ordering");
   AddPermutation(HER, METIS_N, MC64, NULL, ZHERperm_mc64_metis_n, "HERperm_mc64_metis_n", "mwmn", "scaling + MC64 sym. max. weight matching + MeTiS Node ND ordering");
#endif // USE_MC64

//   AddPermutation(, , , , , "", "", "");
}


template<typename T>
void Ilupack<T>::DumpParameters(std::ostream& out)
{
   out << "Ilupack parameters" << std::endl;
   out << "nalu   = " << param_.nalu << std::endl;
   out << "ndaux  = " << param_.ndaux << std::endl;
   out << "ndbuff = " << param_.ndbuff << std::endl;
   out << "niaux  = " << param_.niaux << std::endl;
   out << "nibuff = " << param_.nibuff << std::endl;
   out << "njlu   = " << param_.njlu << std::endl;
   out << "nju    = " << param_.nju << std::endl;
   for(int i = 0; i < ILUPACK_NFPAR; i++)
      out << "fpar[" << i << "] = " << param_.fpar[i] << std::endl;
   for(int i = 0; i < ILUPACK_NIPAR; i++)
      out << "ipar[" << i << "] = " << param_.fpar[i] << std::endl;
}



template<typename T>
void Ilupack<T>::DumpPermutations(std::ostream& out)
{
    out << "Matrix\tOrder\tMatch\tName\tCode\tDescripotion" << std::endl;
    
    for(UInt i = 0; i < permutations.size(); i++)
    {
        Permutation& p = permutations[i];
    
        out << matrix.ToString(p.matrix) << "\t";
        out << ordering.ToString(p.ordering) << "\t";
        out << matching.ToString(p.matching) << "\t";
        out << p.ilupack_name << "\t";
        out << p.code << "\t";
        out << p.description << std::endl;      
    }
       
}


template<typename T>
void Ilupack<T>::DumpSolverStatus(std::ostream& out)
{
    out << "Ilupack solver status:" << std::endl;
    out << " Number of Iterations    : " << param_.ipar[ITERATION_STEPS] << std::endl;
    out << " initial residual norm   : " << param_.fpar[INITIAL_RESIDUAL_NORM] << std::endl;
    out << " target residual norm    : " << param_.fpar[TARGET_RESIDUAL_NORM] << std::endl; 
    out << " current residual norm   : " << param_.fpar[CURRENT_RESIDUAL_NORM] << std::endl;
    out << " convergence rate        : " << param_.fpar[CONVERGENCE_RATE] << std::endl;
    out << " schur drop tolerance    : " << param_.fpar[SCHUR_DROP_TOL] << std::endl;
    out << " relative error tolerance: " << param_.fpar[REL_ERROR_TOL] << std::endl;
    out << " absolute error tolerance: " << param_.fpar[ABS_ERROR_TOL] << std::endl;        
}    


template<typename T>
void Ilupack<T>::DumpFlags(int flag, std::ostream& out)
{
    out << "State\tFlag" << std::endl;
    
    std::map<int, std::string>::iterator iter;
    for(iter = flags.map.begin(); iter != flags.map.end(); iter++)
    {
        int f = iter->first;
        out << ((flag & f) != 0 ? "true" : "false");
        out << "\t" << iter->second;
        out << std::endl;
    }
}

template<typename T>
void Ilupack<T>::ApplyFlagSettings(int& flag)
{
    // the values from xml are stored in ilupack_flag_[on/off]_* strings
    // with the name of the flag as value and * is 0 or greater compact ordered
    std::ostringstream os;
    std::string on_key;
    std::string off_key;

    bool cont = true;

    for(int i = 0; cont; i++)
    {
        os.str("");
        os << "ilupack_flag_on_" << i;
        on_key = os.str();

        os.str("");
        os << "ilupack_flag_off_" << i;
        off_key = os.str();

        cont = false;
        
        if(myParams_->HasKey(on_key))
        {
            cont = true;
            int value = flags.Parse(myParams_->GetStringValue(on_key));
            flag |= value; // or this value in
        }

        if(myParams_->HasKey(off_key))
        {
            cont = true;
            int value = flags.Parse(myParams_->GetStringValue(off_key));
            flag &= ~value; // AND 111110111 
        }
    }
}

template<typename T>
void Ilupack<T>::GetDefaultPermutation(PermutationRole role, Permutation& out)
{
   Matrix m = matrix_;
#ifdef USE_PARDISO
   bool   pard = true;
#else
   bool   pard = false;
#endif

   switch(role)
   {
      case INITIAL:
           // initial and regular with the same settings  
         
      case REGULAR:
           if(m == GNL) { GetPermutation(GNL, MMD, pard ? MWM : PURE, out); return; }
                   else { GetPermutation(SYM, MMD, pard ? MWM : PURE, out); return; }        
           break;        
      
      case FINAL:
           if(m == GNL) { GetPermutation(GNL, PQ, PURE, out); return; }
                   else { GetPermutation(PD, PP, PURE, out); return; }
           break;
                   
      default: throw Exception("PermutationRole case not impemented. Check las-file for available set.", role, __LINE__, __PRETTY_FUNCTION__);
    }
}

template<typename T>
void Ilupack<T>::GetPermutation(Matrix my_matrix, Ordering my_ordering, Matching my_matching, Permutation& out)
{
    // search for the permutation
    for(UInt i = 0; i < permutations.size(); i++)
    {
        Permutation& p = permutations[i]; 
        if(p.matrix   != my_matrix)   continue;
        if(p.ordering != my_ordering) continue;
        if(p.matching != my_matching) continue;        
        
        out = p;
        return; // found!
    }
  
    throw Exception("Permutation (matrix, ordering, matching) not found.", matrix.ToString(my_matrix), ordering.ToString(my_ordering), matching.ToString(my_matching)); 
}

template<typename T>
void Ilupack<T>::GetPermutation(PermutationRole role, Permutation& out)
{
    // to check if we have the value in XML -> therefore use dummy defaults for my_*
    bool in_xml = true; 
    std::ostringstream os;
    std::string type = permutationRole.ToString(role);
    
    os << "ilupack_" << type << "_matrix";
    Matrix my_matrix = GNL;
    if(!myParams_->HasKey(os.str())) in_xml = false;
    else my_matrix = static_cast<Matrix>(matrix.Parse(myParams_->GetStringValue(os.str())));
    
    os.str("");
    os << "ilupack_" << type << "_ordering";
    Ordering my_ordering = NONE;
    if(!myParams_->HasKey(os.str())) in_xml = false;
    else my_ordering = static_cast<Ordering>(ordering.Parse(myParams_->GetStringValue(os.str())));
    
    os.str("");
    os << "ilupack_" << type << "_matching";
    Matching my_matching = PURE;
    if(!myParams_->HasKey(os.str())) in_xml = false;
    else my_matching = static_cast<Matching>(matching.Parse(myParams_->GetStringValue(os.str())));

    // is there stuff in XML? only then the my_* is valid
    if(in_xml) GetPermutation(my_matrix, my_ordering, my_matching, out); 
          else GetDefaultPermutation(role, out);
}


} // end of namespace
