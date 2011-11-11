#include <vector> // no CFS StdVector :(
#include <fstream>

#include "CRS_Matrix.hh"
#include "SCRS_Matrix.hh"
#include "StdMatrix.hh"
#include "Vector.hh"

#include "General/Exception.hh"


using CoupledField::Exception;

namespace CoupledField 
{

template<typename T>
StdMatrix::HarwellBoeing<T>::HarwellBoeing(const StdMatrix* matrix)
{
    this->matrix_ = matrix;
    if(matrix_->GetEntryType() != DOUBLE && matrix_->GetEntryType() != COMPLEX)
        throw Exception("unimplemented entry type");
        
    is_complex_ = matrix->GetEntryType() == COMPLEX;        
}


/** reads the data and writes it to the vector, use only one of the vectors.
 * uses the STL vector as the OLAS StdVector is not templated and doesn't know strings :( 
 * @param doubles access from the 0th element
 * @param ints access from the 0th element 
 * @param length the processed length of data
 * @param the Harwell-Boeing values are written there. */
template<typename T>
void StdMatrix::HarwellBoeing<T>::Fill(const T* values, const UInt* ints, UInt length, std::vector<std::string>& out)
{
    // No label                                                                No key  
    //             77             4            21            41             7
    // RSA                       26            26           163
    // (8I10)          (8I10)          (4E20.13)           (4D20.13)                   
    // F                          1             0                                      
    //          1         9        16        28        39        43        46        54
    //...
    //        161       163       164
    // 1.5954560000000E+11-3.8146972656250E-06 6.1364000000000E+09 0.0000000000000E+00
    // -4.2954600000000E+10-1.3806825000000E+10-3.9886400000000E+10-2.3011375000000E+10
    // ...
    std::ostringstream os;

    // define the settings 
    int width;
    if(ints == NULL) {
        os.setf(std::ios_base::scientific  | std::ios_base::uppercase | std::ios_base::showpoint);
        os.precision(13); // 12 after the dot and one before
        width = 20;
    } else {
        width = 10;
    }

    for(UInt i = 0; i < length; i++)
    {
        // use temporary width to handle comples
        int tw = width; 

        os.width(width);
        if(ints != NULL) os << (ints[i] + 1);
        if(values != NULL && !is_complex_) os << values[i];
        if(values != NULL && is_complex_){
           const Complex* cplx_ptr = reinterpret_cast<const Complex*>(values);
           Complex cplx = cplx_ptr[i];
           os.width(width);
           os << cplx.imag();
           os.width(width);
           os << cplx.real();
           tw = 2 * width;
        }    

        // add line if we cannot add more or are the last element
        if(os.str().length() + tw > 80 || i == length - 1) {
            // fill the line
            os.width(80 - os.str().length());
            os << std::endl;
            
            out.push_back(os.str());
            os.clear();
            os.str("");
        } 
    }
}

template<typename T>
void StdMatrix::HarwellBoeing<T>::Export(const std::string& file, const BaseVector& rhs)
{
   std::vector<std::string> cols;
   std::vector<std::string> rows;
   std::vector<std::string> mat_data;
   std::vector<std::string> rhs_data;      
    
   // we use the integer format (8I10):          1         9        16 ..
   // and the double format (4D20.13) : -1.9073486328125E-06 6.1364000000000E+09-3.8146972656250E-06
   
   // code of system 
   std::string code;
   // the number of elements can be tricky :(
   UInt elements; 

   switch(matrix_->GetStorageType())
   {
      case SPARSE_SYM:
        {
           code = is_complex_ ? "CSA" : "RSA"; 
           const SCRS_Matrix<T>* scrs = dynamic_cast<const SCRS_Matrix<T>*>(matrix_);
           
           // CFS stores the right upper rows, we need the lower columns and we are
           // symmetric -> it's simple babe :)
           Fill(NULL, scrs->GetRowPointer(),scrs->GetNumRows()+1, cols);  // include tail
           
           // Nnz is the total number = Diagonal + 2 * triagonals. But one is skipped!
           elements= scrs->GetNnz() - ((scrs->GetNnz() - scrs->GetNumRows())/2);
                      
           // Nnz is for complete, due to symmetry only half
           Fill(NULL, scrs->GetColPointer(),elements, rows);    
           Fill(scrs->GetDataPointer(),NULL,elements, mat_data);

           break;
        }

     case SPARSE_NONSYM:
        { 
           code = is_complex_ ? "CUA" : "RUA";           
           const CRS_Matrix<T>* crs = dynamic_cast<const CRS_Matrix<T>*>(matrix_);  
           // from compresses row storage transform to column storage (HB)

           elements = crs->GetNnz();
           
           unsigned int* col_ptr = new unsigned int[crs->GetNumCols() + 1];
           unsigned int* row_ptr = new unsigned int[elements];
           T*            val_ptr = new T[elements];
              
           crs->Transpose(col_ptr, row_ptr, val_ptr);
 
           Fill(NULL, col_ptr, crs->GetNumCols() + 1, cols);  // include tail
           Fill(NULL, row_ptr, elements, rows);
           Fill(val_ptr, NULL , elements, mat_data);

           delete[] col_ptr;
           delete[] row_ptr;
           delete[] val_ptr;

           break;
        }
      default: throw "uimplementd storage type in ExportHarwellBoeing()";
   }

   // rhs is common 
   const Vector<T>& vector = dynamic_cast<const Vector<T>&>(rhs);
   Fill(vector.GetPointer(),NULL,vector.GetSize(), rhs_data);

   // construct the 5 lines 80 columns header
   std::ofstream out(file.c_str());

   // the first line is constat
   //   ->|0        1         2         3         4         5         6         7         |<-
   out << "No label                                                                No key  " << std::endl;

   // second line: total lines excl. header + col lines + row lines + mat data lines + rhs lines
   out.width(14); 
   out << (cols.size() + rows.size() + mat_data.size() + rhs_data.size());
   out.width(14);
   out << cols.size();
   out.width(14);
   out << rows.size();
   out.width(14);
   out << mat_data.size();
   out.width(14);
   out << rhs_data.size();
   out << std::endl;

   // third line: matrix code, number or rows (variables), number of columns (elements), number of row indices, 0 (we are assembled)      
   out << code + "           ";
   out.width(14);
   out << matrix_->GetNumRows();        
   out.width(14);
   out << matrix_->GetNumCols();
   out.width(14);
   out << elements; 
   out <<  std::endl;
  
   // forth line is format of the data blocks, hardcoded
   out << "(8I10)          (8I10)          (4E20.13)           (4E20.13)" << std::endl;                   
   
   // fith line is hard coded, we are not sparse (F), have only 1 rhs, and are assembled (0 is ignored) 
   out << "F                          1             0" << std::endl;                                      
  
   // write the data blocks
   for(UInt i = 0; i < cols.size(); i++) out << cols[i];        
   for(UInt i = 0; i < rows.size(); i++) out << rows[i];
   for(UInt i = 0; i < mat_data.size(); i++) out << mat_data[i];
   for(UInt i = 0; i < rhs_data.size(); i++) out << rhs_data[i];      
  
   // finish!
   out.flush();
   out.close();
}

void StdMatrix::ExportHarwellBoeing(const std::string& file, const BaseVector& rhs)
{
    if(GetEntryType() == DOUBLE)  {
        HarwellBoeing<Double> hb(this);
        hb.Export(file, rhs);
    } else {
        HarwellBoeing<Complex> hb(this);
        hb.Export(file, rhs);
    }          
    
}


} // namespace
