#include "pdememento.hh"

namespace CoupledField{

  PDEMemento::PDEMemento()
  {
    ENTER_IFCN( "PDEMemento::PDEMemento");
  
    sol_       = NULL;
    isSet_ = FALSE;
  }

  PDEMemento::~PDEMemento()
  {
    ENTER_IFCN( "PDEMemento::~PDEMemento");
  
    if (sol_)
      delete sol_;

    sol_       = NULL;
  }

  void PDEMemento::Clear()
  {
    ENTER_FCN( "PDEMemento::Clear");
    if (sol_)
      delete sol_;

    sol_       = NULL;
    solDeriv1_.Clear();
    solDeriv2_.Clear();
  }

  std::ostream & operator << ( std::ostream & out, const PDEMemento & mem)
  {
    ENTER_FCN( "operator <<(PDEMemento)" );

    UInt Size;

    Vector<Double>solOut_ = dynamic_cast<Vector<Double>&>(*mem.sol_);

    if(solOut_.GetSize() != mem.solDeriv1_.GetSize() || solOut_.GetSize() != mem.solDeriv2_.GetSize()) {
      Error( "The sizes of sol, solDeriv1 and solDeriv2 are not coincident!", __FILE__, __LINE__  );
      Size=0;
    }
    else {
      Size=solOut_.GetSize();
    }

    for (UInt i=0; i < Size; i++) {
      out << solOut_[i] << " " << mem.solDeriv1_[i] << " " << mem.solDeriv2_[i] << " " << std::endl;
    }

    return out;
  }

  std::ifstream & operator >> ( std::ifstream & in, PDEMemento & mem )
  {
    ENTER_FCN( "operator >>(PDEMemento)" );

    // Buffer to hold the input of one line
    char buffer[64];
    Vector<Double> solIn_;

    solIn_.Resize(mem.size_);
    mem.solDeriv1_.Resize(mem.size_);
    mem.solDeriv2_.Resize(mem.size_);

    // set reading position to the beginning of the file
    in.seekg(0, std::ios_base::beg);

    // not very elegant method to start reading in line number 2
    in.getline( buffer, 64, '\n' );

    for (UInt i=0; i < mem.size_; i++) {
      if ( !in.eof() ) {
        in.getline( buffer, 64, '\n' );
        if ( in.eof() ) {
          Error( "ReadLine: Unexpected end of file! Maybe your resatrt file is not correct",
                 __FILE__, __LINE__  ); 
        }
      }
      
      std::sscanf(buffer,"%lf%lf%lf", &solIn_[i], &mem.solDeriv1_[i], &mem.solDeriv2_[i] );

    }

    mem.sol_ = new Vector<Double>(solIn_);
    mem.isSet_ = TRUE;

    return in;
  }
} // namespace
