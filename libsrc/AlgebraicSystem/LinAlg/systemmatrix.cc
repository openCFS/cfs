#include <fstream>
#include <iostream>
#include <string>

//#include <general_head.hh>
//#include <utils_head.hh>
#include "systemmatrix.hh"

namespace CoupledField
{

template<class T_Matrix>
void SystemMatrix<T_Matrix>::printAb(std::ostream * out, const std::string & title) const
{
 (*out) << " ------- " << title << "------------- " << std::endl;
 (*out) << " ------- System Matrix ------------ " << std::endl;
 (*out) << A;
 (*out) << " ------- Right hand side ----------- " << std::endl;
 (*out) << b;
}

template<class T_Matrix>
void SystemMatrix<T_Matrix>::printx(std::ostream * out, const Double time) const
{
 (*out) << " ------- Solution -----   Step:" << time << std::endl;
 (*out) << x;
}

}
