#include <fstream.h>
#include <iostream.h>
#include <string>

#include <general_head.hh>
#include <utils_head.hh>
#include "systemmatrix.hh"

namespace CoupledField
{

template<class T_Matrix>
void SystemMatrix<T_Matrix>::printAb(ostream * out, const string & title) const
{
 (*out) << " ------- " << title << "------------- " << endl;
 (*out) << " ------- System Matrix ------------ " << endl;
 (*out) << A;
 (*out) << " ------- Right hand side ----------- " << endl;
 (*out) << b;
}

template<class T_Matrix>
void SystemMatrix<T_Matrix>::printx(ostream * out, const Double time) const
{
 (*out) << " ------- Solution -----   Step:" << time << endl;
 (*out) << x;
}

}
