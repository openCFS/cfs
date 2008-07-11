#ifndef FILE_GETDOCUMENTATION_2008
#define FILE_GETDOCUMENTATION_2008

#include <map>
#include <vector>
#include <string>

namespace CoupledField
{
  void GetDocumentation(std::map<std::string, std::vector<unsigned char> >& doc);

  // XML data could also be saved in C++ like this:
  //
  //  #define XML_DATA( data, text )  static char data[] = {  " "#text" " };
  //
  //  int main(int argc, char** argv)
  //  {
  //    XML_DATA(mystr,<myxml myattr="test">\n
  //             </myxml>);
  //    std::string str = mystr;
  //    
  //    std::cout << str << std::endl;
  //    return 0;
  //  }

} // end of namespace

#endif
