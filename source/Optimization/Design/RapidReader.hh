#ifndef RAPIDXML_HH_
#define RAPIDXML_HH_


#include <string>
#include <vector>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "rapidxml.hpp"

namespace CoupledField {


/** This class reads in a density file using rapidxml (http://rapidxml.sourceforge.net). 
 * <p>This is only necessary because Xerces leaks resources and for large 3D problems, this can easily
 * add several hundred MB in Ram usage during the whole optimization process.</p>
 * <p>Unfortunately, while rapidxml uses way less resources and also frees the memory correctly after
 * reading the xml-file, it does not provide schema checking, so we cannot use this as a full replacement
 * for Xerces...</p>
 * @see PtrParamNode*/
class RapidReader 
{
  public:
    /** ctor, parses the xml-file into the buffer-vector */
    RapidReader(const std::string& filename);
    /** dtor, does not need to do anything atm */
    ~RapidReader() {}
    
    /** we only need the header as ParamNode, we read the density information directly */
    PtrParamNode CreateHeaderParamNode();

    /** finds the pointer to the set we need */
    rapidxml::xml_node<> * FindSet(const std::string& setname);
    
    /** recursive implementation of Fill */
    void FillChildren(rapidxml::xml_node<> * node, PtrParamNode out);

  private:
    /** the name of the xml-file */
    const std::string filename_;
    /** the file contents */
    std::vector<char> buffer_;
    /** the parsed document */
    rapidxml::xml_document<> doc_;
};

}

#endif // RAPIDXML_HH_
