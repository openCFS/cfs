#ifndef CONDITION_HH_
#define CONDITION_HH_

#include "Optimization/DesignElement.hh"
#include "General/Enum.hh"

namespace CoupledField
{
   class ParamNode;

   /** our constraint criteria. Can be filled directly from XML */
   class Condition 
   {
     public:
       /** Empty constructor for StdVector */
       Condition() { };
     
       /** we read our own vaues from the ParamNode
        * @param index the unique position in the Constraint array for the 
        *        DesignElement::constraintGradient array indexing */
       Condition(ParamNode* pn, int index);
     
       /** Note the difference to the Type! See Name as equivalent of Kind! */
       typedef enum { VOLUME = 0, GREYNESS = 1, GAUSS_GREYNESS = 2} Name;

       /** Genertal constraint types */             
       typedef enum { EQUAL, LOWER_BOUND, UPPER_BOUND } Type;

       /** Be sure not to mix up with Type! */  
       Name GetName() const { return name_; } 

       /** Be sure not to mix up with Name! */             
       Type GetType() const { return type_; }
       
       int GetIndex() const { return index_; }
       
       /** This is a nice statement for output */
       std::string ToString() const; 
       
       /** This is DEFAULT (= applies always) if not defined */
       DesignElement::Type design;

       /** The value for upper/lower/equal constraint */ 
       double value;
 
       /** The parameter is optional, e.g. "h" for gaussGreyness.
        * Default = 0.0 */
       double parameter;  
       
       /** determine if this is a real constraint or just for logging the information.
        * <b>this has nothing to do with "active" inequalitly constraints in optimizing!</b> */
       bool active;
       
       /** The scaling is evaluated in the Optimizer IPOPT only. Not in OC! */
       double scaling;
       
       static Enum name;              
       static Enum type;

    private:   
       /** this index is the position in the Optimization list and is used to
        * identify the constraint gradient in DesignElement */ 
       int  index_;
       Name name_;
       Type type_;
   };

} // namespace


#endif /*CONDITION_HH_*/
