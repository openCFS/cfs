#ifndef CONDITION_HH_
#define CONDITION_HH_

#include "Optimization/DesignElement.hh"
#include "General/Enum.hh"
#include "General/environment.hh"

namespace CoupledField
{
   class ParamNode;
   class InfoNode;

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
       typedef enum { VOLUME = 0, GREYNESS = 1, GAUSS_GREYNESS = 2, COMPLIANCE = 3, TRACKING = 4} Name;

       /** Genertal constraint types */             
       typedef enum { EQUAL, LOWER_BOUND, UPPER_BOUND } Type;

       /** Be sure not to mix up with Type! */  
       Name GetName() const { return name_; } 

       /** Be sure not to mix up with Name! */             
       Type GetType() const { return type_; }
       
       int GetIndex() const { return index_; }
       
       /** Check whether condition should be calculated for given region */
       bool IsForRegion(RegionIdType regionId);
       
       /** This is a nice statement for output */
       std::string ToString() const; 

       /** log to info.xml */
       void ToInfo(InfoNode* in) const;
       
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

       /**The penalty formulation allows to add constraints via this penalty term to the objective.
        * Actually a penelty method finds iteratively the right value, in practice it is a given
        * parameter. Currently this is *only* implemented for the *Level-Set* method! */
       double penalty;
       
       /** If condition supports constriction to one region */
       RegionIdType region;
       
       /** Used for caching 1.0 / complete_volume per region */
       double volume_fraction_;
       
       static Enum<Name> name;              
       static Enum<Type> type;

    private:   
       /** this index is the position in the Optimization list and is used to
        * identify the constraint gradient in DesignElement */ 
       int  index_;
       Name name_;
       Type type_;
   };

} // namespace


#endif /*CONDITION_HH_*/
