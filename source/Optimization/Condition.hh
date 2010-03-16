#ifndef CONDITION_HH_
#define CONDITION_HH_

#include "Optimization/DesignElement.hh"
#include "General/Enum.hh"
#include "General/environment.hh"
#include "MatVec/matrix.hh"

using std::pair;

namespace CoupledField
{
   class ParamNode;
   class ParamNode;
   class DenseMatrix;

   /** our constraint criteria. Can be filled directly from XML */
   class Condition 
   {
     public:
       /** Empty constructor for StdVector */
       Condition() { };
     
       /** Call this method to append a Condition. This calls the actual (private) constructor.
        * Index is set with position of the relevant list.
        * If it is a homogenization constraint there might be a blow up reCoupledField::DesignElement::DEFAULTsulting in several
        * constraints if either multiple entries are given or the entry position is 'all'.
        * @param pn determines active mode
        * @param constraints stuff is added here if the mode is constraint
        * @param observation stuff is added here in observation mode */
       static void AddCondition(PtrParamNode pn, StdVector<Condition>& constraints, StdVector<Condition>& observation);

       /** Note the difference to the Type! See Name as equivalent of Kind! */
       typedef enum
       {
         VOLUME = 0,
         GREYNESS = 1,
         GAUSS_GREYNESS = 2,
         COMPLIANCE = 3,
         TRACKING = 4,
         REALVOLUME = 5,
         HOMOGENIZATION_TENSOR = 6,   // might blow up to several HOMOGENIZATION_TENSOR if a tensor is given
         HOMOGENIZATION_TRACKING = 7,
         POISSONS_RATIO = 8,          // homogenization only
         YOUNGS_MODULUS = 9,          // homogenization only
         ISOTROPY = 10                // blows up to several HOMOGENIZATION_TENSOR constraints! No ISOTROPY will be left!
       } Name;

       /** Genertal constraint types */             
       typedef enum { EQUAL, LOWER_BOUND, UPPER_BOUND } Type;

       /** Be sure not to mix up with Type! */  
       Name GetName() const { return name_; } 

       /** Be sure not to mix up with Name! */             
       Type GetType() const { return type_; }
       
       /** Has only relevance for type = active! */
       int GetIndex() const { return index_; }
       
       /** Check whether condition should be calculated for given region */
       bool IsForRegion(RegionIdType regionId);
       
       /** Check if this is homogenization constraint */
       bool IsHomogenization() const;

       /** This is a nice statement for output */
       std::string ToString() const; 

       /** log to info.xml */
       void ToInfo(PtrParamNode in) const;
       
       /** the tensor exists only in the homogenization constraint case */
       Matrix<double>& GetTensor();

       /** Read the tensor if it is given, otherwise sets to 1.1
        * @param pn might contain a "tensor" child
        * @param matrix where to store the data
        * @return true if the tensor was read */
       static bool ReadTensor(PtrParamNode pn, Matrix<double>& matrix);

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
       
       /** Shall delta constraints be printed? Is only true if a value is given! */
       bool delta_logging;

       /** If condition supports constriction to one region */
       RegionIdType region;
       
       /** Used for caching 1.0 / complete_volume per region */
       double volume_fraction_;

       /** For the homogenization tensor constraint this gives the actual position within the matrix_.
        * The first entry is for homogenization always set.
        * In the case of a "smart" isotropy constraint also E11-E22 = 0 and
        * E11-E12-2E33 = E11-E12-E33-E33 = 0 are generated. Then coord is 2 or 4 entries.
        * Note, that the entries are 1-based!!! */
       StdVector<pair<unsigned int, unsigned int> > coord;

       static Enum<Name> name;              
       static Enum<Type> type;

    private:   
      /** Helper constructor for AddCondition */
      Condition(PtrParamNode pn);

      /** Reads the coord attribute and sets the coord pair if value is not 'all'
       * @return false if 'all' and the coord pair is not set */
      bool ReadCoord(PtrParamNode pn);

      /** Add a subcondition with only index and value set (to zero) */
      Condition* AppendSubCondition(StdVector<Condition>& list);

      /** Create a new homogenization constraint with the given tensor position
       * @param base the base of cloning. Needs to contain a tensor!
       * @param list where to append the child. Index is set
       * @return the appendend child */
      Condition* AppendSubCondition(StdVector<Condition>& list, PtrParamNode entry_pos);

      /** @see other AppendSubCondition() */
      Condition* AppendSubCondition(StdVector<Condition>& list, unsigned int pos_x, unsigned int pos_y);

       /** this index is the position in the Optimization list and is used to
        * identify the constraint gradient in DesignElement. Only relevant for type = active */
       int  index_;
       Name name_;
       Type type_;

       /** for constraints of type "homogenization" one can give the tensor we want to reach
        * or for multiple constraints the entry sub element 'pos' refers to the constraint value */
       Matrix<double> matrix_;

       bool delta_logging_ignored_;
   };

} // namespace


#endif /*CONDITION_HH_*/
