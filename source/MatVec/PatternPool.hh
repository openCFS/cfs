#ifndef OLAS_PATTERNPOOL
#define OLAS_PATTERNPOOL


#include <vector>
#include <algorithm>
#include <iomanip>

#include "SparsityPatterns.hh"


namespace CoupledField {

  //! This class provides functionality for sharing sparsity pattern
  //! information between different matrix objects.
  class PatternPool {

  public:


    //! Default constructor
    PatternPool() {
    }


    //! Destructor
    ~PatternPool() {


      // Check if there are any undeleted patterns, which we consider to be
      // a sever error, due to the memory problems involved when we go out
      // of scope now!
      bool safe = true;
      for ( UInt i = 0; i < patterns_.size(); i++ ) {
        if ( patterns_[i] != NULL || numPatternUsers_[i] >= 0 ) {
          safe = false;
        }
      }

      // If everything is not okay complain
      if ( safe != true ) {
        EXCEPTION( "PatternPool::~PatternPool: Problem detected!\n\n"
                 << " PatternID | no. users\n ---------------------");
        for ( UInt i = 0; i < patterns_.size(); i++ ) {
          EXCEPTION( std::setw(11) << std::setfill( ' ' ) << i << " | "
                   << std::setw(11) << numPatternUsers_[i]);
        }
      }
    }


    //! Insert a pattern into the pool for common use

    //! This method must be called in order to insert the sparsity pattern
    //! of a matrix into the pattern pool so that it becomes available for
    //! shared use with other matrix objects. The method will automatically
    //! generate and return an identifier for that pattern that allows to
    //! refer to the pattern in future communications with the %PatternPool
    //! object.
    //! \note
    //! - Putting the pattern into the pool does not imply any use of
    //!   that pattern. Thus, if a matrix object puts its pattern into
    //!   the pool in order to share it, it should call RegisterUser()
    //!   right afterwards.
    //! - This basically allows us to also generate a sparsity pattern without
    //!   an associated matrix object and to put it into the pool
    PatternIdType InsertPattern( BaseSparsityPattern *newPattern ) {


      // safety check
      if ( newPattern == NULL ) {
        EXCEPTION( "You do not really want to insert a NULL pointer "
                 << "into the pattern pool, do you!");
      }

      // insert pattern into set and nullify pointer
      patterns_.push_back( newPattern );

      // set usage count to zero
      numPatternUsers_.push_back( 0 );

      // generate and return pattern identifier
      return static_cast<PatternIdType>(numPatternUsers_.size());
    }


    //! Register as user of a certain sparsity pattern

    //! This method must be called to inform the %PatternPool object that
    //! another matrix object requires a specified sparsity pattern. This
    //! call will increase the usage count of the specified pattern. Once
    //! the matrix object does no longer require the sparsity pattern it
    //! should call DeRegisterUser() in order to inform the %PatternPool
    //! about this fact.
    //! \return A pointer to a BaseSparsityPattern object representing the
    //!         pattern belonging to the specified pattern.
    BaseSparsityPattern* RegisterUser( PatternIdType patternID ) {


      // Check that identifier is legal
      if ( patternID > patterns_.size() ) {
        EXCEPTION( "PatternPool::RegisterUser: patternID = '"
                 << patternID << "' is not in range [1,"
                 << patterns_.size() << "]");
      }

      // Check that pattern still exists
      if ( patterns_[patternID - 1] == NULL ) {
        if ( numPatternUsers_[patternID - 1] < 0 ) {
          EXCEPTION( "PatternPool::RegisterUser: Sorry, but I must inform "
                   << "you that the pattern with ID = '" << patternID
                   << "' has been deleted in the mean-time!");
        }
        else {
          EXCEPTION( "PatternPool::RegisterUser: Internal error! "
                   << "Darn! The pattern pointer is NULL, but the usage "
                   << "count is '" << numPatternUsers_[patternID - 1] << "'");
        }
      }

      // Increase usage count of pattern
      numPatternUsers_[patternID - 1]++;

      // Return that pattern
      return patterns_[patternID - 1];
    }


    //! De-register a user of a certain sparsity pattern

    //! This method can be called to inform the %PatternPool object that
    //! a certain user no longer require the specified sparsity pattern.
    //! \note In order to maintain a low memory footprint the current
    //!       implementation of the %PatternPool concept will automatically
    //!       delete a sparsity pattern once the count of register users
    //!       falls by a call to %DeRegisterUser().
    void DeRegisterUser( PatternIdType patternID ) {


      // Check that identifier is legal
      if ( patternID > patterns_.size() ) {
        EXCEPTION( "PatternPool::DeRegisterUser: patternID = '"
                 << patternID << "' is not in range [1,"
                 << patterns_.size() << "]");
      }

      // Check that pattern still exists
      if ( patterns_[patternID - 1] == NULL ) {
        if ( numPatternUsers_[patternID - 1] < 0 ) {
          EXCEPTION( "PatternPool::DeRegisterUser: Sorry, but I must "
                   << "inform you that the pattern with ID = '"
                   << patternID << "' has already been deleted!");
        }
        else {
          EXCEPTION( "PatternPool::DeRegisterUser: Internal error! "
                   << "Darn! The pattern pointer is NULL, but the usage "
                   << "count is '" << numPatternUsers_[patternID - 1] << "'");
        }
      }

      // Check that there are still users left
      if ( numPatternUsers_[patternID - 1] == 0 ) {
        EXCEPTION( "PatternPool::DeRegisterUser: No registered users "
                 << "for pattern '" << patternID << "'. So why do you want "
                 << "to de-register, honey?");
      }

      // Decrease usage count of pattern
      numPatternUsers_[patternID - 1]--;

      // Check if pattern can be deleted
      if ( numPatternUsers_[patternID - 1] == 0 ) {
        delete patterns_[patternID - 1];
        patterns_[patternID - 1] = NULL;
        numPatternUsers_[patternID - 1] = -1;
      }
    }

  private:

    //! Vector of patterns in the pool

    //! This vector stores all the patterns that were put into the pool via
    //! the InsertPattern() method.
    std::vector<BaseSparsityPattern*> patterns_;

    //! Number of registered users of each pattern

    //! This attribute stores for each pattern currently in the pattern pool
    //! the number of current users of that pattern that have registered
    //! themselves via the RegisterPatternUser() method and not yet
    //! de-registed via the DeRegisterPatternUser() method. A negative
    //! usage count indicates that the pattern existed previously, but was
    //! deleted, since the number of registered users dropped to zero in the
    //! mean-time.
    std::vector<Integer> numPatternUsers_;
  };

}

#endif
