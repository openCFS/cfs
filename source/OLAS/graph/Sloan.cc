#include <iostream>
#include <limits>

#include "General/Environment.hh"
#include "OLAS/graph/Sloan.hh"

namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  Sloan::Sloan(StdVector<std::vector<UInt>>& agraph, StdVector<UInt>& order ) :
    BaseOrdering( agraph, order ),
    profOldRem_(0),
    profNewRem_(0),
    profOldMult_(0),
    profNewMult_(0)
  {
    // we do our reordering with vertices starting with one
    // all arrays will have idices starting with one, except
    // order and graph
  }


  // **************
  //   Destructor
  // **************
  Sloan::~Sloan() {
  }


  // ***************
  //   Label Graph
  // ***************
  void Sloan::LabelGraph() {

 
    // label a graph for small profile and rms wavefront

    // working arrays
    Integer *xls, *ls, *hlevel;
    Integer size = order_.GetSize();
    
    // we allocate always for size plus 1; so we can index this
    // arrays starting by 1!!
    NEWARRAY( xls,    Integer, size+2 );
    NEWARRAY( ls,     Integer, size+1 );
    NEWARRAY( hlevel, Integer, size+1 );

    // loop while some nodes remain unnumbered
    Integer startnode = 0;
    Integer nc = 0;
    Integer lstnum = 0;

    while ( lstnum < size ) {

      // find end points of p-diameter for nodes in this component
      // compute distances of nodes from end node
      PseudoDiameter( ls, xls, hlevel, startnode, nc );

      // number nodes in this component
      NumberNodes( nc, startnode, lstnum, ls, xls );
    }

    // Compute old and new profile
    CalcProfile();

    // Clear memory
    delete [] ( xls    );  xls     = NULL;
    delete [] ( ls     );  ls      = NULL;
    delete [] ( hlevel );  hlevel  = NULL;

  }


  // ******************
  //   PseudoDiameter
  // ******************
  void Sloan::PseudoDiameter( Integer *ls, Integer *xls, Integer *hlevel, 
                              Integer &snode, Integer &nc ) {


    // find a pair of nodes with maximum distance within the graph

    // choose first guess for starting node by min degree
    // ( = min number of neighbors)
    // ignore nodes that are invisible (order .ne. 0)
    Integer degree = 0, depth = 0, enode = 0;
    Integer minDegree = 10*order_.GetSize();

    for (UInt n=1, m=order_.GetSize(); n<=m; n++) {
      if (order_[n-1] == 0) {
        degree = graph_[n-1].size();
        if (degree < minDegree) {
          snode = n;
          minDegree = degree;
        }
      }
    }

    // generate level structure for node with min degree
    Integer sdepth, width;
    RootedLevel(snode, order_.GetSize()+1, ls, xls, sdepth, width);

    //store number of nodes in this component
    nc = xls[sdepth+1] -1;

  newStart:

    // store list of nodes that are at max distance from starting node
    // store their degree in xls
    Integer node  = 0;
    Integer hsize = 0;
    Integer istart = xls[sdepth];
    Integer istop  = xls[sdepth+1]-1;
    for (Integer i=istart; i<=istop; i++) {
      node = ls[i];
      hsize += 1;
      hlevel[hsize] = node;
      xls[node] = graph_[node-1].size();
    }


    // sort list of nodes in ascending sequence of their degree
    // use insertion sort algorithm
    if (hsize>1) 
      SortList(hsize, hlevel, xls);

    // remove nodes with duplicate degree
    istop = hsize;
    hsize = 1;
    degree = xls[hlevel[1]];
    for (Integer i=2; i<=istop; i++) {
      node = hlevel[i];
      if ( xls[node] != degree ) {
        degree = xls[node];
        hsize += 1;
        hlevel[hsize] = node;
      }
    }

    // loop over nodes in shrunken level
    Integer ewidth = nc+1;

    for (Integer i=1; i<=hsize; i++) {
      node = hlevel[i];

      //form rooted level structure for each node in shrunken level
      RootedLevel(node, ewidth, ls, xls, depth, width);

      if (width < ewidth) {
        // level structure was not aborted during assembly

        if (depth > sdepth) 
          {
            // level structure of greater depth found
            // store new starting node, new max depth, and begin
            // a new iteration
            snode  = node;
            sdepth = depth;
            goto newStart;
          }
        
        // level structure width for this end node is smallest so far
        // store end node and new min width
        enode = node;
        ewidth = width;
      }
    }

    //  generate level structure rooted at end node if necessary
    if (node != enode)
      RootedLevel(enode, nc+1, ls, xls, depth, width);

    // store distances of each node from end node
    Integer jstart, jstop;
    for (Integer i=1; i<=depth; i++) {
      jstart = xls[i];
      jstop  = xls[i+1] -1;
      for (Integer j=jstart; j<=jstop; j++) 
        order_[ls[j]-1] = i-1;
    }
  }


  void Sloan::RootedLevel( Integer root, Integer maxwidth, Integer *ls,
                           Integer *xls, Integer &depth, Integer &width ) {

 
    // generate rooted level structure using a FORTRAN 77 implementation
    // of the algorithm given by George and Liu

    Integer nc, lstrt, lstop, lwdth, node, nbr;
    StdVector<UInt> mask = order_;
    std::vector<UInt>::iterator liter;

    // initialisation
    mask[root-1] = 1;
    ls[1] = root;
    nc    = 1;
    width = 1;
    depth = 0;
    lstop = 0;
    lwdth = 1;

    while (lwdth > 0 && width < maxwidth) {

      // lwdth is the width of the current level
      // lstrt points to start of current level
      // lstop points to end of current level
      // nc counts the nodes in component

      lstrt = lstop + 1;
      lstop = nc;
      depth += 1;
      xls[depth] = lstrt;

      // generate next level by finding all visible neighbours
      // of nodes in current level
      for (Integer i=lstrt; i<=lstop; i++) {
        node = ls[i];
        for (liter=graph_[node-1].begin(); liter!=graph_[node-1].end(); liter++) {
          nbr = *liter;
          if (mask[nbr-1] == 0) {
            nc  += 1;
            ls[nc] = nbr;
            mask[nbr-1] = 1;
          }
        }
      }

      // compute width of level just assembled and the width of the
      // level structure so far
      lwdth = nc - lstop;
      width = std::max(lwdth, width);
    }

    xls[depth+1] = lstop + 1;

    // reset mask = 0 for nodes in the level structure
    for (Integer i =1; i<=nc; i++)
      mask[ls[i]-1] = 0;

  }


  void Sloan::NumberNodes( Integer nc, Integer snode, Integer &lstnum,
                           Integer* q, Integer* p ) {


    // number nodes in component of graph for small profile and rms
    // wavefront priority is assigned so that nodes with a large distance
    // from the end node and
    // having a low degree (number of neighbors) get a high priority
 
    Integer node;
    StdVector<Integer> s;
    Integer w1 = 1;
    Integer w2 = 2;

    s.Resize(order_.GetSize());
    for(UInt i=0, n = order_.GetSize(); i<n; i++) 
    {
      s[i] = static_cast<Integer>(order_[i]);
    }
    
    // initialise priorities and status for each node of this component
    // initial priority = w1*dist - w2*degree  where:
    // w1     = a positive weight
    // w2     = a positive weight
    // dist   = distance of node from end node
    // degree = initial current degree for node

    for (Integer i=1; i<=nc; i++) {
      node = q[i];
      p[node] = w1*s[node-1] - w2*(graph_[node-1].size()+1);
      s[node-1] = -2;
    }

    // Insert starting node in queue and assign it a preactive status
    // nn is the size of queue
    Integer nn = 1;
    q[nn]    = snode;
    s[snode-1] = -1;

    Integer addres, maxprt, prty, next, nabor, nbr;
    std::vector<UInt>::iterator liter, liter1;

    // loop while queue is not empty
    while ( nn>0 ) {
      //scan queue for node with max priority
      addres = 1;
      maxprt = p[q[1]];
      for (Integer i=2; i<=nn; i++) {
        prty = p[q[i]];
        if (prty>maxprt) {
          addres = i;
          maxprt = prty;
        }
      }

      // next is the node to be numbered next
      next = q[addres];

      // delete node next from queue
      q[addres] = q[nn];
      nn -= 1;
  
      if (s[next-1] == -1) {
        // node next is preactive, examine its neighbours

        for (liter=graph_[next-1].begin(); liter!=graph_[next-1].end(); liter++) {
          // decrease current degree of neighbour by -1
          nbr = *liter;
          p[nbr] += w2;

          // add neighbour to queue if it is inactive
          // assign it a preactive status
          if (s[nbr-1] == -2) {
            nn += 1;
            q[nn] = nbr;
            s[nbr-1] = -1;
          }
        }
      }

      // store new node number for node next
      // status for node next is now postactive
      lstnum += 1;
      s[next-1] = lstnum;

      // search for preactive neighbours of node next 
      for (liter=graph_[next-1].begin(); liter!=graph_[next-1].end(); liter++) {
        nbr = *liter;
        if (s[nbr-1] == -1) {
          // decrease current degree of preactive neighbour by -1
          // assign neighbour an active status
          p[nbr] += w2;
          s[nbr-1] = 0;

          // loop over nodes adjacent ro preactive neighbour
          for ( liter1 = graph_[nbr-1].begin(); liter1 != graph_[nbr-1].end();
                liter1++ ) {
            nabor = *liter1;

            // decrease current degree of adjacent node by -1
            p[nabor] += w2;
            if (s[nabor-1] == -2) {
              //insert inactive node in queue with a preactive status
              nn       += 1;
              q[nn]    = nabor;
              s[nabor-1] = -1;
            }
          }
        }
      }
    }

    for(UInt i=0, n = order_.GetSize(); i<n; i++) 
    {
      order_[i] = static_cast<UInt>(s[i]);
    }
  }


  // ************
  //   SortList
  // ************
  void Sloan::SortList( Integer nl, Integer* list, Integer *key ) {


    Integer t, value;
    bool found;

    // order a list of integers in asending sequence of their keys
    // using insertion sort
    for (Integer i=2; i<=nl; i++) {
      found = false;
      t     = list[i];
      value = key[t];
      for (Integer j=i-1; j>=1; --j) {
        if (value >= key[list[j]]) {
          list[j+1] = t;
          found     = true;
          continue;
        }
        list[j+1] = list[j];
      }
      if (!found) 
        list[1] = t;
    }

  }


  // ***************
  //   CalcProfile
  // ***************
  void Sloan::CalcProfile() {


    std::vector<UInt>::iterator liter;
    UInt node, oldMin, newMin;
    UInt profOld = 0;
    UInt profNew = 0;
    UInt aux = 0;

    for ( UInt i = 1, n=order_.GetSize(); i <= n; i++ ) {

      oldMin = n;
      newMin = n;

      for ( liter = graph_[i-1].begin(); liter != graph_[i-1].end(); liter++ ) {
        node = *liter;
        oldMin = std::min( oldMin, node );
        newMin = std::min( newMin, order_[node-1] );
      }

      if ( i > oldMin ) {
        aux = profOld;
        profOld += i - oldMin;
        if ( profOld < aux ) {
          profOldMult_++;
        }
      }

      if ( order_[i-1] > newMin ) {
        aux = profNew;
        profNew += order_[i-1] - newMin;
        if ( profNew < aux ) {
          profNewMult_++;
        }
      }
    }

    // add diagonal terms
    aux = profOld;
    profOld += order_.GetSize();
    if ( profOld < aux ) {
      profOldMult_++;
    }

    aux = profNew;
    profNew += order_.GetSize();
    if ( profNew < aux ) {
      profNewMult_++;
    }

    // assign remainders
    profOldRem_ = profOld;
    profNewRem_ = profNew;

  }


  // **************
  //   GetProfile
  // **************
  void Sloan::GetProfile( Double &oldprof, Double &newprof ) {


    // Check if we must compute the profile
    if ( profOldRem_ == 0 ) {
      CalcProfile();
    }

    UInt max = std::numeric_limits<UInt>::max();
    
    // Convert profile information to floating point representation
    oldprof = profOldMult_ * max + profOldRem_;
    newprof = profNewMult_ * max + profNewRem_;

  }


}//namespace
