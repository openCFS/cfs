// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
|
| Revision 1.5  2001/06/15 08:36:56  prkipfer
| fixed range overflow bug
|
| Revision 1.4  2000/07/20 13:32:38  prkipfer
| complete rework: now KCC is the standard compiler for Linux
|
| Revision 1.3  2000/07/17 10:48:36  prkipfer
| changed to static memory block handling because g++ doesn't like partially specialized templates
|
| Revision 1.2  2000/07/03 16:15:02  prkipfer
| extended IO subsystem, created class registry, added tetrahedron
|
| Revision 1.1  2000/06/14 15:39:13  prkipfer
| improved base classes and added funcstruct processing for mesh
|
|
|
+---------------------------------------------------------------------*/

#include <algorithm>
#include <iomanip>

// init static members
template<class T>
GbMemPool<T> *
GbMemPool<T>::headOfFreeList_;

template <class T>
unsigned long
GbMemPool<T>::poolSize_ = 0;

template<class T>
const int
GbMemPool<T>::poolIncrement_ = MEMPOOLINC;

template <class T>
GbMemPool<T> *
GbMemPool<T>::poolMem_[MEMPOOLBLOCKS];

template <class T>
unsigned int
GbMemPool<T>::poolMemSize_ = 0;

template <class T>
GbMemPool<T>::GbMemPool()
{
  debugmsg("constructor: item size "<<sizeof(GbMemPool<T>)<<" bytes");
  debugmsg("             ptr 0x"<<std::hex<<sizeof(GbMemPool<T> *)<<" struct 0x"<<sizeof(T)<<std::dec);
  debugmsg("             block size "<<poolIncrement_ * sizeof(T)<<" bytes");
  unsigned long a = MEMPOOLBLOCKS * poolIncrement_ * sizeof(T);
  char s[8];
  formatNumber(&a,s);
  debugmsg("             max linear capacity "<<MEMPOOLBLOCKS<<" blocks ("<<a<<" "<<s<<"B)");
  // nothing else to do here, alloc is done in operator new
}

template <class T>
GbMemPool<T>::~GbMemPool()
{
  if (poolSize_==0) {
    freePoolMem();
  }
  else {
    warningmsg("destructing memory pool that contains "<<poolSize_<< " items -> not releasing memory");
  }
}

template<class T>
void *
GbMemPool<T>::alloc(size_t n)
{
  // if not this object type redirect to global new
  if (n != sizeof(T)) {
    warningmsg("size of "<<n<<" does not match pool's item size ("<<sizeof(T)<<") -> calling global new");
    return ::operator new(n);
  }

  // get start of free space
  GbMemPool<T> *p = headOfFreeList_;

  // take first free element, if available
  if (p) {
    headOfFreeList_ = p->next_;
  }
  else {
    // capacity limit ?
    if (MEMPOOLBLOCKS == poolMemSize_) {
      warningmsg("memory pool capacity limit reached");
      return NULL;
    }

    // need more space -> get a whole block
    GbMemPool<T> *newBlock = static_cast<GbMemPool<T>*>(::operator new(poolIncrement_ * sizeof(T)));
    poolMem_[poolMemSize_] = newBlock; poolMemSize_++;
    debugmsg("allocating new memory block of "<<poolIncrement_ * sizeof(T)<<" bytes");

    // chain the free space block
    for (int i=1; i<poolIncrement_-1; ++i)
      newBlock[i].next_ = &newBlock[i+1];

    // last element has no successor
    newBlock[poolIncrement_-1].next_ = 0;

    // adjust state, return first free element adress
    p = newBlock;
    headOfFreeList_ = &newBlock[1];
  }

  poolSize_++;
  return p;
}

template<class T>
void 
GbMemPool<T>::free(void *p, size_t n)
{
  // its legal to delete NULL pointers
  if (p == 0) return;

  // if not this object type redirect to global delete
  if (n != sizeof(T)) {
    warningmsg("size of "<<n<<" does not match pool's item size ("<<sizeof(T)<<") -> calling global delete");
    ::operator delete(p);
    return;
  }

  // put space back into free list
  GbMemPool<T> *dead = static_cast<GbMemPool<T> *>(p);
  dead->next_ = headOfFreeList_;
  headOfFreeList_ = dead;
  poolSize_--;
#ifdef MINIMEM
  if(!poolSize_) freePoolMem();
#endif
}

template<class T>
void 
GbMemPool<T>::shrink()
{
  debugmsg("shrinking pool to minimum size...");

  GbMemPool<T> *slot = headOfFreeList_;
//    std::vector<GbMemPool<T> *>::iterator deathobj;
  std::vector<GbMemPool<T> *> deathrow;
//  std::vector<char *>::iterator dead;
//    std::vector<char *> alreadydead;

  if (slot == NULL) {
    debugmsg("no free space in pool -> no memory cleared");
    return;
  }

  // build lookup tables of dead objects and start of memblocks
//    while (slot != NULL) {
//      alreadydead.push_back((char*)slot);
//    deathobj = find(poolMem_.begin(), poolMem_.end(), slot);
//    if (deathobj!=poolMem_.end())
//      deathrow.push_back(slot);
//      slot = slot->next_;
//    }
//    debugmsg(alreadydead.size()<<" dead objects in "<<deathrow.size()<<" blocks");

//  int killed=0;
  if (deathrow.size()>0) {

    warningmsg("implementation not finished yet, aborting routine");
//!! find ueber block grenzen hinweg nicht sicher
#if 0
    // check the memblocks found
    for (std::vector<GbMemPool<T> *>::iterator iter = deathrow.begin();
	 iter != deathrow.end();
	 iter++) {
      
      int i;//,counter=0;
      // search through dead object list
      for (i=0; i<poolIncrement_; i++) {
	dead = find(alreadydead.begin(), alreadydead.end(), (char *)((*iter)+i));
	if (dead==alreadydead.end())
	  break;
//	for(dead = alreadydead.begin(); dead != alreadydead.end(); dead++)
//	  if ((*dead)==((char *)((*iter)+i))) counter++;
      }
//      std::cerr<<i<<"found ";

      // we can delete the memblock if all objects are dead
      if (i==poolIncrement_) {
	// all dead -> kill block
//	deathobj = find(poolMem_.begin(), poolMem_.end(), (*iter));
//	poolMem_.erase(deathobj);
	::operator delete(*iter);
//	debugmsg("memory block deleted");
	killed++;
      }
      else {
	// not all dead -> remove from slot chain fix list
	debugmsg(i<<" slots will remain free because of memory block alignment");
	for (int j=0; j<i; j++) {
	  dead = find(alreadydead.begin(), alreadydead.end(), (char *)((*iter)+j));
	  alreadydead.erase(dead);
	}
      }

    } // for deathrow

    // fix free slot chain -> process dead object list
    // all objects in list have no memory block anymore
    slot = headOfFreeList_;
    while (slot != NULL) {
      dead = find(alreadydead.begin(), alreadydead.end(), (char *)slot);
      // if we find one the chain is pointing to nowhere -> set NULL
      if (dead!=alreadydead.end()) {
	slot->next_ = NULL;
	break;
      }
      slot = slot->next_;
    }

    debugmsg(killed*poolIncrement_*sizeof(T)<<" bytes of memory cleared");
#endif
  }
  else {
    debugmsg("no potentially empty memory blocks -> nothing cleared");
  }
//    alreadydead.clear();
//    deathrow.clear();
}

template<class T>
void 
GbMemPool<T>::freePoolMem()
{
  // no elements in pool -> be smart and free memory
  debugmsg("pool is empty -> freeing memory");

  for(unsigned int i=0; i<poolMemSize_; ++i) {
    ::operator delete(poolMem_[i]);
  }

  debugmsg(poolMemSize_*poolIncrement_*sizeof(T)<<" bytes of memory cleared");
  poolMemSize_=0;
  headOfFreeList_ = NULL;
}

template<class T>
void 
GbMemPool<T>::getStatistic()
{
  unsigned long b = poolSize_*sizeof(T);
  unsigned long c = poolMemSize_*poolIncrement_*sizeof(T);
  unsigned long d = MEMPOOLBLOCKS*poolIncrement_*sizeof(T);
  double bd = (double)(b*100)/(double)d;
  double bc;
  char s[8];

  if(c>0)
    bc = (double)(b*100)/(double)c;

  infomsg("statistics:");

  unsigned long a = b;
  formatNumber(&a,s);
  infomsg("elements: "<<poolSize_<<" elements x "<<sizeof(T)<<" bytes = "<<b<<" ("<<a<<" "<<s<<"B)");

  a = c;
  formatNumber(&a,s);
  infomsg("memory:   "<<poolMemSize_<<" blocks x "<<poolIncrement_<<" elements = "<<c<<" ("<<a<<" "<<s<<"B)");

  infomsg("\t--> usage "<<std::setprecision(3)<<bd<<"%");
  if(c>0)
    infomsg("\t--> efficiency "<<std::setprecision(3)<<bc<<"%");
}

template<class T>
void
GbMemPool<T>::formatNumber(unsigned long* a, char* s)
{
  if(*a>1024) { *a/=1024; strcpy(s,"k"); }
  if(*a>1024) { *a/=1024; strcpy(s,"M"); }
  if(*a>1024) { *a/=1024; strcpy(s,"G"); }
}
