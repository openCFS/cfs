/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GBDEFINES_HH
#define GBDEFINES_HH

#ifndef INLINE
#ifdef OUTLINE
#define INLINE
#else
#define INLINE inline
#endif
#endif

#include <typeinfo>
#include <iostream>

#ifdef DEBUG
#define infomsg(s) std::cout<<typeid(*this).name()<<": \033[32mINFO:\033[0m "<<s<<std::endl
#define staticinfomsg(n,s) std::cout<<n<<": \033[32mINFO:\033[0m "<<s<<std::endl
#else
#define infomsg(s) std::cout<<s<<std::endl
#endif

#ifdef DEBUG
#define debugmsg(s) std::cerr<<typeid(*this).name()<<": \033[36mDEBUG:\033[0m "<<s<<std::endl
#define staticdebugmsg(n,s) std::cerr<<n<<": \033[36mDEBUG:\033[0m "<<s<<std::endl
#else
#define debugmsg(s)
#endif

#ifdef DEBUG
#define warningmsg(s) std::cerr<<typeid(*this).name()<<": \033[33mWARNING:\033[0m "<<s<<std::endl
#define staticwarningmsg(n,s) std::cerr<<n<<": \033[33mWARNING:\033[0m "<<s<<std::endl
#else
#define warningmsg(s) std::cerr<<"\033[33mWARNING:\033[0m "<<s<<std::endl
#endif

#ifdef DEBUG
#define errormsg(s) std::cerr<<typeid(*this).name()<<": \033[31mERROR:\033[0m "<<s<<std::endl
#define staticerrormsg(n,s) std::cerr<<n<<": \033[31mERROR:\033[0m "<<s<<std::endl
#else
#define errormsg(s) std::cerr<<"\033[31mERROR:\033[0m "<<s<<std::endl
#endif

// the following is here to steer memory pooling
#define MEMPOOLINC 64*1024
//1048576
//5242880
//2048
// this defines the paging size
#define MEMPOOLBLOCKS 4096

// some macros for better code readability used in rendering subsystem
#define GR_DECLARE_DEFAULT_STATE(classname) \
    public: \
        static classname* getDefault () { return kDefault_; } \
    protected: \
        static classname* kDefault_; \
        friend class _##classname##InitTermDS

#define GR_IMPLEMENT_DEFAULT_STATE(classname) \
    classname* classname::kDefault_ = NULL; \
    class _##classname##InitTermDS { \
    public: \
        _##classname##InitTermDS () { \
            classname::kDefault_ = new classname; \
            GrRenderState::Type eType = classname::kDefault_->getType(); \
            classname::defaultStates_[eType] = classname::kDefault_; \
        } \
        ~_##classname##InitTermDS () { \
            delete classname::kDefault_; \
        } \
    }; \
    static _##classname##InitTermDS _force##classname##InitTermDS

// control endianess
#define swabi2(i2) (((i2) >> 8) + (((i2) & 255) << 8))
#define swabi4(i4) (((i4) >> 24) + (((i4) >> 8) & 65280) + (((i4) & 65280) << 8) + (((i4) & 255) << 24))

#endif // GBDEFINES_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:56  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.10  2001/12/11 12:41:05  prkipfer
| fixes for KCC compiler on new PCs
|
| Revision 1.9  2001/06/11 10:33:11  prkipfer
| moved defines to GbDefines
|
| Revision 1.8  2001/01/02 15:21:33  prkipfer
| introduced cloning and support for new base classes
|
| Revision 1.7  2000/10/23 09:36:12  prkipfer
| small modif
|
| Revision 1.6  2000/09/07 16:57:17  prkipfer
| added subdivision
|
| Revision 1.5  2000/08/22 17:15:22  prkipfer
| added possibility to query mesh contents
|
| Revision 1.4  2000/07/17 11:13:28  prkipfer
| moved compiler commandline defines here
|
| Revision 1.3  2000/07/03 16:15:01  prkipfer
| extended IO subsystem, created class registry, added tetrahedron
|
| Revision 1.2  2000/06/16 13:33:43  prkipfer
| further hierarchy fixes
|
| Revision 1.1  2000/06/14 15:39:10  prkipfer
| improved base classes and added funcstruct processing for mesh
|
|
+---------------------------------------------------------------------*/
