/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/

#ifndef GSI_ENDIAN_CONVERTER
#define GSI_ENDIAN_CONVERTER

#include <vector>

#include "GSITypeDefs.hh"

namespace GridlibSocketInterface
{
  
extern bool BigEndianSystem;

void InitEndian( void );		//makes use of a clever trick in Quake 2

int16 INT16_Swap( int16 i );
int16 INT16_NoSwap( int16 i );

uint16 UINT16_Swap( uint16 i );
uint16 UINT16_NoSwap( uint16 i );

int32 INT32_Swap( int32 i );
int32 INT32_NoSwap( int32 i );

uint32 UINT32_Swap( uint32 i );
uint32 UINT32_NoSwap( uint32 i );

real32 REAL32_Swap( real32 f );
real32 REAL32_NoSwap( real32 f );

real64 REAL64_Swap( real64 f );
real64 REAL64_NoSwap( real64 f );

//Use these functions
extern int16 (*BigINT16) ( int16 i );
extern int16 (*LittleINT16) ( int16 i );
extern uint16 (*BigUINT16) ( uint16 i );
extern uint16 (*LittleUINT16) ( uint16 i );
extern int32 (*BigINT32) ( int32 i );
extern int32 (*LittleINT32) ( int32 i );
extern uint32 (*BigUINT32) ( uint32 i );
extern uint32 (*LittleUINT32) ( uint32 i );

extern real32 (*BigREAL32) ( real32 f );
extern real32 (*LittleREAL32) ( real32 f );
extern real64 (*BigREAL64) ( real64 f );
extern real64 (*LittleREAL64) ( real64 f );

bool BigINT16Array(int16 *dst, int16 *src, int size);
bool LittleINT16Array(int16 *dst, int16 *src, int size);
bool BigINT16Vector(std::vector<int16> &dst, std::vector<int16> &src);
bool LittleINT16Vector(std::vector<int16> &dst, std::vector<int16> &src);

bool BigUINT16Array(uint16 *dst, uint16 *src, int size);
bool LittleUINT16Array(uint16 *dst, uint16 *src, int size);
bool BigUINT16Vector(std::vector<uint16> &dst, std::vector<uint16> &src);
bool LittleUINT16Vector(std::vector<uint16> &dst, std::vector<uint16> &src);

bool BigINT32Array(int32 *dst, int32 *src, int size);
bool LittleINT32Array(int32 *dst, int32 *src, int size);
bool BigINT32Vector(std::vector<int32> &dst, std::vector<int32> &src);
bool LittleINT32Vector(std::vector<int32> &dst, std::vector<int32> &src);

bool BigUINT32Array(uint32 *dst, uint32 *src, int size);
bool LittleUINT32Array(uint32 *dst, uint32 *src, int size);
bool BigUINT32Vector(std::vector<uint32> &dst, std::vector<uint32> &src);
bool LittleUINT32Vector(std::vector<uint32> &dst, std::vector<uint32> &src);

bool BigREAL32Array(real32 *dst, real32 *src, int size);
bool LittleREAL32Array(real32 *dst, real32 *src, int size);
bool BigREAL32Vector(std::vector<real32> &dst, std::vector<real32> &src);
bool LittleREAL32Vector(std::vector<real32> &dst, std::vector<real32> &src);

bool BigREAL64Array(real64 *dst, real64 *src, int size);
bool LittleREAL64Array(real64 *dst, real64 *src, int size);
bool BigREAL64Vector(std::vector<real64> &dst, std::vector<real64> &src);
bool LittleREAL64Vector(std::vector<real64> &dst, std::vector<real64> &src);

}


#endif // ENDIANCONVERTER
