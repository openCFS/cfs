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

float32 FLOAT32_Swap( float32 f );
float32 FLOAT32_NoSwap( float32 f );

float32 FLOAT32_Swap( float64 f );
float32 FLOAT32_NoSwap( float64 f );

//Use these functions
extern int16 (*BigINT16) ( int16 i );
extern int16 (*LittleINT16) ( int16 i );
extern uint16 (*BigUINT16) ( uint16 i );
extern uint16 (*LittleUINT16) ( uint16 i );
extern int32 (*BigINT32) ( int32 i );
extern int32 (*LittleINT32) ( int32 i );
extern uint32 (*BigUINT32) ( uint32 i );
extern uint32 (*LittleUINT32) ( uint32 i );

extern float32 (*BigFLOAT32) ( float32 f );
extern float32 (*LittleFLOAT32) ( float32 f );
extern float64 (*BigFLOAT64) ( float64 f );
extern float64 (*LittleFLOAT64) ( float64 f );

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

bool BigFLOAT32Array(float32 *dst, float32 *src, int size);
bool LittleFLOAT32Array(float32 *dst, float32 *src, int size);
bool BigFLOAT32Vector(std::vector<float32> &dst, std::vector<float32> &src);
bool LittleFLOAT32Vector(std::vector<float32> &dst, std::vector<float32> &src);

bool BigFLOAT64Array(float64 *dst, float64 *src, int size);
bool LittleFLOAT64Array(float64 *dst, float64 *src, int size);
bool BigFLOAT64Vector(std::vector<float64> &dst, std::vector<float64> &src);
bool LittleFLOAT64Vector(std::vector<float64> &dst, std::vector<float64> &src);

}


#endif // ENDIANCONVERTER
