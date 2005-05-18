#ifdef USE_RCSID
static const char RCSid_GSIEndianConverter[] = "$Id$";
#endif

/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.2  2005/05/18 19:26:02  strieben
| Upgraded GSI library to newest available version.
|
| Revision 1.2  2004/09/01 15:24:58  simon
| Added support for writing int16 and uint16
|
| Revision 1.1.1.1  2004/08/31 15:53:00  simon
| Initial GSI import
|
|
|
| adapted from Writing Endian Independent Code in C++
| by Promit Roy
| http://www.gamedev.net/reference/articles/article2091.asp
|
+---------------------------------------------------------------------*/

#include "GSIEndianConverter.hh"

namespace GridlibSocketInterface
{

bool BigEndianSystem;

int16 (*BigINT16) ( int16 i );
int16 (*LittleINT16) ( int16 i );
uint16 (*BigUINT16) ( uint16 i );
uint16 (*LittleUINT16) ( uint16 i );
int32 (*BigINT32) ( int32 i );
int32 (*LittleINT32) ( int32 i );
uint32 (*BigUINT32) ( uint32 i );
uint32 (*LittleUINT32) ( uint32 i );
real32 (*BigREAL32) ( real32 f );
real32 (*LittleREAL32) ( real32 f );
real64 (*BigREAL64) ( real64 f );
real64 (*LittleREAL64) ( real64 f );


//adapted from Quake 2 source

int16 INT16_Swap( int16 s )
{
	uint8 b1, b2;

	b1 = s & 255;
	b2 = (s >> 8) & 255;

	return (b1 << 8) + b2;
}

int16 INT16_NoSwap( int16 s )
{
	return s;
}

uint16 UINT16_Swap( uint16 s )
{
	uint8 b1, b2;

	b1 = s & 255;
	b2 = (s >> 8) & 255;

	return (b1 << 8) + b2;
}

uint16 UINT16_NoSwap( uint16 s )
{
	return s;
}


int32 INT32_Swap (int32 i)
{
	uint8 b1, b2, b3, b4;

	b1 = i & 255;
	b2 = ( i >> 8 ) & 255;
	b3 = ( i>>16 ) & 255;
	b4 = ( i>>24 ) & 255;

	return ((int32)b1 << 24) + ((int32)b2 << 16) + ((int32)b3 << 8) + b4;
}

int32 INT32_NoSwap( int32 i )
{
	return i;
}


uint32 UINT32_Swap (uint32 i)
{
	uint8 b1, b2, b3, b4;

	b1 = i & 255;
	b2 = ( i >> 8 ) & 255;
	b3 = ( i>>16 ) & 255;
	b4 = ( i>>24 ) & 255;

	return ((uint32)b1 << 24) + ((uint32)b2 << 16) + ((uint32)b3 << 8) + b4;
}

uint32 UINT32_NoSwap( uint32 i )
{
	return i;
}


real32 REAL32_Swap( real32 f )
{
	union
	{
		real32 f;
		uint8 b[4];
	} dat1, dat2;

	dat1.f = f;
	dat2.b[0] = dat1.b[3];
	dat2.b[1] = dat1.b[2];
	dat2.b[2] = dat1.b[1];
	dat2.b[3] = dat1.b[0];
	return dat2.f;
}

real32 REAL32_NoSwap( real32 f )
{
	return f;
}

real64 REAL64_Swap( real64 f )
{
	union
	{
		real64 f;
		uint8 b[8];
	} dat1, dat2;

	dat1.f = f;
	dat2.b[0] = dat1.b[7];
	dat2.b[1] = dat1.b[6];
	dat2.b[2] = dat1.b[5];
	dat2.b[3] = dat1.b[4];
	dat2.b[4] = dat1.b[3];
	dat2.b[5] = dat1.b[2];
	dat2.b[6] = dat1.b[1];
	dat2.b[7] = dat1.b[0];

	return dat2.f;
}

real64 REAL64_NoSwap( real64 f )
{
	return f;
}

bool BigINT16Array(int16 *dst, int16 *src, int size)
{
	if((dst == NULL) || (src == NULL))
		return false;

	for(int i = 0; i<size; i++)
		dst[i] = BigINT16(src[i]);
	return true;
}

bool LittleINT16Array(int16 *dst, int16 *src, int size)
{
	if((dst == NULL) || (src == NULL))
		return false;

	for(int i = 0; i<size; i++)
		dst[i] = LittleINT16(src[i]);
	return true;
}

bool BigINT16Vector(std::vector<int16> &dst, std::vector<int16> &src)
{
	if(dst.size() != src.size())
		return false;

	int size = dst.size();
	for(int i = 0; i<size; i++)
		dst[i] = BigINT16(src[i]);
	return true;
}

bool LittleINT16Vector(std::vector<int16> &dst, std::vector<int16> &src)
{
	if(dst.size() != src.size())
		return false;

	int size = dst.size();
	for(int i = 0; i<size; i++)
		dst[i] = LittleINT16(src[i]);
	return true;
}

bool BigUINT16Array(uint16 *dst, uint16 *src, int size)
{
	if((dst == NULL) || (src == NULL))
		return false;

	for(int i = 0; i<size; i++)
		dst[i] = BigUINT16(src[i]);
	return true;
}

bool LittleUINT16Array(uint16 *dst, uint16 *src, int size)
{
	if((dst == NULL) || (src == NULL))
		return false;

	for(int i = 0; i<size; i++)
		dst[i] = LittleUINT16(src[i]);
	return true;
}

bool BigUINT16Vector(std::vector<uint16> &dst, std::vector<uint16> &src)
{
	if(dst.size() != src.size())
		return false;

	int size = dst.size();
	for(int i = 0; i<size; i++)
		dst[i] = BigUINT16(src[i]);
	return true;
}

bool LittleUINT16Vector(std::vector<uint16> &dst, std::vector<uint16> &src)
{
	if(dst.size() != src.size())
		return false;

	int size = dst.size();
	for(int i = 0; i<size; i++)
		dst[i] = LittleUINT16(src[i]);
	return true;
}

bool BigINT32Array(int32 *dst, int32 *src, int size)
{
	if((dst == NULL) || (src == NULL))
		return false;

	for(int i = 0; i<size; i++)
		dst[i] = BigINT32(src[i]);
	return true;
}

bool LittleINT32Array(int32 *dst, int32 *src, int size)
{
	if((dst == NULL) || (src == NULL))
		return false;

	for(int i = 0; i<size; i++)
		dst[i] = LittleINT32(src[i]);
	return true;
}

bool BigINT32Vector(std::vector<int32> &dst, std::vector<int32> &src)
{
	if(dst.size() != src.size())
		return false;

	int size = dst.size();
	for(int i = 0; i<size; i++)
		dst[i] = BigINT32(src[i]);
	return true;
}

bool LittleINT32Vector(std::vector<int32> &dst, std::vector<int32> &src)
{
	if(dst.size() != src.size())
		return false;

	int size = dst.size();
	for(int i = 0; i<size; i++)
		dst[i] = LittleINT32(src[i]);
	return true;
}

bool BigUINT32Array(uint32 *dst, uint32 *src, int size)
{
	if((dst == NULL) || (src == NULL))
		return false;

	for(int i = 0; i<size; i++)
		dst[i] = BigUINT32(src[i]);
	return true;
}

bool LittleUINT32Array(uint32 *dst, uint32 *src, int size)
{
	if((dst == NULL) || (src == NULL))
		return false;

	for(int i = 0; i<size; i++)
		dst[i] = LittleUINT32(src[i]);
	return true;
}

bool BigUINT32Vector(std::vector<uint32> &dst, std::vector<uint32> &src)
{
	if(dst.size() != src.size())
		return false;

	int size = dst.size();
	for(int i = 0; i<size; i++)
		dst[i] = BigUINT32(src[i]);
	return true;
}

bool LittleUINT32Vector(std::vector<uint32> &dst, std::vector<uint32> &src)
{
	if(dst.size() != src.size())
		return false;

	int size = dst.size();
	for(int i = 0; i<size; i++)
		dst[i] = LittleUINT32(src[i]);
	return true;
}

bool BigREAL32Array(real32 *dst, real32 *src, int size)
{
	if((dst == NULL) || (src == NULL))
		return false;

	for(int i = 0; i<size; i++)
		dst[i] = BigREAL32(src[i]);
	return true;
}

bool LittleREAL32Array(real32 *dst, real32 *src, int size)
{
	if((dst == NULL) || (src == NULL))
		return false;

	for(int i = 0; i<size; i++)
		dst[i] = LittleREAL32(src[i]);
	return true;
}

bool BigREAL32Vector(std::vector<real32> &dst, std::vector<real32> &src)
{
	if(dst.size() != src.size())
		return false;

	int size = dst.size();
	for(int i = 0; i<size; i++)
		dst[i] = BigREAL32(src[i]);
	return true;
}

bool LittleREAL32Vector(std::vector<real32> &dst, std::vector<real32> &src)
{
	if(dst.size() != src.size())
		return false;

	int size = dst.size();
	for(int i = 0; i<size; i++)
		dst[i] = LittleREAL32(src[i]);
	return true;
}

bool BigREAL64Array(real64 *dst, real64 *src, int size)
{
	if((dst == NULL) || (src == NULL))
		return false;

	for(int i = 0; i<size; i++)
		dst[i] = BigREAL64(src[i]);
	return true;
}

bool LittleREAL64Array(real64 *dst, real64 *src, int size)
{
	if((dst == NULL) || (src == NULL))
		return false;

	for(int i = 0; i<size; i++)
		dst[i] = LittleREAL64(src[i]);
	return true;
}

bool BigREAL64Vector(std::vector<real64> &dst, std::vector<real64> &src)
{
	if(dst.size() != src.size())
		return false;

	int size = dst.size();
	for(int i = 0; i<size; i++)
		dst[i] = BigREAL64(src[i]);
	return true;
}

bool LittleREAL64Vector(std::vector<real64> &dst, std::vector<real64> &src)
{
	if(dst.size() != src.size())
		return false;

	int size = dst.size();
	for(int i = 0; i<size; i++)
		dst[i] = LittleREAL64(src[i]);
	return true;
}


void InitEndian( void )
{
	//clever little trick from Quake 2 to determine the endian
	//of the current system without depending on a preprocessor define

	uint8 SwapTest[2] = { 1, 0 };

	if( *(int16 *) SwapTest == 1 )
	{
		//little endian
		BigEndianSystem = false;

		//set func point32ers to correct funcs
		BigINT16 = INT16_Swap;
		LittleINT16 = INT16_NoSwap;
		BigUINT16 = UINT16_Swap;
		LittleUINT16 = UINT16_NoSwap;
		BigINT32 = INT32_Swap;
		LittleINT32 = INT32_NoSwap;
		BigUINT32 = UINT32_Swap;
		LittleUINT32 = UINT32_NoSwap;
		BigREAL32 = REAL32_Swap;
		LittleREAL32 = REAL32_NoSwap;
		BigREAL64 = REAL64_Swap;
		LittleREAL64 = REAL64_NoSwap;
	}
	else
	{
		//big endian
		BigEndianSystem = true;

		BigINT16 = INT16_NoSwap;
		LittleINT16 = INT16_Swap;
		BigUINT16 = UINT16_NoSwap;
		LittleUINT16 = UINT16_Swap;
		BigINT32 = INT32_NoSwap;
		LittleINT32 = INT32_Swap;
		BigUINT32 = UINT32_NoSwap;
		LittleUINT32 = UINT32_Swap;
		BigREAL32 = REAL32_NoSwap;
		LittleREAL32 = REAL32_Swap;
		BigREAL64 = REAL64_NoSwap;
		LittleREAL64 = REAL64_Swap;
	}
}


/*
int32 main(int32 argc, char** argv)
{
	InitEndian();

	int16 array16[4] = {0x0102, 0x0203, -0x0304, -0x0405};
	int16 test16[4];
	printf("Array INT16: %d %d %d %d\n", array16[0], array16[1], array16[2], array16[3]);

	LittleINT16Array(test16, array16, 4);
	printf("Little INT16: %d %d %d %d\n", test16[0], test16[1], test16[2], test16[3]);
	BigINT16Array(test16, array16, 4);
	printf("Big INT16: %d %d %d %d\n", test16[0], test16[1], test16[2], test16[3]);

	uint16 uarray16[4] = {0x0102, 0x0203, 0x0304, 0x0405};
	uint16 utest16[4];
	printf("Array UINT16: %d %d %d %d\n", uarray16[0], uarray16[1], uarray16[2], uarray16[3]);

	LittleUINT16Array(utest16, uarray16, 4);
	printf("Little UINT16: %d %d %d %d\n", utest16[0], utest16[1], utest16[2], utest16[3]);
	BigUINT16Array(utest16, uarray16, 4);
	printf("Big UINT16: %d %d %d %d\n", utest16[0], utest16[1], utest16[2], utest16[3]);

	int32 array32[4] = {0x01020304, 0x02030405, -0x03040506, -0x04050607};
	int32 test32[4];
	printf("Array INT32: %d %d %d %d\n", array32[0], array32[1], array32[2], array32[3]);

	LittleINT32Array(test32, array32, 4);
	printf("Little INT32: %d %d %d %d\n", test32[0], test32[1], test32[2], test32[3]);
	BigINT32Array(test32, array32, 4);
	printf("Big INT32: %d %d %d %d\n", test32[0], test32[1], test32[2], test32[3]);

	uint32 uarray32[4] = {0x01020304, 0x02030405, 0x03040506, 0x04050607};
	uint32 utest32[4];
	printf("Array INT32: %d %d %d %d\n", uarray32[0], uarray32[1], uarray32[2], uarray32[3]);

	LittleUINT32Array(utest32, uarray32, 4);
	printf("Little UINT32: %d %d %d %d\n", utest32[0], utest32[1], utest32[2], utest32[3]);
	BigUINT32Array(utest32, uarray32, 4);
	printf("Big UINT32: %d %d %d %d\n", utest32[0], utest32[1], utest32[2], utest32[3]);


	real32 farray32[4] = {-1, 1, 0.5589693, 7834.99};
	real32 ftest32[4];
	printf("Array FLOAT32: %f %f %f %f\n", farray32[0], farray32[1], farray32[2], farray32[3]);

	LittleFLOAT32Array(ftest32, farray32, 4);
	printf("Little FLOAT32: %f %f %f %f\n", ftest32[0], ftest32[1], ftest32[2], ftest32[3]);
	BigFLOAT32Array(ftest32, farray32, 4);
	printf("Big FLOAT32: %f %f %f %f\n", ftest32[0], ftest32[1], ftest32[2], ftest32[3]);

	float64 farray64[4] = {-1, 1, 0.5589693, 7834.99};
	float64 ftest64[4];
	printf("Array FLOAT64: %f %f %f %f\n", farray64[0], farray64[1], farray64[2], farray64[3]);

	LittleFLOAT64Array(ftest64, farray64, 4);
	printf("Little FLOAT64: %f %f %f %f\n", ftest64[0], ftest64[1], ftest64[2], ftest64[3]);
	BigFLOAT64Array(ftest64, farray64, 4);
	printf("Big FLOAT64: %f %f %f %f\n", ftest64[0], ftest64[1], ftest64[2], ftest64[3]);

	return 0;
}
*/

}
