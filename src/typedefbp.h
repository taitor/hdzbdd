#include <stdint.h>

#ifndef _TYPEDEFBP_H_
#define _TYPEDEFBP_H_

#ifndef DEF_byte
#define DEF_byte
typedef unsigned char byte;
#endif
#ifndef DEF_word
#define DEF_word
typedef unsigned short word;
#endif
#ifndef DEF_dword
#define DEF_dword
typedef unsigned int dword;
#endif

#ifndef DEF_i64
#define DEF_i64
typedef long int i64;
#endif
#ifndef DEF_u64
#define DEF_u64
typedef unsigned long int u64;
#endif

typedef dword pb;
#define logD 5
//typedef u64 pb;
//#define logD 6
#define D (1<<logD)


//typedef byte sb_type;
typedef word sb_type;

#ifdef INDEX64
  typedef i64 mb_type;
#else
//  typedef dword mb_type;
  typedef int mb_type;
#endif

#endif
