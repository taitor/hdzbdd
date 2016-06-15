/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#ifndef _bitvector_H_
#define _bitvector_H_

#include "typedefbv.h"
#include "mman.h"

#ifdef __SSE4_2__
#include <smmintrin.h>
#define POPCOUNT(x) _mm_popcnt_u64(x)
#else
#define POPCOUNT(x) _mm_popcnt_u64(x)
//#define POPCOUNT(x) popcount(x)
#endif

#define SDARRAY_RANK1 (1<<0)
//#define SDARRAY_RANK0 (1<<1)
//#define SDARRAY_RR (1<<1)
#define SDARRAY_SPARSE (1<<1)
#define SDARRAY_SELECT1 (1<<2)
#define SDARRAY_SELECT0 (1<<3)
//#define SDARRAY_LENHUF (1<<4)
//#define SDARRAY_NOBUF (1<<5)
//#define SDARRAY_COMPRANK (1<<6)
//#define SDARRAY_COMPPTR (1<<7)
//#define SDARRAY_SUC (1<<8)
#define SDARRAY_RRR (1<<4)
#define SDARRAY_HUFFMAN (1<<5)


#ifndef _BITVEC_T_
#define _BITVEC_T_
typedef u64 bitvec_t;
#endif

typedef struct bitvector {
  i64 n; // ベクトルの長さ
  i64 m; // 1の数
  i64 k; // 整数のバイト数
  i64 size; // 索引サイズ (ベクトルは含まない)
  bitvec_t *buf; // ベクトル
  int opt; // サポートする操作
  i64 ml[2],ms[2]; // selectの表の要素数
  i64 rrr; // rankの中ブロックのサイズ

  i64 N; // 確保した領域の長さ

// for select 1,0
  uchar *lp[2];
  uchar *sl[2];
  word *ss[2];
  i64 *p[2];

// for rank
  uchar *rl;
  word *rs;

// for pointers to compressed blocks
  uchar *pl;
  word *ps;

// for sparsearray
  int low_width;
  bitvec_t *low;

  int (*getbit)(struct bitvector *da, i64 i);
  void (*setbit)(struct bitvector *da, i64 i, int x);
  u64 (*getbits)(struct bitvector *da, i64 i, int d);
  void (*setbits)(struct bitvector *da, i64 i, int d, u64 x);
  void (*push)(struct bitvector *da, int x);

  i64 (*rank)(struct bitvector *da, i64 i, int c);
  i64 (*select)(struct bitvector *da, i64 i, int c);
  i64 (*succ)(struct bitvector *da, i64 i, int c);
  i64 (*pred)(struct bitvector *da, i64 i, int c);


} bitvector;

#ifdef __cplusplus
extern "C" {
#endif
void bitvector_make_selecttbl(void);
bitvector *bitvector_new(i64 n);
i64 bitvector_makeindex(bitvector *da, int blocksize, int opt);
i64 bitvector_write(bitvector *da, FILE *f);
bitvector *bitvector_read(uchar **map);
bitvector *bitvector_read_from_file(uchar *filename);
void bitvector_free(bitvector *da);
#ifdef __cplusplus
}
#endif
#endif // _DENSEARRAY_H_
