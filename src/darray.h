#ifndef _DARRAY_H_
#define _DARRAY_H_

#include "typedefbp.h"

#define OPT_NO_RANK (1<<30)
#define OPT_DARRAY_COMP_RLE (1<<29)
#define USE_CACHE

#define logR 16
#define R1 (1<<logR)
#define logRR 10
//#define logRR 8
#define RR (1<<logRR)
#define logRRR 7
#define RRR (1<<logRRR)

typedef struct {
  i64 n,m;
  i64 size;
  pb *buf;
  pb *bufc;
#ifdef INDEX64
  i64 *lp;
  i64 *sl;
#else
  dword *lp;
  dword *sl;
#endif
  word *ss;
  dword *p;  // if (n < 2^{32} * L(=2^{10})), dword is enough.

  //dword *rl;
  i64 *rl;
  word *rm;
  byte *rs;

#ifdef INDEX64
  i64 *pcl; // rle_compressed_vector の大ブロックへのポインタ
#else
  dword *pcl;
#endif
  word *pcm; //  rle_compressed_vector の中ブロックへのポインタ


  int opt;
  i64 idx_size;
} darray;

#ifdef __cplusplus
extern "C" {
#endif
int setbit(pb *B, i64 i,i64 x);
int setbits(pb *B, i64 i, i64 d, i64 x);
int getbit(pb *B, i64 i);
dword getbits(pb *B, i64 i, i64 d);
unsigned int popcount(pb x);
int decodegamma(pb *B,i64 p,i64 *ans);

void darray_make_selecttbl(void);
int darray_construct(darray *da, i64 n, pb *buf,int opt);
int darray_comp_rle(darray *da);
void darray_decode_block_rle(darray *da, i64 s, i64 t, pb *out);
i64 darray_decode_block_rle2(darray *da, i64 s, i64 t, i64 **out);

i64 darray_select(darray *da, i64 i,int f);
i64 darray_rank(darray *da, i64 i);
i64 darray_getbit(darray *da, i64 i);
i64 darray_pat_construct(darray *da, i64 n, pb *buf, i64 k, pb pat, int opt);
i64 darray_pat_select(darray *da, i64 i, pb (*getpat)(pb *));
i64 darray_pat_rank(darray *da, i64 i, pb (*getpat)(pb *));

i64 darray_select_bsearch(darray *da, i64 i, pb (*getpat)(pb *));
    
//    added by Lee
void darray_free(darray *da);

void block_cache_init(int size);
void block_cache_free(void);


extern i64 block_cache_hit, block_cache_miss;

#ifdef __cplusplus
}
#endif

#endif
