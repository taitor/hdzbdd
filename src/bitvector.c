/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/

// supports select1 and select0

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/timeb.h>
#include "bitvector.h"

#define ID_BITVECTOR 0x12

#define logQ 6
#define Q (1<<logQ) // 小ブロック(1ワード)のサイズ
//#define logM 5
//#define M (1<<logM)
//#define logP 8
//#define P (1<<logP)
#define logLL 16    // size of word
#define LL (1<<logLL) // selectの索引を変えるときの閾値
#define logLLL 5
//#define logLLL 2
#define LLL (1<<logLLL) // selectの中ブロックのサイズ
//#define logL 10
//#define logL (logLL-3)
#define logL (logLL-1-5)
#define L (1<<logL) // selectの大ブロックのサイズ

#define logRR 16
#define RR (1<<logRR)  // rankの大ブロックのサイズ

#ifndef min
 #define min(x,y) ((x)<(y)?(x):(y))
#endif


static int msize=0;
#define mymalloc(p,n,f) {p = malloc((n)*sizeof(*p)); msize += (f)*(n)*sizeof(*p); /* if (f) printf("malloc %d bytes at line %d total %d\n",(n)*sizeof(*p),__LINE__,msize);  */ if ((p)==NULL) {printf("not enough memory (%d bytes) in line %d\n",msize,__LINE__); exit(1);};}
#define myrealloc(p,n,t) {p = (t *)realloc((p),(n)*sizeof(*p)); if ((p)==NULL) {printf("not enough memory at line %d\n",__LINE__); exit(1);};}

static int blog(i64 x)
{
i64 l;
  l = -1;
  while (x>0) {
    x>>=1;
    l++;
  }
  return l;
}


static void writeuint(int k,u64 x,FILE *f)
{
  int i;
  for (i=k-1; i>=0; i--) {
    fputc(x & 0xff,f);
    x >>= 8;
  }
}

static u64 getuint(uchar *s, i64 i, i64 w)
{
  u64 x;
  i64 j;
  s += i*w;
  x = 0;
  for (j=0; j<w; j++) {
    x += ((u64)(*s++)) << (j*8);
  }
  return x;
}

static void putuint(uchar *s, i64 i, i64 x, i64 w)
{
  i64 j;
  s += i*w;
  for (j=0; j<w; j++) {
    *s++ = x & 0xff;
    x >>= 8;
  }
}


static int setbit(bitvec_t *B, i64 i,int x)
{
  i64 j,l;

  j = i / Q;
  l = i & (Q-1);
  if (x==0) B[j] &= (~(1L<<(Q-1-l)));
  else if (x==1) B[j] |= (1L<<(Q-1-l));
  else {
    printf("error setbit x=%d\n",x);
    exit(1);
  }
  return x;
}

static void bitvector_setbit(bitvector *da, i64 i,int x)
{
  i64 j,l;
  bitvec_t *B;

  B = da->buf;
  j = i / Q;
  l = i & (Q-1);
  if (x==0) B[j] &= (~(1L<<(Q-1-l)));
  else if (x==1) B[j] |= (1L<<(Q-1-l));
  else {
    printf("error setbit x=%d\n",x);
    exit(1);
  }
}

static void bitvector_push(bitvector *da, int x)
{
  i64 i, m;
  i64 new_N;
  if (da->n >= da->N) {
    new_N = da->N + 1024;
    m = (new_N+Q-1)/Q - (da->N+Q-1)/Q;
    myrealloc(da->buf, (new_N+Q-1)/Q, bitvec_t);
    for (i=(da->N+Q-1)/Q; i<(new_N+Q-1)/Q; i++) da->buf[i] = 0;
    da->size += m * sizeof(da->buf[0]);
    da->N = new_N;
  }
  bitvector_setbit(da, da->n, x);
  da->n++;
}


static void setbits(bitvec_t *B, i64 i, int d, u64 x)
{
  u64 y,m;
  int d2;
#if 0
  int j;

  for (j=0; j<d; j++) {
    setbit(B,i+j,(x>>(d-j-1))&1);
  }
#else
  B += (i>>logQ);
  i &= (Q-1);

  while (i+d > Q) {
    d2 = Q-i; // x の上位 d2 ビットを格納
    y = x >> (d-d2);
    m = (1LL<<d2)-1;
    *B = (*B & (~m)) | y;
    B++;  i=0;
    d -= d2;
    x &= (1LL<<d)-1; // x の上位ビットを消去
  }
  m = (1LL<<d)-1;
  y = x << (Q-i-d);
  m <<= (Q-i-d);
  *B = (*B & (~m)) | y;
#endif
}

static void bitvector_setbits(bitvector *da, i64 i, int d, u64 x)
{
  u64 y,m;
  int d2;
  bitvec_t *B;
  
  B = da->buf;
  B += (i>>logQ);
  i &= (Q-1);

  while (i+d > Q) {
    d2 = Q-i; // x の上位 d2 ビットを格納
    y = x >> (d-d2);
    m = (1<<d2)-1;
    *B = (*B & (~m)) | y;
    B++;  i=0;
    d -= d2;
    x &= (1<<d)-1; // x の上位ビットを消去
  }
  m = (1<<d)-1;
  y = x << (Q-i-d);
  m <<= (Q-i-d);
  *B = (*B & (~m)) | y;

}

static int getbit(bitvec_t *B, i64 i)
{
  i64 j,l;

  j = i >> logQ;
  l = i & (Q-1);
  return (B[j] >> (Q-1-l)) & 1;
}

static int bitvector_getbit(bitvector *da, i64 i)
{
  i64 j,l;

  j = i >> logQ;
  l = i & (Q-1);
  return (da->buf[j] >> (Q-1-l)) & 1;
}

static bitvec_t getbits(bitvec_t *B, i64 i, int d)
{
  bitvec_t x,z;

  if (d == 0) return 0;
  B += (i >>logQ);
  i &= (Q-1);
  if (i+d <= Q) {
    x = B[0];
    x <<= i;
    x >>= (Q-d);  // Q==64, d==0 だと動かない
  } else {
    x = B[0] << i;
    x >>= Q-d;
    z = B[1] >> (Q-(i+d-Q));
    x += z;
  }
#if 0
  if (d < Q) x &= ((1LL << d) -1);
#endif
  return x;
}

static bitvec_t bitvector_getbits(bitvector *da, i64 i, int d)
{
  bitvec_t x,z;
  bitvec_t *B;

  if (d == 0) return 0;
  B = da->buf;
  B += (i >>logQ);
  i &= (Q-1);
  if (i+d <= Q) {
    x = B[0];
    x <<= i;
    x >>= (Q-d);  // Q==64, d==0 だと動かない
  } else {
    x = B[0] << i;
    x >>= Q-d;
    z = B[1] >> (Q-(i+d-Q));
    x += z;
  }
#if 0
  if (d < Q) x &= (1LL << d) -1;
#endif
  return x;
}

static unsigned int selecttbl[8*256];

static unsigned int popcount(bitvec_t x)
{
  bitvec_t r;
//  __uint128_t rr;
  r = x;
  r = ((r & 0xaaaaaaaaaaaaaaaa)>>1) + (r & 0x5555555555555555);
  r = ((r & 0xcccccccccccccccc)>>2) + (r & 0x3333333333333333);
  r = ((r>>4) + r) & 0x0f0f0f0f0f0f0f0f;
//  r = (r>>8) + r;
//  r = (r>>16) + r;
//  r = ((r>>32) + r) & 127;

  r *= 0x0101010101010101;
  r >>= 64-8;
//  printf("r1 %016lx\n",r);
//  rr = r;
//  rr *= 0x0101010101010101;
//  r = rr >> 56;
//  printf("r2 %016lx\n",r);
//  r &= 0xff;
  return r;
}

#if 0
static int select_sub(bitvec_t x, int s)
{
  bitvec_t r, w;
  __uint128_t rr;
  int p;

  p = 0;

  r = x;
  r = ((r & 0xaaaaaaaaaaaaaaaa)>>1) + (r & 0x5555555555555555);
  r = ((r & 0xcccccccccccccccc)>>2) + (r & 0x3333333333333333);
  r = ((r>>4) + r) & 0x0f0f0f0f0f0f0f0f;

  rr = r;
  rr *= 0x0101010101010101;
  w = rr >> 56;
  w += 0x8080808080808080;
  w -= s * 0x0101010101010101;
  w &= 0x8080808080808080;
  if ((w & 0xffffffff00000000) == 0) {
    p += 32;
    x <<= 32;
  }
  if ((w & 0xffff000000000000) == 0) {
    p += 16;
    x <<= 16;
  }
  if ((w & 0xff00000000000000) == 0) {
    p += 8;
    x <<= 8;
  }
  return p;
}
#endif

///////////////////////////////////////////////////////////
// rank(i, c)
//   returns the number of c's in [0,i]
//   if i < 0, it returns 0
//   if i >= n, it returns total number of ones
///////////////////////////////////////////////////////////
static i64 bitvector_rank(bitvector *da, i64 i, int c)
{
  i64 r,j;
  bitvec_t *p;
  i64 rrr;

#ifdef DEBUG
   if (c < 0 || c > 1) {
    printf("bitvector_rank: error c = %d\n", c);
    exit(1);
  }
#endif
#if 1
  if (i < 0) return 0;
  if (i >= da->n) i = da->n-1;
#else
  if (i < 0 || i >= da->n) {
    printf("bitvector_rank: i=%ld\n", i);
  }
  if (i < 0) return 0;
  if (i >= da->n) i = da->n-1;
#endif
  
  rrr = da->rrr;
  r = getuint(da->rl,i>>logRR,da->k) + da->rs[i/rrr];
  p = da->buf + ((i & (~(rrr-1))) >> logQ);
  j = i & (rrr-1);
  while (j >= Q) {
    r += POPCOUNT(*p++);
    j -= Q;
  }
  r += POPCOUNT(*p >> (Q-1-j));

  if (c == 0) {
    r = (i+1) - r;
  }

  return r;
}

static i64 bitvector_rank_and_bit(bitvector *da, i64 i, int *c)
{
  i64 r,j;
  bitvec_t *p,x;
  i64 rrr;

  if (i < 0) return 0;
  if (i >= da->n) i = da->n-1;

  rrr = da->rrr;
  r = getuint(da->rl,i>>logRR,da->k) + da->rs[i/rrr];
  p = da->buf + ((i & (~(rrr-1))) >> logQ);
  j = i & (rrr-1);
  while (j >= Q) {
    r += POPCOUNT(*p++);
    j -= Q;
  }
  x = *p >> (Q-1-j);
  r += POPCOUNT(x);
  *c = (int)(x & 1);

  return r;
}

#if 0
static i64 bitvector_rank0(bitvector *da, i64 i)
{
  if (i < 0) return 0;
  if (i >= da->n) i = da->n-1;
  return i+1 - bitvector_rank(da,i, 1);
}

static i64 bitvector_select0_naive(bitvector *da, i64 i)
{
  i64 l,r,m;

  l = 0; r = da->n-1;
  while (l <= r) {
    m = (l+r)/2;
    if (i <= bitvector_rank0(da,m)) r = m-1;  else l = m+1;
  }
  return l;
}
#endif

static i64 bitvector_select_by_rank(bitvector *da, i64 ith, int c)
{
  bitvec_t x, *buf;
  
  i64 j,k;
  i64 ofs, ofs2;
  i64 r, r2;
  int rr;
  i64 d,runlen;

  i64 ll, rl, ml, pl; // 大ブロックの2分探索用
  i64 p; // 答え
  i64 ii;
  i64 rrr;
  bitvec_t *q;
  i64 p0;
  i64 l0;

  rrr = da->rrr;

  ii = ith;

  ll = 0;  rl = (da->n-1) >> logRR;
  pl = ll;  r2 = 0;
  while (ll <= rl) {
    ml = (ll+rl) >> 1;
    r = getuint(da->rl,ml,da->k);
    if (c == 0) {r = (ml<<logRR) - r;} // select0 の場合
    if (r < ith) {
      pl = ml;
      ll = ml+1;
      r2 = r;
    } else {
      rl = ml-1;
    }
  }
  ith -= r2; // 大ブロック内のランク
  ll = (pl<<logRR)/rrr;
  rl = min(((pl+1)<<logRR)/rrr-1, (da->n-1)/rrr);

  pl = ll;  r2 = 0;
  l0 = ll;
  while (ll <= rl) {
    ml = (ll+rl) >> 1;
    r = da->rs[ml];
    if (c == 0) {r = ((ml-l0)*rrr) - r;} // select0 の場合
    if (r < ith) {
      pl = ml;
      ll = ml+1;
      r2 = r;
    } else {
      rl = ml-1;
    }
  }
  ith -= r2; // 小ブロック内のランク

  p = pl * rrr; // ith を含む小ブロックの最初のビット位置
  rr = 0;

  p0 = p;

  q = &(da->buf[p>>logQ]);

  while (1) {
    x = *q;
    if (c == 0) x = ~x;
    rr = POPCOUNT(x);
    if (rr >= ith) break;
    ith -= rr;
    p += Q;
    q++;
  }
      
  x = *q;
  if (c == 0) x = ~x;
  while (1) {
    rr = POPCOUNT(x >> (Q-8));
    if (rr >= ith) break;
    ith -= rr;
    p += 8;
    x <<= 8;
  }
  p += selecttbl[((ith-1)<<8)+(x>>(Q-8))];

  if (p >= p0 + rrr) {
    printf("i = %ld p = %ld p0 = %ld\n", ii, p, p0);
  }

#if 0
  if (p != bitvector_select0_naive(da, ii)) {
    printf("select0: i = %ld p = %ld naive = %ld\n", 
      ii, p, bitvector_select0_naive(da, ii));
    exit(1);
  }
#endif
  return p;
  
}

///////////////////////////////////////////////////////////
// select(i, f)
//   returns the position of i-th f. The position is in [0,n-1].
//   if i == 0, it returns -1
//   if f is larger than the number of f's in the vector, returns n
///////////////////////////////////////////////////////////
static i64 bitvector_select(bitvector *da, i64 i,int f)
{
  
  //  dword *s;
  i64 p,r;
  i64 il;
  i64 rr;
  bitvec_t x;
  bitvec_t *q;
  i64 m;
  i64 stmp;

  if (i == 0) return -1;
  m = (f == 1) ? da->m : da->n - da->m;
  if (i > m) return da->n;

  if (m == da->n) return i-1; // sada 2015/1/16

  if (f == 0 && (da->opt & SDARRAY_SELECT0) == 0) {
    return bitvector_select_by_rank(da,i,0);
  }
  if (f == 1 && (da->opt & SDARRAY_SELECT1) == 0) {
    return bitvector_select_by_rank(da,i,1);
  }

  i--;

  il = da->p[f][i>>logL];
  if (il < 0) {
    il = -il-1;
    p = getuint(da->sl[f],(il<<logL)+(i & (L-1)),da->k);
  } else {
    p = getuint(da->lp[f],i>>logL,da->k);
    p += da->ss[f][(il<<(logL-logLLL))+(i & (L-1))/LLL];
    r = i - (i & (LLL-1));

    q = &(da->buf[p>>logQ]);

    if (f == 1) {
      rr = p & (Q-1);
      r -= POPCOUNT(*q >> (Q-1-rr));
      p = p - rr;
      
      while (1) {
        rr = POPCOUNT(*q);
        if (r + rr >= i) break;
        r += rr;
        p += Q;
        q++;
      }
      
      x = *q;
      while (1) {
        rr = POPCOUNT(x >> (Q-8));
        if (r + rr >= i) break;
        r += rr;
        p += 8;
        x <<= 8;
      }
      p += selecttbl[((i-r-1)<<8)+(x>>(Q-8))];
    } else {
      rr = p & (Q-1);
      r -= POPCOUNT((~(*q))  >> (Q-1-rr));
      p = p - rr;
      
      while (1) {
        rr = POPCOUNT(~(*q));
        if (r + rr >= i) break;
        r += rr;
        p += Q;
        q++;
      }
      
      x = ~(*q);
      while (1) {
        rr = POPCOUNT(x >> (Q-8));
        if (r + rr >= i) break;
        r += rr;
        p += 8;
        x <<= 8;
      }
      p += selecttbl[((i-r-1)<<8)+(x>>(Q-8))];
    }
  }
  
  return p;
}

static i64 bitvector_pred_naive(bitvector *da, i64 x, int c)
{
  i64 s0;
  i64 r;
  i64 m;
  i64 s;

//  r = bitvector_rank(da, x, c);
//  s = bitvector_select(da, r, c);
  r = da->rank(da, x, c);
  s = da->select(da, r, c);
  
  return s;
}

///////////////////////////////////////////////////////////
// pred(x, c)
//   returns max{i | B[i]=c and i <= x}
//   if there is no such i, returns -1
//   (if x >= n, error, or the position of last c?)
///////////////////////////////////////////////////////////
static i64 bitvector_pred(bitvector *da, i64 x, int c)
{
  i64 s0;
  i64 r;
  i64 m;
  i64 s;

  bitvec_t *q, z;
  int i,d, rr;

  if (x < 0) return -1;
  if (x >= da->n) x = da->n-1;

  s = x;

  q = &(da->buf[x>>logQ]);
  z = *q--;
  if (c == 0) z = ~z;

  d = x & (Q-1);
  z >>= (Q-1-d);
  z <<= (Q-1-d); // 下位 Q-1-d ビットを 0 にする
  s &= ~(Q-1);

  for (i=0; i < 4; i++) {
    rr = POPCOUNT(z);
    if (rr > 0) break;
    s -= Q;
    if (s < 0) return -1;
    z = *q--;
    if (c == 0) z = ~z;
  }

  if (i < 4) { // scan して見つかった
    // 一番右の 1 のbitを見つける
    if (z & 0x00000000ffffffffULL) {
      z <<= 32;
      s += 32;
    }
    if (z & 0x0000ffff00000000ULL) {
      z <<= 16;
      s += 16;
    }
    if (z & 0x00ff000000000000ULL) {
      z <<= 8;
      s += 8;
    }
    if (z & 0x0f00000000000000ULL) {
      z <<= 4;
      s += 4;
    }
    if (z & 0x3000000000000000ULL) {
      z <<= 2;
      s += 2;
    }
    if (z & 0x4000000000000000ULL) {
      z <<= 1;
      s += 1;
    }
#if 0
    if (s != bitvector_pred_naive(da, x, c)) {
      printf("x = %ld s = %ld ans = %ld\n", x, s, bitvector_pred_naive(da, x, c));
    }
#endif
  } else {
    s = bitvector_pred_naive(da, x, c);
  }
  return s;
}


// rank(min{i | B[i]=c and i >= x})
static i64 bitvector_succ_rank(bitvector *da, i64 x)
{
  i64 r;

  r = (x > 0) ? bitvector_rank(da, x-1, 1) + 1 : 1;

  return r;
}

static i64 bitvector_succ_naive(bitvector *da, i64 x, int c)
{
  i64 r;
  i64 m;
  i64 s;

//  r = (x > 0) ? bitvector_rank(da, x-1, c) + 1 : 1;
  r = (x > 0) ? da->rank(da, x-1, c) + 1 : 1;
  m = (c == 1) ? da->m : da->n - da->m;
//  s = (r <= m) ? bitvector_select(da, r, c) : da->n;
  s = (r <= m) ? da->select(da, r, c) : da->n;
  return s;
}

///////////////////////////////////////////////////////////
// succ(x, c)
//   returns min{i | B[i]=c and i >= x}
//   if there is no such i, returns n
///////////////////////////////////////////////////////////
static i64 bitvector_succ(bitvector *da, i64 x, int c)
{
  i64 s0;
  i64 r;
  i64 m;
  i64 s;

  bitvec_t *q, z;
  int i,d, rr;

  if (x >= da->n) return da->n;
  if (x < 0) x = 0; // debug

  s = x;

  q = &(da->buf[x>>logQ]);
  z = *q++;
  if (c == 0) z = ~z;

  d = x & (Q-1);
  z <<= d;
  z >>= d; // 上位 d ビットを 0 にする
  s &= ~(Q-1);

  for (i=0; i < 4; i++) {
    rr = POPCOUNT(z);
    if (rr > 0) break;
    s += Q;
    if (s >= da->n) return da->n;
    z = *q++;
    if (c == 0) z = ~z;
  }

  if (i < 4) { // scan して見つかった
    // 一番左の 1 のbitを見つける
    while (1) {
      rr = POPCOUNT(z >> (Q-8));
      if (rr > 0) break;
      s += 8;
      z <<= 8;
    }
    s += selecttbl[((1-1)<<8)+(z>>(Q-8))];
    if (s > da->n) s = da->n;

  } else {
    s = bitvector_succ_naive(da, x, c);
  }
  return s;
}


static i64 sparsearray_rank(bitvector *sa, i64 i, int f)
{
  i64 d,x,w,y;
  i64 j;

  if (i < 0) return 0;
  if (i >= sa->n) i = sa->n-1;

  if (sa->m == 0) return 0; // sada 4/25

  d = sa->low_width;

  y = bitvector_select(sa,i>>d,0)+1;
  x = y - (i>>d);

  j = i - ((i>>d)<<d);

  while (1) {
    if (getbit(sa->buf,y)==0) break;
    w = getbits(sa->low,x*d,d);
    if (w >= j) {
      if (w == j) x++;
      break;
    }
    x++;
    y++;
  }
  if (f == 0) { // rank0
    x = (i+1) - x;
  }
  return x;
}

static int sparsearray_getbit(bitvector *sa, i64 i)
{
  i64 d,x,w,y;
  i64 j;
  int b;

  if (i < 0) return 0;
  if (i >= sa->n) i = sa->n-1;

  if (sa->m == 0) return 0;

  d = sa->low_width;

  y = bitvector_select(sa,i>>d,0)+1;
  x = y - (i>>d);
  j = i - ((i>>d)<<d);
  b = 0;

  while (1) {
    if (getbit(sa->buf,y)==0) break;
    w = getbits(sa->low,x*d,d);
    if (w >= j) {
      if (w == j) b = 1;
      break;
    }
    x++;
    y++;
  }
  return b;
}



static i64 sparsearray_select(bitvector *sa, i64 i, int f)
{
  i64 d,x;

  if (i == 0) return -1;
  if (i < 0) {
    printf("sparsearray_select: i=%ld\n",i);
    exit(1);
  }

  if (f == 1) { // select1
    if (i > sa->m) return sa->n;

    if (sa->m == sa->n) return i-1; // sada 2015/1/16

    d = sa->low_width;

    x = bitvector_select(sa,i,1) - (i-1);
    x <<= d;
    x += getbits(sa->low,(i-1)*d,d);
  } else { // select0 (slow!)
    i64 l,r,m;
    l = 0; r = sa->n-1;
    while (l <= r) {
      m = (l+r)/2;
      if (i <= sparsearray_rank(sa,m,0)) r = m-1;  else l = m+1;
    }
    x = l;
  }
  return x;

}

static i64 sparsearray_pred_rank(bitvector *sa, i64 i, i64 *rank, int f)
{
  i64 d,x,w,y;
  i64 j;
  i64 w0, y0;

  if (i < 0) return 0;
  if (i >= sa->n) i = sa->n-1;

  d = sa->low_width;

  y = bitvector_select(sa,i>>d,0)+1;
  y0 = y;
  x = y - (i>>d);

  j = i - ((i>>d)<<d);

  w0 = -1;
  while (1) {
    if (getbit(sa->buf,y)==0) break;
    w = getbits(sa->low,x*d,d);
    if (w >= j) {
      if (w == j) {
        w0 = w;
        x++;
      }
      break;
    }
    x++;
    y++;
    w0 = w;
  }
  *rank = x;

  if (y == y0) {
    return sparsearray_select(sa, x, 1);
  }

  y0 = bitvector_select(sa,x,1) - (x-1);

  return (y0<<d) + w0;
}

void bitvector_make_selecttbl(void)
{
  i64 i,x,r;
  bitvec_t buf[1];

  for (x = 0; x < 256; x++) {
    setbits(buf,0,8,x);
    for (r=0; r<8; r++) selecttbl[(r<<8)+x] = -1;
    r = 0;
    for (i=0; i<8; i++) {
      if (getbit(buf,i)) {
        selecttbl[(r<<8)+x] = i;
        r++;
      }
    }
  }

}

static void bitvector_notsupported(void)
{
  printf("bitvector: not supported.\n");
  exit(1);
}

bitvector *bitvector_new(i64 n)
{
  bitvector *da;
  i64 i;

  mymalloc(da, 1, 0);

  da->ntmp = da->mtmp = 0;

  da->n = n;
  da->size = 0;
  mymalloc(da->buf, (n+Q-1)/Q, 0);
  da->size += (n+Q-1)/Q * sizeof(da->buf[0]);
  for (i=0; i<(n+Q-1)/Q; i++) da->buf[i] = 0;

  da->N = n;


  da->rl = NULL;
  da->rs = NULL;
  da->lp[0] = da->lp[1] = NULL;
  da->sl[0] = da->sl[1] = NULL;
  da->ss[0] = da->ss[1] = NULL;
  da->p[0] = da->p[1] = NULL;

  da->low = NULL;
  da->low_width = 0;

  da->getbit = bitvector_getbit;
  da->setbit = bitvector_setbit;
  da->getbits = bitvector_getbits;
  da->setbits = bitvector_setbits;
  da->push = bitvector_push;


  da->rank = (i64 (*)(struct bitvector *, i64, int))bitvector_notsupported;
  da->select = (i64 (*)(struct bitvector *, i64, int))bitvector_notsupported;
  da->succ = (i64 (*)(struct bitvector *, i64, int))bitvector_notsupported;
  da->pred = (i64 (*)(struct bitvector *, i64, int))bitvector_notsupported;

  return da;
}

static void change_nm(bitvector *da, i64 new_n, i64 new_m)
{
i64 mm;
int new_d, d;
i64 current_size, new_size;
i64 i,p;
bitvec_t *buf, *low;
bitvec_t *new_buf, *new_low;
//  printf("change_nm: current n=%ld m=%ld d=%ld\n", da->ntmp, da->mtmp, da->low_width);
  mm = new_m;
  new_d = 0;
  while (mm < new_n) {
    mm <<= 1;
    new_d++;
  }
  da->size = sizeof(*da);
//  current_size = (2*da->m+Q-1)/Q+1;
  new_size = (2*new_m+Q-1)/Q+1;
  mymalloc(new_buf, new_size,1);
  for (i=0; i<new_size; i++) new_buf[i] = 0;
  da->size += new_size;

//  current_size = (da->low_width*da->m+Q-1)/Q+1;
  new_size = (new_d*new_m+Q-1)/Q+1;
  mymalloc(new_low, new_size,1);
  for (i=0; i<new_size; i++) new_low[i] = 0;
  da->size += new_size;

  buf = da->buf;
  low = da->low;
  d = da->low_width;
  p = 0; // position
  for (i=0; i<da->m; i++) {
    int c;
    i64 x;

    // obtain current value
    while (getbit(buf, p) == 0) {
      p++;
      if (p >= da->mtmp*2) {
        printf("change_nm: ??? %ld %ld\n", p, da->mtmp);
        exit(1);
      }
    }
    x = (p-i) << d;
    x += getbits(low, i*d, d);
    p++;
//    printf("%ld: %ld\n", i, x);

    // rewrite in the new buffer
    setbit(new_buf,(x>>new_d)+i,1);
    setbits(new_low,i*new_d,new_d,x & ((1<<new_d)-1));
  }

  free(da->buf);
  free(da->low);
  da->buf = new_buf;
  da->low = new_low;
  da->low_width = new_d;
  da->ntmp = new_n;
  da->mtmp = new_m;
//  printf("change_nm: new     n=%ld m=%ld d=%ld\n", da->ntmp, da->mtmp, da->low_width);
}

////////////////////////////////////////////////////////////////////////
// void bitvector_setbit_sparse(bitvector *da, i64 r, i64 x)
//   set a bit at position x
//   x >= 0
//   r >= 0 (rank of x)
////////////////////////////////////////////////////////////////////////
void bitvector_setbit_sparse(bitvector *da, i64 r, i64 x)
{
i64 new_n, new_m;
int d;

  new_m = da->mtmp;
  new_n = da->ntmp;

  if (r >= da->mtmp) {
    new_m = da->mtmp * 3/2;  if (new_m < r) new_m = r;
  }
#if 0
  if (x >= da->ntmp) {
    new_n = da->ntmp * 5/4;  if (new_n < x) new_n = x;
  }
#else
  if (x >= da->ntmp) {
    printf("x = %ld must be < n = %ld\n", x, da->ntmp);
    exit(1);
  }
#endif
  if (new_m > da->mtmp || new_n > da->ntmp) {
    change_nm(da, new_n, new_m);
    da->mtmp = new_m;
    da->ntmp = new_n;
  }
  d = da->low_width;
  setbit(da->buf,(x>>d)+r,1);
  setbits(da->low,r*d,d,x & ((1<<d)-1));
  da->m++; // #stored elements
#if 0
  if (x > da->n) da->n = x; // maximum element
#endif
}

////////////////////////////////////////////////////////////////////////
// mtmp is an estimated value
////////////////////////////////////////////////////////////////////////
bitvector *bitvector_new_sparse(i64 ntmp, i64 mtmp)
{
  bitvector *da;
  i64 i;

  mymalloc(da, 1, 0);
  //da->n = n;
  da->ntmp = ntmp;
  da->mtmp = mtmp;
  da->n = ntmp;
  da->m = 0;

  da->size = 0;
//  mymalloc(da->buf, (n+Q-1)/Q, 0);
//  da->size += (n+Q-1)/Q * sizeof(da->buf[0]);
//  for (i=0; i<(n+Q-1)/Q; i++) da->buf[i] = 0;
//  da->low = NULL;
//  da->low_width = 0;

  da->low_width = 0;
  da->buf = NULL;
  da->low = NULL;
  change_nm(da, ntmp, mtmp);

  da->rl = NULL;
  da->rs = NULL;
  da->lp[0] = da->lp[1] = NULL;
  da->sl[0] = da->sl[1] = NULL;
  da->ss[0] = da->ss[1] = NULL;
  da->p[0] = da->p[1] = NULL;


  da->getbit = (int (*)(struct bitvector *da, i64 i))bitvector_notsupported;
  da->setbit = (void (*)(struct bitvector *da, i64 i, int x))bitvector_notsupported;
  da->getbits = (u64 (*)(struct bitvector *da, i64 i, int d))bitvector_notsupported;
  da->setbits = (void (*)(struct bitvector *da, i64 i, int d, u64 x))bitvector_notsupported;

  da->rank = (i64 (*)(struct bitvector *, i64, int))bitvector_notsupported;
  da->select = (i64 (*)(struct bitvector *, i64, int))bitvector_notsupported;
  da->succ = (i64 (*)(struct bitvector *, i64, int))bitvector_notsupported;
  da->pred = (i64 (*)(struct bitvector *, i64, int))bitvector_notsupported;

  return da;
}
void bitvector_free(bitvector *da)
{
  int i;
  free(da->buf);
  if (da->rl) free(da->rl);
  if (da->rs) free(da->rs);
  for (i=0; i<2; i++) {
    if (da->lp[i]) free(da->lp[i]);
    if (da->sl[i]) free(da->sl[i]);
    if (da->ss[i]) free(da->ss[i]);
    if (da->p[i]) free(da->p[i]);
  }
  if (da->low) free(da->low);
  free(da);
}

#if 0
int bitvector_construct_set(bitvector *da, i64 i, int x)
{
  if (x > 0) setbit(da->buf,i,x);
  return 0;
}

i64 bitvector_construct_end(bitvector *da, ushort l, int opt)
{
  bitvector_construct(da, da->n, da->buf, opt);
//  return da->size;
  return da->size + sizeof(*da->buf) * ((da->n+Q-1)/Q);
}
#endif

static void bitvector_construct(bitvector *da, i64 n, i64 m, bitvec_t *buf, int blocksize, int opt)
{
  i64 i,j,k;
  i64 nl;
  i64 p,pp;
  i64 il,is,ml[2],ms[2];
  i64 r;
  i64 size;
//  dword *s;
  i64 p1,m1,p2,m2;
  i64 rrr;

i64 s1,s2,s3,s4;
  
#if 0
  i64 b;
  i64 freq[RRR+1];
#endif

  size = sizeof(bitvector);

  da->rrr = rrr = blocksize;

  if (L/LLL == 0) {
    printf("ERROR: L=%d LLL=%d\n",L,LLL);
    exit(1);
  }

  da->k = k = (blog(n+1)+1+8-1)/8;

  da->buf = buf;

  if ((opt & SDARRAY_SELECT1) && (m>0)) { // sada 2015/1/16
    nl = (m-1) / L + 1;
    mymalloc(da->lp[1],(nl+1)*k,1);  size += sizeof(*(da->lp[1]))*(nl+1)*k;
    for (il = 0; il < nl+1; il++) putuint(da->lp[1],il,0,k);

    mymalloc(da->p[1],nl+1,1);  size += sizeof(*(da->p[1]))*(nl+1);
    for (il = 0; il < nl+1; il++) da->p[1][il] = 0;

    for (r = 0; r < 2; r++) {
      ml[1] = ms[1] = 0;
      p1 = p2 = -1;
      m1 = m2 = -1;
      for (il = 0; il < nl; il++) {
        i64 p1tmp, m1tmp;
#if 0
        while (m1 < il*L) m1 += getbit(buf,++p1);
#else
        while (m1 < il*L - Q) {
          m1 += POPCOUNT(getbits(buf,p1+1, Q));
          p1 += Q;
        }
        while (m1 < il*L) m1 += getbit(buf,++p1);
#endif
        if (p1 >= n) printf("???1 p1 = %ld\n",p1);
        pp = p1;
        putuint(da->lp[1],il,pp,k);
        i = min((il+1)*L-1,m-1);
#if 0
        while (m2 < i) m2 += getbit(buf,++p2);
#else
        while (m2 < i - Q) {
          m2 += POPCOUNT(getbits(buf,p2+1, Q));
          p2 += Q;
        }
        while (m2 < i) m2 += getbit(buf,++p2);
#endif
        if (p2 >= n) printf("???2 p2 = %ld\n",p2);
        p = p2;
        if (p - pp >= LL) {
          if (r == 1) {
            for (is = 0; is < L; is++) {
              if (il*L+is >= m) break;
#if 0
              while (m1 < il*L+is) m1 += getbit(buf,++p1);
#else
              while (m1 < il*L+is - Q) {
                m1 += POPCOUNT(getbits(buf,p1+1, Q));
                p1 += Q;
              }
              while (m1 < il*L+is) m1 += getbit(buf,++p1);
#endif
              if (p1 >= n) printf("???3 p1 = %ld\n",p1);
              putuint(da->sl[1],ml[1]*L+is,p1,k);
            }
          }
          da->p[1][il] = -(ml[1]+1);
          ml[1]++;
        } else {
          if (r == 1) {
            for (is = 0; is < L/LLL; is++) {
              if (il*L+is*LLL >= m) break;
#if 0
              while (m1 < il*L+is*LLL) m1 += getbit(buf,++p1);
#else
              while (m1 < il*L+is*LLL - Q) {
                m1 += POPCOUNT(getbits(buf,p1+1, Q));
                p1 += Q;
              }
              while (m1 < il*L+is*LLL) m1 += getbit(buf,++p1);
#endif
              if (p1 >= n) printf("???4 p1 = %ld\n",p1);
              da->ss[1][ms[1]*(L/LLL)+is] = p1 - pp;
            }
          }
          da->p[1][il] = ms[1];
          ms[1]++;
        }
      }
      if (r == 0) {
        da->ml[1] = ml[1];  da->ms[1] = ms[1];
        mymalloc(da->sl[1],(ml[1]*L+1)*k,1);  size += sizeof(*(da->sl[1]))*(ml[1]*L+1)*k;
        for (il = 0; il < ml[1]*L+1; il++) putuint(da->sl[1],il,0,k);
        mymalloc(da->ss[1],ms[1]*(L/LLL)+1,1);
        for (il = 0; il < ms[1]*(L/LLL)+1; il++) da->ss[1][il] = 0;
        size += sizeof(*(da->ss[1]))*(ms[1]*(L/LLL)+1);
      }
    }
  } else {
    da->lp[1] = NULL;  da->p[1] = NULL;
    da->sl[1] = NULL;  da->ss[1] = NULL;
  }

  if (opt & SDARRAY_SELECT0) {
    i64 m0;
    m0 = n - m;
    nl = (m0-1) / L + 1;
    mymalloc(da->lp[0],(nl+1)*k,1);  size += sizeof(*(da->lp[0]))*(nl+1)*k;
    for (il = 0; il < nl+1; il++) putuint(da->lp[0],il,0,k);
    mymalloc(da->p[0],nl+1,1);  size += sizeof(*(da->p[0]))*(nl+1);
    for (il = 0; il < nl+1; il++) da->p[0][il] = 0;

    for (r = 0; r < 2; r++) {
      ml[0] = ms[0] = 0;
      p1 = p2 = -1;
      m1 = m2 = -1;
      for (il = 0; il < nl; il++) {
        while (m1 < il*L) m1 += 1-getbit(buf,++p1);
        if (p1 >= n) printf("???1 p1 = %ld\n",p1);
        pp = p1;
        putuint(da->lp[0],il,pp,k);
        i = min((il+1)*L-1,m0-1);
        while (m2 < i) m2 += 1-getbit(buf,++p2);
        if (p2 >= n) printf("???2 p2 = %ld\n",p2);
        p = p2;
        if (p - pp >= LL) {
          if (r == 1) {
            for (is = 0; is < L; is++) {
              if (il*L+is >= m0) break;
              while (m1 < il*L+is) m1 += 1-getbit(buf,++p1);
              if (p1 >= n) printf("???3 p1 = %ld\n",p1);
              putuint(da->sl[0],ml[0]*L+is,p1,k);
            }
          }
          da->p[0][il] = -(ml[0]+1);
          ml[0]++;
        } else {
          if (r == 1) {
            for (is = 0; is < L/LLL; is++) {
              if (il*L+is*LLL >= m0) break;
              while (m1 < il*L+is*LLL) m1 += 1-getbit(buf,++p1);
              if (p1 >= n) printf("???4 p1 = %ld\n",p1);
              da->ss[0][ms[0]*(L/LLL)+is] = p1 - pp;
            }
          }
          da->p[0][il] = ms[0];
          ms[0]++;
        }
      }
      if (r == 0) {
        da->ml[0] = ml[0];  da->ms[0] = ms[0];
        mymalloc(da->sl[0],(ml[0]*L+1)*k,1);  size += sizeof(*(da->sl[0]))*(ml[0]*L+1)*k;
        for (il = 0; il < ml[0]*L+1; il++) putuint(da->sl[0],il,0,k);
        mymalloc(da->ss[0],ms[0]*(L/LLL)+1,1);
        for (il = 0; il < ms[0]*(L/LLL)+1; il++) da->ss[0][il] = 0;
        size += sizeof(*(da->ss[0]))*(ms[0]*(L/LLL)+1);
      }
    }
  } else {
    da->lp[0] = NULL;  da->p[0] = NULL;
    da->sl[0] = NULL;  da->ss[0] = NULL;
  }


// rank index
  if (opt & SDARRAY_RANK1) {
    mymalloc(da->rl,(n+RR-1)/RR*k,1);  
      size += sizeof(*(da->rl))*((n+RR-1)/RR)*k;
    mymalloc(da->rs,(n+rrr-1)/rrr,1);
      size += sizeof(*(da->rs))*((n+rrr-1)/rrr);
    r = 0;
    for (i=0; i<n; i+=RR) {
      putuint(da->rl,i/RR,r,k);
      m = 0;
#if 0
      for (j=0; j<RR; j++) {
        if (j % rrr == 0 && i+j < n) {
          da->rs[(i+j)/rrr] = m;
        }
        if (i+j < n && getbit(buf,i+j)==1) m++;
      }
#else
      for (j=0; j<RR; j+=Q) {
        if (j % rrr == 0 && i+j < n) {
          da->rs[(i+j)/rrr] = m;
        }
        if (i+j < n) m += POPCOUNT(buf[(i+j)/Q]);
      }
#endif
      r += m;
    }
  } else {
    da->rl = NULL;  da->rs = NULL;
  }
//
  da->opt = opt;
  da->size = size;


  return;
}

#if 0
static int encode_enum(bitvec_t *out, i64 i, int n, bitvec_t *in, i64 j)
{
  int k,d,w,m;
  bitvec_t x;
  int len;
  
  len = 0;

  m = 0;
  for (k=0; k<n; k++) m += getbit(in,j+k);

  if (m == Q || m == 0) {
    setbits(out, i, blog(Q-1)+1, 0);
    len += blog(Q-1)+1;
    setbits(out, i+len, 1, m/Q);
    len++;
  } else {
    setbits(out, i, blog(Q-1)+1, m);
    len += blog(Q-1)+1;

    d = 0;
    x = 0;
    for (k = 0; k < n; k++) {
      if (getbit(in,j+k)==1) {
        x += nCk[m][n-1-k];
        d++;
        m--;
      }
    }

    w = blog(nCk[d][n]-1)+1;
    setbits(out,i+len,w,x);
    len += w;
  }
  return len;
}
#endif

i64 bitvector_makeindex(bitvector *da, int blocksize, int opt)
{
  i64 i,n,m;
  bitvec_t *buftmp;

//  bitvector_make_selecttbl(); // 最初に1回だけでいいので消す?

  if (da->mtmp == 0) {
    n = da->n;  m = 0;
    buftmp = da->buf;
#if 0
    for (i=0; i<n; i++) {
      if (getbit(buftmp,i)) {
        m++;
      }
    }
#else
    for (i=0; i<n; i+=Q) {
      m += POPCOUNT(getbits(buftmp,i,Q));
    }
#endif
    da->m = m;
  }

  if (opt & SDARRAY_SPARSE) { // sparse array
    i64 mm,r;
    int d;
    int opt2;

    opt2 = SDARRAY_RANK1 | SDARRAY_SPARSE;
    if (opt & SDARRAY_SELECT1) {
      opt2 |= SDARRAY_SELECT1;
    }
    if (opt & SDARRAY_RANK1) {
      opt2 |= SDARRAY_SELECT0;
    }

    if (da->mtmp == 0) {
      mm = m;
      d = 0;
      while (mm < n) {
        mm <<= 1;
        d++;
      }
      da->low_width = d;
      mymalloc(da->buf,(2*m+Q-1)/Q+1,1); // high
      mymalloc(da->low,(d*m+Q-1)/Q+1,1); // low
      for (i=0; i<(2*m+Q-1)/Q+1; i++) da->buf[i] = 0;
      for (i=0; i<(d*m+Q-1)/Q+1; i++) da->low[i] = 0;

      r = 0;
      for (i=0; i<n; i++) {
        if (getbit(buftmp,i)==1) {
          setbit(da->buf,(i>>d)+r,1);
          setbits(da->low,r*d,d,i & ((1<<d)-1));
          r++;
        }
      }
      free(buftmp);
	} else {
#if 0
      da->n++;
#endif
      n = da->n;
      m = da->m;
      change_nm(da, n, m);
      //if (m>0) change_nm(da, n, m); // sada 4/25
      d = da->low_width;
      da->size = 0;
	}

    if (m>0) bitvector_construct(da, 2*m, m, da->buf, blocksize, opt2);
    da->getbit = sparsearray_getbit;
    da->setbit = (void (*)(struct bitvector *da, i64 i, int x))bitvector_notsupported;
    da->rank = sparsearray_rank;
    da->select = sparsearray_select;
    da->succ = bitvector_succ_naive;
    da->pred = bitvector_pred_naive;
    da->getbits = (u64 (*)(struct bitvector *, i64, int))bitvector_notsupported;
    da->setbits = (void (*)(struct bitvector *, i64, int, u64))bitvector_notsupported;
    da->size += sizeof(da->buf[0]) * ((2*m+Q-1)/Q+1);
    da->size += sizeof(da->low[0]) * ((d*m+Q-1)/Q+1);
#if 0
  } else if (opt & SDARRAY_RRR) {
    i64 jl,jm,js;
    i64 p; // pointer to blocks
    i64 p2; // pointer to blocks inside large block
    int pass;
    int w; // #ones in a block

    mymalloc(da->pl,(n+RR-1)/RR*k,1);  
      size += sizeof(*(da->pl))*((n+RR-1)/RR)*k;
    mymalloc(da->ps,(n+blocksize-1)/blocksize,1);
      size += sizeof(*(da->ps))*((n+blocksize-1)/blocksize);
    p = 0;
    for (jl=0; jl<n; jl+=RR) {
      putuint(da->pl,jl/RR,p,k);
      p2 = 0;
      for (jm=0; jm<RR; jm+=blocksize) {
        if (jl+jm < n) {
          da->ps[(jl+jm)/blocksize] = p2;
        }
        for (js=0; js<blocksize; js+=Q) {
          w = 0;
          for (j=0; j<Q; j++) {
            if (jl+jm+js+j < n && getbit(buftmp,jl+jm+js+j)==1) w++;
          }
          p2 += 0000;
        }
      }
      p += p2;
    }

    for (pass = 1; pass <= 2; pass++) {
      if (pass == 1) {
      }
    }
#endif
  } else {
    bitvector_construct(da, n, m, da->buf, blocksize, opt);
    da->size += (da->n+Q-1)/Q * sizeof(da->buf[0]);
    da->rank = bitvector_rank;
    da->select = bitvector_select;
    da->succ = bitvector_succ;
    da->pred = bitvector_pred;
    da->setbits = (void (*)(struct bitvector *, i64, int, u64))bitvector_notsupported;

  }

  return 0;
}

i64 bitvector_write(bitvector *da, FILE *f)
{
  i64 i,k;
  i64 nl;
  i64 rrr;
  i64 n, m;

  i64 size;

  k = da->k;
  rrr = da->rrr;
  size = 0;
  writeuint(1,ID_BITVECTOR,f);
  writeuint(1, da->k, f);
  writeuint(sizeof(da->n), da->n, f);
  writeuint(sizeof(da->m), da->m, f);
  writeuint(4, da->opt, f);
  size += 1 + sizeof(da->n) + sizeof(da->m) + 1;
  if (da->opt & SDARRAY_RANK1) {
    writeuint(sizeof(int), RR, f);
    writeuint(sizeof(int), rrr, f);
    size += 2*sizeof(int);
  }
  if (da->opt & SDARRAY_SELECT1) {
    writeuint(sizeof(da->ml[1]), da->ml[1], f);  size += sizeof(da->ml[1]);
    writeuint(sizeof(da->ms[1]), da->ms[1], f);  size += sizeof(da->ms[1]);
  }
  if (da->opt & SDARRAY_SELECT0) {
    writeuint(sizeof(da->ml[0]), da->ml[0], f);  size += sizeof(da->ml[0]);
    writeuint(sizeof(da->ms[0]), da->ms[0], f);  size += sizeof(da->ms[0]);
  }
  if (da->opt & SDARRAY_SPARSE) {
    writeuint(sizeof(int), da->low_width, f);
  }

  if (da->opt & SDARRAY_SPARSE) {
    n = 2*da->m;  m = da->m;
  } else {
    n = da->n;  m = da->m;
  }

//  writeuint(sizeof(i64), size, f);

  if (da->opt & SDARRAY_RANK1) {
    for (i=0; i<(n+RR-1)/RR; i++) {
      writeuint(k, getuint(da->rl,i,k), f); size += k;
    }
    for (i=0; i<(n+rrr-1)/rrr; i++) {
      writeuint(sizeof(*da->rs), da->rs[i], f);
      size += sizeof(*da->rs);
    }
  }

  if (da->opt & SDARRAY_SELECT1) {
    nl = (m-1) / L + 1;
    for (i=0; i<nl+1; i++) {
//      writeuint(sizeof(*da->lp[1]), da->lp[1][i], f);  size += sizeof(*da->lp[1]);
      writeuint(k, getuint(da->lp[1],i,k), f);  size += k;
    }
    for (i=0; i<nl+1; i++) {
      writeuint(sizeof(*da->p[1]), da->p[1][i], f);  size += sizeof(*da->p[1]);
    }
    for (i=0; i<da->ml[1]*L+1; i++) {
//      writeuint(sizeof(*da->sl[1]), da->sl[1][i], f);  size += sizeof(*da->sl[1]);
      writeuint(k, getuint(da->sl[1],i,k), f);  size += k;
    }
    for (i=0; i<da->ms[1]*(L/LLL)+1; i++) {
      writeuint(sizeof(*da->ss[1]), da->ss[1][i], f);  size += sizeof(*da->ss[1]);
    }
  }
  if (da->opt & SDARRAY_SELECT0) {
    nl = (n-m-1) / L + 1;
    for (i=0; i<nl+1; i++) {
//      writeuint(sizeof(*da->lp[0]), da->lp[0][i], f);  size += sizeof(*da->lp[0]);
      writeuint(k, getuint(da->lp[0],i,k), f);  size += k;
    }
    for (i=0; i<nl+1; i++) {
      writeuint(sizeof(*da->p[0]), da->p[0][i], f);  size += sizeof(*da->p[0]);
    }
    for (i=0; i<da->ml[0]*L+1; i++) {
//      writeuint(sizeof(*da->sl[0]), da->sl[0][i], f);  size += sizeof(*da->sl[0]);
      writeuint(k, getuint(da->sl[0],i,k), f);  size += k;
    }
    for (i=0; i<da->ms[0]*(L/LLL)+1; i++) {
      writeuint(sizeof(*da->ss[0]), da->ss[0][i], f);  size += sizeof(*da->ss[0]);
    }
  }

  if (da->opt & SDARRAY_SPARSE) {
    nl = (2*da->m+Q-1)/Q+1;
  } else {
    nl = (da->n+Q-1)/Q;
  }
  for (i=0; i<nl; i++) {
    writeuint(sizeof(*da->buf), da->buf[i], f);
    size += sizeof(*da->buf);
  }

  if (da->opt & SDARRAY_SPARSE) {
    nl = (da->low_width*da->m+Q-1)/Q+1;
    for (i=0; i<nl; i++) {
      writeuint(sizeof(*da->low), da->low[i], f);
      size += sizeof(*da->low);
    }
  }

  return size;
}

bitvector *bitvector_read(uchar **map)
{
  bitvector *da;
  i64 nl;
  uchar *p;
  i64 n,m;
  i64 rr, rrr;

  mymalloc(da, 1, 0);

  p = *map;

  rr = getuint(p,0,1);  p += 1;
  if (rr != ID_BITVECTOR) {
    printf("bitvector_read: id = %ld is not supported.\n",rr);
    exit(1);
  }

  da->k = getuint(p,0,1);  p += 1;
  da->n = getuint(p,0,sizeof(da->n));  p += sizeof(da->n);
  da->m = getuint(p,0,sizeof(da->m));  p += sizeof(da->m);
  da->opt = getuint(p,0,4);  p += 4;
//  printf("densearray_read: n=%ld m=%ld opt=%d\n",da->n,da->m,da->opt);
  if (da->opt & SDARRAY_RANK1) {
    rr = getuint(p,0,sizeof(int));  p += sizeof(int);
    if (rr != RR) {
      printf("error2 RR=%ld must be %d\n",rr,RR);
    }
    rrr = getuint(p,0,sizeof(int));  p += sizeof(int);
    da->rrr = rrr;
//    printf("RRR = %ld\n",rrr);
#if 0
    if (rrr != RRR) {
      printf("error RRR=%ld must be %d\n",rrr,RRR);
    }
#endif
  }
  if (da->opt & SDARRAY_SELECT1) {
    da->ml[1] = getuint(p,0,sizeof(da->ml[1]));  p += sizeof(da->ml[1]);
    da->ms[1] = getuint(p,0,sizeof(da->ms[1]));  p += sizeof(da->ms[1]);
  }
  if (da->opt & SDARRAY_SELECT0) {
    da->ml[0] = getuint(p,0,sizeof(da->ml[0]));  p += sizeof(da->ml[0]);
    da->ms[0] = getuint(p,0,sizeof(da->ms[0]));  p += sizeof(da->ms[0]);
  }
  if (da->opt & SDARRAY_SPARSE) {
    da->low_width = getuint(p,0,sizeof(int));  p += sizeof(int);
  }
  if (da->opt & SDARRAY_SPARSE) {
    n = 2*da->m;  m = da->m;
  } else {
    n = da->n;  m = da->m;
  }
//  size = getuint(p,0,sizeof(size));  p += sizeof(size);
//  printf("size %ld\n",size);

  if (da->opt & SDARRAY_RANK1) {
    da->rl = p;
    p += sizeof(*da->rl) * ((n+rr-1)/rr) * da->k;
    da->rs = (word *)p;
    p += sizeof(*da->rs) * ((n+rrr-1)/rrr);
  }

  if (da->opt & SDARRAY_SELECT1) {
    nl = (m-1) / L + 1;
//    da->lp[1] = (dword *)p;
//    p += sizeof(*da->lp[1]) * (nl+1);
    da->lp[1] = p;
    p += sizeof(*da->lp[1]) * (nl+1) * da->k;
    da->p[1] = (i64 *)p;
    p += sizeof(*da->p[1]) * (nl+1);
//    da->sl[1] = (dword *)p;
//    p += sizeof(*da->sl[1]) * (da->ml[1]*L+1);
    da->sl[1] = p;
    p += sizeof(*da->sl[1]) * (da->ml[1]*L+1) * da->k;
    da->ss[1] = (word *)p;
    p += sizeof(*da->ss[1]) * (da->ms[1]*(L/LLL)+1);
  }
  if (da->opt & SDARRAY_SELECT0) {
    nl = (n - m-1) / L + 1;
//    da->lp[0] = (dword *)p;
//    p += sizeof(*da->lp[0]) * (nl+1);
    da->lp[0] = p;
    p += sizeof(*da->lp[0]) * (nl+1) * da->k;
    da->p[0] = (i64 *)p;
    p += sizeof(*da->p[0]) * (nl+1);
//    da->sl[0] = (dword *)p;
//    p += sizeof(*da->sl[0]) * (da->ml[0]*L+1);
    da->sl[0] = p;
    p += sizeof(*da->sl[0]) * (da->ml[0]*L+1) * da->k;
    da->ss[0] = (word *)p;
    p += sizeof(*da->ss[0]) * (da->ms[0]*(L/LLL)+1);
  }
  
  if (da->opt & SDARRAY_SPARSE) {
    nl = (2*da->m+Q-1)/Q+1;
  } else {
    nl = (da->n+Q-1)/Q;
  }
  da->buf = (bitvec_t *)p;
  p += sizeof(*da->buf) * nl;

  if (da->opt & SDARRAY_SPARSE) {
    nl = (da->low_width*da->m+Q-1)/Q+1;
    da->low = (bitvec_t *)p;
    p += sizeof(*da->buf) * nl;
  } else {
    da->low = NULL;
  }

  *map = p;


  if (da->opt & SDARRAY_SPARSE) {
    da->getbit = sparsearray_getbit;
    da->rank = sparsearray_rank;
    da->select = sparsearray_select;
    da->succ = bitvector_succ_naive;
    da->pred = bitvector_pred_naive;
  } else {
    da->getbit = bitvector_getbit;
    da->setbit = bitvector_setbit;
    da->getbits = bitvector_getbits;
    da->setbits = bitvector_setbits;

    da->rank = bitvector_rank;
    da->select = bitvector_select;
    da->succ = bitvector_succ;
    da->pred = bitvector_pred;
  }

  bitvector_make_selecttbl();

  return da;
}

bitvector *bitvector_read_from_file(uchar *filename)
{
  MMAP *map;
  uchar *mapp;

  map = mymmap(filename);
  if (map->addr==NULL) {
    perror("mmap2\n");
    exit(1);
  }
  mapp = (uchar *)map->addr;

  return bitvector_read(&mapp);
}


#ifdef MAIN
typedef struct timeb mytimestruct;

void mygettime(mytimestruct *t)
{
  ftime(t);
}

double mylaptime(mytimestruct *before,mytimestruct *after)
{
  double t;
  t = after->time - before->time;
  t += (double)(after->millitm - before->millitm)/1000;
  return t;
}


#define N 10485760
#define CHECK
//#define RANDOM


int main(int argc, char *argv[])
{
  bitvector *s1, *s2;

  i64 i,r,n,m,rr;
  u64 hoge,sum;
  FILE *infp = NULL;
  FILE *out;
  int b;
  i64 m2;

//  bitvec_t *B;
//  byte *B2;
  dword *S,*R;
//  MMAP *map;
//  uchar *mapp;

  double t;
  mytimestruct before,after;

#ifdef __SSE4_2__
  printf("SSE4.2 is available.\n");
#else
  printf("no SSE4.2\n");
#endif

  srand(2);

  n = N; // length of bit vector
  r = 2; // ratio of ones
  m = 0; // number of ones
  if (argc >= 2){
    r = atoi(argv[1]);
    if (r == 0) {
      infp = fopen(argv[2],"rb");
      if (infp == NULL){
        printf("cannot open %s\n",argv[2]);
        return -1;
      }
      fseek(infp,0,SEEK_END);
      n = ftell(infp);
      rewind(infp);
      //printf("n: %d\n",n);
    } else if (argc >= 3) {
      n = atoi(argv[2]);
    }
  }
  rr = r;

  s1 = bitvector_new(n);
  s2 = bitvector_new(n);

  if (!infp) {
    m = 0;
    for (i = 0; i < n; i++) {
      if (rand() % 100 < r) {
//        setbit(B,i,1);
        s1->setbit(s1,i,1);
        s2->setbit(s2,i,1);
        m++;
      } else {
//        setbit(B,i,0);
        s1->setbit(s1,i,0);
        s2->setbit(s2,i,0);
      }
    }
  } else {
    m = 0;
    for (i = 0; i < n; i++){
      int c = fgetc(infp);
      if (c == EOF){
        printf("unexpected error at reading from %s\n",argv[2]);
        return -1;
      }
      if (c == '1' || c == '(') {
//        setbit(B,i,1);
        s1->setbit(s1,i,1);
        s2->setbit(s2,i,1);
        m++;
      } else if (c == '0' || c == ')') {
//        setbit(B,i,0);
        s1->setbit(s1,i,0);
        s2->setbit(s2,i,0);
      } else {
        printf("unexpected error (2) at reading from %s\n",argv[2]);
        return -1;
      }
    }
  }


//  bitvector_makeindex(s1, 512, SDARRAY_RANK1);
  bitvector_makeindex(s1, 512, SDARRAY_RANK1 | SDARRAY_SELECT1);
//  bitvector_makeindex(s1, 512, SDARRAY_RANK1 | SDARRAY_SELECT1 | SDARRAY_SELECT0);
//  bitvector_makeindex(s2, 512, SDARRAY_RANK1 | SDARRAY_SELECT1 | SDARRAY_SELECT0 | SDARRAY_SPARSE);
  bitvector_makeindex(s2, 512, SDARRAY_SELECT1 | SDARRAY_SPARSE);
//  bitvector_makeindex(s1, 512, 0);
//  bitvector_makeindex(s2, 512, SDARRAY_SPARSE);

  printf("da: used memory: %d bytes (%lf bpc)\n",s1->size,(double)s1->size*8/n);
  printf("sa: used memory: %d bytes (%lf bpc)\n",s2->size,(double)s2->size*8/n);

#if 1
  out = fopen("bitvectortmp.dat","w");
  bitvector_write(s1, out);
  fclose(out);
  s1 = bitvector_read_from_file("bitvectortmp.dat");

  out = fopen("bitvectortmp2.dat","w");
  bitvector_write(s2, out);
  fclose(out);
  s2 = bitvector_read_from_file("bitvectortmp2.dat");
#endif

  s2->select(s2, 8010, 1);
  s1->select(s1, 8010, 1);



#ifdef CHECK
  mymalloc(S,n+1,0);
  mymalloc(R,n+1,0);
#endif

  mygettime(&before);
  sum = 0;
  for (i = 0; i < n; i++) {
    if (s1->getbit(s1, i) != s2->getbit(s2, i)) {
      printf("i = %ld getbit(dense) %d getbit(sparse) %d\n", s1->getbit(s1, i), s2->getbit(s2, i));
    }
  }
  mygettime(&after);
  t = mylaptime(&before,&after);
  printf("%ld %f bitvector getbit time\n",rr,t);


for (b = 0; b <= 1; b++) {

#if 0
  mygettime(&before);
  sum = 0;
  for (i = 0; i < n; i++) {
    sum += s1->pred(s1, i, b);
  }
  mygettime(&after);
  t = mylaptime(&before,&after);
  printf("%ld %f bitvector(dense) pred(%d) time sum=%ld\n",rr,t,b,sum);
#endif

#ifdef CHECK
  r = 0;
  S[r] = -1;
  for (i=0; i<n; i++) {
    if (s1->getbit(s1,i) == b) {
      r++;
      S[r] = i;
    }
    R[i] = r;
  }
  for (i = r+1; i <= n; i++) S[i] = n;
#endif

#if 1
  srand(4);

  mygettime(&before);
  hoge = rand();
  sum = 0;
  for (i = 0; i < 100000000; i++) {
    int j;
    //j = (rand() % n);
#ifdef RANDOM
    j = hoge % n;
#else
    j = i % n;
#endif
#ifdef CHECK
    if (s1->rank(s1,j,b) != R[j]) {
      printf("ERROR: (%d) R[%d] = %d, r = %d\n", s1->getbit(s1,i),j, R[j],
                                                 s1->rank(s1,j,b));
    }
    sum += s1->rank(s1,j,b);
#else
    sum += s1->rank(s1,j,b);
#endif
    hoge = hoge * 1566083941U + 1;
  }
  mygettime(&after);
  t = mylaptime(&before,&after);
  printf("%ld %f bitvector(dense) rank(%d) time sum=%ld\n",rr,t,b,sum);
#endif

#if 1
  srand(4);

  mygettime(&before);
  hoge = rand();
  sum = 0;
  for (i = 0; i < 100000000; i++) {
    int j;
    //j = (rand() % n);
#ifdef RANDOM
    j = hoge % n;
#else
    j = i % n;
#endif
#ifdef CHECK
    if (s2->rank(s2,j,b) != R[j]) {
      printf("ERROR: (%d) R[%d] = %d, r = %d\n", s1->getbit(s1,i),j, R[j],
                                                 s2->rank(s2,j,b));
    }
    sum += s2->rank(s2,j,b);
#else
    sum += s2->rank(s2,j,b);
#endif
    hoge = hoge * 1566083941U + 1;
  }
  mygettime(&after);
  t = mylaptime(&before,&after);
  printf("%ld %f bitvector(sparse) rank(%d) time sum=%ld\n",rr,t,b,sum);
#endif

  if (b == 1) m2 = m; else m2 = n-m;

#if 0
  srand(3);

  mygettime(&before);
  hoge = rand();
  sum = 0;
  for (i = 0; i < 100000000; i++) {
    i64 j;
    //j = (rand() % r)+1;
    j = i % m2 + 1;
    if (s1->select(s1,j,b) != s2->select(s2,j,b)) {
      printf("ERROR: s1 = %d, s2 = %d\n",j,s1->select(s1,j,b),s2->select(s2,j,b));
    }
  }
  mygettime(&after);
  t = mylaptime(&before,&after);
  printf("select(%d) time %f\n",b,t);
#endif


#if 1
  srand(3);

  mygettime(&before);
  hoge = rand();
  sum = 0;
  for (i = 0; i < 100000000; i++) {
    i64 j;
    //j = (rand() % r)+1;
#ifdef RANDOM
    j = hoge % m + 1;
#else
    j = i % m2 + 1;
#endif
#ifdef CHECK
    if (s1->select(s1,j,b) != S[j]) {
      printf("ERROR: S[%d] = %d, s = %d\n",j,S[j],s1->select(s1,j,b));
    }
    sum += s1->select(s1,j,b);
#else
    sum += s1->select(s1,j,b);
#endif
    hoge = hoge * 1566083941U + 1;
  }
  mygettime(&after);
  t = mylaptime(&before,&after);
  printf("%ld %f bitvector(dense) select(%d) time sum=%ld\n",rr,t,b,sum);
#endif

#if 1
  srand(3);

  mygettime(&before);
  hoge = rand();
  sum = 0;
  for (i = 0; i < 100000000; i++) {
    i64 j;
    //j = (rand() % r)+1;
#ifdef RANDOM
    j = hoge % m + 1;
#else
    j = i % m2 + 1;
#endif
#ifdef CHECK
    if (s2->select(s2,j,b) != S[j]) {
      printf("ERROR: S[%d] = %d, s = %d\n",j,S[j],s2->select(s2,j,b));
    }
    sum += s2->select(s2,j,b);
#else
    sum += s2->select(s2,j,b);
#endif
    hoge = hoge * 1566083941U + 1;
  }
  mygettime(&after);
  t = mylaptime(&before,&after);
  printf("%ld %f bitvector(sparse) select(%d) time sum=%ld\n",rr,t,b,sum);
#endif

}

  return 0;
}

#endif //  MAIN
