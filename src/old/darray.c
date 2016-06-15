#include <stdio.h>
#include <stdlib.h>
#include "typedefbp.h"
#include "darray.h"

#ifdef __SSE4_2__
#include <smmintrin.h>
#define POPCOUNT(x) _mm_popcnt_u64(x)
#else
#define POPCOUNT(x) _mm_popcnt_u64(x)
//#define POPCOUNT(x) popcount(x)
#endif

#define GAMMATBL 10

#define PBS (sizeof(pb)*8)
#define logM 5
#define M (1<<logM)
#define logP 8
#define P (1<<logP)
#define logLL 16    // size of word
#define LL (1<<logLL)
#define logLLL 7
//#define logLLL 5
#define LLL (1<<logLLL)
//#define logL 10
//#define logL (logLL-3)
#define logL (logLL-1-5)
#define L (1<<logL)

#ifndef min
 #define min(x,y) ((x)<(y)?(x):(y))
#endif

#define mymalloc(p,n,f) {p = malloc((n)*sizeof(*p)); if ((p)==NULL) {printf("not enough memory\n"); exit(1);}}

int setbit(pb *B, i64 i,i64 x)
{
  i64 j,l;

  j = i / D; //何番目の要素か
  l = i % D; //何ビット目か
  if (x==0) B[j] &= (~(1<<(D-1-l)));
  else if (x==1) B[j] |= (1<<(D-1-l));
  else {
    printf("error setbit x=%d\n",x);
    exit(1);
  }
  return x;
}

static int setbits_fast(pb *B, i64 i, int d, u64 x)
{
  u64 y,m;
  int d2;

  B += (i>>logD);
  i &= (D-1);

  while (i+d > D) {
    d2 = D-i; // x の上位 d2 ビットを格納
    y = x >> (d-d2);
    m = (1<<d2)-1;
    *B = (*B & (~m)) | y;
    B++;  i=0;
    d -= d2;
    x &= (1<<d)-1; // x の上位ビットを消去
  }
  m = (1<<d)-1;
  y = x << (D-i-d);
  m <<= (D-i-d);
  *B = (*B & (~m)) | y;
    
    return 0; //edited by Lee
}

int setbits(pb *B, i64 i, i64 d, i64 x)
{
  i64 j;

  for (j=0; j<d; j++) {
    setbit(B,i+j,(x>>(d-j-1))&1);
  }
  return x;
}


int getbit(pb *B, i64 i)
{
  i64 j,l;

  //j = i / D;
  //l = i % D;
  j = i >> logD;
  l = i & (D-1);
  return (B[j] >> (D-1-l)) & 1;
}

dword getbits0(pb *B, i64 i, i64 d)
{
  dword x;
  i64 j;

  x = 0;
  for (j=0; j<d; j++) {
    x <<= 1;
    x += getbit(B,i+j);
  }
  return x;
}

dword getbits(pb *B, i64 i, i64 d)
{
  dword x,z;

  if (d == 0) return 0;

  B += (i >>logD);
  i &= (D-1);
  if (i+d <= D) {
    x = B[0];
    x <<= i;
    x >>= (D-d);  // D==32, d==0 だと動かない
  } else {
    x = B[0] << i;
    x >>= D-d;
    z = B[1] >> (D-(i+d-D));
    x += z;
  }
  return x;
}


i64 getpattern(pb *B, i64 i, i64 k, pb pat)
{
  i64 j;
  i64 x;
  x = 1;
  for (j=0; j<k; j++) {
    x &= getbit(B,i+j) ^ (~(pat>>(k-1-j)));
  }
  //printf("getpattern(%d) = %d\n",i,x);
  return x;
}


static const unsigned int popCount[] = {
0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8
};

static int selecttbl[8*256];

unsigned int popcount(pb x)
{
  pb r;
#if 0
  r = x;
  r = r - ((r>>1) & 0x77777777) - ((r>>2) & 0x33333333) - ((r>>3) & 0x11111111);
  r = ((r + (r>>4)) & 0x0f0f0f0f) % 0xff;
#elif 1
  r = x;
  r = ((r & 0xaaaaaaaa)>>1) + (r & 0x55555555);
  r = ((r & 0xcccccccc)>>2) + (r & 0x33333333);
  //r = ((r & 0xf0f0f0f0)>>4) + (r & 0x0f0f0f0f);
  r = ((r>>4) + r) & 0x0f0f0f0f;
  //r = ((r & 0xff00ff00)>>8) + (r & 0x00ff00ff);
  r = (r>>8) + r;
  //r = ((r & 0xffff0000)>>16) + (r & 0x0000ffff);
  r = ((r>>16) + r) & 63;
#else
  r = popCount[x & 0xff];
  x >>= 8;
  r += popCount[x & 0xff];
  x >>= 8;
  r += popCount[x & 0xff];
  x >>= 8;
  r += popCount[x & 0xff];
#endif
  return r;
}

unsigned int popcount8(pb x)
{
  dword r;
#if 1
  r = x;
  r = ((r & 0xaa)>>1) + (r & 0x55);
  r = ((r & 0xcc)>>2) + (r & 0x33);
  r = ((r>>4) + r) & 0x0f;
#else
  r = popCount[x & 0xff];
#endif
  return r;
}

void darray_make_selecttbl(void)
{
  i64 i,x,r;
  pb buf[1];

  buf[0] = 0;
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


int darray_construct(darray *da, i64 n, pb *buf, int opt)
{
  i64 i,j,k,m;
  i64 nl;
  i64 p,pp;
  i64 il,is,ml,ms;
  i64 r,m2;
  i64 p1,p2,p3,p4,s1,s2,s3,s4;

  da->idx_size = 0;

//  darray_make_selecttbl();

  da->buf = da->bufc = NULL;


  if (L/LLL == 0) {
    printf("ERROR: L=%d LLL=%d\n",L,LLL);
    exit(1);
  }

  m = 0;
#if 0
  for (i=0; i<n; i++) m += getbit(buf,i);
#else
  for (i=0; i<n-D; i+=D) m += POPCOUNT(getbits(buf,i, D));
  for (; i<n; i++) m += getbit(buf,i);
#endif
  da->n = n;
  da->m = m;
  //printf("n=%d m=%d\n",n,m);

  da->buf = buf;
  da->opt = opt;
  da->pcm = NULL;
  da->pcl = NULL;

  if (opt & (~OPT_NO_RANK)) {  // construct select table
#if 0
    mymalloc(s,m,0);
    m = 0;
    for (i=0; i<n; i++) {
      if (getbit(buf,i)) {
        m++;
        s[m-1] = i;
      }
    }
#endif    
    nl = (m-1) / L + 1;
    mymalloc(da->lp,nl+1,1);  da->idx_size += (nl+1)*sizeof(*da->lp);
    mymalloc(da->p,nl+1,1);  da->idx_size += (nl+1)*sizeof(*da->p);
#if 0
    printf("lp table: %d bytes (%1.2f bpc)\n",(nl+1)*sizeof(*da->lp), (double)(nl+1)*sizeof(*da->lp) * 8/n);
    printf("p table: %d bytes (%1.2f bpc)\n",(nl+1)*sizeof(*da->p), (double)(nl+1)*sizeof(*da->p) * 8/n);
#endif    

    for (r = 0; r < 2; r++) {
      s1 = s2 = s3 = s4 = 0;
      p1 = p2 = p3 = p4 = -1;
    
      ml = ms = 0;
      for (il = 0; il < nl; il++) {
        //pp = s[il*L];
#if 0
        while (s1 <= il*L) {
          if (getbit(buf,p1+1)) s1++;
          p1++;
        }
#else
        while (s1 <= il*L - D) {
          s1 += POPCOUNT(getbits(buf,p1+1, D));
          p1 += D;
        }
        while (s1 <= il*L) {
          if (getbit(buf,p1+1)) s1++;
          p1++;
        }
#endif
        pp = p1;
        da->lp[il] = pp;
        i = min((il+1)*L-1,m-1);
        //p = s[i];
#if 0
        while (s2 <= i) {
          if (getbit(buf,p2+1)) s2++;
          p2++;
        }
#else
        while (s2 <= i - D) {
          s2 += POPCOUNT(getbits(buf,p2+1, D));
          p2 += D;
        }
        while (s2 <= i) {
          if (getbit(buf,p2+1)) s2++;
          p2++;
        }
#endif
        p = p2;
        //printf("%d ",p-pp);
        if (p - pp >= LL) {
          if (r == 1) {
            for (is = 0; is < L; is++) {
              if (il*L+is >= m) break;
              //da->sl[ml*L+is] = s[il*L+is];
#if 0
              while (s3 <= il*L+is) {
                if (getbit(buf,p3+1)) s3++;
                p3++;
              }
#else
              while (s3 <= il*L+is - D) {
                s3 += POPCOUNT(getbits(buf,p3+1, D));
                p3 += D;
              }
              while (s3 <= il*L+is) {
                if (getbit(buf,p3+1)) s3++;
                p3++;
              }
#endif
              da->sl[ml*L+is] = p3;
            }
          }
          da->p[il] = -(ml+1);
          ml++;
        } else {
          if (r == 1) {
            for (is = 0; is < L/LLL; is++) {
              if (il*L+is*LLL >= m) break;
#if 0
              while (s4 <= il*L+is*LLL) {
                if (getbit(buf,p4+1)) s4++;
                p4++;
              }
#else
              while (s4 <= il*L+is*LLL - D) {
                s4 += POPCOUNT(getbits(buf,p4+1, D));
                p4 += D;
              }
              while (s4 <= il*L+is*LLL) {
                if (getbit(buf,p4+1)) s4++;
                p4++;
              }
#endif
              //da->ss[ms*(L/LLL)+is] = s[il*L+is*LLL] - pp;
              da->ss[ms*(L/LLL)+is] = (word)(p4 - pp);
            }
          }
          da->p[il] = ms;
          ms++;
        }
      }
      if (r == 0) {
        mymalloc(da->sl,ml*L+1,1);  da->idx_size += (ml*L+1)*sizeof(*da->sl);
        mymalloc(da->ss,ms*(L/LLL)+1,1);  da->idx_size += (ms*(L/LLL)+1)*sizeof(*da->ss);
#if 0
        printf("sl table: %d bytes (%1.2f bpc)\n",(ml*L+1)*sizeof(*da->sl), (double)(ml*L+1)*sizeof(*da->sl) * 8/n);
        printf("ss table: %d bytes (%1.2f bpc)\n",(ms*(L/LLL)+1)*sizeof(*da->ss), (double)(ms*(L/LLL)+1)*sizeof(*da->ss) * 8/n);
#endif
      }
    }
    //free(s);
  } else { // no select table
    da->lp = NULL;
    da->p = NULL;
    da->sl = NULL;
    da->ss = NULL;
  }

  // construct rank table

  if ((opt & OPT_NO_RANK) == 0) {
    mymalloc(da->rl,n/R1+2,1);  da->idx_size += (n/R1+2)*sizeof(*da->rl);
    mymalloc(da->rm,n/RR+2,1);  da->idx_size += (n/RR+2)*sizeof(*da->rm);
    mymalloc(da->rs,n/RRR+2,1);  da->idx_size += (n/RRR+2)*sizeof(*da->rs);
#if 0
    printf("rl table: %d bytes (%1.2f bpc)\n",(n/R1+2)*sizeof(*da->rl), (double)(n/R1+2)*sizeof(*da->rl) * 8/n);
    printf("rm table: %d bytes (%1.2f bpc)\n",(n/RR+2)*sizeof(*da->rm), (double)(n/RR+2)*sizeof(*da->rm) * 8/n);
    printf("rs table: %d bytes (%1.2f bpc)\n",(n/RRR+2)*sizeof(*da->rs), (double)(n/RRR+2)*sizeof(*da->rs) * 8/n);
#endif
    r = 0;
    for (k=0; k<=n+R1; k+=R1) {
      da->rl[k/R1] = r;
      m2 = 0;
      for (i=0; i<R1; i+=RR) {
        if (k+i <= n) da->rm[(k+i)/RR] = (word)m2;
        m = 0;
#if 1
        for (j=0; j<RR; j++) {
          if (k+i+j < n && getbit(buf,k+i+j)==1) m++;
          if (j % RRR == RRR-1) {
            if (k+i+j <= n) da->rs[(k+i+j)/RRR] = (byte)m;
            m2 += m;
            m = 0;
          }
#else
        for (j=0; j<RR; j+=D) {
          if (k+i+j < n) m += POPCOUNT(buf[(k+i+j)/D]);
          if (j % RRR == RRR-1) {
            if (k+i+j <= n) da->rs[(k+i+j)/RRR] = (byte)m;
            m2 += m;
            m = 0;
          }
#endif
        }
        if (m != 0) {
          printf("???\n");
        }
        //m2 += m;
      }
      r += m2;
    }
  }

  return 0;
}

static int blog(i64 x) // [0,n] の数を格納するには blog(n)+1 ビット必要
{
int l;
  l = -1;
  while (x>0) {
    x>>=1;
    l++;
  }
  return l;
}

static int decodegamma_naive(pb *B,i64 p,i64 *ans)
{
i64 w,w2;
i64 x,y;
int i;

//  w = getzerorun(B,p);
  w = 0;
  while (getbit(B,p+w)==0) w++;

  x = 1;
  x <<= w;
  x += getbits(B, p+w+1, w);

  *ans = x;
  return 2*w+1;
}


//static int *decodegamma_tbl_w;
//static int *decodegamma_tbl_x;
static int decodegamma_tbl_w[1<<GAMMATBL];
static int decodegamma_tbl_x[1<<GAMMATBL];

static void decodegamma_maketable(int d)
{
  i64 w, ans;
  i64 i,j;
  pb x[4];
  int d2;
  w = 1 << d; // テーブルのサイズ
//  mymalloc(decodegamma_tbl_w, w, 0);
//  mymalloc(decodegamma_tbl_x, w, 0);
  for (i = 0; i < w; i++) {
    x[0] = i << (D-d);
    d2 = decodegamma_naive(x, 0, &ans);
    if (d2 <= d) {
      decodegamma_tbl_w[i] = d2;
      decodegamma_tbl_x[i] = ans;
    } else {
      decodegamma_tbl_w[i] = -1;
    }
  }
}

static int encodegamma(pb *B,i64 p,i64 x) /* x >= 1 */
{
i64 j,w;
  if (x==0) {
    fprintf(stderr,"encodegamma %ld\n",x);  exit(1);
  }
  w = blog(x)+1;
  for (j=0;j<w-1;j++) setbit(B,p+j,0);
  for (j=w-1;j>=0;j--) setbit(B,p+(w-1)+(w-1)-j,(x >> j)&1);
  return 2*w-1;
}

static int getzerorun(pb *B,i64 p)
{
i64 w,w2;
#if 0
  w = 0;
  while (getbit(B,p+w)==0) w++;
#endif

#if 1
  w = 0;
  w2 = getbits(B,p,D);
  while ((w2 & 0x80000000) == 0) {
    w2 <<= 1;
    w++;
  }
#endif

#if 0
  w = 0;
  w2 = getbits(B,p,D);
  if ((w2 & 0xffff0000) == 0) {
    w += 16;
    w2 <<= 16;
  }
  if ((w2 & 0xff000000) == 0) {
    w += 8;
    w2 <<= 8;
  }
  if ((w2 & 0xf0000000) == 0) {
    w += 4;
    w2 <<= 4;
  }
  if ((w2 & 0xc0000000) == 0) {
    w += 2;
    w2 <<= 2;
  }
  if ((w2 & 0x80000000) == 0) {
    w += 1;
//    w2 <<= 1;
  }
#endif
  return w;
}

static int decodegamma(pb *B,i64 p,i64 *ans)
{
i64 w;
u64 w2;
i64 x,y;
int i;

  w2 = getbits(B, p, D);
  i = decodegamma_tbl_w[w2 >> (D-GAMMATBL)];
  if (i >= 0) {
    *ans = decodegamma_tbl_x[w2 >> (D-GAMMATBL)];
    return i;
  }

//  w = getzerorun(B,p);
  w = 0;
  while ((w2 & 0x80000000) == 0) {
    w2 <<= 1;
    w++;
  }


  x = 1;
#if 0
  y = 0;
  for (i=0;i<w;i++) {
    y <<= 1;
    y += getbit(B,p+w+1+i);
  }
  x <<= w;
  x += y;
#else
  x <<= w;
  x += getbits(B, p+w+1, w);
#endif

  *ans = x;
  return 2*w+1;
}

int darray_comp_rle(darray *da)
{
i64 n;
i64 i1, i2, i3;
int loop;
pb *B, *compB;
i64 p1, p2;
int r;
int c, c2;
pb tmpbuf[4];

  decodegamma_maketable(GAMMATBL);

  n = da->n;
  B = da->buf;
  // ポインタの配列の確保
  mymalloc(da->pcl, (n+RR-1)/RR, 0);
  da->idx_size += (n+RR-1)/RR * sizeof(*da->pcl);
  mymalloc(da->pcm, (n+RRR-1)/RRR, 0);
  da->idx_size += (n+RRR-1)/RRR * sizeof(*da->pcm);

  for (loop = 1; loop <= 2; loop++) {
    p1 = 0; // 圧縮されたベクトルへのポインタ
    for (i1 = 0; i1 < n; i1 += RR) {
      da->pcl[i1/RR] = p1;
      p2 = 0; // 圧縮されたベクトルへのポインタ (大ブロック内)
      for (i2 = 0; i2 < RR; i2 += RRR) {
        if (i1 + i2 >= n) break;
        da->pcm[(i1+i2)/RRR] = p2;
        c = getbit(B, i1 + i2); // 中ブロックの最初のビット
        if (loop == 2) {
          setbit(compB, p1 + p2, c); // 最初のビットはそのまま格納
        }
        p2++;
        r = 1; // run length
        for (i3 = 1; i3 < RRR; i3 += 1) {
          if (i1 + i2 + i3 >= n) break;
          c2 = getbit(B, i1 + i2 + i3);
          if (c2 == c) {
            r++;
          } else {
            if (loop == 1) {
              p2 += encodegamma(tmpbuf, 0, r);
            } else {
              p2 += encodegamma(compB, p1 + p2, r);
            }
            c = c2;
            r = 1;
          }
        }
        if (r > 0) { // 必ず成立
          if (loop == 1) {
            p2 += encodegamma(tmpbuf, 0, r);
          } else {
            p2 += encodegamma(compB, p1 + p2, r);
          }
        }
      }
      p1 += p2;
    }
    if (loop == 1) {
      printf("original length = %ld, RLE = %ld\n", n, p1);
      mymalloc(compB, (p1+D-1)/D, 0);
      da->idx_size += (p1+D-1)/D * sizeof(pb);
    }
  }

  free(B); // sada 2015/1/16
  da->bufc = compB;
  da->buf = NULL;
  da->opt |= OPT_DARRAY_COMP_RLE;
  return 0;
}

void darray_decode_block_rle(darray *da, i64 s, i64 t, pb *out)
{
pb *buf;
i64 p, q, r;
i64 i;
int c;
  if (t >= da->n) t = da->n-1;
  p = da->pcl[s >> logRR] + da->pcm[s >> logRRR];
  buf = da->bufc;

  q = s & ~(RRR-1); // decode開始位置

  while (q <= t) {
    if ((q & (RRR-1)) == 0) { // 中ブロック境界
      c = getbit(buf, p++);
    }
    p += decodegamma(buf, p, &r);
    while (r > 0) {
      if (s <= q && q <= t) setbit(out, q-s, c);
      r--;
      q++;
    }
    c = 1-c;
  }
#if 0 // for debug
  if ((t-s) % D != D-1) {
    i64 t2;
    t2 = t + D - ((t-s) % D) - 1;
    for (i = t+1;  i <= t2; i++) {
      setbit(out, i-s, 0);
    }
  }
#endif
}

i64 darray_decode_block_rle2(darray *da, i64 s, i64 t, i64 *out)
{
pb *buf;
i64 p, q, r;
i64 i;
int c;
i64 num;

  if (s % RRR != 0) {
    printf("darray_decode_block_rle2: s = %ld must be a multiple of %d\n", s, RRR);
    exit(1);
  }

  if (t >= da->n) t = da->n-1;
  p = da->pcl[s >> logRR] + da->pcm[s >> logRRR];
  buf = da->bufc;

//  q = s & ~(RRR-1); // decode開始位置
  q = s;

  num = 0;
  while (q <= t) {
    if ((q & (RRR-1)) == 0) { // 中ブロック境界
      c = getbit(buf, p++);
    }
    p += decodegamma(buf, p, &r);
    if (c == 1) {
      out[num] = r;
    } else {
      out[num] = -r;
    }
    q += r;
    num++;
    c = 1-c;
  }
  return num;
}

i64 darray_rank0(darray *da, i64 i)
{
  i64 r,j;
  pb *p;

#if (RRR == D*2)
  r = da->rl[i>>logR] + da->rm[i>>logRR] + da->rs[i>>logRRR];
  p = da->buf + ((i>>logRRR)<<(logRRR-logD));
  j = i & (RRR-1);
  if (j < D) r += POPCOUNT(*p >> (D-1-j));
  else r += POPCOUNT(*p) + POPCOUNT(p[1] >> (D-1-(j-D)));
#else

  j = i & (RRR-1);
  if (j < RRR/2) {
    r = da->rl[i>>logR] + da->rm[i>>logRR] + da->rs[i>>logRRR];
    p = da->buf + ((i>>logRRR)<<(logRRR-logD));
    while (j >= D) {
      r += POPCOUNT(*p++);
      j -= D;
    }
    r += POPCOUNT(*p >> (D-1-j));
  } else {
    j = RRR-1 - (i & (RRR-1));
    i += j+1;
    r = da->rl[i>>logR] + da->rm[i>>logRR] + da->rs[i>>logRRR];
    p = da->buf + ((i>>logRRR)<<(logRRR-logD));
    while (j >= D) {
      r -= POPCOUNT(*--p);
      j -= D;
    }
    if (j > 0) r -= POPCOUNT(*--p << (D-j));
  }

#endif

  return r;
}

i64 darray_rank(darray *da, i64 i)
{
  i64 r,j;
  pb *p;
  pb *p2;
  pb blocktmp[RRR/D];

  r = da->rl[i>>logR] + da->rm[i>>logRR];
  j = (i>>logRRR) & (RR/RRR-1);
  while (j > 0) {
    r += da->rs[((i>>logRR)<<(logRR-logRRR))+j-1];
    j--;
  }

//////////////////////////////// *p は圧縮されていないデータへのポインタ
  if (da->opt & OPT_DARRAY_COMP_RLE) {
    i64 s;
    s = (i>>logRRR)<<logRRR;
    darray_decode_block_rle(da, s, s+RRR-1, blocktmp);
    p = &blocktmp[0];
    p2 = da->buf + ((i>>logRRR)<<(logRRR-logD));
  } else {
    p = da->buf + ((i>>logRRR)<<(logRRR-logD));
  }
  j = i & (RRR-1);
  while (j >= D) {
//////////////////////////////// *p は圧縮されていないデータへのポインタ
//    if (da->opt & OPT_DARRAY_COMP_RLE) {
//      if (*p != *p2) {
//        printf("darray_rank: p %x p2 %x\n", p, p2);
//      }
//    }
    r += POPCOUNT(*p++);
    p2++;
    j -= D;
  }
//////////////////////////////// *p は圧縮されていないデータへのポインタ
//  if (da->opt & OPT_DARRAY_COMP_RLE) {
//    if (*p != *p2) {
//    i64 s;
//    s = (i>>logRRR)<<logRRR;
//      printf("darray_rank: p %x p2 %x\n", p, p2);
//    darray_decode_block_rle(da, s, s+RRR-1, blocktmp);
//    }
//  }
  r += POPCOUNT(*p >> (D-1-j));

  return r;
}

i64 darray_select_bsearch(darray *da, i64 i, pb (*getpat)(pb *))
{
  i64 l,r,m,n;
  i64 t,x,rr;
  pb *p,w;

  // for debug
  //s = darray_select(da,i,1);
  //
  //printf("select(%d)=%d\n",i,s);



  if (i == 0) return -1;

  //printf("i=%d da->m=%d\n",i,da->m);
  if (i > da->m) {
    return -1;
  }

  i--;



  //n = da->m;
  n = da->n;

  t = i;

  l = 0;  r = (n-1)>>logR;
  // find the smallest index x s.t. rl[x] >= t
  while (l < r) {
    m = (l+r)/2;
    //printf("[%d,%d] m=%d rl[m+1]=%d t=%d\n",l,r,m,da->rl[m+1],t);
    if ((i64)da->rl[m+1] > t) { // m+1 is out of range
      r = m;  // new r = m >= l
    } else {
      l = m+1; // new l = m+1 <= r
    }
  }
  x = l;
  t -= da->rl[x];

  x <<= logR;

  l = x >> logRR;  r = (min(x+R1-1,n))>>logRR;
  while (l < r) {
    m = (l+r)/2;
    //printf("[%d,%d] m=%d rm[m+1]=%d t=%d\n",l,r,m,da->rm[m+1],t);
    if (da->rm[m+1] > t) { // m+1 is out of range
      r = m;
    } else {
      l = m+1; // new l = m+1 <= r
    }
  }
  x = l;
  t -= da->rm[x];

  x <<= logRR;

#if 0
  l = x >> logRRR;  r = (min(x+RR-1,n))>>logRRR;
  while (l < r) {
    m = (l+r)/2;
    //printf("m=%d rs[m+1]=%d t=%d\n",m,da->rs[m+1],t);
    if (da->rs[m+1] > t) { // m+1 is out of range
      r = m;
    } else {
      l = m+1; // new l = m+1 <= r
    }
  }
  x = l;
  t -= da->rs[x];
#else
  l = x >> logRRR;
  while (t > da->rs[l]) {
    //printf("l=%d rs[l]=%d t=%d\n",l,da->rs[l],t);
    t -= da->rs[l];
    l++;
  }
  x = l;
#endif

  x <<= logRRR;

//////////////////////////////// *p は圧縮されていないデータへのポインタ
  p = &da->buf[x >> logD];
  while (1) {
    //printf("x=%d m=%d t=%d\n",x,m,t);
//////////////////////////////// *p は圧縮されていないデータへのポインタ
    m = POPCOUNT(getpat(p));
    if (m > t) break;
    t -= m;
    x += D;
    p++;
  }

//////////////////////////////// *p は圧縮されていないデータへのポインタ
  w = getpat(p);
  while (1) {
    rr = popCount[w >> (D-8)];
    if (rr > t) break;
    t -= rr;
    x += 8;
    w <<= 8;
  }
  x += selecttbl[((t-0)<<8)+(w>>(D-8))];

#if 0
  if (x != s) {
    printf("error x=%d s=%d\n",x,s);
  }
#endif
  return x;
}

i64 darray_pat_rank(darray *da, i64 i, pb (*getpat)(pb *))
{
  i64 r,j;
  pb *p;

  r = da->rl[i>>logR] + da->rm[i>>logRR];
  j = (i>>logRRR) & (RR/RRR-1);
  while (j > 0) {
    r += da->rs[((i>>logRR)<<(logRR-logRRR))+j-1];
    j--;
  }

//////////////////////////////// *p は圧縮されていないデータへのポインタ
  p = da->buf + ((i>>logRRR)<<(logRRR-logD));
  j = i & (RRR-1);
  while (j >= D) {
//////////////////////////////// *p は圧縮されていないデータへのポインタ
    r += POPCOUNT(getpat(p));
    p++;
    j -= D;
  }
//////////////////////////////// *p は圧縮されていないデータへのポインタ
  r += POPCOUNT(getpat(p) >> (D-1-j));

  return r;
}


i64 darray_select(darray *da, i64 i,int f)
{
  i64 p,r;
  i64 il;
  i64 rr;
  pb x;
  pb *q;

  if (i == 0) return -1;

  if (i > da->m) {
    return -1;
    //printf("ERROR: m=%d i=%d\n",da->m,i);
    //exit(1);
  }

  i--;

  il = da->p[i>>logL];
  if (il < 0) {
    il = -il-1;
    p = da->sl[(il<<logL)+(i & (L-1))];
  } else {
    p = da->lp[i>>logL];
    p += da->ss[(il<<(logL-logLLL))+(i & (L-1))/LLL];
    r = i - (i & (LLL-1));

//////////////////////////////// *q は圧縮されていないデータへのポインタ
    q = &(da->buf[p>>logD]);

    if (f == 1) {
      rr = p & (D-1);
//////////////////////////////// *q は圧縮されていないデータへのポインタ
      r -= POPCOUNT(*q >> (D-1-rr));
      p = p - rr;
      
      while (1) {
//////////////////////////////// *q は圧縮されていないデータへのポインタ
        rr = POPCOUNT(*q);
        if (r + rr >= i) break;
        r += rr;
        p += D;
        q++;
      }
      
//////////////////////////////// *q は圧縮されていないデータへのポインタ
      x = *q;
      while (1) {
        rr = POPCOUNT(x >> (D-8));
        //rr = popCount[x >> (D-8)];
        //rr = popcount8(x >> (D-8));
        if (r + rr >= i) break;
        r += rr;
        p += 8;
        x <<= 8;
      }
      p += selecttbl[((i-r-1)<<8)+(x>>(D-8))];
    } else {
      rr = p & (D-1);
//////////////////////////////// *q は圧縮されていないデータへのポインタ
      r -= POPCOUNT((~(*q))  >> (D-1-rr));
      p = p - rr;
      
      while (1) {
//////////////////////////////// *q は圧縮されていないデータへのポインタ
        rr = POPCOUNT(~(*q));
        if (r + rr >= i) break;
        r += rr;
        p += D;
        q++;
      }
      
//////////////////////////////// *q は圧縮されていないデータへのポインタ
      x = ~(*q);

      while (1) {
        rr = POPCOUNT(x >> (D-8));
        //rr = popCount[x >> (D-8)];
        //rr = popcount8(x >> (D-8));
        if (r + rr >= i) break;
        r += rr;
        p += 8;
        x <<= 8;
      }
      p += selecttbl[((i-r-1)<<8)+(x>>(D-8))];
    }
  }
  return p;
}

i64 darray_pat_select(darray *da, i64 i, pb (*getpat)(pb *))
{
  i64 p,r;
  i64 il;
  i64 rr;
  pb x;
  pb *q;

  if (i == 0) return -1;

  if (i > da->m) {
    return -1;
    //printf("ERROR: m=%d i=%d\n",da->m,i);
    //exit(1);
  }

  i--;

  il = da->p[i>>logL];
  if (il < 0) {
    il = -il-1;
    p = da->sl[(il<<logL)+(i & (L-1))];
  } else {
    p = da->lp[i>>logL];
    p += da->ss[(il<<(logL-logLLL))+(i & (L-1))/LLL];
    r = i - (i & (LLL-1));

//////////////////////////////// *q は圧縮されていないデータへのポインタ
    q = &(da->buf[p>>logD]);

    rr = p & (D-1);
//////////////////////////////// *q は圧縮されていないデータへのポインタ
    r -= POPCOUNT(getpat(q) >> (D-1-rr));
    p = p - rr;
    
    while (1) {
//////////////////////////////// *q は圧縮されていないデータへのポインタ
      rr = POPCOUNT(getpat(q));
      if (r + rr >= i) break;
      r += rr;
      p += D;
      q++;
    }
    
//////////////////////////////// *q は圧縮されていないデータへのポインタ
    x = getpat(q);
    while (1) {
      rr = POPCOUNT(x >> (D-8));
      //rr = popCount[x >> (D-8)];
      //rr = popcount8(x >> (D-8));
      if (r + rr >= i) break;
      r += rr;
      p += 8;
      x <<= 8;
    }
    p += selecttbl[((i-r-1)<<8)+(x>>(D-8))];
  }
  return p;
}

i64 darray_pat_construct(darray *da, i64 n, pb *buf, i64 k, pb pat, int opt)
{
  i64 i;
  pb *b;
  mymalloc(b,(n+D-1)/D,0);

  for (i=0; i<n-k+1; i++) {
    setbit(b,i,getpattern(buf,i,k,pat));
  }
  for (i=n-k+1; i<n; i++) {
    setbit(b,i,0);
  }

  darray_construct(da,n,b,opt);
  da->buf = buf;

  free(b);
  
  return 0;
}

//added by Lee
void darray_free(darray *da) {
    if (da->buf) free(da->buf);
    if (da->bufc) free(da->bufc);
#ifdef INDEX64
    if (da->lp) free(da->lp);
    if (da->sl) free(da->sl);
#else
    if (da->lp) free(da->lp);
    if (da->sl) free(da->sl);
#endif
    if (da->ss) free(da->ss);
    if (da->p) free(da->p);
    if (da->rl) free(da->rl);
    if (da->rm) free(da->rm);
    if (da->rs) free(da->rs);
#ifdef INDEX64
    if (da->pcl) free(da->pcl);
#else
    if (da->pcl) free(da->pcl);
#endif
    if (da->pcm) free(da->pcm);
}
