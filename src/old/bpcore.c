#include <stdio.h>
#include <stdlib.h>
#include "typedefbp.h"
#include "darray.h"
#include "bp.h"

#define SBid(i) ((i)>>logSB)
#define SBfirst(i) ((i) & (~(SB-1)))
#define SBlast(i) ((i) | (SB-1))

#define MBid(i) ((i)>>logMB)
#define MBfirst(i) ((i) & (~(MB-1)))
#define MBlast(i) ((i) | (MB-1))


#define NOTFOUND -2
#define CONTINUE -3
#define END -4
#define FOUND -5

i64 blog(i64 x)
{
i64 l;
  l = 0;
  while (x>0) {
    x>>=1;
    l++;
  }
  return l;
}

//////////////////////////////// *b ‚Íˆ³k‚³‚ê‚Ä‚¢‚È‚¢ƒf[ƒ^‚Ö‚Ìƒ|ƒCƒ“ƒ^
pb getpat_preorder(pb *b)
{
  printf("error 1\n");  exit(1);
  return *b;
}

//////////////////////////////// *b ‚Íˆ³k‚³‚ê‚Ä‚¢‚È‚¢ƒf[ƒ^‚Ö‚Ìƒ|ƒCƒ“ƒ^
pb getpat_postorder(pb *b)
{
  printf("error 2\n");  exit(1);
  return ~(*b);
}

//////////////////////////////// *b ‚Íˆ³k‚³‚ê‚Ä‚¢‚È‚¢ƒf[ƒ^‚Ö‚Ìƒ|ƒCƒ“ƒ^
pb getpat_leaf(pb *b)
{
  pb w1,w2,w;
  printf("error 3\n");  exit(1);
  w1 = b[0];
  w2 = (w1 << 1) + (b[1] >> (D-1));
  w = w1 & (~w2);
  return w;
}

//////////////////////////////// *b ‚Íˆ³k‚³‚ê‚Ä‚¢‚È‚¢ƒf[ƒ^‚Ö‚Ìƒ|ƒCƒ“ƒ^
pb getpat_inorder(pb *b)
{
  pb w1,w2,w;
  printf("error 4\n");  exit(1);
  w1 = b[0];
  w2 = (w1 << 1) + (b[1] >> (D-1));
  w = (~w1) & w2;
  return w;
}

//////////////////////////////// *b ‚Íˆ³k‚³‚ê‚Ä‚¢‚È‚¢ƒf[ƒ^‚Ö‚Ìƒ|ƒCƒ“ƒ^
pb getpat_dfuds_leaf(pb *b)
{
  pb w1,w2,w;
  printf("error 5\n");  exit(1);
  w1 = b[0];
  w2 = (w1 << 1) + (b[1] >> (D-1));
  w = (~w1) & (~w2);
  return w;
}



///////////////////////////////////////////
//  depth(bp *b, i64 s)
//    returns the depth of s
//  The root node has depth 1
///////////////////////////////////////////
i64 depth(bp *b, i64 s)
{
  i64 d;
  if (s < 0) return 0;
  d = 2 * darray_rank(b->da,s) - (s+1);
#if 0
  if (d != naive_depth(b,s)) {
    d = naive_depth(b,s);
    darray_rank(b->da,s);
  }
  //printf("depth(%d)=%d\n",s,d);
#endif
  return d;
}


////////////////////////////////////////////////
// i64 fast_depth(bp *b, i64 s)
// ƒuƒƒbƒN‹«ŠE‚Ìdepth‚È‚Ì‚ÅŠ‡ŒÊ—ñ‚ÍŽg‚í‚È‚¢
///////////////////////////////////////////////
i64 fast_depth(bp *b, i64 s)
{
  i64 d;
  darray *da;
  i64 r,j;

  s++;
#if 0
  if ((s & (RRR-1)) != 0) {
    printf("fast_depth:warning s=%d\n",s);
    return depth(b,s);
  }
#endif
  da = b->da;
  //d = 2 * (da->rl[s>>logR] + da->rm[s>>logRR] + da->rs[s>>logRRR]) - s;

  r = da->rl[s>>logR] + da->rm[s>>logRR];
  j = (s>>logRRR) & (RR/RRR-1);
  while (j > 0) {
    r += da->rs[((s>>logRR)<<(logRR-logRRR))+j-1];
    j--;
  }
  d = 2 * r - s;

  return d;
}

i64 search_SB_r_rle(bp *b, i64 i, i64 rel) // b[i..] ‚ðŒŸõ (i ‚ðŠÜ‚Þ)
{
  i64 j,r,n,il;
  pb *p,x,w;
  pb *p2;
  i64 blocktmp2[SB];
  i64 num;
  i64 s,t,u,d;
  i64 z;
  i64 ii;

//  if (i == 979456) {
//    printf("hoge\n");
//  }
//  ii = i;
//  z = search_SB_r(b, i, rel);

  n = b->n;
  il = min(SBlast(i),n-1); // i ‚ðŠÜ‚ÞƒuƒƒbƒN“à‚ÌÅŒã‚Ì“_

//    darray_decode_block_rle(b->da, SBfirst(i), SBlast(i), blocktmp);
  num = darray_decode_block_rle2(b->da, SBfirst(i), SBlast(i), blocktmp2);

  s = SBfirst(i)-1;

// i ‚æ‚è¶‘¤‚ð“Ç‚Ý”ò‚Î‚·
  for (u = 0; u < num; u++) {
    r = blocktmp2[u];
    if (s + abs(r) >= i) break;
    s += abs(r);
  }
  u++;

// rel ‚ÌŒŸõŠJŽn
// i-1‚ÌŽŸ‚©‚çRLE‚ªŽn‚Ü‚Á‚Ä‚¢‚é‚Æ‚Ý‚È‚·
  if (r > 0) {
    r -= i-s-1;
  } else {
    r += i-s-1;
  }
  if (r == 0) {
    printf("??? r == 0\n");
  }
  d = 0; // i-1‚Ì‘Š‘Î[‚³

  while (1) {
    if (r > 0) {
      if (d+1 <= rel && rel <= d+r) { // found
        t = i-1 + rel; // Œ©‚Â‚©‚Á‚½êŠ
        if (t >= n) {
//          if (z != NOTFOUND) {
//            printf("???1\n");
//          }
          return NOTFOUND;
        }
//        if (z != t) {
//          printf("???2\n");
//        }
        return t;
      }
      i += r;
      rel -= r;
    } else {
      if (d-1 >= rel && rel >= d+r) { // found
        t = i-1 - rel; // Œ©‚Â‚©‚Á‚½êŠ
        if (t >= n) {
//          if (z != NOTFOUND) {
//            printf("???3\n");
//          }
          return NOTFOUND;
        }
//        if (z != t) {
//          printf("???4\n");
//        }
        return t;
      }
      i -= r;
      rel -= r;
    }
    if (i > il) break;
    if (u >= num) {
      printf("??? u %ld num %ld\n", u, num);
    }
    r = blocktmp2[u++];
  }

//  if (z != CONTINUE) {
//    printf("???5\n");
//  }
  return CONTINUE;
}

#if 0
i64 search_SB_r(bp *b, i64 i, i64 rel)
{
  i64 j,r,n,il;
  pb *p,x,w;

  n = b->n;
  il = min((SBid(i) + 1) << logSB,n);
  p = &b->B[i>>logD];
  while (i<il) {
    x = *p++;
    j = i & (D-1);
    x <<= j;
    j = D-j;
    while (j>0) {
      w = (x >> (D-ETW)) & ((1<<ETW)-1);
      if (rel >= -ETW && rel <= ETW) {
        r = fwdtbl[((rel+ETW)<<ETW)+w];
        if (r<ETW && r<j) {
          if (i+r >= n) return NOTFOUND;
          return i+r;
        }
      }
      r = min(j,ETW);
      rel -= 2*popCount[w]-r;
      x <<= r;
      i += r;
      j -= r;
    }
  }
  return CONTINUE;
}
#else
i64 search_SB_r(bp *b, i64 i, i64 rel) // b[i..] ‚ðŒŸõ (i ‚ðŠÜ‚Þ)
{
  i64 j,r,n,il;
  pb *p,x,w;
  pb blocktmp[SB/D];
  pb *p2;

  if (b->opt & OPT_COMP_RLE) {
    return search_SB_r_rle(b, i, rel);
  }

  n = b->n;
  il = min(SBlast(i),n-1); // i ‚ðŠÜ‚ÞƒuƒƒbƒN“à‚ÌÅŒã‚Ì“_

//////////////////////////////// *p ‚Íˆ³k‚³‚ê‚Ä‚¢‚È‚¢ƒf[ƒ^‚Ö‚Ìƒ|ƒCƒ“ƒ^
  if (b->opt & OPT_COMP_RLE) {
    darray_decode_block_rle(b->da, SBfirst(i), SBlast(i), blocktmp);
    p = &blocktmp[(i-SBfirst(i))>>logD];
    p2 = &b->B[i>>logD];
  } else {
    p = &b->B[i>>logD];
  }

  while (i<=il) {
  // i ‚ðŠÜ‚ÞƒZƒOƒƒ“ƒg‚Ìƒrƒbƒg‚ðŽæ‚Á‚Ä‚­‚é (ÅãˆÊ‚©‚ç j ƒrƒbƒg)
//////////////////////////////// *p ‚Íˆ³k‚³‚ê‚Ä‚¢‚È‚¢ƒf[ƒ^‚Ö‚Ìƒ|ƒCƒ“ƒ^
//    if (b->opt & OPT_COMP_RLE) {
//      if (*p != *p2) {
//        printf("search_SB_r: p %x p2 %x\n", p, p2);
//    darray_decode_block_rle(b->da, SBfirst(i), SBlast(i), blocktmp);
//      }
//    }
    x = *p++;  p2++;
    j = i & (D-1);
    x <<= j;
    j = D-j;

    while (j >= ETW) { // •\‚ÌƒTƒCƒYˆÈãƒrƒbƒg‚ª‚ ‚é‚Æ‚«
      w = (x >> (D-ETW)) & ((1<<ETW)-1);
      if (rel >= -ETW && rel <= ETW) { // “š‚¦‚ª‚ ‚é‰Â”\«‚ª‚ ‚é‚È‚ç
        r = fwdtbl[((rel+ETW)<<ETW)+w];
        if (r < ETW) { // “š‚¦‚ªŒ©‚Â‚©‚Á‚½
          if (i+r >= n) return NOTFOUND; // ƒrƒbƒg—ñ‚Ì”ÍˆÍŠO‚È‚ç“š‚¦‚Í‚È‚¢
          return i+r;
        }
      }
      r = ETW;
      rel -= 2*popCount[w]-r;
      x <<= r;
      i += r;
      j -= r;
    }

    if (j>0) { // ‚Ü‚¾ x ‚Éƒrƒbƒg‚ªŽc‚Á‚Ä‚¢‚é‚Æ‚«
      w = (x >> (D-ETW)) & ((1<<ETW)-1);
      if (rel >= -ETW && rel <= ETW) {
        r = fwdtbl[((rel+ETW)<<ETW)+w];
        if (r < j) { // “š‚¦‚ª—LŒø‚Èƒrƒbƒg“à‚É‚ ‚é‚È‚ç
          if (i+r >= n) return NOTFOUND;
          return i+r;
        }
      }
      rel -= 2*popCount[w]-j;
      i += j;
    }
  }
  return CONTINUE;
}
#endif


i64 search_MB_r(bp *b, i64 i, i64 td) // i ‚ðŠÜ‚ÞSB‚©‚ç‰E‚ÉŒŸõ.  i ‚ÍSB‚Ìæ“ª. td ‚Í[‚³‚Ìâ‘Î’l
{
  i64 il,d;
  i64 m,M,n;
  pb *B;

  B = b->B;
  n = b->n;

  il = min((MBid(i) + 1) << logMB,n); // i ‚ðŠÜ‚ÞMB‚ÌÅŒã‚ÌˆÊ’u +1
  for (  ;  i < il;  i+=SB) {
#if (SB % RRR != 0)
    d = depth(b,i-1);
#else
    d = fast_depth(b,i-1); // i ‚Ì’¼‘O‚Ì[‚³
#endif
    m = d + b->sm[SBid(i)] - SB; // i ‚ðŠÜ‚ÞSB‚Å‚ÌÅ¬’l
    M = d + b->sM[SBid(i)] - 1;
    if (m <= td && td <= M) {
      return search_SB_r(b,i,td-d); // i ‚ðŠÜ‚ÞSB‚ð’T‚·
    }
  }
  if (i >= n) return NOTFOUND; // ƒrƒbƒg—ñ‚ÌÅŒã‚Ü‚Å’T‚µ‚½
  return CONTINUE; // ‚±‚ÌMB“à‚É‚Í“š‚¦‚ª‚È‚©‚Á‚½
}

///////////////////////////////////////////
//  fwd_excess(bp *b,i64 s, i64 rel)
//    find the leftmost value depth(s)+rel to the right of s (exclusive)
///////////////////////////////////////////
i64 fwd_excess(bp *b,i64 s, i64 rel)
{
  i64 i,n;
  i64 d,td;
  i64 m,M;
  i64 m_ofs;
//  pb *B;
  n = b->n;
//  B = b->B;

  i = s+1;
#if 0
  d = search_SB_r(b,i,rel);
  if (d >= NOTFOUND) return d;

  i = min((SBid(i) + 1) << logSB,n);
  td = depth(b,s) + rel;
  d = search_MB_r(b,i,td);
  if (d >= NOTFOUND) return d;
#else
  if (i != SBfirst(i)) { // i ‚ªƒuƒƒbƒN‚Ì“r’†‚©‚çŽn‚Ü‚éê‡
    d = search_SB_r(b,i,rel); // i ‚ðŠÜ‚ÞƒuƒƒbƒN‚ðŒŸõ
    if (d >= NOTFOUND) return d; // ƒrƒbƒg—ñ‚ÌÅŒã‚Ü‚Ås‚Á‚ÄŒ©‚Â‚©‚ç‚È‚¯‚ê‚Î‚¨‚µ‚Ü‚¢
  }

  td = depth(b,s) + rel; // ’T‚µ‚Ä‚¢‚é[‚³

  i = SBid(i+SB-1) << logSB; // i ‚ðŠÜ‚ÞSB‚ÌŽŸ‚ÌSB‚Ìæ“ª

  if (i != MBfirst(i)) { // MB ‚Ì“r’†‚ÌSB‚©‚çŽn‚Ü‚éê‡
    d = search_MB_r(b,i,td);
    if (d >= NOTFOUND) return d;
  }
#endif

  m_ofs = b->m_ofs;
  i = MBid(s) + m_ofs;
  while (i > 0) {
    if ((i&1) == 0) {
      i++;
      m = b->mm[i];
      M = b->mM[i];
      if (m <= td && td <= M) break;
    }
//    i = PARENT(i);
    i >>= 1;
  }
  if (i < 0) return NOTFOUND;
#if 0
  while (i < m_ofs) {
    i64 k;
    for (k=0; k<W1; k++) {
      m = b->mm[CHILD(i,k)];
      M = b->mM[CHILD(i,k)];
      if (m <= td && td <= M) break;
    }
    i = CHILD(i,k);
  }
#else
  while (i < m_ofs) {
    i <<= 1;
    m = b->mm[i];
    M = b->mM[i];
    if (!(m <= td && td <= M)) i++;
  }
#endif
  i -= m_ofs;
  i <<= logMB;

  d = search_MB_r(b,i,td);
  if (d >= NOTFOUND) return d;
  
  // unexpected (bug)
  printf("fwd_excess: ???\n");
  return -99;

}

#if 0
i64 degree_SB_slow(bp *b, i64 i, i64 t, i64 rel, i64 *ans, i64 ith)
{
  i64 j,r,n,il;
  pb *p,x,w,w2;
  i64 d, deg, v;

  n = t;
  il = min((SBid(i) + 1) << logSB,n);
  d = deg = 0;

  while (i < il) {
    if (getbit(b->B,i)==OP) {
      d++;
    } else {
      d--;
    }
    if (d < rel) {  // reached the end
      if (ith > 0) {
        return NOTFOUND;
      } else {
        *ans = deg;
        return END;
      }
    }
    if (d == rel) {  // found the same depth
      deg++;
      if (deg == ith) {
        *ans = i;
        return FOUND;
      }
    }
    i++;
  }
  *ans = deg;
  return CONTINUE;
}
#endif

i64 degree_SB(bp *b, i64 i, i64 t, i64 rel, i64 *ans, i64 ith)
{
  i64 j,r,n,il;
  pb *p,x,w,w2;
  i64 d, deg, v;
  pb *p2;

  n = t;
  il = min((SBid(i) + 1) << logSB,n);
  d = deg = 0;

//////////////////////////////// *p ‚Íˆ³k‚³‚ê‚Ä‚¢‚È‚¢ƒf[ƒ^‚Ö‚Ìƒ|ƒCƒ“ƒ^
  if (b->opt & OPT_COMP_RLE) {
    pb blocktmp[SB/D];
    darray_decode_block_rle(b->da, SBfirst(i), SBlast(i), blocktmp);
    p = &blocktmp[(i-SBfirst(i))>>logD];
    p2 = &b->B[i>>logD];
  } else {
    p = &b->B[i>>logD];
  }

  while (i < il) {
//////////////////////////////// *p ‚Íˆ³k‚³‚ê‚Ä‚¢‚È‚¢ƒf[ƒ^‚Ö‚Ìƒ|ƒCƒ“ƒ^
//    if (b->opt & OPT_COMP_RLE) {
//      if (*p != *p2) {
//        printf("degree_SB: p %x p2 %x\n", p, p2);
//      }
//    }
    x = *p++;  p2++;
    j = i & (D-1);
    x <<= j;
    j = min(D-j,il-i);
    while (j>0) {
      w = (x >> (D-ETW)) & ((1<<ETW)-1);
      w2 = 0;
      if (j < ETW || il-i < ETW) {
        r = max(ETW-j,ETW-(il-i));
        w2 = (1<<r)-1;
      }
      v = minmaxtbl_v[0][w | w2];
      if (d + v < rel) {
        if (ith > 0) {
#if 0
          for (r = 0; r < ETW; r++) {
            if (w & 0x80) {
              d++;
            } else {
              d--;
              if (d < rel) break;
              if (d == rel) {
                ith--;
                if (ith == 0) {
                  *ans = i + r;
                  return FOUND;
                }
              }
            }
          w <<= 1;
        }
        return NOTFOUND;
#else
          r = childtbl2[rel-d+ETW][ith-1][w];
          if (r >= 0) {
            *ans = i + r;
            return FOUND;
          }
          return NOTFOUND;
#endif
        }
        r = ETW-1-minmaxtbl_i[0][w | w2];
        w2 = (1<<r)-1;
        deg += degtbl2[((rel-d+ETW)<<ETW) + (w & (~w2))];
        *ans = deg;
        return END;
      }
      if (d + v == rel) {
        r = degtbl[w | w2];
        deg += r;
        if (ith > 0) {
          if (ith <= r) {
            *ans = i + childtbl[((ith-1)<<ETW) + (w | w2)];
            return FOUND;
          }
          ith -= r;
        }
      }

      r = min(j,ETW);
      d += 2*popCount[w]-r;
      x <<= r;
      i += r;
      j -= r;
    }
  }

  *ans = deg;
  return CONTINUE;
}

i64 degree_MB(bp *b, i64 i, i64 t, i64 td, i64 *ans,i64 ith)
{
  i64 il,d;
  i64 m,n,r;
//  pb *B;
  i64 deg,degtmp;

  d = 0;
//  B = b->B;
  n = t;

  il = min((MBid(i) + 1) << logMB,n);
  deg = 0;
  for (  ;  i+SB-1 < il;  i+=SB) {
#if (SB % RRR != 0)
    d = depth(b,i-1);
#else
    d = fast_depth(b,i-1);
#endif
    m = d + b->sm[SBid(i)] - SB;
    if (m < td) {
      r = degree_SB(b,i,n,td-d,&degtmp,ith);
      if (ith > 0) {
        if (r == NOTFOUND) return NOTFOUND;
        *ans = degtmp;
        return FOUND;
      } else {
        *ans = deg + degtmp;
        return END;
      }
    }
    if (m == td) {
      if (ith > 0) {
        if (ith <= b->sd[SBid(i)]) break;
        ith -= b->sd[SBid(i)];
      }
      deg += b->sd[SBid(i)];
    }
  }
  if (i < il) {
#if (SB % RRR != 0)
    d = depth(b,i-1);
#else
    d = fast_depth(b,i-1);
#endif
    r = degree_SB(b,i,n,td-d,&degtmp,ith);
    if (ith > 0) {
      if (r == NOTFOUND) return NOTFOUND;
      if (r == FOUND) {
        *ans = degtmp;
        return FOUND;
      }
    } else {
      deg += degtmp;
    }
  }
  *ans = deg;
  if (i >= n) return END;
  return CONTINUE;
}

static i64 partition_range(i64 s,i64 t)
{
  i64 h;

  printf("partition [%d,%d] => ",s,t);
  h = 1;
  while (s <= t) {
    if (s & h) {
      if (s+h-1 <= t) {
        printf("[%d,%d] ",s,s+h-1);
        s += h;
      }
    } else {
      if (s+h > t) break;
    }
    h <<= 1;
  }
  while (h > 0) {
    if (s+h-1 <= t) {
      printf("[%d,%d] ",s,s+h-1);
      s += h;
    }
    h >>= 1;
  }
  printf("\n");
    return 0; //editted by Lee
}




///////////////////////////////////////////
//  fast_degree(bp *b,i64 s, i64 t, i64 ith)
//    returns the number of children of s, to the left of t
//    returns the position of (ith)-th child of s if (ith > 0)
///////////////////////////////////////////
i64 fast_degree(bp *b,i64 s, i64 t, i64 ith)
{
  i64 i,j,n;
  i64 d,td;
  i64 m_ofs;
//  pb *B;
  i64 deg,degtmp;
  i64 sm,tm,ss,h;

  n = t;  
//  B = b->B;

  deg = 0;

  i = s+1;


#if 1
  d = degree_SB(b,i,n,0,&degtmp,ith);
  if (ith > 0) {
    if (d == NOTFOUND) return -1;
    if (d == FOUND) return degtmp;
    ith -= degtmp;
  }
  if (d == END) return degtmp;
  deg += degtmp;
  //i = min((SBid(i) + 1) << logSB,n);
  i = (SBid(i) + 1) << logSB;
  if (i >= n) return deg;
#else
  if (i != SBfirst(i)) {
    d = degree_SB(b,i,n,0,&degtmp,ith);
    if (ith > 0) {
      if (d == NOTFOUND) return -1;
      if (d == FOUND) return degtmp;
      ith -= degtmp;
    }
    if (d == END) return degtmp;
      deg += degtmp;
  }
  i = SBid(i+SB-1) << logSB;
#endif

  td = depth(b,s);


  if (i != MBfirst(i)) {
    d = degree_MB(b,i,n,td,&degtmp,ith);
    if (ith > 0) {
      if (d == NOTFOUND) return -1;
      if (d == FOUND) return degtmp;
      ith -= degtmp;
      deg += degtmp;
    } else {
      deg += degtmp;
      if (d == END) {
        return deg;
      }
    }
  }

#if 0
  // sequential search

  sm = MBid(i+MB-1);
  tm = MBid((n-1)+1)-1; // the rightmost MB fully contained in [0,n-1]

  m_ofs = b->m_ofs;
  sm += m_ofs;  tm += m_ofs;
  for (i=sm; i<=tm; i++) {
    if (b->mm[i] < td) {
      break;
    }
    if (b->mm[i] == td) {
      if (ith > 0) {
        if (ith <= b->md[i]) break;
        ith -= b->md[i];
      }
      deg += b->md[i];
    }
  }
  ss = i - m_ofs;
#else
  sm = MBid(i+MB-1);
  tm = MBid((n-1)+1)-1; // the rightmost MB fully contained in [0,n-1]

  m_ofs = b->m_ofs;
  sm += m_ofs;  tm += m_ofs;
  ss = sm;

  //partition_range(sm,tm);

  //printf("partition [%d,%d] => ",sm,tm);
  h = 1;
  while (sm <= tm) {
    if (sm & h) {
      if (sm+h-1 <= tm) {
        //printf("[%d,%d] ",sm,sm+h-1);
        j = sm / h;
        if (b->mm[j] < td) {
          h >>= 1;
          break;
        }
        if (b->mm[j] == td) {
          if (ith > 0) {
            if (ith <= b->md[j]) {
              h >>= 1;
              break;
            }
            ith -= b->md[j];
          }
          deg += b->md[j];
        }
        sm += h;
      }
    } else {
      if (sm+h > tm) break;
    }
    h <<= 1;
  }
  while (h > 0) {
    if (sm+h-1 <= tm) {
      //printf("[%d,%d] ",sm,sm+h-1);
      j = sm / h;
      if (ith > 0) {
        if (b->mm[j] >= td) {
          if (b->mm[j] == td) {
            if (ith > b->md[j]) {
              ith -= b->md[j];
              sm += h;
            } else {
              deg += b->md[j];
            }
          } else {
            sm += h;
          }
        }
      } else {
        if (b->mm[j] >= td) {
          if (b->mm[j] == td) {
            deg += b->md[j];
          }
          sm += h;
        }
      }
    }
    h >>= 1;
  }
  //printf("\n");
  ss = sm;

  ss -= m_ofs;

#endif

  ss <<= logMB;

  d = degree_MB(b,ss,n,td,&degtmp,ith);
  if (ith > 0) {
    if (d == NOTFOUND) return -1;
    if (d == FOUND) return degtmp;
  }
  deg += degtmp;
  if (d == END) return deg;
  return deg;
  
  // unexpected (bug)
  printf("degree: ???\n");
  return -99;

}

i64 search_SB_l_rle(bp *b, i64 i, i64 rel) // b[..i] ‚ðŒŸõ (i ‚ðŠÜ‚Þ)
{
  i64 j,r,n,il;
  pb *p,x,w;
  pb *p2;
  i64 blocktmp2[SB];
  i64 num;
  i64 s,t,u,d;
  i64 z;
  i64 ii;

//  if (i == 979456) {
//    printf("hoge\n");
//  }
//  ii = i;
//  z = search_SB_r(b, i, rel);

  n = b->n;
  il = min(SBlast(i),n-1); // i ‚ðŠÜ‚ÞƒuƒƒbƒN“à‚ÌÅŒã‚Ì“_

//  num = darray_decode_block_rle2(b->da, SBfirst(i), SBlast(i), blocktmp2);
  num = darray_decode_block_rle2(b->da, SBfirst(i), i, blocktmp2);

  s = SBfirst(i)-1;

// i ‚æ‚è¶‘¤‚ð“Ç‚Ý”ò‚Î‚·
  for (u = 0; u < num; u++) {
    r = blocktmp2[u];
    if (s + abs(r) >= i) break;
    s += abs(r);
  }
//  u++;

// rel ‚ÌŒŸõŠJŽn
// i‚ÅRLE‚ªI‚í‚Á‚Ä‚¢‚é‚Æ‚Ý‚È‚·
  if (r > 0) {
    r = i-s;
  } else {
    r = -(i-s);
  }
  if (r == 0) {
    printf("??? r == 0\n");
  }
  d = 0; // i+1‚Ì‘Š‘Î[‚³

  il = SBfirst(i);

  while (1) {
    if (r < 0) {
      if (d+1 <= rel && rel <= d-r) { // found
        t = i - rel; // Œ©‚Â‚©‚Á‚½êŠ
        if (t < 0) {
//          if (z != NOTFOUND) {
//            printf("???1\n");
//          }
          return NOTFOUND;
        }
//        if (z != t) {
//          printf("???2\n");
//        }
        return t;
      }
      i += r;
      rel += r;
    } else {
      if (d-1 >= rel && rel >= d-r) { // found
        t = i + rel; // Œ©‚Â‚©‚Á‚½êŠ
        if (t < -1) {
//          if (z != NOTFOUND) {
//            printf("???3\n");
//          }
          return NOTFOUND;
        }
//        if (z != t) {
//          printf("???4\n");
//        }
        return t;
      }
      i -= r;
      rel += r;
    }
    if (i < il) break;
    if (u <= 0) {
      printf("??? u %ld\n", u);
    }
    r = blocktmp2[--u];
  }

//  if (z != CONTINUE) {
//    printf("???5\n");
//  }
  if (i < 0) return NOTFOUND;
  return CONTINUE;
}


i64 search_SB_l(bp *b, i64 i, i64 rel)
{
  i64 j,r,il;
  pb *p,x,w;
  pb *p2;
  pb blocktmp[SB/D];

  if (b->opt & OPT_COMP_RLE) {
    return search_SB_l_rle(b, i, rel);
  }

  il = SBfirst(i);

//////////////////////////////// *p ‚Íˆ³k‚³‚ê‚Ä‚¢‚È‚¢ƒf[ƒ^‚Ö‚Ìƒ|ƒCƒ“ƒ^
  if (b->opt & OPT_COMP_RLE) {
    darray_decode_block_rle(b->da, SBfirst(i), SBlast(i), blocktmp);
    p = &blocktmp[(i-SBfirst(i))>>logD];
    p2 = &b->B[i>>logD];
  } else {
    p = &b->B[i>>logD];
  }

  while (i>=il) {
//////////////////////////////// *p ‚Íˆ³k‚³‚ê‚Ä‚¢‚È‚¢ƒf[ƒ^‚Ö‚Ìƒ|ƒCƒ“ƒ^
//    if (b->opt & OPT_COMP_RLE) {
//      if (*p != *p2) {
//        printf("search_SB_l: p %x p2 %x\n", p, p2);
//      }
//    }
    x = *p--;  p2--;
    j = (i & (D-1))+1;
    x >>= D-j;
    while (j>0) {
      w = x & ((1<<ETW)-1);
      if (rel >= -ETW && rel <= ETW) {
        r = bwdtbl[((rel+ETW)<<ETW)+w];
        if (r<ETW && r<j) {
          if (i-r < 0) return NOTFOUND;
          return i-r-1;
        }
      }
      r = min(j,ETW);
      rel += 2*popCount[w]-r;
      x >>= r;
      i -= r;
      j -= r;
    }

  }
  if (i < 0) return NOTFOUND;
  return CONTINUE;
}

i64 search_MB_l(bp *b, i64 i, i64 td)
{
  i64 il,d;
  i64 m,M;
//  pb *B;

#if 0
  if (i % SB != SB-1) {
    printf("search_MB_l:error!!! i=%d SB=%d\n",i,i%SB);
  }
#endif
//  B = b->B;

  il = MBfirst(i);
  for (  ;  i >= il;  i-=SB) {
#if (SB % RRR != 0)
    d = depth(b,i-SB);
#else
    d = fast_depth(b,i-SB);
#endif
    m = d + b->sm[SBid(i)] - SB;
    M = d + b->sM[SBid(i)] - 1;
    if (m <= td && td <= M) {
#if (SB % RRR != 0)
      d = depth(b,i);
#else
      d = fast_depth(b,i);
#endif
      if (d == td) return i;
#if 0
{
  i64 z, z2;
  z = search_SB_l(b,i,td-d);
  z2 = search_SB_l_rle(b,i,td-d);
  if (z != z2) {
    printf("z %ld z2 %ld\n", z, z2);
    z = search_SB_l(b,i,td-d);
    z2 = search_SB_l_rle(b,i,td-d);
  }
}
#endif
      return search_SB_l(b,i,td-d);
    }
  }
  return CONTINUE;
}

///////////////////////////////////////////
//  bwd_excess(bp *b,i64 s, i64 rel)
//    find the rightmost value depth(s)+rel to the left of s (exclusive)
///////////////////////////////////////////
i64 bwd_excess(bp *b,i64 s, i64 rel)
{
  i64 i,n;
  i64 d,td;
  i64 m,M;
  i64 m_ofs;
  int f;
//  pb *B;

  n = b->n;
//  B = b->B;

  i = s;
#if 0
{
  i64 z, z2;
  z = search_SB_l(b,i,rel);
  z2 = search_SB_l_rle(b,i,rel);
  if (z != z2) {
    printf("z %ld z2 %ld\n", z, z2);
    z = search_SB_l(b,i,rel);
    z2 = search_SB_l_rle(b,i,rel);
  }
}
#endif
  d = search_SB_l(b,i,rel);
  if (d >= NOTFOUND) return d;

  i = SBfirst(i) -1;

  td = depth(b,s) + rel;

  d = search_MB_l(b,i,td);
  if (d >= NOTFOUND) return d;

  f = 0;
  m_ofs = b->m_ofs;
  i = (s>>logMB) + m_ofs;
  while (i > 0) {
    if ((i&1) == 1) {
      i--;
      m = b->mm[i];
      M = b->mM[i];
      if (m <= td && td <= M) {
//        f = 1;
//        i >>= 1;
        break;
      }
    }
    i >>= 1;
  }
//  if (f == 0 && i == 0) { // ‚È‚ñ‚©ƒoƒO‚Á‚Ä‚é?
  if (i == 0) { // ‚È‚ñ‚©ƒoƒO‚Á‚Ä‚é?
    if (td == 0) return -1;
    else return NOTFOUND;
  }

  while (i < m_ofs) {
    i = i*2 + 1;
    m = b->mm[i];
    M = b->mM[i];
    if (!(m <= td && td <= M)) i--;
  }
  i -= m_ofs;
  i = ((i+1)<<logMB)-1;

  d = search_MB_l(b,i,td);
  if (d >= NOTFOUND) return d;

  // unexpected (bug)
  printf("bwd_excess: ???\n");
  return -99;

}

i64 rmq_SB(bp *b, i64 s, i64 t, i64 opt, i64 *dm)
{
  i64 i,d;
  i64 is,ds;
  pb *p,x;
  i64 lr;
  i64 op;
  i64 j,r,v;
  pb w,w2;

  lr = 0;  if (opt & OPT_RIGHT) lr = 1;
  op = opt & (OPT_RIGHT | OPT_MAX);

  is = s;  ds = d = 0;
  i = s+1;

#if SB >= ETW
//////////////////////////////// *p ‚Íˆ³k‚³‚ê‚Ä‚¢‚È‚¢ƒf[ƒ^‚Ö‚Ìƒ|ƒCƒ“ƒ^
  if (b->opt & OPT_COMP_RLE) {
    pb blocktmp[SB/D];
    darray_decode_block_rle(b->da, SBfirst(i), SBlast(i), blocktmp);
    p = &blocktmp[(i-SBfirst(i))>>logD];
  } else {
    p = &b->B[i>>logD];
  }

  while (i <= t) {
//////////////////////////////// *p ‚Íˆ³k‚³‚ê‚Ä‚¢‚È‚¢ƒf[ƒ^‚Ö‚Ìƒ|ƒCƒ“ƒ^
    x = *p++;
    j = i & (D-1);
    x <<= j;
    j = min(D-j,t-i+1);
    while (j>0) {
      w = (x >> (D-ETW)) & ((1<<ETW)-1);
      w2 = 0;
      if (j < ETW || t-i < ETW-1) {
        r = max(ETW-j,ETW-1-(t-i));
        w2 = (1<<r)-1;
      }

      if (op & OPT_MAX) {
        v = minmaxtbl_v[op][w & (~w2)];
        if (d + v + lr > ds) {
          ds = d + v;
          is = i + minmaxtbl_i[op][w & (~w2)];
        }
      } else {
        v = minmaxtbl_v[op][w | w2];
        if (d + v < ds +lr) {
          ds = d + v;
          is = i + minmaxtbl_i[op][w | w2];
        }
      }

      r = min(j,ETW);
      d += 2*popCount[w]-r;
      x <<= r;
      i += r;
      j -= r;
    }
  }
#else
  while (i <= t) {
//////////////////////////////// B ‚Íˆ³k‚³‚ê‚Ä‚¢‚È‚¢ƒf[ƒ^‚Ö‚Ìƒ|ƒCƒ“ƒ^
    if (getbit(b->B,i)==OP) {
      d++;
      if (op & OPT_MAX) {
        if (d + lr > ds) {
          ds = d;  is = i;
        }
      }
    } else {
      d--;
      if (!(op & OPT_MAX)) {
        if (d < ds + lr) {
          ds = d;  is = i;
        }
      }
    }
    i++;
  }
#endif
  *dm = ds;
  return is;
}

i64 rmq_MB(bp *b, i64 s, i64 t, i64 opt, i64 *dm)
{
  i64 i,d,m;
  i64 mi,md;
  i64 lr;

  lr = 0;  if (opt & OPT_RIGHT) lr = 1;

  md = *dm;  mi = -1;
  for (i = s;  i <= t;  i++) {
#if (SB % RRR != 0)
    d = depth(b,(i<<logSB)-1);
#else
    d = fast_depth(b,(i<<logSB)-1);
#endif
    if (opt & OPT_MAX) {
      m = d + b->sM[i] - 1;
      if (m + lr > md) {
        md = m;  mi = i;
      }
    } else {
      m = d + b->sm[i] - SB;
      if (m < md + lr) {
        md = m;  mi = i;
      }
    }
  }
  *dm = md;
  return mi;
}

#if 0
i64 fast_rmq(bp *b, i64 s, i64 t,int opt)
{
  i64 h,w;
  i64 i,j,m;
  i64 sl,tl,ss,tt;
  i64 ds,is;

  sl = (s >> logFLB)+1;
  tl = (t >> logFLB)-1;

  if (tl <= sl) {
    return std_rmq(b,s,t,0);
  }

  ss = s;  tt = (sl << logFLB) -1;
  is = std_rmq(b,ss,tt,0);
  ds = depth(b,is);

  m = (b->n+FLB-1)/FLB;
  h = blog(tl-sl);
  w = 1 << (h-1);

  i = b->rmq_tbl[(h-1) * m + sl];
  j = depth(b,i);
  if (j < ds) {
    ds = j;  is = i;
  }

  i = b->rmq_tbl[(h-1) * m + tl - w + 1];
  j = depth(b,i);
  if (j < ds) {
    ds = j;  is = i;
  }

  ss = (tl+1) << logFLB;  tt = t;
  i = std_rmq(b,ss,tt,0);
  j = depth(b,i);
  if (j < ds) {
    ds = j;  is = i;
  }

  return is;
}
#endif


i64 rmq(bp *b, i64 s, i64 t, i64 opt)
{
  i64 ss, ts;  // SB index of s and t
  i64 sm, tm;  // MB index of s and t
  i64 ds;   // current min value
  i64 is;   // current min index
  i64 ys;   // type of current min index
               // 0: is is the index of min
               // 1: is is the SB index
               // 2: is is the MB index
  i64 m_ofs;
  i64 i,j,il,d,n;
  i64 dm;
  i64 lr;

  lr = 0;  if (opt & OPT_RIGHT) lr = 1;

  n = b->n;
  if (s < 0 || t >= n || s > t) {
    printf("rmq: error s=%d t=%d n=%d\n",s,t,n);
    return -1;
  }
  if (s == t) return s;


  ////////////////////////////////////////////////////////////
  // search the SB of s
  ////////////////////////////////////////////////////////////

  il = min(SBlast(s),t);
  is = rmq_SB(b,s,il,opt,&dm);
  if (il == t) {  // scan reached the end of the range
    return is;
  }
  ds = depth(b,s) + dm;  ys = 0;

  ////////////////////////////////////////////////////////////
  // search the MB of s
  ////////////////////////////////////////////////////////////

  ss = SBid(s) + 1;
  il = min(SBid(MBlast(s)),SBid(t)-1);
  dm = ds;
  j = rmq_MB(b,ss,il,opt,&dm);
  //if (dm < ds + lr) {
  if (j >= 0) {
    ds = dm;  is = j;  ys = 1;
  }

  ////////////////////////////////////////////////////////////
  // search the tree
  ////////////////////////////////////////////////////////////

  sm = MBid(s) + 1;
  tm = MBid(t) - 1;

  if ((b->opt & OPT_FAST_LCA) && (opt == 0)) {
    i64 m,h,w;
    
    if (sm == tm) {
      j = b->mm[sm+b->m_ofs];
      if (j < ds) {
        ds = j;  
        is = b->rmq_tbl[sm];
        ys = 0;
      }
    } else if (sm+1 == tm) {
      j = b->mm[sm+b->m_ofs];
      if (j < ds) {
        ds = j;  is = b->rmq_tbl[sm];  ys = 0;
      }
      
      j = b->mm[tm+b->m_ofs];
      if (j < ds) {
        ds = j;  is = b->rmq_tbl[tm];  ys = 0;
      }
    } else if (sm < tm) {
      m = (n+MB-1) >> logMB;
      h = blog(tm-sm);
      w = 1 << (h-1);
      
      i = b->rmq_tbl[(h-1) * m + sm];
      j = depth(b,i);
      if (j < ds) {
        ds = j;  is = i;  ys = 0;
      }
      
      i = b->rmq_tbl[(h-1) * m + tm - w + 1];
      j = depth(b,i);
      if (j < ds) {
        ds = j;  is = i;  ys = 0;
      }
    }

  } else {
#if 0
  // sequential search
  m_ofs = b->m_ofs;
  sm += m_ofs;  tm += m_ofs;
  for (i=sm; i<=tm; i++) {
    if (opt & OPT_MAX) {
      if (b->mM[i] + lr > ds) {
        ds = b->mM[i];  is = i - m_ofs;  ys = 2;
      }
    } else {
      if (b->mm[i] < ds + lr) {
        ds = b->mm[i];  is = i - m_ofs;  ys = 2;
      }
    }
  }

#else
  if (sm <= tm) {
    i64 h;
    h = blog(sm ^ tm);

    m_ofs = b->m_ofs;
    sm += m_ofs;  tm += m_ofs;

    if (opt & OPT_MAX) {
      if (b->mM[sm] + lr > ds) {
        ds = b->mM[sm];  is = sm;  ys = 2;
      }
      for (i=0; i<=h-2; i++) {
        j = sm>>i;
        if ((j&1) == 0) {
          if (b->mM[j+1] + lr > ds) {
            ds = b->mM[j+1];  is = j+1;  ys = 2;
          }
        }
      }
      for (i=h-2; i>=0; i--) {
        j = tm>>i;
        if ((j&1)==1) {
          if (b->mM[j-1] + lr > ds) {
            ds = b->mM[j-1];  is = j-1;  ys = 2;
          }
        }
      }
      if (b->mM[tm] + lr > ds) {
        ds = b->mM[tm];  is = tm;  ys = 2;
      }
      if (ys == 2) {
        while (is < m_ofs) {
          is <<= 1;
          if (b->mM[is+1] + lr > b->mM[is]) is++;
        }
        is -= m_ofs;
      }
    } else { // MIN
      if (b->mm[sm] < ds + lr) {
        ds = b->mm[sm];  is = sm;  ys = 2;
      }
      for (i=0; i<=h-2; i++) {
        j = sm>>i;
        if ((j&1) == 0) {
          if (b->mm[j+1] < ds + lr) {
            ds = b->mm[j+1];  is = j+1;  ys = 2;
          }
        }
      }
      for (i=h-2; i>=0; i--) {
        j = tm>>i;
        if ((j&1)==1) {
          if (b->mm[j-1] < ds + lr) {
            ds = b->mm[j-1];  is = j-1;  ys = 2;
          }
        }
      }
      if (b->mm[tm] < ds + lr) {
        ds = b->mm[tm];  is = tm;  ys = 2;
      }
      if (ys == 2) {
        while (is < m_ofs) {
          is <<= 1;
          if (b->mm[is+1] < b->mm[is] + lr) is++;
        }
        is -= m_ofs;
      }
    }
  }

#endif
  }

  ////////////////////////////////////////////////////////////
  // search the MB of t
  ////////////////////////////////////////////////////////////


  ts = max(SBid(MBfirst(t)),SBid(s)+1);
  il = SBid(SBfirst(t)-1);
  dm = ds;
  j = rmq_MB(b,ts,il,opt,&dm);
  //if (dm < ds + lr) {
  if (j >= 0) {
    ds = dm;  is = j;  ys = 1;
  }

  ////////////////////////////////////////////////////////////
  // search the SB of t
  ////////////////////////////////////////////////////////////

  i = SBfirst(t);
  j = rmq_SB(b,i,t,opt,&dm);
  d = depth(b,i) + dm;
  if (opt & OPT_MAX) {
    if (d + lr > ds) {
      ds = d;  is = j;  ys = 0;
    }
  } else {
    if (d < ds + lr) {
      ds = d;  is = j;  ys = 0;
    }
  }

  ////////////////////////////////////////////////////////////
  // search the rest
  ////////////////////////////////////////////////////////////

  if (ys == 2) {
    ss = SBid(is << logMB);
    il = SBid(MBlast(is << logMB));
    if (opt & OPT_MAX) {
      dm = -n-1;
    } else {
      dm = n+1;
    }
    j = rmq_MB(b,ss,il,opt,&dm);
    ds = dm;  is = j;  ys = 1;
  }

  if (ys == 1) {
    ss = is << logSB;
    il = SBlast(is << logSB);
    is = rmq_SB(b,ss,il,opt,&dm);
    //ds = depth(b,ss) + dm;  ys = 0;
  }

  return is;
}

