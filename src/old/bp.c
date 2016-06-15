#include <stdio.h>
#include <stdlib.h>
#include "typedefbp.h"
#include "darray.h"
#include "bp.h"

//#define CHECK
#define RANDOM

i64 msize=0;
#define mymalloc(p,n,f) {p = malloc((n)*sizeof(*p)); msize += (f)*(n)*sizeof(*p); /* if (f) printf("malloc %d bytes at line %d total %d\n",(n)*sizeof(*p),__LINE__,msize);  */ if ((p)==NULL) {printf("not enough memory (%d bytes) in line %d\n",msize,__LINE__); exit(1);};}



i64 naive_depth(bp *b, i64 s)
{
  i64 i,d;
  if (s < 0) return 0;
  d = 0;
  for (i=0; i<=s; i++) {
//////////////////////////////// B は圧縮されていないデータへのポインタ
    if (getbit(b->B,i)==OP) {
      d++;
    } else {
      d--;
    }
  }
  return d;
}

void printbp(bp *b, i64 s, i64 t)
{
  i64 i,c,d;
  d = 0;
  for (i=s; i<=t; i++) {
    if (getbit(b->B,i)==OP) {
      c = '(';
      d++;
    } else {
      c = ')';
      d--;
    }
    putchar(c);
  }
}

i64 *matchtbl,*parenttbl;
void make_naivetbl(pb *B,i64 n)
{
  i64 i,v;
  i64 *stack,s;
  i64 maxdepth;

  maxdepth = 1 << 20;

  mymalloc(matchtbl,n,0);
  mymalloc(parenttbl,n,0);
  mymalloc(stack,maxdepth,0);

  for (i=0; i<n; i++) matchtbl[i] = -1;
  //for (i=0; i<n; i++) parenttbl[i] = -1;

  s = 0;
  v = 0;
  stack[0] = -1;
  for (i=0; i<n; i++) {
    if (getbit(B,i)==OP) {
      v++;
      if (v >= maxdepth) {
        printf("exceed maxdepth=%d\n",maxdepth);
        exit(1);
      }
      if (v > 0) {
        parenttbl[i] = stack[v-1];
        stack[v] = i;
      }
    } else {
      if (v > 0) {
        matchtbl[stack[v]] = i;  // close
        matchtbl[i] = stack[v];   // open
      }
      v--;
    }
  }

  free(stack);
}

int popCount[1<<ETW];
int fwdtbl[(2*ETW+1)*(1<<ETW)];
int bwdtbl[(2*ETW+1)*(1<<ETW)];
#if 0
int mi64bl_li[1<<ETW], mi64bl_lv[1<<ETW];
int mi64bl_ri[1<<ETW], mi64bl_rv[1<<ETW];
int maxtbl_li[1<<ETW], maxtbl_lv[1<<ETW];
int maxtbl_ri[1<<ETW], maxtbl_rv[1<<ETW];
#endif
int minmaxtbl_i[4][1<<ETW], minmaxtbl_v[4][1<<ETW];
int degtbl[1<<ETW];
int degtbl2[(2*ETW+1)*(1<<ETW)];
int childtbl[(ETW)*(1<<ETW)];
int childtbl2[2*ETW+1][ETW][(1<<ETW)];
int depthtbl[(2*ETW+1)*(1<<ETW)];

void bp_make_matchtbl(void)
{
  i64 i,j,x,r;
  i64 m,M;
  pb buf[1];
  i64 deg;

  for (x = 0; x < (1<<ETW); x++) {
    setbits(buf,0,ETW,x);
    for (r=-ETW; r<=ETW; r++) fwdtbl[((r+ETW)<<ETW)+x] = ETW;
    for (r=-ETW; r<=ETW; r++) bwdtbl[((r+ETW)<<ETW)+x] = ETW;
    for (r=-ETW; r<=ETW; r++) degtbl2[((r+ETW)<<ETW)+x] = 0;
    for (r=-ETW; r<=ETW; r++) depthtbl[((r+ETW)<<ETW)+x] = 0;

    r = 0;
    for (i=0; i<ETW; i++) {
      if (getbit(buf,i)==OP) {
        r++;
      } else {
        r--;
      }
      if (fwdtbl[((r+ETW)<<ETW)+x] == ETW) fwdtbl[((r+ETW)<<ETW)+x] = i;
    }

    r = 0;
    for (i=ETW-1; i>=0; i--) {
      if (getbit(buf,i)==OP) {
        r--;
      } else {
        r++;
      }
      if (bwdtbl[((r+ETW)<<ETW)+x] == ETW) bwdtbl[((r+ETW)<<ETW)+x] = ETW-1-i;
    }

    r = 0;
    for (i=0; i<ETW; i++) {
      if (getbit(buf,i)==OP) {
        r++;
      } else {
        r--;
      }
      depthtbl[((r+ETW)<<ETW)+x] += (1<<(ETW-1));
    }

    r = 0;
    for (i=0; i<ETW; i++) {
      if (getbit(buf,i)==OP) r++;
    }
    popCount[x] = r;

    r = 0;  m = 0;  M = 0;
    m = ETW+1;  M = -ETW-1;
    //maxtbl_lv[x] = -ETW-1;
    //mi64bl_lv[x] = ETW+1;
    minmaxtbl_v[OPT_MAX | OPT_LEFT][x] = -ETW-1;
    minmaxtbl_v[OPT_MIN | OPT_LEFT][x] = ETW+1;
    deg = 0;
    for (i=0; i<ETW; i++) {
      if (getbit(buf,i)==OP) {
        r++;
        if (r > M) {
          M = r;  
          //maxtbl_li[x] = i;  maxtbl_lv[x] = r;
          minmaxtbl_i[OPT_MAX | OPT_LEFT][x] = i;
          minmaxtbl_v[OPT_MAX | OPT_LEFT][x] = r;
        }
      } else {
        r--;
        if (r == m) {
          deg++;
          childtbl[((deg-1)<<ETW) + x] = i;
        }
        if (r < m) {
          m = r;
          //mi64bl_li[x] = i;  mi64bl_lv[x] = r;
          minmaxtbl_i[OPT_MIN | OPT_LEFT][x] = i;
          minmaxtbl_v[OPT_MIN | OPT_LEFT][x] = r;
          deg = 1;
          childtbl[((deg-1)<<ETW) + x] = i;
        }
      }
      if (r <= m) degtbl2[((r+ETW)<<ETW)+x]++;
    }
    degtbl[x] = deg;

    r = 0;  m = 0;  M = 0;
    //maxtbl_rv[x] = -ETW-1;
    //mi64bl_rv[x] = ETW+1;
    minmaxtbl_v[OPT_MAX | OPT_RIGHT][x] = -ETW-1;
    minmaxtbl_v[OPT_MIN | OPT_RIGHT][x] = ETW+1;
    for (i=0; i<ETW; i++) {
      if (getbit(buf,i)==OP) {
        r++;
        if (r >= M) {
          M = r;  
          //maxtbl_ri[x] = i;  maxtbl_rv[x] = r;
          minmaxtbl_i[OPT_MAX | OPT_RIGHT][x] = i;
          minmaxtbl_v[OPT_MAX | OPT_RIGHT][x] = r;
        }
      } else {
        r--;
        if (r <= m) {
          m = r;
          //mi64bl_ri[x] = i;  mi64bl_rv[x] = r;
          minmaxtbl_i[OPT_MIN | OPT_RIGHT][x] = i;
          minmaxtbl_v[OPT_MIN | OPT_RIGHT][x] = r;
        }
      }
    }

    for (i = 0; i < ETW; i++) {
      for (j = -ETW; j <= ETW; j++) {
        childtbl2[j+ETW][i][x] = -1;
      }
    }

    for (j=-ETW; j<=ETW; j++) {
      i64 ith;
      ith = 0;
      r = 0;
      for (i = 0; i < ETW; i++) {
        if (getbit(buf,i)==OP) {
          r++;
        } else {
          r--;
          if (r < j) break;
          if (r == j) {
            ith++;
            childtbl2[j+ETW][ith-1][x] = i;
          }
        }
      }
    }
  }

}

i64 bp_construct(bp *b,i64 n, pb *B, i64 opt)
{
  i64 i,j,d;
  i64 m,M,ds, m2, M2;
  i64 ns,nm;
  sb_type *sm, *sM;
  sb_type *sd;
  mb_type *mm, *mM;
  mb_type *md;

  i64 m_ofs;
  i64 r, r2; // # of minimum values
  int opt2;

#if 0
  if (SB % D != 0) {
    printf("warning: SB=%d should be a multiple of D=%d\n",SB,D);
    // not necessarily?
  }
  if (SB % RRR != 0) {
    printf("warning: SB=%d should be a multiple of RRR=%d\n",SB,RRR);
  }
#endif

  opt2 = opt & ~OPT_COMP_RLE; // まだ圧縮されていない

  b->B = B;
  b->n = n;
  b->opt = opt2;
  b->idx_size = 0;
  mymalloc(b->da,1,0);
  darray_construct(b->da,n,B, opt & OPT_FAST_PREORDER_SELECT);
//  b->idx_size += b->da->idx_size;
//  printf("preorder rank/select table: %d bytes (%1.2f bpc)\n",
//         b->da->idx_size,(double)b->da->idx_size*8/n);

//  bp_make_matchtbl();

  ns = (n+SB-1)/SB;
  mymalloc(sm, ns, 0);  b->idx_size += ns * sizeof(*sm);
  mymalloc(sM, ns, 0);  b->idx_size += ns * sizeof(*sM);
  b->sm = sm;
  b->sM = sM;
  if (opt & OPT_DEGREE) {
    mymalloc(sd, ns, 0);  b->idx_size += ns * sizeof(*sd);
    b->sd = sd;
    //printf("SB degree table: %d bytes (%1.2f bpc)\n",ns * sizeof(*sd), (double)ns * sizeof(*sd) * 8/n);
  }
//  printf("SB table: %d bytes (%1.2f bpc)\n",ns * sizeof(*sm) * 2, (double)ns * sizeof(*sm)*2 * 8/n);


  nm = (n+MB-1)/MB;
#if 1
  m_ofs = 1 << blog(nm-1);
  b->m_ofs = m_ofs;
#else
  i = 1;  m_ofs = -1;
  while (i < nm) {
    i *= W1;
    m_ofs = (m_ofs+1)*W1;
  }
  b->m_ofs = m_ofs;
#endif

  mymalloc(mm, nm + m_ofs, 0);  b->idx_size += (nm+m_ofs) * sizeof(*mm);
  mymalloc(mM, nm + m_ofs, 0);  b->idx_size += (nm+m_ofs) * sizeof(*mM);
  b->mm = mm;
  b->mM = mM;
  if (opt & OPT_DEGREE) {
    mymalloc(md, nm + m_ofs, 0);  b->idx_size += (nm+m_ofs) * sizeof(*md);
    b->md = md;
    //printf("MB degree table: %d bytes (%1.2f bpc)\n",(nm+m_ofs) * sizeof(*md), (double)(nm+m_ofs) * sizeof(*md) * 8/n);
  }
//  printf("MB table: %d bytes (%1.2f bpc)\n",(nm+m_ofs) * sizeof(*mm) * 2, (double)(nm+m_ofs) * sizeof(*mm)*2 * 8/n);


  d = 0;
  for (i=0; i<n; i++) {
    pb w;
#if 0
    if (i % SB == 0) {
      ds = depth(b,i);
      m = M = ds;
      r = 1;
    } else {
      d = depth(b,i);
      if (d == m) r++;
      if (d < m) {
        m = d;
        r = 1;
      }
      if (d > M) M = d;
    }
#else
    if (i % D == 0) {
      w = B[i/D];
    }
    d += (int)(((w>>(D-1))&1))*2-1; // d = depth(b, i);
    w <<= 1;
//    if (d != depth(b, i)) {
//      printf("i=%ld d %ld depth %ld\n", i, d, depth(b, i));
//    }
    if (i % SB == 0) {
      ds = d;
      m = M = ds;
      r = 1;
    } else {
      if (d == m) r++;
      if (d < m) {
        m = d;
        r = 1;
      }
      if (d > M) M = d;
    }
#endif
    if (i % SB == SB-1 || i==n-1) {
      ds = depth(b,(i/SB)*SB-1);
      if (m - ds + SB < 0 || m - ds + SB >= (1<<(8*sizeof(sb_type)))) {
        printf("error m=%d ds=%d\n",m,ds);
      }
      if (M - ds + 1 < 0 || M - ds + 1 >= (1<<(8*sizeof(sb_type)))) {
        printf("error M=%d ds=%d\n",M,ds);
      }
      sm[i/SB] = (sb_type)(m - ds + SB);
      sM[i/SB] = (sb_type)(M - ds + 1);
      if (opt & OPT_DEGREE) sd[i/SB] = (sb_type)r;
    }


    if (i % MB == 0) {
      m2 = M2 = d;
      r2 = 1;
    } else {
      if (d == m2) r2++;
      if (d < m2) {
        m2 = d;
        r2 = 1;
      }
      if (d > M2) M2 = d;
    }
    if (i % MB == MB-1 || i==n-1) {
      mm[m_ofs+ i/MB] = m2;
      mM[m_ofs+ i/MB] = M2;
      if (opt & OPT_DEGREE) md[m_ofs+ i/MB] = r2;
    }

  }

#if 0
  printf("sd: ");
  for (i=0;i<n/SB;i++) printf("%d ",sd[i]);
  printf("\n");
#endif


#if 0
  for (i=0; i<n; i++) {
    d = depth(b,i);
    if (i % MB == 0) {
      m = M = d;
      r = 1;
    } else {
      if (d == m) r++;
      if (d < m) {
        m = d;
        r = 1;
      }
      if (d > M) M = d;
    }
    if (i % MB == MB-1 || i==n-1) {
      mm[m_ofs+ i/MB] = m;
      mM[m_ofs+ i/MB] = M;
      if (opt & OPT_DEGREE) md[m_ofs+ i/MB] = r;
    }
  }
#endif


#if 0
  for (j=m_ofs-1; j >= 0; j--) {
    i64 k;
    m = n+1;
    for (k=0; k<W1; k++) {
      if (CHILD(j,k) < nm + m_ofs) m = min(m,mm[CHILD(j,k)]);
    }
    M = -n-1;
    for (k=0; k<W1; k++) {
      if (CHILD(j,k) < nm + m_ofs) M = max(M,mM[CHILD(j,k)]);
    }
    mm[j] = m;  mM[j] = M;
    if (opt & OPT_DEGREE) {
      i64 min_m;
      d = 0;  min_m = n+1;
      for (k=0; k<W1; k++) {
        if (CHILD(j,k) < nm + m_ofs) {
          if (min_m == mm[CHILD(j,k)]) d += md[CHILD(j,k)];
          if (min_m > mm[CHILD(j,k)]) d = md[CHILD(j,k)];
        }
      md[j] = d;
      }
    }
  }
#else
  for (j=m_ofs-1; j > 0; j--) {
    m = 0;
    if (j*2 < nm + m_ofs) m = mm[j*2];
    if (j*2+1 < nm + m_ofs) m = min(m,mm[j*2+1]);
    M = 0;
    if (j*2 < nm + m_ofs) M = mM[j*2];
    if (j*2+1 < nm + m_ofs) M = max(M,mM[j*2+1]);
    mm[j] = m;  mM[j] = M;
    if (opt & OPT_DEGREE) {
      d = 0;
      if (j*2 < nm + m_ofs) d = md[j*2];
      if (j*2+1 < nm + m_ofs) {
        if (mm[j*2] == mm[j*2+1]) d += md[j*2+1];
        if (mm[j*2] > mm[j*2+1]) d = md[j*2+1];
      }
      md[j] = d;
    }
  }
  mm[0] = -1;
  mM[0] = mM[1];
  if (opt & OPT_DEGREE) {
    md[0] = -1;
  }
#endif
  
#if 0
  mm[0] = -1;
  mM[0] = mM[1];
  if (opt & OPT_DEGREE) {
    md[0] = -1;
  }
#endif

#if 0
  printf("md: ");
  for (i=0;i<m_ofs + n/MB;i++) printf("%d ",md[i]);
  printf("\n");
#endif

  if (opt & OPT_LEAF) {
    mymalloc(b->da_leaf,1,0);
    darray_pat_construct(b->da_leaf, n, B, 2, 0x2, opt & OPT_FAST_LEAF_SELECT);
    printf("leaf rank/select table: %d bytes (%1.2f bpc)\n", 
           b->da_leaf->idx_size, (double)b->da_leaf->idx_size*8/n);
    b->idx_size += b->da_leaf->idx_size;  
  } else {
    b->da_leaf = NULL;
  }

  if (opt & OPT_INORDER) {
    mymalloc(b->da_inorder,1,0);
    darray_pat_construct(b->da_inorder, n, B, 2, 0x1, opt & OPT_FAST_INORDER_SELECT);
    printf("inorder rank/select table: %d bytes (%1.2f bpc)\n", 
           b->da_inorder->idx_size,(double)b->da_inorder->idx_size*8/n);
    b->idx_size += b->da_inorder->idx_size;
  } else {
    b->da_inorder = NULL;
  }

  if (opt & OPT_FAST_POSTORDER_SELECT) {
    mymalloc(b->da_postorder,1,0);
    darray_pat_construct(b->da_postorder, n, B, 1, 0x0, (opt & OPT_FAST_POSTORDER_SELECT) | OPT_NO_RANK);
    printf("postorder rank/select table: %d bytes (%1.2f bpc)\n",
           b->da_postorder->idx_size,(double)b->da_postorder->idx_size*8/n);
    b->idx_size += b->da_postorder->idx_size;
  } else {
    b->da_postorder = NULL;
  }

  if (opt & OPT_DFUDS_LEAF) {
    mymalloc(b->da_dfuds_leaf,1,0);
    darray_pat_construct(b->da_dfuds_leaf, n, B, 2, 0x0, opt & OPT_FAST_DFUDS_LEAF_SELECT);
    printf("dfuds leaf rank/select table: %d bytes (%1.2f bpc)\n",
           b->da_dfuds_leaf->idx_size,(double)b->da_dfuds_leaf->idx_size*8/n);
    b->idx_size += b->da_dfuds_leaf->idx_size;
  } else {
    b->da_dfuds_leaf = NULL;
  }

  if (opt & OPT_FAST_LCA) {
    i64 h,s,t;
    mb_type *rmq_tbl;

    b->opt &= ~OPT_FAST_LCA;
    h = blog(nm-1);
    mymalloc(rmq_tbl,nm*h,0);
    b->rmq_tbl = rmq_tbl;
    for (i=0; i<nm; i++) {
      for (j=0; j<h; j++) {
        s = i*MB;
        t = min((i+(1<<j))*MB-1,n-1);
        rmq_tbl[j*nm + i] = rmq(b,s,t,0);
      }
    }
//    printf("fast_lca table: %d bytes (%1.2f bpc)\n",
//           nm*h*sizeof(*rmq_tbl),(double)nm*h*sizeof(*rmq_tbl)*8/n);
    b->opt |= OPT_FAST_LCA;
  }

  if (opt & OPT_COMP_RLE) { // 圧縮する
    darray_comp_rle(b->da);
    b->idx_size += b->da->idx_size;
//    printf("preorder rank/select table + compressed BP: %d bytes (%1.2f bpc)\n",
//           b->da->idx_size,(double)b->da->idx_size*8/n);
  } else {
    i64 bp_size;
    b->idx_size += b->da->idx_size;
    bp_size = (n*2+D-1)/D * sizeof(pb);
    b->idx_size += bp_size;
//    printf("preorder rank/select table + uncompressed BP: %d bytes (%1.2f bpc)\n",
//           b->da->idx_size + bp_size,(double)(b->da->idx_size + bp_size)*8/n);
  }

  b->opt = opt;

  return 0;
}

i64 naive_fwd_excess(bp *b,i64 s, i64 rel)
{
  i64 i,v,n;
  pb *B;
  n = b->n;
//////////////////////////////// *B は圧縮されていないデータへのポインタ
  B = b->B;
  v = 0;
  for (i=s+1; i<n; i++) {
    if (getbit(B,i)==OP) {
      v++;
    } else {
      v--;
    }
    if (v == rel) return i;
  }
  return -1;
}

i64 naive_bwd_excess(bp *b,i64 s, i64 rel)
{
  i64 i,v;
  pb *B;
//////////////////////////////// *B は圧縮されていないデータへのポインタ
  B = b->B;
  v = 0;
  for (i=s; i>=0; i--) {
    if (getbit(B,i)==OP) {
      v--;
    } else {
      v++;
    }
    if (v == rel) return i-1;
  }
  return -2;
}

i64 naive_search_SB_l(bp *b, i64 i, i64 rel)
{
  i64 il;

  il = (i / SB) * SB;
  for (; i>=il; i--) {
    if (getbit(b->B,i)==OP) {
      rel++;
    } else {
      rel--;
    }
    if (rel == 0) return i-1;
  }
  if (i < 0) return -2;
  return -3;
}

i64 naive_rmq(bp *b, i64 s, i64 t,i64 opt)
{
  i64 d,i,dm,im;

  if (opt & OPT_RIGHT) {
    d = dm = depth(b,t);  im = t;
    i = t-1;
    while (i >= s) {
      if (getbit(b->B,i+1)==CP) {
        d++;
        if (opt & OPT_MAX) {
          if (d > dm) {
            dm = d;  im = i;
          }
        }
      } else {
        d--;
        if (!(opt & OPT_MAX)) {
          if (d < dm) {
            dm = d;  im = i;
          }
        }
      }
      i--;
    }
  } else {
    d = dm = depth(b,s);  im = s;
    i = s+1;
    while (i <= t) {
      if (getbit(b->B,i)==OP) {
        d++;
        if (opt & OPT_MAX) {
          if (d > dm) {
            dm = d;  im = i;
          }
        }
      } else {
        d--;
        if (!(opt & OPT_MAX)) {
          if (d < dm) {
            dm = d;  im = i;
          }
        }
      }
      i++;
    }
  }
  return im;
}

i64 root_node(bp *b)
{
  return 0;
}


i64 rank_open(bp *b, i64 s)
{
  return darray_rank(b->da,s);
}

i64 rank_close(bp *b, i64 s)
{
  return s+1 - darray_rank(b->da,s);
}

i64 select_open(bp *b, i64 s)
{
  if (b->opt & OPT_FAST_PREORDER_SELECT) {
    return darray_select(b->da,s,1);
  } else {
    return darray_select_bsearch(b->da,s,getpat_preorder);
  }
}

i64 select_close(bp *b, i64 s)
{
  if (b->opt & OPT_FAST_POSTORDER_SELECT) {
    return darray_pat_select(b->da_postorder,s,getpat_postorder);
  } else {
    return postorder_select_bsearch(b,s);
  }
}

///////////////////////////////////////////
//  find_close(bp *b,i64 s)
//    returns the matching close parenthesis of s
///////////////////////////////////////////
i64 find_close(bp *b,i64 s)
{
  return fwd_excess(b,s,-1);
}

///////////////////////////////////////////
//  find_open(bp *b,i64 s)
//    returns the matching open parenthesis of s
///////////////////////////////////////////
i64 find_open(bp *b,i64 s)
{
  i64 r;
  r = bwd_excess(b,s,0);
  if (r >= -1) return r+1;
  return -1;
}

///////////////////////////////////////////
//  parent(bp *b,i64 s)
//    returns the parent of s
//            -1 if s is the root
///////////////////////////////////////////
i64 parent(bp *b,i64 s)
{
  i64 r;
  r = bwd_excess(b,s,-2);
  if (r >= -1) return r+1;
  return -1;
}

i64 enclose(bp *b,i64 s)
{
  return parent(b,s);
}

///////////////////////////////////////////
//  level_ancestor(bp *b,i64 s,i64 d)
//    returns the ancestor of s with relative depth d (d < 0)
//            -1 if no such node
///////////////////////////////////////////
i64 level_ancestor(bp *b,i64 s,i64 d)
{
  i64 r;
  r = bwd_excess(b,s,d-1);
  if (r >= -1) return r+1;
  return -1;
}


///////////////////////////////////////////
//  rmq(bp *b, i64 s, i64 t, i64 opt)
//    returns the index of leftmost/rightmost minimum/maximum value
//                 in the range [s,t] (inclusive)
//    returns the leftmost minimum if opt == 0
//            the rightmost one if (opt & OPT_RIGHT)
//            the maximum if (opt & OPT_MAX)
///////////////////////////////////////////
#if 0
i64 rmq(bp *b, i64 s, i64 t, i64 opt)
{
  if ((b->opt & OPT_FAST_LCA) && (opt == 0)) {
    return fast_rmq(b,s,t,0);
  }
  return std_rmq(b,s,t,opt);
}
#endif


///////////////////////////////////////////
//  lca(bp *b, i64 s, i64 t)
//    returns the lowest common ancestor of s and t
///////////////////////////////////////////
i64 lca(bp *b, i64 s, i64 t)
{
  return parent(b,rmq(b,s,t,0)+1);
}


///////////////////////////////////////////
//  preorder_rank(bp *b,i64 s)
//    returns the preorder (>= 1) of node s (s >= 0)
///////////////////////////////////////////
i64 preorder_rank(bp *b,i64 s)
{
  return darray_rank(b->da,s);
}

///////////////////////////////////////////
//  preorder_select(bp *b,i64 s)
//    returns the node with preorder s (s >= 1)
//            -1 if no such node
///////////////////////////////////////////
i64 preorder_select(bp *b,i64 s)
{
  //  no error handling
  if (b->opt & OPT_FAST_PREORDER_SELECT) {
    return darray_select(b->da,s,1);
  } else {
    return darray_select_bsearch(b->da,s,getpat_preorder);
  }
}

///////////////////////////////////////////
//  postorder_rank(bp *b,i64 s)
//    returns the postorder (>= 1) of node s (s >= 0)
//            -1 if s-th bit is not OP
///////////////////////////////////////////
i64 postorder_rank(bp *b,i64 s)
{
  i64 t;
  if (inspect(b,s) == CP) return -1;
  t = find_close(b,s);
  //  return t+1 - darray_rank(b->da,t);
  return rank_close(b,t);
}

i64 postorder_select_bsearch(bp *b,i64 s)
{
  i64 l,r,m;

  if (s == 0) return -1;

  if (s > b->da->n - b->da->m) {
    return -1;
  }
  l = 0;  r = b->da->n - 1;

  while (l < r) {
    m = (l+r)/2;
    //printf("m=%d rank=%d s=%d\n",m,m+1 - darray_rank(b->da,m),s);
    if (m+1 - darray_rank(b->da,m) >= s) {
      r = m;
    } else {
      l = m+1;
    }
  }
  return l;
}

///////////////////////////////////////////
//  postorder_select(bp *b,i64 s)
//    returns the position of CP of the node with postorder s (>= 1)
///////////////////////////////////////////
i64 postorder_select(bp *b,i64 s)
{
#if 0
  if (b->opt & OPT_FAST_POSTORDER_SELECT) {
    return darray_pat_select(b->da_postorder,s,getpat_postorder);
  } else {
    return postorder_select_bsearch(b->da,s);
  }
#else
  return select_close(b,s);
#endif
}

///////////////////////////////////////////
//  leaf_rank(bp *b,i64 s)
//    returns the number of leaves to the left of s
///////////////////////////////////////////
i64 leaf_rank(bp *b,i64 s)
{
  if ((b->opt & OPT_LEAF) == 0) {
    printf("leaf_rank: error!!!   not supported\n");
    return -1;
  }
  if (s >= b->n-1) {
    s = b->n-2;
  }
  return darray_pat_rank(b->da_leaf,s,getpat_leaf);
}

///////////////////////////////////////////
//  leaf_select(bp *b,i64 s)
//    returns the position of s-th leaf
///////////////////////////////////////////
i64 leaf_select(bp *b,i64 s)
{
  if ((b->opt & OPT_LEAF) == 0) {
    printf("leaf_select: error!!!   not supported\n");
    return -1;
  }
  if (s > b->da_leaf->m) return -1;
  if (b->opt & OPT_FAST_LEAF_SELECT) {
    return darray_pat_select(b->da_leaf,s,getpat_leaf);
  } else {
    return darray_select_bsearch(b->da_leaf,s,getpat_leaf);
  }
}


///////////////////////////////////////////
//  inorder_rank(bp *b,i64 s)
//    returns the number of ")("  (s >= 0)
///////////////////////////////////////////
i64 inorder_rank(bp *b,i64 s)
{
  if ((b->opt & OPT_INORDER) == 0) {
    printf("inorder_rank: error!!!   not supported\n");
    return -1;
  }
  if (s >= b->n-1) {
    s = b->n-2;
  }
  return darray_pat_rank(b->da_inorder,s,getpat_inorder);
}

///////////////////////////////////////////
//  inorder_select(bp *b,i64 s)
//    returns the s-th position of ")("  (s >= 1)
///////////////////////////////////////////
i64 inorder_select(bp *b,i64 s)
{
  if ((b->opt & OPT_INORDER) == 0) {
    printf("inorder_select: error!!!   not supported\n");
    return -1;
  }
  if (b->opt & OPT_FAST_INORDER_SELECT) {
    return darray_pat_select(b->da_inorder,s,getpat_inorder);
  } else {
    return darray_select_bsearch(b->da_inorder,s,getpat_inorder);
  }
}

///////////////////////////////////////////
//  leftmost_leaf(bp *b, i64 s)
///////////////////////////////////////////
i64 leftmost_leaf(bp *b, i64 s)
{
  if ((b->opt & OPT_LEAF) == 0) {
    printf("leftmost_leaf: error!!!   not supported\n");
    return -1;
  }
  return leaf_select(b,leaf_rank(b,s)+1);
}

///////////////////////////////////////////
//  rightmost_leaf(bp *b, i64 s)
///////////////////////////////////////////
i64 rightmost_leaf(bp *b, i64 s)
{
  i64 t;
  if ((b->opt & OPT_LEAF) == 0) {
    printf("leftmost_leaf: error!!!   not supported\n");
    return -1;
  }
  t = find_close(b,s);
  return leaf_select(b,leaf_rank(b,t));
}



///////////////////////////////////////////
//  inspect(bp *b, i64 s)
//    returns OP (==1) or CP (==0) at s-th bit (0 <= s < n)
///////////////////////////////////////////
i64 inspect(bp *b, i64 s)
{
  if (s < 0 || s >= b->n) {
    printf("inspect: error s=%d is out of [%d,%d]\n",s,0,b->n-1);
  }
  if (b->opt & OPT_COMP_RLE) {
    pb tmp[1];
    darray_decode_block_rle(b->da, s, s, tmp);
    return getbit(tmp,0);
  } else {
    return getbit(b->B,s);
  }
}

i64 isleaf(bp *b, i64 s)
{
  if (inspect(b,s) != OP) {
    printf("isleaf: error!!! B[%d] = OP\n",s);
  }
  if (inspect(b,s+1) == CP) return 1;
  else return 0;
}


///////////////////////////////////////////
//  subtree_size(bp *b, i64 s)
//    returns the number of nodes in the subtree of s
///////////////////////////////////////////
i64 subtree_size(bp *b, i64 s)
{
  return (find_close(b,s) - s + 1) / 2;
}

///////////////////////////////////////////
//  first_child(bp *b, i64 s)
//    returns the first child
//            -1 if s is a leaf
///////////////////////////////////////////
i64 first_child(bp *b, i64 s)
{
  if (inspect(b,s+1) == CP) return -1;
  return s+1;
}

i64 last_child(bp *b, i64 s)
{
    if (inspect(b, s+1) == CP) return -1;
    return find_open(b, find_close(b, s)-1);
}

///////////////////////////////////////////
//  next_sibling(bp *b,i64 s)
//    returns the next sibling of parent(s)
//            -1 if s is the last child
//////////////////////////////////////////
i64 next_sibling(bp *b, i64 s)
{
  i64 t;
  t = find_close(b,s)+1;
  if (t >= b->n) {
    printf("next_sibling: error s=%ld t=%ld\n",s,t);
  }
  if (inspect(b,t) == CP) return -1;
  return t;
}

///////////////////////////////////////////
//  prev_sibling(bp *b,i64 s)
//    returns the previous sibling of parent(s)
//            -1 if s is the first child
//////////////////////////////////////////
i64 prev_sibling(bp *b, i64 s)
{
  i64 t;
  if (s < 0) {
    printf("prev_sibling: error s=%ld\n",s);
  }
  if (s == 0) return -1;
  if (inspect(b,s-1) == OP) return -1;
  t = find_open(b,s-1);
  return t;
}

///////////////////////////////////////////
//  deepest_node(bp *b,i64 s)
//    returns the first node with the largest depth in the subtree of s
///////////////////////////////////////////
i64 deepest_node(bp *b,i64 s)
{
  i64 t,m;
  t = find_close(b,s);
  m = rmq(b,s,t, OPT_MAX);
  return m;
}

///////////////////////////////////////////
//  subtree_height(bp *b,i64 s)
//    returns the height of the subtree of s
//            0 if s is a leaf
///////////////////////////////////////////
i64 subtree_height(bp *b,i64 s)
{
  i64 t;
  t = deepest_node(b,s);
  return depth(b,t) - depth(b,s);
}

i64 naive_degree(bp *b, i64 s)
{
  i64 t,d;
  d = 0;
  t = first_child(b,s);
  while (t >= 0) {
    d++;
    t = next_sibling(b,t);
  }
  return d;
}

///////////////////////////////////////////
//  degree(bp *b, i64 s)
//    returns the number of children of s
//            0 if s is a leaf
///////////////////////////////////////////
i64 degree(bp *b, i64 s)
{
  if (b->opt & OPT_DEGREE) {
    return fast_degree(b,s,b->n,0);
  } else {
    return naive_degree(b,s);
  }
}

i64 naive_child(bp *b, i64 s, i64 d)
{
  i64 t,i;
  t = first_child(b,s);
  for (i = 1; i < d; i++) {
    if (t == -1) break;
    t = next_sibling(b,t);
  }
  return t;
}

///////////////////////////////////////////
//  child(bp *b, i64 s, i64 d)
//    returns the d-th child of s  (1 <= d <= degree(s))
//            -1 if no such node
///////////////////////////////////////////
i64 child(bp *b, i64 s, i64 d)
{
  i64 r;
  if (b->opt & OPT_DEGREE) {
    //return find_open(b,fast_degree(b,s,b->n,d));
    if (d==1) return first_child(b,s);
    r = fast_degree(b,s,b->n,d-1)+1;
    if (inspect(b,r) == CP) return -1;
    return r;
  } else {
    return naive_child(b,s,d);
  }
    
}


i64 naive_child_rank(bp *b, i64 t)
{
  i64 d;
  d = 0;
  while (t != -1) {
    d++;
    t = prev_sibling(b,t);
  }
  return d;
}

///////////////////////////////////////////
//  child_rank(bp *b, i64 t)
//    returns d if t is the d-th child of the parent of t (d >= 1)
//            1 if t is the root
///////////////////////////////////////////
i64 child_rank(bp *b, i64 t)
{
  i64 r;
  if (t == root_node(b)) return 1;
  if (b->opt & OPT_DEGREE) {
    r = parent(b,t);
    return fast_degree(b,r,t,0)+1;
  } else {
    return naive_child_rank(b,t);
  }
}



///////////////////////////////////////////
//  is_ancestor(bp *b, i64 s, i64 t)
//     returns 1 if s is an ancestor of t
//             0 otherwise
///////////////////////////////////////////
i64 is_ancestor(bp *b, i64 s, i64 t)
{
  i64 v;
  v = find_close(b,s);
  if (s <= t && t <= v) return 1;
  return 0;
}

///////////////////////////////////////////
//  distance(bp *b, i64 s, i64 t)
//    returns the length of the shortest path from s to t in the tree
///////////////////////////////////////////
i64 distance(bp *b, i64 s, i64 t)
{
  i64 v,d;
  v = lca(b,s,t);
  d = depth(b,v);
  return (depth(b,s) - d) + (depth(b,t) - d);
}

///////////////////////////////////////////
//  level_next(bp *b, i64 d)
///////////////////////////////////////////
i64 level_next(bp *b,i64 s)
{
  i64 t;
  t = fwd_excess(b,s,0);
  return t;
}

///////////////////////////////////////////
//  level_prev(bp *b, i64 d)
///////////////////////////////////////////
i64 level_prev(bp *b,i64 s)
{
  i64 t;
  t = find_open(b,bwd_excess(b,s,0)+1);
  return t;
}

///////////////////////////////////////////
//  level_leftmost(bp *b, i64 d)
///////////////////////////////////////////
i64 level_leftmost(bp *b, i64 d)
{
  i64 t;
  if (d < 1) return -1;
  if (d == 1) return 0;
  t = fwd_excess(b,0,d);
  return t;
}

///////////////////////////////////////////
//  level_rigthmost(bp *b, i64 d)
///////////////////////////////////////////
i64 level_rigthmost(bp *b, i64 d)
{
  i64 t;
  if (d < 1) return -1;
  if (d == 1) return 0;
  t = bwd_excess(b,b->n-1,d)+1;
  return find_open(b,t);
}

///////////////////////////////////////////
//  leaf_size(bp *b, i64 s)
///////////////////////////////////////////
i64 leaf_size(bp *b, i64 s)
{
  i64 t;
  if ((b->opt & OPT_LEAF) == 0) {
    printf("leaf_size: error!!!   not supported\n");
    return -1;
  }
  t = find_close(b,s);
  return leaf_rank(b,t) - leaf_rank(b,s);
}


#if 0
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
// 新しい実装

typedef struct
{
  i64 sum;
  i64 min;
  i64 max;
} smm;

// vec[s..t] に対する MB 構造を作成
// t <= s+MB-1;
void construct_mblock(mblock *mb, pb *vec, i64 s, i64 t, smm *ret)
{
  i64 sb_sum,sb_min,sb_max;
  i64 mb_sum,mb_min,mb_max;
  i64 i,j,k;
  int b;

  if (s % MB != 0) {
    printf("construct_mblock: s = %ld is not a multiple of MB = %ld\n",s,MB);
//    exit(1);
  }

  mb_sum = 0;
  mb_min = MB+1;  mb_max = -MB-1;
  for (i=s; i<s+MB; i+=SB) {
    sb_sum = 0;
    sb_min = SB+1;  sb_max = -SB-1;
    for (j=0; j<SB; j++) {
//////////////////////////////// *vec は圧縮されていないデータへのポインタ
      b = (i+j <= t) ? getbit(vec,i+j) : CP;  // 括弧列の最後にダミーの )))... を書く
      if (b == OP) sb_sum++; else sb_sum--;
      sb_max = max(sb_sum,sb_max);
      sb_min = min(sb_sum,sb_min);
    }
    mb->pi_sum[(i-s)/SB] = sb_sum;
    mb->pi_min[(i-s)/SB] = sb_min;
    mb->pi_max[(i-s)/SB] = sb_max;
    mb_min = min(mb_min,mb_sum+sb_min);
    mb_max = min(mb_max,mb_sum+sb_max);
    mb_sum += sb_sum;
  }
  ret->sum = mb_sum;
  ret->min = mb_min;
  ret->max = mb_max;
}

// mb[0..m-1] に対する LB 構造を作成
void construct_lblock_leaf(lblock *lb, mblock *mb, i64 m)
{
  i64 i;
  i64 lb_sum,lb_min,lb_max;
  
  lb_sum = 0;
  lb_min = MB+1;  lb_max = -MB-1;
  for (i=0; i<m; i++) {
    lb_min = min(lb_min,lb_sum+mb[i]->pi_min);
    lb_max = max(lb_max,lb_sum+mb[i]->pi_max);
    lb_sum += mb[i]->pi_sum;
    lb->pi_min[i] = lb_min;
    lb->pi_max[i] = lb_max;
    lb->pi_sum[i] = lb_sum;
  }
}

bp2 *construct_bp2(i64 n, pb *vec)
{
  i64 i,j,l;
  bp2 *b;
  i64 ns; // SBの数
  i64 nm; // MBの数
  i64 depth; //  rmm-treeの深さ
  i64 m, m_ofs;
  smm ret,lb;

  b = (bp2 *)mymalloc(sizeof(bp2));
  ns = (n+SB-1)/SB;
  nm = (n+MB-1)/MB;

  depth = 0;
  m = 1;  m_ofs = -1;
  while (m < nm) {
    m *= BF;
    m_ofs = (m_ofs+1)*BF;
    depth++;
  }
  b->m_ofs = m_ofs;

  b->lb = (lblock *)mymalloc(m_ofs * sizeof(lblock));
  b->mb = (mblock *)mymalloc(nm * sizeof(mblock));

  for (l=0; l<n; l+=MB*BF) {
    lb.sum = 0;  lb.min = n+1;  lb.max = -n-1;
    for (i=l; i<n; i+=MB) {
      j = min(n-1,MBlast(i));
      construct_mblock(&b->mb[MBid(i)],vec,i,j,&ret);
      j = MM_P(MBid(i)+m_ofs);
      
    }
    
    for (i=0; i<ns; i+=BF) {
    j = min(ns-1,i+BF-1);
    construct_mblock(&b->mb[m_ofs + 00000000],vec,i,j);
  }
}
#endif

void bp_free(bp *b) {
//    if (b->B) free(b->B);
    if (b->da) {
      darray_free(b->da);
      free(b->da); // sada 2015/1/16
    }
    if (b->sm) free(b->sm);
    if (b->sM) free(b->sM);
//    if (b->sd) free(b->sd);
    if (b->mm) free(b->mm);
    if (b->mM) free(b->mM);
//    if (b->md) free(b->md);
    if (b->rmq_tbl) free(b->rmq_tbl);
    if (b->da_leaf) darray_free(b->da_leaf);
    if (b->da_inorder) darray_free(b->da_inorder);
    if (b->da_postorder) darray_free(b->da_postorder);
    if (b->da_dfuds_leaf) darray_free(b->da_dfuds_leaf);
}
