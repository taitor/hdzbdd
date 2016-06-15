#ifndef _BP_H_
#define _BP_H_
#include "typedefbp.h"
#include "darray.h"

#define OP 1
#define CP 0

#define OPT_MIN 0
#define OPT_MAX 1
#define OPT_LEFT 0
#define OPT_RIGHT 2

#define OPT_LEAF (1<<0)
#define OPT_INORDER (1<<1)
#define OPT_DEGREE (1<<2)
#define OPT_FAST_PREORDER_SELECT (1<<3)
#define OPT_FAST_LEAF_SELECT (1<<4)
#define OPT_FAST_INORDER_SELECT (1<<5)
#define OPT_FAST_POSTORDER_SELECT (1<<6)
#define OPT_DFUDS_LEAF (1<<7)
#define OPT_FAST_DFUDS_LEAF_SELECT (1<<8)
#define OPT_FAST_LCA (1<<9)
#define OPT_COMP_RLE (1<<10)
#define OPT_MIN_POS (1<<11)

//#define logSB 9
#define logSB 8
//#define logSB 2
#define SB (1<<logSB)
//#define logMB 15
//#define logMB 8
#define logMB 10
//#define logMB 12
//#define logMB 3
#define MB (1<<logMB)
#define BF (MB/SB) // branching factor

#define ETW 8  // width of excess lookup table
#define W1 2    // branching factor

#ifndef min
 #define min(x,y) ((x)<(y)?(x):(y))
#endif
#ifndef max
 #define max(x,y) ((x)>(y)?(x):(y))
#endif

#define CHILD(i,k) (((i)+1)*W1+(k))
#define PARENT(i) ((i)/W1-1)

#if 0
typedef struct {
//  sb_type len; // SB中のビット列の長さ (static な時は不要 == SB)
//  pb *vec; // ビット列へのポインタ (static な時は不要)
  sb_type pi_sum; // SBでの pi の和
  sb_type pi_min; // SBでの pi の最小値 (SBの直前からの相対値)
  sb_type pi_max; // SBでの pi の最大値 (SBの直前からの相対値)
} sblock;
#endif

typedef struct {
//  sb_type len[BF];  // SB中のビット列の長さ (static な時は不要 == SB)
//  pb *vec[BF];      // ビット列へのポインタ (static な時は不要)
  sb_type pi_sum[BF]; // 各SBの最後での pi の和 (MBの直前からの相対値)
  sb_type pi_min[BF]; // 各SBでの pi の最小値 (MBの直前らの相対値)
  sb_type pi_max[BF]; // 各SBでの pi の最大値 (MBの直前からの相対値)
} mblock;

typedef struct {
//  mb_type len[BF];  // MB中のビット列の長さ (static な時は不要 == MB)
// lblock *son[BF];   // 子のlblock (またはmblock) へのポインタ (static な時は不要)
// pb lm_flag[(BF+D-1)/D]; // 子がlblockかmblockかを表すフラグ
  mb_type pi_sum[BF]; // 各MBの最後での pi の和 (LBの直前からの相対値)
  mb_type pi_min[BF]; // 各MBでの pi の最小値 (LBの直前らの相対値)
  mb_type pi_max[BF]; // 各MBでの pi の最大値 (LBの直前らの相対値)
} lblock;

#define MM_P(x) ((x)/BF-1)
#define MM_SON(x) (((x)+1)*BF)

typedef struct {
  i64 n;

  mblock *mb;
  lblock *lb;
  i64 m_ofs;

  i64 idx_size;
  i64 opt;
} bp2;


typedef struct {
  i64 n;
  pb *B;
  darray *da;
  sb_type *sm, *sM;
  sb_type *sd;
  sb_type *smi;
  mb_type *mm, *mM;
  mb_type *md;
//  mb_type *mmi;

  mb_type *rmq_tbl;

  i64 m_ofs;

  darray *da_leaf;
  darray *da_inorder;
  darray *da_postorder;
  darray *da_dfuds_leaf;

  

  i64 idx_size;
  int opt;
} bp;

#ifdef __cplusplus
extern "C" {
#endif
i64 bp_construct(bp *b,i64 n, pb *B, i64 opt);
void printbp(bp *b, i64 s, i64 t);

typedef struct {
  void *env;
  void (*rewind)(void *env);
  int (*next)(void *env);
} bpsub;

i64 bp_construct2(bp *b,i64 n, pb *B, i64 opt, bpsub *e);


i64 rank_open(bp *b, i64 s);
i64 rank_close(bp *b, i64 s);
i64 select_open(bp *b, i64 s);
i64 select_close(bp *b, i64 s);


i64 root_node(bp *b);
i64 find_close(bp *b,i64 s);
i64 find_open(bp *b,i64 s);
i64 enclose(bp *b,i64 s);
i64 parent(bp *b,i64 s);
i64 level_ancestor(bp *b,i64 s,i64 d);
i64 depth(bp *b, i64 s);
i64 fast_depth(bp *b, i64 s);
i64 preorder_rank(bp *b,i64 s);
i64 postorder_rank(bp *b,i64 s);
i64 inspect(bp *b, i64 s);
i64 isleaf(bp *b, i64 s);
i64 rmq(bp *b, i64 s, i64 t, i64 opt);
i64 std_rmq(bp *b, i64 s, i64 t, i64 opt);
i64 fast_rmq(bp *b, i64 s, i64 t, int opt);
i64 subtree_size(bp *b, i64 s);
i64 first_child(bp *b, i64 s);
i64 last_child(bp *b, i64 s);
i64 next_sibling(bp *b, i64 s);
i64 prev_sibling(bp *b, i64 s);
i64 deepest_node(bp *b,i64 s);
i64 subtree_height(bp *b,i64 s);
i64 is_ancestor(bp *b, i64 s, i64 t);
i64 distance(bp *b, i64 s, i64 t);
i64 level_lefthmost(bp *b, i64 d);
i64 level_rigthmost(bp *b, i64 d);
i64 degree(bp *b,i64 s);

// not efficient
i64 naive_degree(bp *b, i64 s);
i64 naive_child(bp *b, i64 s, i64 d);
i64 naive_child_rank(bp *b, i64 t);
i64 naive_rmq(bp *b, i64 s, i64 t,i64 opt);
i64 postorder_select(bp *b,i64 s);
i64 postorder_select_bsearch(bp *b,i64 s);
i64 naive_bwd_excess(bp *b,i64 s, i64 rel);

// using preorder select index
i64 preorder_select(bp *b,i64 s);

// using leaf index
i64 leaf_rank(bp *b,i64 s);
i64 leaf_size(bp *b, i64 s);
i64 leftmost_leaf(bp *b, i64 s);
i64 rightmost_leaf(bp *b, i64 s);

// using leaf select index
i64 leaf_select(bp *b,i64 s);

// using inorder index
i64 inorder_rank(bp *b,i64 s);

// using inorder select index
i64 inorder_select(bp *b,i64 s);

// using degree index
i64 fast_degree(bp *b,i64 s, i64 t, i64 ith);
i64 child_rank(bp *b, i64 t);
i64 child(bp *b, i64 s, i64 d);


i64 blog(i64 x);
pb getpat_preorder(pb *b);
pb getpat_leaf(pb *b);
pb getpat_inorder(pb *b);
pb getpat_postorder(pb *b);


i64 fwd_excess(bp *b,i64 s, i64 rel);
i64 bwd_excess(bp *b,i64 s, i64 rel);

extern i64 *matchtbl,*parenttbl;
void make_naivetbl(pb *B,i64 n);

void bp_make_matchtbl(void);

extern int popCount[1<<ETW];
extern int fwdtbl[(2*ETW+1)*(1<<ETW)];
extern int bwdtbl[(2*ETW+1)*(1<<ETW)];
extern int mintbl_li[1<<ETW], mintbl_lv[1<<ETW];
extern int mintbl_ri[1<<ETW], mintbl_rv[1<<ETW];
extern int maxtbl_li[1<<ETW], maxtbl_lv[1<<ETW];
extern int maxtbl_ri[1<<ETW], maxtbl_rv[1<<ETW];

extern int minmaxtbl_i[4][1<<ETW], minmaxtbl_v[4][1<<ETW];
extern int degtbl[1<<ETW];
extern int degtbl2[(2*ETW+1)*(1<<ETW)];
extern int childtbl[(ETW)*(1<<ETW)];
extern int depthtbl[(2*ETW+1)*(1<<ETW)];
extern int childtbl2[2*ETW+1][ETW][(1<<ETW)];
    
//    added by Lee
void bp_free(bp *b);
#ifdef __cplusplus
}
#endif
#endif
