/*****************************************
 *  BDD Package (SAPPORO-1.56)   - Body   *
 *  (C) Shin-ichi MINATO  (Dec. 13, 2012) *
 ******************************************/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "hdbdd.h"

//include for Hybrid DenseBDD
#include "bp.h"
#include "bitvector.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include <limits.h>
//include for debug
#include <sys/resource.h>
#include <time.h>

#include "bitio.h"

using namespace std;

unsigned long compress_count = 0;
unsigned long long int trans_time = 0;
unsigned long long int calcPrerank_time = 0;
unsigned long long int convertToDense_time = 0;
unsigned long long int getnewdense_time = 0;

double THETA = 0.1f;

#define Z_LIMIT 300000

/* ----------------- MACRO Definitions ---------------- */
/* Operation IDs in Cache */
#define BC_NULL        0
#define BC_AND         1
#define BC_XOR         2
#define BC_AT0         3
#define BC_AT1         4
#define BC_LSHIFT      5
#define BC_RSHIFT      6
#define BC_COFACTOR    7
#define BC_UNIV        8
#define BC_SUPPORT     9
#define BC_INTERSEC   10
#define BC_UNION      11
#define BC_SUBTRACT   12
#define BC_OFFSET     13
#define BC_ONSET      14
#define BC_CHANGE     15
#define BC_CARD       16
#define BC_LIT        17
#define BC_LEN        18

//added by Lee
#define BC_MULT       20
#define BC_DIV        21
#define BC_RSTR       22
#define BC_PERMIT     23
#define BC_PERMITSYM  24
#define BC_SYMCHK     25
#define BC_ALWAYS     26
#define BC_SYMSET     27
#define BC_COIMPSET   28
#define BC_MEET       29

/* Macros for malloc, realloc */
#define B_MALLOC(type, size) \
(type *)malloc(sizeof(type) * size)
#define B_REALLOC(ptr, type, size) \
(type *)realloc(ptr, sizeof(type) * size)

/* Printf format of bddp */
#ifdef B_64
#  define B_BDDP_FD "%lld"
#  define B_BDDP_FX "0x%llX"
#else
#  define B_BDDP_FD "%d"
#  define B_BDDP_FX "0x%X"
#endif

/* strtol or strtoll */
#ifdef B_64
#  define B_STRTOI strtoll
#else
#  define B_STRTOI strtol
#endif

/* Table spaces */
#define B_NODE_MAX (B_VAL_MASK>>1U) /* Max number of BDD nodes */
#define B_NODE_SPC0 256 /* Default initial node size */
#define B_VAR_SPC0   16 /* Initial var table size */
#define B_HASH_SPC0   4 /* Initial hash size */
#define B_RFCT_SPC0   4 /* Initial RFCT size */

/* Negative edge manipulation */
#define B_NEG(f)  ((f) & B_INV_MASK)
#define B_NOT(f)  ((f) ^ B_INV_MASK)
#define B_ABS(f)  ((f) & ~B_INV_MASK)

/* Constant node manipulation */
#define B_CST(f)  ((f) & B_CST_MASK)
#define B_VAL(f)  ((f) & B_VAL_MASK)

/* Conversion of bddp and node index/pointer  */
//#define B_NP(f)       (Node+(B_ABS(f)>>1U))
//node ID of Node[0] is Dense_Spc+1
#define B_NP(f)       (Node+((B_ABS(f)>>1U)-FirstIdx)) //convert bddp to pointer in array of B_NodeTable
#define B_NDX(f)      ((B_ABS(f)>>1U)-FirstIdx) //convert bddp to index in array of B_NodeTable
#define B_BDDP_NP(p)  ((bddp)((((p)-Node)+FirstIdx)<<1U)) //convert pointer of B_NodeTable to bddp

/* Read & Write of bddp field in the tables */
#ifdef B_64
#  define B_LOW32(f) ((bddp_32)((f)&((1ULL<<32U)-1U)))
#  define B_HIGH8(f) ((bddp_h8)((f)>>32U))
#  define B_SET_NXP(p, f, i) \
(p ## _h8 = f ## _h8 + i, p ## _32 = f ## _32 + i)
#  define B_GET_BDDP(f) \
((bddp) f ## _32 | ((bddp) f ## _h8 << 32U))
#  define B_SET_BDDP(f, g) \
(f ## _h8 = B_HIGH8(g), f ## _32 = B_LOW32(g))
#  define B_CPY_BDDP(f, g) \
(f ## _h8 = g ## _h8, f ## _32 = g ## _32)
#else
#  define B_SET_NXP(p, f, i) (p ## _32 = f ## _32 + i)
#  define B_GET_BDDP(f) (f ## _32)
#  define B_SET_BDDP(f, g) (f ## _32 = g)
#  define B_CPY_BDDP(f, g) (f ## _32 = g ## _32)
#endif /* B_64 */

/* var & rfc manipulation */
#define B_VAR_NP(p)    ((p)->varrfc & B_VAR_MASK)
#define B_RFC_MASK  (~B_VAR_MASK)
#define B_RFC_UNIT  (1U << B_VAR_WIDTH)
#define B_RFC_NP(p)    ((p)->varrfc >> B_VAR_WIDTH)
#define B_RFC_ZERO_NP(p) ((p)->varrfc < B_RFC_UNIT)
#define B_RFC_ONE_NP(p) (((p)->varrfc & B_RFC_MASK) == B_RFC_UNIT)
#if 1 // sada 4/14
#define B_RFC_INC_NP(p) \
(((p)->varrfc < B_RFC_MASK - B_RFC_UNIT)? \
((p)->varrfc += B_RFC_UNIT, 0) : rfc_inc_ovf(p))
#define B_RFC_DEC_NP(p) \
(((p)->varrfc >= B_RFC_MASK)? rfc_dec_ovf(p): \
(B_RFC_ZERO_NP(p))? \
err("B_RFC_DEC_NP: rfc under flow", p-Node): \
((p)->varrfc -= B_RFC_UNIT, 0))
#else
#define B_RFC_INC_NP(p)
#define B_RFC_DEC_NP(p)
#endif
/* ----------- Stack overflow limitter ------------ */
//const int BDD_RecurLimit = 8192;
const int BDD_RecurLimit = 8192*2;
int BDD_RecurCount = 0;
#if 1 // sada 4/14
#define BDD_RECUR_INC \
{if(++BDD_RecurCount >= BDD_RecurLimit) \
err("BDD_RECUR_INC: Recursion Limit", BDD_RecurCount);}
#define BDD_RECUR_DEC BDD_RecurCount--
#else
#define BDD_RECUR_INC
#define BDD_RECUR_DEC
#endif

/* Conversion of ZBDD node flag */
#define B_Z_NP(p) ((p)->f0_32 & (bddp_32)B_INV_MASK)

/* Hash Functions */
#define B_HASHKEY(f0, f1, hashSpc) \
(((B_CST(f0)? (f0): (f0)+2U) \
^(B_NEG(f0)? ~((f0)>>1U): ((f0)>>1U)) \
^(B_CST(f1)? (f1)<<1U: ((f1)+2U)<<1U) \
^(B_NEG(f1)? ~((f1)>>1U): ((f1)>>1U)) )\
& (hashSpc-1U))
/*  (((f0)^((f0)>>10)^((f0)>>31)^(f1)^((f1)>>8)^((f1)>>31)) \*/
#define B_CACHEKEY(op, f, g) \
((((bddp)(op)<<2U) \
^(B_CST(f)? (f): (f)+2U) \
^(B_NEG(f)? ~((f)>>1U): ((f)>>1U)) \
^(B_CST(g)? (g)<<3U: ((g)+2U)<<3U) \
^(B_NEG(g)? ~((g)>>1U): ((g)>>1U)) )\
& (CacheSpc-1U))

/* ------- Declaration of static (internal) data ------- */
/* typedef of bddp field in the tables */
typedef unsigned int bddp_32;
#ifdef B_64
typedef unsigned char bddp_h8;
#endif

/* Declaration of Node table */
struct B_NodeTable
{
  bddp_32      f0_32;  /* 0-edge */
  bddp_32      f1_32;  /* 1-edge */
  bddp_32      nx_32;  /* Node index */
  unsigned int varrfc; /* VarID & Reference counter */
#ifdef B_64
  bddp_h8      f0_h8;  /* Extention of 0-edge */
  bddp_h8      f1_h8;  /* Extention of 1-edge */
  bddp_h8      nx_h8;  /* Extention of node index */
#endif /* B_64 */
};
static struct B_NodeTable *Node = 0; /* Node Table */
static bddp NodeLimit = 0;    /* Final limit size */
static bddp NodeUsed = 0;     /* Number of used node */
static bddp Avail = bddnull;  /* Head of available node */
static bddp NodeSpc = 0;      /* Current Node-Table size */
static bddp FirstIdx = 1; //the index of Node[0], added by Lee

/* Declaration of Hash-table per Var */
struct B_VarTable
{
  bddp    hashSpc;  /* Current hash-table size */
  bddp    hashUsed;  /* Current used entries */
  bddvar  lev;      /* Level of the variable */
  bddp_32 *hash_32; /* Hash-table */
#ifdef B_64
  bddp_h8 *hash_h8; /* Extension of hash-table */
#endif /* B_64 */
};
static struct B_VarTable *Var = 0; /* Var-tables */
static bddvar *VarID = 0;     /* VarID reverse table */
static bddvar VarUsed = 0;    /* Number of used Var */
static bddvar VarSpc = 0;     /* Current Var-table size */

/* Declaration of Operation Cache */
struct B_CacheTable
{
  bddp_32       f_32; /* an operand BDD */
  bddp_32       g_32; /* an operand BDD */
  bddp_32       h_32; /* Result BDD */
  unsigned char op;   /* Operation code */
#ifdef B_64
  bddp_h8       f_h8; /* Extention of an operand BDD */
  bddp_h8       g_h8; /* Extention of an operand BDD */
  bddp_h8       h_h8; /* Extention of result BDD */
#endif /* B_64 */
};
static struct B_CacheTable *Cache = 0; /* Opeartion cache */
static bddp CacheSpc = 0;           /* Current cache size */

/* Declaration of RFC-table */
struct B_RFC_Table
{
  bddp_32 nx_32;   /* Node index */
  bddp_32 rfc_32;  /* RFC */
#ifdef B_64
  bddp_h8 nx_h8;   /* Extension of Node index */
  bddp_h8 rfc_h8;  /* Extension of RFC */
#endif /* B_64 */
};
static struct B_RFC_Table *RFCT = 0; /* RFC-Table */
static bddp RFCT_Spc;   /* Current RFC-table size */
static bddp RFCT_Used;  /* Current RFC-table used entries */

/* Declaration of DenseBDD, added by Lee */
struct DenseBDD {
  bp *U; //zero edge tree
  bitvector *M; //dummy node vector
  bitvector *I; //one child array
  bddp width; //width of one child array
  bitvector *visit; //visit flag array
};
static struct DenseBDD *DB = 0;
#define IN_DENSE(f) ((B_ABS(f)>>1U)<FirstIdx) //return 1 if node f is in DenseBDD

#define BDD_CACHE_CHK_RETURN(op, f, g) \
{ bddp h = bddrcache(op, f, g); \
if (h!=bddnull) return h; \
BDD_RECUR_INC; }

#define BDD_CACHE_ENT_RETURN(op, f, g, h) \
{ BDD_RECUR_DEC; \
if (h!=bddnull) bddwcache(op, f, g, h); \
return h; }


/* ----- Declaration of static (internal) functions ------ */
/* Private procedure */
static int  err(const char msg[], bddp num);
static int  rfc_inc_ovf(struct B_NodeTable *np);
static int  rfc_dec_ovf(struct B_NodeTable *np);
static void var_enlarge();
static int  node_enlarge();
static int  hash_enlarge(bddvar v);
static bddp getnode(bddvar v, bddp f0, bddp f1);
static bddp getdnode(bddvar v, bddp f0, bddp f1); //added by Lee
static bddp getbddp(bddvar v, bddp f0, bddp f1);
static bddp getzbddp(bddvar v, bddp f0, bddp f1);
static bddp dzero(bddp f); //added by Lee
static bddp done(bddp f); // added by Lee
static bddvar dvar(bddp f); // added by Lee
static bddvar dlev(bddp f); // added by Lee
static bddp dvisit(bddp f); // added by Lee
static void setdvisit(bddp f, int flag); // added by Lee
static bddp apply(bddp f, bddp g, unsigned char op, unsigned char skip);
static void gc1(struct B_NodeTable *np);
static bddp count(bddp f);
static void dump(bddp f);
static void reset(bddp f);

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
  
  
  /* ------------------ Body of program -------------------- */
  /* ----------------- External functions ------------------ */
  int bddinit(bddp initsize, bddp limitsize)
  /* Returns 1 if not enough memory (usually 0) */
  {
    bddp   ix;
    bddvar i;
    
    /* Check dupulicate initialization */
    if(Node) {free(Node);  Node=NULL;}
    if(Var)
    {
      for(i=0; i<VarSpc; i++)
      {
        if(Var[i].hash_32) free(Var[i].hash_32);
#ifdef B_64
        if(Var[i].hash_h8) free(Var[i].hash_h8);
#endif
      }
      free(Var);  Var=NULL;
    }
    if(VarID) {free(VarID);  VarID=NULL;}
    if(Cache) {free(Cache);  Cache=NULL;}
    if(RFCT) {free(RFCT);  RFCT=NULL;}
    
    /* Set NodeLimit */
    if(limitsize < B_NODE_SPC0) NodeLimit = B_NODE_SPC0;
    else if(limitsize > B_NODE_MAX) NodeLimit = B_NODE_MAX;
    else NodeLimit = limitsize;
    
    /* Set NodeSpc */
    if(initsize < B_NODE_SPC0) NodeSpc = B_NODE_SPC0;
    else if(initsize > NodeLimit) NodeSpc = NodeLimit;
    else NodeSpc = initsize;
    
    /* Set CacheSpc */
    for(CacheSpc=B_NODE_SPC0; CacheSpc<NodeSpc>>1; CacheSpc<<=1U)
      ; /* empty */
    
    /* Set VarSpc */
    VarSpc = B_VAR_SPC0;
    
    /* Memory allocation */
    Node = B_MALLOC(struct B_NodeTable, NodeSpc);
    Var = B_MALLOC(struct B_VarTable, VarSpc);
    VarID = B_MALLOC(bddvar, VarSpc);
    Cache = B_MALLOC(struct B_CacheTable, CacheSpc);
    
    /* Check overflow */
    if(Node == 0 || Var == 0 || VarID == 0 || Cache == 0)
    {
      if(Cache){ free(Cache); Cache = 0; }
      if(VarID){ free(VarID); VarID = 0; }
      if(Var){ free(Var); Var = 0; }
      if(Node){ free(Node); Node = 0; }
      NodeLimit = 0;
      return 1;
    }
    
    /* Initialize */
    NodeUsed = 0;
    Node[NodeSpc-1U].varrfc = 0;
    B_SET_BDDP(Node[NodeSpc-1U].nx, bddnull);
    for(ix=0; ix<NodeSpc-1U; ix++)
    {
      Node[ix].varrfc = 0;
      B_SET_BDDP(Node[ix].nx, ix+1U);
    }
    Avail = 0;
    
    VarUsed = 0;
    for(i=0; i<VarSpc; i++)
    {
      Var[i].hashSpc = 0;
      Var[i].hashUsed = 0;
      Var[i].lev = i;
      VarID[i] = i;
      Var[i].hash_32 = 0;
#ifdef B_64
      Var[i].hash_h8 = 0;
#endif
    }
    
    for(ix=0; ix<CacheSpc; ix++) Cache[ix].op = BC_NULL;
    
    RFCT_Spc = 0;
    RFCT_Used = 0;
    
    return 0;
  }
  
  bddp bddcopy(bddp f)
  {
    struct B_NodeTable *fp;
    
    if(f == bddnull) return bddnull;
    if(B_CST(f)||IN_DENSE(f)) return f; /* Constant or in Dense */
    fp = B_NP(f);
    if(fp >= Node+NodeSpc || fp->varrfc == 0)
      err("bddcopy: Invalid bddp", f);
    B_RFC_INC_NP(fp);
    return f;
  }
  
  void bddfree(bddp f)
  {
    struct B_NodeTable *fp;
    
    if(f == bddnull) return;
    if(B_CST(f)||IN_DENSE(f)) return; /* Constant or in Dense */
    fp = B_NP(f);
    if(fp >= Node+NodeSpc || fp->varrfc == 0)
      err("bddfree: Invalid bddp", f);
    B_RFC_DEC_NP(fp);
  }
  
  int bddgc()
  /* Returns 1 if there are no free node (usually 0) */
  {
    bddp i, n, f;
    struct B_NodeTable *fp;
    struct B_CacheTable *cachep;
    struct B_NodeTable *np;
    struct B_VarTable *varp;
    bddvar v;
    bddp oldSpc, newSpc, nx, key, f0, f1;
    bddp_32 *newhash_32, *p_32, *p2_32;
#ifdef B_64
    bddp_h8 *newhash_h8, *p_h8, *p2_h8;
#endif
    
    
    n = NodeUsed;
    for(fp=Node; fp<Node+NodeSpc; fp++)
      if(fp->varrfc != 0 && B_RFC_ZERO_NP(fp))
        gc1(fp);
    if(n == NodeUsed) return 1; /* No free node */
    
    /* Cache clear */
    for(cachep=Cache; cachep<Cache+CacheSpc; cachep++)
    {
      switch(cachep->op)
      {
        case BC_NULL:
          break;
        case BC_AND:
        case BC_XOR:
        case BC_INTERSEC:
        case BC_UNION:
        case BC_SUBTRACT:
        case BC_CHANGE:
          f = B_GET_BDDP(cachep->f);
          //            added by Lee
          if(!B_CST(f) && !IN_DENSE(f) && (fp=B_NP(f))<Node+NodeSpc && fp->varrfc == 0)
          {
            cachep->op = BC_NULL;
            break;
          }
          //            end of addition
          //      if(!B_CST(f) && (fp=B_NP(f))<Node+NodeSpc && fp->varrfc == 0)
          //      {
          //        cachep->op = BC_NULL;
          //        break;
          //      }
          f = B_GET_BDDP(cachep->g);
          //            added by Lee
          if(!B_CST(f) && !IN_DENSE(f) && (fp=B_NP(f))<Node+NodeSpc && fp->varrfc == 0)
          {
            cachep->op = BC_NULL;
            break;
          }
          //            end of addition
          //      if(!B_CST(f) && (fp=B_NP(f))<Node+NodeSpc && fp->varrfc == 0)
          //      {
          //        cachep->op = BC_NULL;
          //        break;
          //      }
          f = B_GET_BDDP(cachep->h);
          //            added by Lee
          if(!B_CST(f) && !IN_DENSE(f) && (fp=B_NP(f))<Node+NodeSpc && fp->varrfc == 0)
          {
            cachep->op = BC_NULL;
            break;
          }
          //            end of addition
          //      if(!B_CST(f) && (fp=B_NP(f))<Node+NodeSpc && fp->varrfc == 0)
          //      {
          //        cachep->op = BC_NULL;
          //        break;
          //      }
          break;
        case BC_AT0:
        case BC_AT1:
        case BC_OFFSET:
        case BC_ONSET:
          f = B_GET_BDDP(cachep->f);
          //            added by Lee
          if(!B_CST(f) && !IN_DENSE(f) && (fp=B_NP(f))<Node+NodeSpc && fp->varrfc == 0)
          {
            cachep->op = BC_NULL;
            break;
          }
          //            end of addition
          //      if(!B_CST(f) && (fp=B_NP(f))<Node+NodeSpc && fp->varrfc == 0)
          //      {
          //        cachep->op = BC_NULL;
          //        break;
          //      }
          f = B_GET_BDDP(cachep->h);
          //            added by Lee
          if(!B_CST(f) && !IN_DENSE(f) && (fp=B_NP(f))<Node+NodeSpc && fp->varrfc == 0)
          {
            cachep->op = BC_NULL;
            break;
          }
          //            end of addition
          //      if(!B_CST(f) && (fp=B_NP(f))<Node+NodeSpc && fp->varrfc == 0)
          //      {
          //        cachep->op = BC_NULL;
          //        break;
          //      }
          break;
        case BC_CARD:
        case BC_LIT:
        case BC_LEN:
          f = B_GET_BDDP(cachep->f);
          //            added by Lee
          if(!B_CST(f) && !IN_DENSE(f) && (fp=B_NP(f))<Node+NodeSpc && fp->varrfc == 0)
          {
            cachep->op = BC_NULL;
            break;
          }
          //            end of addition
          //      if(!B_CST(f) && (fp=B_NP(f))<Node+NodeSpc && fp->varrfc == 0)
          //      {
          //        cachep->op = BC_NULL;
          //        break;
          //      }
          break;
        default:
          cachep->op = BC_NULL;
          break;
      }
    }
    
    /* Hash-table packing */
    for(v=1; v<=VarUsed; v++)
    {
      varp = &Var[v];
      
      /* Get new size */
      oldSpc = varp->hashSpc;
      newSpc = oldSpc;
      while(newSpc > B_HASH_SPC0)
      {
        if(newSpc>>2 < varp->hashUsed) break;
        newSpc >>= 1;
      }
      if(newSpc == oldSpc) continue;
      
      /* Reduce space */
#ifdef B_64
      newhash_32 = B_MALLOC(bddp_32, newSpc);
      newhash_h8 = B_MALLOC(bddp_h8, newSpc);
      if(!newhash_32 || !newhash_h8)
      {
        if(newhash_32) {free(newhash_32);  newhash_32=NULL;}
        if(newhash_h8) {free(newhash_h8);  newhash_h8=NULL;}
        break; /* Not enough memory */
      }
#else
      newhash_32 = B_MALLOC(bddp_32, newSpc);
      if(!newhash_32) break; /* Not enough memory */
#endif
      
      /* Initialize new hash entry */
      for(i=0; i<newSpc; i++)
      {
        B_SET_NXP(p, newhash, i);
        B_SET_BDDP(*p, bddnull);
      }
      
      /* restore hash entry */
      for(i=0; i<oldSpc; i++)
      {
        key = i & (newSpc-1U);
        np = 0;
        B_SET_NXP(p, newhash, key);
        nx = B_GET_BDDP(*p);
        while(nx != bddnull)
        {
          np = Node + nx;
          nx = B_GET_BDDP(np->nx);
        }
        if(np) { B_SET_NXP(p2, varp->hash, i); B_CPY_BDDP(np->nx, *p2); }
        else
        {
          B_SET_NXP(p, newhash, key);
          B_SET_NXP(p2, varp->hash, i);
          B_CPY_BDDP(*p, *p2);
        }
      }
      varp->hashSpc = newSpc;
      free(varp->hash_32);
      varp->hash_32 = newhash_32;
#ifdef B_64
      free(varp->hash_h8);
      varp->hash_h8 = newhash_h8;
#endif
    }
    return 0;
  }
  
  bddp bddused() { return NodeUsed; }
  
  bddp bddsize(bddp f)
  /* Returns 0 for bddnull */
  {
    bddp num;
    struct B_NodeTable *fp;
    
    if(f == bddnull) return 0;
    if(B_CST(f)) return 0; /* Constant */
    //    added by Lee
    if (!IN_DENSE(f)) {
      if((fp=B_NP(f))>=Node+NodeSpc || fp->varrfc == 0)
        err("bddsize: Invalid bddp", f);
    }
    //    end of addition
    //  if((fp=B_NP(f))>=Node+NodeSpc || fp->varrfc == 0)
    //    err("bddsize: Invalid bddp", f);
    
    num = count(f);
    reset(f);
    return num;
  }
  
  bddp bddvsize(bddp *p, int lim)
  /* Returns 0 for bddnull */
  {
    bddp num;
    struct B_NodeTable *fp;
    int n, i;
    
    /* Check operand */
    n = lim;
    for(i=0; i<n; i++)
    {
      if(p[i] == bddnull)
      {
        n = i;
        break;
      }
      //      added by Lee
      if (!IN_DENSE(p[i])) {
        if(!B_CST(p[i])&&
           ((fp=B_NP(p[i]))>=Node+NodeSpc || fp->varrfc==0))
          err("bddvsize: Invalid bddp", p[i]);
      }
      //      end of addition
      //    if(!B_CST(p[i])&&
      //       ((fp=B_NP(p[i]))>=Node+NodeSpc || fp->varrfc==0))
      //      err("bddvsize: Invalid bddp", p[i]);
    }
    num = 0;
    for(i=0; i<n; i++) num += count(p[i]);
    for(i=0; i<n; i++) reset(p[i]);
    return num;
  }
  
  void bdddump(bddp f)
  {
    struct B_NodeTable *fp;
    
    /* Check indexes */
    if(f == bddnull) { printf("RT = NULL\n\n"); return; }
    //    added by Lee
    if (!IN_DENSE(f)) {
      if(!B_CST(f)&&
         ((fp=B_NP(f))>=Node+NodeSpc || fp->varrfc==0))
        err("bdddump: Invalid bddp", f);
    }
    //    end of addition
    //  if(!B_CST(f)&&
    //     ((fp=B_NP(f))>=Node+NodeSpc || fp->varrfc==0))
    //      err("bdddump: Invalid bddp", f);
    
    /* Dump nodes */
    dump(f);
    reset(f);
    
    /* Dump top node */
    printf("RT = ");
    if(B_NEG(f)) putchar('~');
    if(B_CST(f)) printf(B_BDDP_FD, B_ABS(B_VAL(f)));
    else { printf("N"); printf(B_BDDP_FD, B_NDX(f));
    }
    printf("\n\n");
  }
  
  void bddvdump(bddp *p, int n)
  {
    struct B_NodeTable *fp;
    int i;
    
    /* Check operands */
    for(i=0; i<n; i++)
    {
      if(p[i] == bddnull) return;
      //      added by Lee
      if (!IN_DENSE(p[i])) {
        if(!B_CST(p[i])&&
           ((fp=B_NP(p[i]))>=Node+NodeSpc || fp->varrfc==0))
          err("bddvdump: Invalid bddp", p[i]);
      }
      //      end of addition
      //    if(!B_CST(p[i])&&
      //       ((fp=B_NP(p[i]))>=Node+NodeSpc || fp->varrfc==0))
      //      err("bddvdump: Invalid bddp", p[i]);
    }
    
    /* Dump nodes */
    for(i=0; i<n; i++) if(p[i] != bddnull) dump(p[i]);
    for(i=0; i<n; i++) if(p[i] != bddnull) reset(p[i]);
    
    /* Dump top node */
    for(i=0; i<n; i++)
    {
      printf("RT%d = ", i);
      if(p[i] == bddnull) printf("NULL");
      else
      {
        if(B_NEG(p[i])) putchar('~');
        if(B_CST(p[i])) printf(B_BDDP_FD, B_ABS(B_VAL(p[i])));
        else { printf("N"); printf(B_BDDP_FD, B_NDX(p[i])); }
      }
      putchar('\n');
    }
    printf("\n");
  }
  
  bddp bddrcache(unsigned char op, bddp f, bddp g)
  {
    struct B_CacheTable *cachep;
    
    cachep = Cache + B_CACHEKEY(op, f, g);
    if(op == cachep->op &&
       f == B_GET_BDDP(cachep->f) &&
       g == B_GET_BDDP(cachep->g))
      return B_GET_BDDP(cachep->h); /* Hit */
    return bddnull;
  }
  
  void bddwcache(unsigned char op, bddp f, bddp g, bddp h)
  {
    struct B_CacheTable *cachep;
    
    if(op < 20) err("bddwcache: op < 20", op);
    if(h == bddnull) return;
    cachep = Cache + B_CACHEKEY(op, f, g);
    cachep->op = op;
    B_SET_BDDP(cachep->f, f);
    B_SET_BDDP(cachep->g, g);
    B_SET_BDDP(cachep->h, h);
  }
  
  bddp bddnot(bddp f)
  {
    if(f == bddnull) return bddnull;
    return B_NOT(bddcopy(f));
  }
  
  bddvar bddlevofvar(bddvar v)
  {
    if(v > VarUsed)
      err("bddlevofvar: Invalid VarID", v);
    return Var[v].lev;
  }
  
  bddvar bddvaroflev(bddvar lev)
  {
    if(lev > VarUsed)
      err("bddvaroflev: Invalid level", lev);
    return VarID[lev];
  }
  
  bddvar bddvarused()
  {
    return VarUsed;
  }
  
  bddvar bddnewvar()
  {
    if(++VarUsed == VarSpc) var_enlarge();
    return VarUsed;
  }
  
  /* not supported yet */
  bddvar bddnewvaroflev(bddvar lev)
  {
    bddvar i, v;
    
    if(lev == 0 || lev > ++VarUsed)
      err("bddnewvaroflev: Invalid level", lev);
    if(VarUsed == VarSpc) var_enlarge();
    for(i=VarUsed; i>lev; i--) Var[ VarID[i] = VarID[i-1U] ].lev = i;
    Var[ VarID[lev] = VarUsed ].lev = lev;
    return VarUsed;
  }
  
  bddvar bddtop(bddp f)
  {
    struct B_NodeTable *fp;
    
    if(f == bddnull) return 0;
    if(B_CST(f)) return 0; /* Constant */
    //    added by Lee
    if (IN_DENSE(f)) {
      return dvar(f);
    }
    //    end of addition
    fp = B_NP(f);
    if(fp >= Node+NodeSpc || fp->varrfc == 0)
      err("bddtop: Invalid bddp", f);
    return B_VAR_NP(fp);
  }
  
  bddp bddsupport(bddp f)
  /* Returns bddnull if not enough memory */
  {
    struct B_NodeTable *fp;
    
    /* Check operands */
    if(f == bddnull) return bddnull;
    if(B_CST(f)) return bddfalse;
    //    added by Lee
    if (!IN_DENSE(f)) {
      fp = B_NP(f);
      if(fp>=Node+NodeSpc || !fp->varrfc)
        err("bddsupport: Invalid bddp", f);
    }
    //    end of addition
    //  fp = B_NP(f);
    //  if(fp>=Node+NodeSpc || !fp->varrfc)
    //    err("bddsupport: Invalid bddp", f);
    
    return apply(f, bddfalse, BC_SUPPORT, 0);
  }
  
  bddp bddlshift(bddp f, bddvar shift)
  /* Returns bddnull if not enough memory */
  {
    struct B_NodeTable *fp;
    bddvar flev;
    
    /* Check operands */
    if(shift >= VarUsed)
      err("bddlshift: Invalid shift", shift);
    if(f == bddnull) return bddnull;
    if(B_CST(f)) return f;
    if(shift == 0) return bddcopy(f);
    //      added by Lee
    if (!IN_DENSE(f)) {
      if((fp=B_NP(f))>=Node+NodeSpc || !fp->varrfc)
        err("bddlshift: Invalid bddp", f);
    }
    //    end of addition
    //  if((fp=B_NP(f))>=Node+NodeSpc || !fp->varrfc)
    //    err("bddlshift: Invalid bddp", f);
    
    return apply(f, (bddp)shift, BC_LSHIFT, 0);
  }
  
  bddp bddrshift(bddp f, bddvar shift)
  /* Returns bddnull if not enough memory */
  {
    struct B_NodeTable *fp;
    bddvar flev;
    
    /* Check operands */
    if(shift >= VarUsed)
      err("bddrshift: Invalid shift", shift);
    if(f == bddnull) return bddnull;
    if(B_CST(f)) return f;
    if(shift == 0) return bddcopy(f);
    //      added by Lee
    if (!IN_DENSE(f)) {
      if((fp=B_NP(f))>=Node+NodeSpc || !fp->varrfc)
        err("bddrshift: Invalid bddp", f);
    }
    //    end of addition
    //  if((fp=B_NP(f))>=Node+NodeSpc || !fp->varrfc)
    //    err("bddrshift: Invalid bddp", f);
    
    return apply(f, (bddp)shift, BC_RSHIFT, 0);
  }
  
  bddp    bddoffset(bddp f, bddvar v)
  /* Returns bddnull if not enough memory */
  {
    struct B_NodeTable *fp;
    
    /* Check operands */
    if(v > VarUsed || v == 0) err("bddoffset: Invalid VarID", v);
    if(f == bddnull) return bddnull;
    if(B_CST(f)) return f;
    //    added by Lee
    if (IN_DENSE(f)) {
      if (!(dzero(f)&1)) err("bddoffset: applying non-ZBDD node", f);
    }
    else {
      fp = B_NP(f);
      if(fp>=Node+NodeSpc || !fp->varrfc)
        err("bddoffset: Invalid bddp", f);
      if(!B_Z_NP(fp)) err("bddoffset: applying non-ZBDD node", f);
    }
    //    end of addition
    //  fp = B_NP(f);
    //  if(fp>=Node+NodeSpc || !fp->varrfc)
    //    err("bddoffset: Invalid bddp", f);
    //  if(!B_Z_NP(fp)) err("bddoffset: applying non-ZBDD node", f);
    
    return apply(f, (bddp)v, BC_OFFSET, 0);
  }
  
  bddp    bddonset0(bddp f, bddvar v)
  /* Returns bddnull if not enough memory */
  {
    struct B_NodeTable *fp;
    
    /* Check operands */
    if(v > VarUsed || v == 0) err("bddonset0: Invalid VarID", v);
    if(f == bddnull) return bddnull;
    if(B_CST(f)) return bddfalse;
    //    added by Lee
    if (IN_DENSE(f)) {
      if (!(dzero(f)&1)) err("bddonset0: applying non-ZBDD node", f);
    }
    else {
      fp = B_NP(f);
      if(fp>=Node+NodeSpc || !fp->varrfc)
        err("bddonset0: Invalid bddp", f);
      if(!B_Z_NP(fp)) err("bddonset0: applying non-ZBDD node", f);
    }
    //    end of addition
    //  fp = B_NP(f);
    //  if(fp>=Node+NodeSpc || !fp->varrfc)
    //    err("bddonset0: Invalid bddp", f);
    //  if(!B_Z_NP(fp)) err("bddonset0: applying non-ZBDD node", f);
    
    return apply(f, (bddp)v, BC_ONSET, 0);
  }
  
  bddp    bddonset(bddp f, bddvar v)
  /* Returns bddnull if not enough memory */
  {
    bddp g, h;
    
    g = bddonset0(f, v);
    h = bddchange(g, v);
    bddfree(g);
    return h;
  }
  
  bddp    bddchange(bddp f, bddvar v)
  /* Returns bddnull if not enough memory */
  {
    struct B_NodeTable *fp;
    
    /* Check operands */
    if(v > VarUsed || v == 0) err("bddchange: Invalid VarID", v);
    if(f == bddnull) return bddnull;
    if(!B_CST(f))
    {
      //      added by Lee
      if (IN_DENSE(f)) {
        if (!(dzero(f)&1)) err("bddchange: applying non-ZBDD node", f);
      }
      else {
        fp = B_NP(f);
        if(fp>=Node+NodeSpc || !fp->varrfc)
          err("bddchange: Invalid bddp", f);
        if(!B_Z_NP(fp)) err("bddchange: applying non-ZBDD node", f);
      }
      //      end of addition
      //    fp = B_NP(f);
      //    if(fp>=Node+NodeSpc || !fp->varrfc)
      //      err("bddchange: Invalid bddp", f);
      //    if(!B_Z_NP(fp)) err("bddchange: applying non-ZBDD node", f);
    }
    
    return apply(f, (bddp)v, BC_CHANGE, 0);
  }
  
  bddp bddintersec(bddp f, bddp g)
  /* Returns bddnull if not enough memory */
  {
    struct B_NodeTable *fp;
    
    /* Check operands */
    if(f == bddnull) return bddnull;
    if(g == bddnull) return bddnull;
    if(B_CST(f))
    { if(B_ABS(f) != bddfalse) err("bddintersec: Invalid bddp", f); }
    else
    {
      //      added by Lee
      if (IN_DENSE(f)) {
        if (!(dzero(f)&1)) err("bddintersec: applying non-ZBDD node", f);
      }
      else {
        fp = B_NP(f);
        if(fp>=Node+NodeSpc || !fp->varrfc)
          err("bddintersec: Invalid bddp", f);
        if(!B_Z_NP(fp)) err("bddintersec: applying non-ZBDD node", f);
      }
      //      end of addition
      //    fp = B_NP(f);
      //    if(fp>=Node+NodeSpc || !fp->varrfc)
      //      err("bddintersec: Invalid bddp", f);
      //    if(!B_Z_NP(fp)) err("bddintersec: applying non-ZBDD node", f);
    }
    if(B_CST(g))
    { if(B_ABS(g) != bddfalse) err("bddintersec: Invalid bddp", g); }
    else
    {
      //      added by Lee
      if (IN_DENSE(g)) {
        if (!(dzero(g)&1)) err("bddintersec: applying non-ZBDD node", g);
      }
      else {
        fp = B_NP(g);
        if(fp>=Node+NodeSpc || !fp->varrfc)
          err("bddintersec: Invalid bddp", g);
        if(!B_Z_NP(fp)) err("bddintersec: applying non-ZBDD node", g);
      }
      //      end of addition
      //    fp = B_NP(g);
      //    if(fp>=Node+NodeSpc || !fp->varrfc)
      //      err("bddintersec: Invalid bddp", g);
      //    if(!B_Z_NP(fp)) err("bddintersec: applying non-ZBDD node", g);
    }
    
    return apply(f, g, BC_INTERSEC, 0);
  }
  
  bddp bddunion(bddp f, bddp g)
  /* Returns bddnull if not enough memory */
  {
    struct B_NodeTable *fp;
    
    /* Check operands */
    if(f == bddnull) return bddnull;
    if(g == bddnull) return bddnull;
    if(B_CST(f))
    { if(B_ABS(f) != bddfalse) err("bddunion: Invalid bddp", f); }
    else
    {
      //      added by Lee
      if (IN_DENSE(f)) {
        if (!(dzero(f)&1)) err("bddunion: applying non-ZBDD node", f);
      }
      else {
        fp = B_NP(f);
        if(fp>=Node+NodeSpc || !fp->varrfc)
          err("bddunion: Invalid bddp", f);
        if(!B_Z_NP(fp)) err("bddunion: applying non-ZBDD node", f);
      }
      //      end of addition
      //    fp = B_NP(f);
      //    if(fp>=Node+NodeSpc || !fp->varrfc)
      //      err("bddunion: Invalid bddp", f);
      //    if(!B_Z_NP(fp)) err("bddunion: applying non-ZBDD node", f);
    }
    if(B_CST(g))
    { if(B_ABS(g) != bddfalse) err("bddunion: Invalid bddp", g); }
    else
    {
      //      added by Lee
      if (IN_DENSE(g)) {
        if (!(dzero(g)&1)) err("bddunion: applying non-ZBDD node", f);
      }
      else {
        fp = B_NP(g);
        if(fp>=Node+NodeSpc || !fp->varrfc)
          err("bddunion: Invalid bddp", g);
        if(!B_Z_NP(fp)) err("bddunion: applying non-ZBDD node", g);
      }
      //      end of addition
      //    fp = B_NP(g);
      //    if(fp>=Node+NodeSpc || !fp->varrfc)
      //      err("bddunion: Invalid bddp", g);
      //    if(!B_Z_NP(fp)) err("bddunion: applying non-ZBDD node", g);
    }
    
    return apply(f, g, BC_UNION, 0);
  }
  
  bddp bddsubtract(bddp f, bddp g)
  /* Returns bddnull if not enough memory */
  {
    struct B_NodeTable *fp;
    
    /* Check operands */
    if(f == bddnull) return bddnull;
    if(g == bddnull) return bddnull;
    if(B_CST(f))
    { if(B_ABS(f) != bddfalse) err("bddsubtract: Invalid bddp", f); }
    else
    {
      //      added by Lee
      if (IN_DENSE(f)) {
        if (!(dzero(f)&1)) err("bddsubtarct: applying non-ZBDD node", f);
      }
      else {
        fp = B_NP(f);
        if(fp>=Node+NodeSpc || !fp->varrfc)
          err("bddsubtarct: Invalid bddp", f);
        if(!B_Z_NP(fp)) err("bddsubtarct: applying non-ZBDD node", f);
      }
      //      end of addition
      //    fp = B_NP(f);
      //    if(fp>=Node+NodeSpc || !fp->varrfc)
      //      err("bddsubtarct: Invalid bddp", f);
      //    if(!B_Z_NP(fp)) err("bddsubtarct: applying non-ZBDD node", f);
    }
    if(B_CST(g))
    { if(B_ABS(g) != bddfalse) err("bddsubtarct: Invalid bddp", g); }
    else
    {
      //      added by Lee
      if (IN_DENSE(g)) {
        if (!(dzero(g)&1)) err("bddsubtarct: applying non-ZBDD node", g);
      }
      else {
        fp = B_NP(g);
        if(fp>=Node+NodeSpc || !fp->varrfc)
          err("bddsubtarct: Invalid bddp", g);
        if(!B_Z_NP(fp)) err("bddsubtarct: applying non-ZBDD node", g);
      }
      //      end of addition
      //    fp = B_NP(g);
      //    if(fp>=Node+NodeSpc || !fp->varrfc)
      //      err("bddsubtarct: Invalid bddp", g);
      //    if(!B_Z_NP(fp)) err("bddsubtarct: applying non-ZBDD node", g);
    }
    
    return apply(f, g, BC_SUBTRACT, 0);
  }
  
  bddp bddcard(bddp f)
  {
    struct B_NodeTable *fp;
    
    if(f == bddnull) return 0;
    if(B_CST(f)) return (f == bddempty)? 0: 1;
    //    added by Lee
    if (IN_DENSE(f)) {
      if (!(dzero(f)&1)) err("bddcard: applying non-ZBDD node", f);
    }
    else {
      fp = B_NP(f);
      if(fp>=Node+NodeSpc || !fp->varrfc)
        err("bddcard: Invalid bddp", f);
      if(!B_Z_NP(fp)) err("bddcard: applying non-ZBDD node", f);
    }
    //    end of addition
    //  fp = B_NP(f);
    //  if(fp>=Node+NodeSpc || !fp->varrfc)
    //    err("bddcard: Invalid bddp", f);
    //  if(!B_Z_NP(fp)) err("bddcard: applying non-ZBDD node", f);
    
    return apply(f, bddempty, BC_CARD, 0);
  }
  
  bddp bddlit(bddp f)
  {
    struct B_NodeTable *fp;
    
    if(f == bddnull) return 0;
    if(B_CST(f)) return 0;
    //    added by Lee
    if (IN_DENSE(f)) {
      if (!(dzero(f)&1)) {
        err("bddlit: applying non-ZBDD node", f);
      }
    }
    else {
      fp = B_NP(f);
      if(fp>=Node+NodeSpc || !fp->varrfc)
        err("bddlit: Invalid bddp", f);
      if(!B_Z_NP(fp)) err("bddlit: applying non-ZBDD node", f);
    }
    //    end of addition
    //  fp = B_NP(f);
    //  if(fp>=Node+NodeSpc || !fp->varrfc)
    //    err("bddlit: Invalid bddp", f);
    //  if(!B_Z_NP(fp)) err("bddlit: applying non-ZBDD node", f);
    
    return apply(f, bddempty, BC_LIT, 0);
  }
  
  bddp bddlen(bddp f)
  {
    struct B_NodeTable *fp;
    
    if(f == bddnull) return 0;
    if(B_CST(f)) return 0;
    //    added by Lee
    if (IN_DENSE(f)) {
      if (! dzero(f)&1) err("bddlen: applying non-ZBDD node", f);
    }
    else {
      fp = B_NP(f);
      if(fp>=Node+NodeSpc || !fp->varrfc)
        err("bddlen: Invalid bddp", f);
      if(!B_Z_NP(fp)) err("bddlen: applying non-ZBDD node", f);
    }
    //  fp = B_NP(f);
    //  if(fp>=Node+NodeSpc || !fp->varrfc)
    //    err("bddlen: Invalid bddp", f);
    //  if(!B_Z_NP(fp)) err("bddlen: applying non-ZBDD node", f);
    
    return apply(f, bddempty, BC_LEN, 0);
  }
  
  //added by Lee
  bddp bddmult(bddp f, bddp g)
  {
    struct B_NodeTable *fp;
    
    /* Check operands */
    if(f == bddnull) return bddnull;
    if(g == bddnull) return bddnull;
    if(B_CST(f))
    { if(B_ABS(f) != bddfalse) err("bddmult: Invalid bddp", f); }
    else
    {
      if (IN_DENSE(f)) {
        if (!(dzero(f)&1)) err("bddmult: applying non-ZBDD node", f);
      }
      else {
        fp = B_NP(f);
        if(fp>=Node+NodeSpc || !fp->varrfc)
          err("bddmult: Invalid bddp", f);
        if(!B_Z_NP(fp)) err("bddmult: applying non-ZBDD node", f);
      }
    }
    if(B_CST(g))
    { if(B_ABS(g) != bddfalse) err("bddmult: Invalid bddp", g); }
    else
    {
      //      added by Lee
      if (IN_DENSE(g)) {
        if (!(dzero(g)&1)) err("bddmult: applying non-ZBDD node", g);
      }
      else {
        fp = B_NP(g);
        if(fp>=Node+NodeSpc || !fp->varrfc)
          err("bddmult: Invalid bddp", g);
        if(!B_Z_NP(fp)) err("bddmult: applying non-ZBDD node", g);
      }
    }
    
    //    check trivial case
    if (f==bddempty||g==bddempty) {
      return bddempty;
    }
    if (f==bddsingle) {
      return g;
    }
    if (g==bddsingle) {
      return f;
    }
    
    bddvar ftop = bddtop(f);
    bddvar gtop = bddtop(g);
    if (bddlevofvar(ftop)<bddlevofvar(gtop)) {
      //        swap operands
      bddp tempp = f;
      f = g;
      g = tempp;
      bddvar tempv = ftop;
      ftop = gtop;
      gtop = tempv;
    }
    if (ftop==gtop && f<g) {
      //        swap operands
      bddp tempp = f;
      f = g;
      g = tempp;
    }
    
    BDD_CACHE_CHK_RETURN(BC_MULT, f, g);
    
    bddp f1 = bddonset0(f, ftop);
    bddp f0 = bddoffset(f, ftop);
    bddp h;
    if (ftop!=gtop) {
      h = bddmult(f1, g);
      h = bddunion(bddchange(h, ftop), bddmult(f0, g));
    }
    else {
      bddp g1 = bddonset0(g, ftop);
      bddp g0 = bddoffset(g, ftop);
      h = bddunion(bddunion(bddmult(f1, g1), bddmult(f1, g0)), bddmult(f0, g0));
      h = bddunion(bddchange(h, ftop), bddmult(f0, g0));
    }
    BDD_CACHE_ENT_RETURN(BC_MULT, f, g, h);
  }
  
  //added by Lee
  bddp bdddiv(bddp f, bddp p)
  {
    struct B_NodeTable *fp;
    
    /* Check operands */
    if(f == bddnull) return bddnull;
    if(p == bddnull) return bddnull;
    if(B_CST(f))
    { if(B_ABS(f) != bddfalse) err("bdddiv: Invalid bddp", f); }
    else
    {
      if (IN_DENSE(f)) {
        if (!(dzero(f)&1)) err("bdddiv: applying non-ZBDD node", f);
      }
      else {
        fp = B_NP(f);
        if(fp>=Node+NodeSpc || !fp->varrfc)
          err("bdddiv: Invalid bddp", f);
        if(!B_Z_NP(fp)) err("bdddiv: applying non-ZBDD node", f);
      }
    }
    if(B_CST(p))
    { if(B_ABS(p) != bddfalse) err("bdddiv: Invalid bddp", p); }
    else
    {
      //      added by Lee
      if (IN_DENSE(p)) {
        if (!(dzero(p)&1)) err("bdddiv: applying non-ZBDD node", p);
      }
      else {
        fp = B_NP(p);
        if(fp>=Node+NodeSpc || !fp->varrfc)
          err("bdddiv: Invalid bddp", p);
        if(!B_Z_NP(fp)) err("bdddiv: applying non-ZBDD node", p);
      }
    }
    
    //    check trivial case
    if (p==bddsingle) {
      return f;
    }
    if (f==p) {
      return bddsingle;
    }
    if (p==bddempty) {
      err("bdddiv: divided by zero", p);
    }
    bddvar top = bddtop(p);
    if (bddlevofvar(bddtop(f))<bddlevofvar(top)) {
      return bddempty;
    }
    
    BDD_CACHE_CHK_RETURN(BC_DIV, f, p);
    
    bddp q = bdddiv(bddonset0(f, top), bddonset0(p, top));
    if (q!=bddempty) {
      bddp p0 = bddoffset(p, top);
      if (p0!=bddempty) {
        q = bddintersec(q, bdddiv(bddoffset(f, top), p0));
      }
    }
    BDD_CACHE_ENT_RETURN(BC_DIV, f, p, q);
  }
  
  //added by Lee
  bddp bddquot(bddp f, bddp p)
  {
    struct B_NodeTable *fp;
    
    /* Check operands */
    if(f == bddnull) return bddnull;
    if(p == bddnull) return bddnull;
    if(B_CST(f))
    { if(B_ABS(f) != bddfalse) err("bddquot: Invalid bddp", f); }
    else
    {
      if (IN_DENSE(f)) {
        if (!(dzero(f)&1)) err("bddquot: applying non-ZBDD node", f);
      }
      else {
        fp = B_NP(f);
        if(fp>=Node+NodeSpc || !fp->varrfc)
          err("bdddiv: Invalid bddp", f);
        if(!B_Z_NP(fp)) err("bddquot: applying non-ZBDD node", f);
      }
    }
    if(B_CST(p))
    { if(B_ABS(p) != bddfalse) err("bddquot: Invalid bddp", p); }
    else
    {
      //      added by Lee
      if (IN_DENSE(p)) {
        if (!(dzero(p)&1)) err("bddquot: applying non-ZBDD node", p);
      }
      else {
        fp = B_NP(p);
        if(fp>=Node+NodeSpc || !fp->varrfc)
          err("bddquot: Invalid bddp", p);
        if(!B_Z_NP(fp)) err("bddquot: applying non-ZBDD node", p);
      }
    }
    
    return bddsubtract(f, bddmult(bdddiv(f, p), p));
  }
  
  void bddexport(FILE *strm, bddp *p, int lim) {printf("error\n"); exit(1);}
  int bddimport(FILE *strm, bddp *p, int lim) {printf("error\n"); exit(1);}
  int bddimportz(FILE *strm, bddp *p, int lim) {printf("error\n"); exit(1);}
  
  bddp   bddprime(bddvar v){printf("error\n"); exit(1);}
  bddp   bddand(bddp f, bddp g){printf("error\n"); exit(1);}
  bddp   bddor(bddp f, bddp g){printf("error\n"); exit(1);}
  bddp   bddxor(bddp f, bddp g){printf("error\n"); exit(1);}
  bddp   bddnand(bddp f, bddp g){printf("error\n"); exit(1);}
  bddp   bddnor(bddp f, bddp g){printf("error\n"); exit(1);}
  bddp   bddxnor(bddp f, bddp g){printf("error\n"); exit(1);}
  bddp   bddat0(bddp f, bddvar v){printf("error\n"); exit(1);}
  bddp   bddat1(bddp f, bddvar v){printf("error\n"); exit(1);}
  bddp   bddcofactor(bddp f, bddp g){printf("error\n"); exit(1);}
  
  
  
#ifdef __cplusplus
}
#endif /* __cplusplus */

/* ----------------- Internal functions ------------------ */
static void var_enlarge()
{
  bddvar i, newSpc;
  struct B_VarTable *newVar;
  unsigned int *newVarID;
  
  /* Get new size */
  if(VarSpc == bddvarmax+1U)
    err("var_enlarge: var index range full", VarSpc);
  newSpc = VarSpc << 2U;
  if(newSpc > bddvarmax+1) newSpc = bddvarmax+1U;
  
  /* Enlarge space */
  newVar = B_MALLOC(struct B_VarTable, newSpc);
  newVarID = B_MALLOC(unsigned int, newSpc);
  if(newVar && newVarID)
  {
    for(i=0; i<VarSpc; i++)
    {
      newVar[i].hashSpc = Var[i].hashSpc;
      newVar[i].hashUsed = Var[i].hashUsed;
      newVar[i].lev = Var[i].lev;
      newVar[i].hash_32 = Var[i].hash_32;
      newVarID[i] = VarID[i];
#ifdef B_64
      newVar[i].hash_h8 = Var[i].hash_h8;
#endif
    }
    free(Var);
    free(VarID);
    Var = newVar;
    VarID = newVarID;
  }
  else
  {
    if(newVar) {free(newVar);  newVar=NULL;}
    if(newVarID) {free(newVarID);  newVarID=NULL;}
    err("var_enlarge: memory allocation failed", VarSpc);
  }
  
  /* Initialize new space */
  for(i=VarSpc; i<newSpc; i++)
  {
    Var[i].hashSpc = 0;
    Var[i].hashUsed = 0;
    Var[i].lev = i;
    Var[i].hash_32 = 0;
    VarID[i] = i;
#ifdef B_64
    Var[i].hash_h8 = 0;
#endif
  }
  VarSpc = newSpc;
}

static int node_enlarge()
/* Returns 1 if not enough memory */
{
  bddp i, newSpc;
  struct B_NodeTable *newNode;
  struct B_CacheTable *newCache, *cp, *cp1;
  
  /* Get new size */
  if(NodeSpc == NodeLimit) return 1; /* Cannot enlarge */
  newSpc = NodeSpc << 1U;
  if(newSpc > NodeLimit) newSpc = NodeLimit;
  
  /* Enlarge space */
  newNode = B_MALLOC(struct B_NodeTable, newSpc);
  if(newNode)
  {
    for(i=0; i<NodeSpc; i++)
    {
      newNode[i].varrfc = Node[i].varrfc;
      newNode[i].f0_32 = Node[i].f0_32;
      newNode[i].f1_32 = Node[i].f1_32;
      newNode[i].nx_32 = Node[i].nx_32;
#ifdef B_64
      newNode[i].f0_h8 = Node[i].f0_h8;
      newNode[i].f1_h8 = Node[i].f1_h8;
      newNode[i].nx_h8 = Node[i].nx_h8;
#endif /* B_64 */
    }
    free(Node);
    Node = newNode;
  }
  else return 1; /* Not enough memory */
  
  /* Initialize new space */
  Node[newSpc-1U].varrfc = 0;
  B_SET_BDDP(Node[newSpc-1U].nx, Avail);
  for(i=NodeSpc; i<newSpc-1U; i++)
  {
    Node[i].varrfc = 0;
    B_SET_BDDP(Node[i].nx, i+1U);
  }
  Avail = NodeSpc;
  NodeSpc = newSpc;
  //    added by Lee
#ifdef MY_DEBUG
  cout << "node enlarge:" << newSpc << endl;
#endif
  /* Realloc Cache */
  for(newSpc=CacheSpc; newSpc<NodeSpc>>1U; newSpc<<=1U)
    ; /* empty */
  newCache = B_MALLOC(struct B_CacheTable, newSpc);
  if(newCache)
  {
    for(i=0; i<CacheSpc; i++)
    {
      cp = newCache + i;
      cp1 = Cache + i;
      cp->op = cp1->op;
      B_CPY_BDDP(cp->f, cp1->f);
      B_CPY_BDDP(cp->g, cp1->g);
      B_CPY_BDDP(cp->h, cp1->h);
    }
    free(Cache);
    Cache = newCache;
  }
  else return 0; /* Only NodeTable enlarged */
  
  /* Reconstruct Cache */
  for(i=CacheSpc; i<newSpc; i++)
  {
    cp = Cache + i;
    cp1 = Cache + (i & (CacheSpc - 1));
    cp->op = cp1->op;
    B_CPY_BDDP(cp->f, cp1->f);
    B_CPY_BDDP(cp->g, cp1->g);
    B_CPY_BDDP(cp->h, cp1->h);
  }
  CacheSpc = newSpc;
  
  return 0;
}

static int hash_enlarge(bddvar v)
/* Returns 1 if not enough memory */
{
  struct B_NodeTable *np, *np0;
  struct B_VarTable *varp;
  bddp i, oldSpc, newSpc, nx, key, f0, f1;
  bddp_32 *newhash_32, *p_32;
#ifdef B_64
  bddp_h8 *newhash_h8, *p_h8;
#endif
  
  varp = &Var[v];
  /* Get new size */
  oldSpc = varp->hashSpc;
  if(oldSpc == B_NODE_MAX + 1U)
    return 0; /*  Cancel enlarging */
  newSpc = oldSpc << 1U;
  
  /* Enlarge space */
#ifdef B_64
  newhash_32 = B_MALLOC(bddp_32, newSpc);
  newhash_h8 = B_MALLOC(bddp_h8, newSpc);
  if(newhash_32 && newhash_h8)
  {
    for(i=0; i<varp->hashSpc; i++)
    {
      newhash_32[i] = varp->hash_32[i];
      newhash_h8[i] = varp->hash_h8[i];
    }
    free(varp->hash_32);
    free(varp->hash_h8);
    varp->hash_32 = newhash_32;
    varp->hash_h8 = newhash_h8;
  }
  else
  {
    if(newhash_32) {free(newhash_32);  newhash_32=NULL;}
    if(newhash_h8) {free(newhash_h8);  newhash_h8=NULL;}
    return 1;
  }
#else
  newhash_32 = B_MALLOC(bddp_32, newSpc);
  if(newhash_32)
  {
    for(i=0; i<varp->hashSpc; i++) newhash_32[i] = varp->hash_32[i];
    free(varp->hash_32);
    varp->hash_32 = newhash_32;
  }
  else return 1; /* Not enough memory */
#endif
  varp->hashSpc = newSpc;
  
  /* Initialize new hash entry */
  for(i=oldSpc; i<newSpc; i++)
  {
    B_SET_NXP(p, varp->hash, i);
    B_SET_BDDP(*p, bddnull);
  }
  
  /* restore hash entry */
  for(i=0; i<oldSpc; i++)
  {
    np0 = 0;
    B_SET_NXP(p, varp->hash, i);
    nx = B_GET_BDDP(*p);
    while(nx != bddnull)
    {
      np = Node + nx;
      f0 = B_GET_BDDP(np->f0);
      f1 = B_GET_BDDP(np->f1);
      key = B_HASHKEY(f0, f1, newSpc);
      if(key == i) np0 = np;
      else
      {
        if(np0) B_CPY_BDDP(np0->nx, np->nx);
        else { B_SET_NXP(p, varp->hash, i); B_CPY_BDDP(*p, np->nx); }
        B_SET_NXP(p, varp->hash, key);
        B_CPY_BDDP(np->nx, *p);
        B_SET_BDDP(*p, nx);
      }
      if(np0) nx = B_GET_BDDP(np0->nx);
      else { B_SET_NXP(p, varp->hash, i); nx = B_GET_BDDP(*p); }
    }
  }
  return 0;
}

static bddp getnode(bddvar v, bddp f0, bddp f1)
/* Returns bddnull if not enough memory */
{
  /* After checking elimination rule & negative edge rule */
#ifdef MEASURE_TIME
  getnode_count++;
  clock_t t0 = clock();
#endif
  struct B_NodeTable *np, *fp;
  struct B_VarTable *varp;
  bddp ix, nx, key;
  bddp_32 *p_32;
#ifdef B_64
  bddp_h8 *p_h8;
#endif
  
  varp = &Var[v];
  if(varp->hashSpc == 0)
  /* Create hash-table */
  {
    varp->hash_32 = B_MALLOC(bddp_32, B_HASH_SPC0);
    if(!varp->hash_32) return bddnull;
#ifdef B_64
    varp->hash_h8 = B_MALLOC(bddp_h8, B_HASH_SPC0);
    if(!varp->hash_h8)
    {
      free(varp->hash_32);
      varp->hash_32=NULL;
      return bddnull;
    }
#endif
    for(ix=0; ix<B_HASH_SPC0; ix++)
    {
      B_SET_NXP(p, varp->hash, ix);
      B_SET_BDDP(*p, bddnull);
    }
    varp->hashSpc = B_HASH_SPC0;
    key = B_HASHKEY(f0, f1, varp->hashSpc);
  }
  else
  /* Looking for equivalent existing node */
  {
    key = B_HASHKEY(f0, f1, varp->hashSpc);
    B_SET_NXP(p, varp->hash, key);
    nx = B_GET_BDDP(*p);
    while(nx != bddnull)
    {
      np = Node + nx;
      if(f0 == B_GET_BDDP(np->f0) &&
         f1 == B_GET_BDDP(np->f1) )
      {
        /* Sharing equivalent node */
        //          added by Lee
        if(!(B_CST(f0)||IN_DENSE(f0))) { fp = B_NP(f0); B_RFC_DEC_NP(fp); }
        if(!(B_CST(f1)||IN_DENSE(f1))) { fp = B_NP(f1); B_RFC_DEC_NP(fp); }
        //          end of addition
        //        if(!B_CST(f0)) { fp = B_NP(f0); B_RFC_DEC_NP(fp); }
        //        if(!B_CST(f1)) { fp = B_NP(f1); B_RFC_DEC_NP(fp); }
        B_RFC_INC_NP(np);
        return B_BDDP_NP(np);
      }
      nx = B_GET_BDDP(np->nx);
    }
  }
  
  //  added by Lee
  if ((B_CST(f0)||IN_DENSE(f0))&&(B_CST(f1)||IN_DENSE(f1))) {
    bddp f = getdnode(v, f0, f1);
    if (f!=bddnull) {
#ifdef MEASURE_TIME
      getnode_time += clock()-t0;
#endif
      return f;
    }
  }
  
  /* Check hash-table overflow */
  if(++ varp->hashUsed >= varp->hashSpc)
  {
    if(hash_enlarge(v)) return bddnull; /* Hash-table overflow */
    key = B_HASHKEY(f0, f1, varp->hashSpc); /* Enlarge success */
  }
  
  /* Check node-table overflow */
  if(NodeUsed >= NodeSpc-1U)
  {
    if(node_enlarge())
    {
      if(bddgc()) return bddnull; /* Node-table overflow */
      key = B_HASHKEY(f0, f1, varp->hashSpc);
    }
    /* Node-table enlarged or GC succeeded */
  }
  NodeUsed++;
  
  /* Creating a new node */
  nx = Avail;
  np = Node + nx;
  Avail = B_GET_BDDP(np->nx);
  B_SET_NXP(p, varp->hash, key);
  B_CPY_BDDP(np->nx, *p);
  B_SET_BDDP(*p, nx);
  B_SET_BDDP(np->f0, f0);
  B_SET_BDDP(np->f1, f1);
  np->varrfc = v;
  B_RFC_INC_NP(np);
#ifdef MEASURE_TIME
  getnode_time += clock() - t0;
#endif
  return B_BDDP_NP(np);
}

//static long long int getdnode_count = 0;
//static long long int miss_count = 0;

static bddp getdnode(bddvar v, bddp f0, bddp f1)
/* Returns bddnull if cannot find such a node, added by Lee */
{
  if (DB) {
    bddp f0p = (B_CST(f0))?0:f0>>1;
    i64 next = DB->M->succ(DB->M, f0p+1, 1);
    if (!(next==-1||next==DB->M->n)) {
      i64 rd = v - depth(DB->U, next);
      if (rd<0) {
        i64 w = level_ancestor(DB->U, next, rd);
        if (w!=-1) {
//          getdnode_count++;
          // wのchildに目的のノードがあるはず
          f1 = B_VAL_MASK&f1;
          i64 d = degree(DB->U, w);
          i64 left = 1;
          i64 right = d;
          i64 cp;
          bddp c1;
          {
            // 枝刈り
            cp = last_child(DB->U, w);
            if (!DB->M->getbit(DB->M, cp)) {
//              miss_count++;
              return bddnull;
            }
            c1 = B_VAL_MASK&done(bddp(cp<<1));
            if (c1>f1) {
//              miss_count++;
              return bddnull;
            }
            else if (c1==f1) {
              return bddp(cp<<1);
            }
            right--;
            if (left>right) {
//              miss_count++;
              return bddnull;
            }
            cp = first_child(DB->U, w);
            if (!DB->M->getbit(DB->M, cp)) {
              left++;
              cp = next_sibling(DB->U, cp);
            }
            if (left>right) {
//              miss_count++;
              return bddnull;
            }
            c1 = B_VAL_MASK&done(bddp(cp<<1));
            if (c1<f1) {
//              miss_count++;
              return bddnull;
            }
            else if (c1==f1) {
              return bddp(cp<<1);
            }
            left++;
          }
          
          // binary searchを高速化
          // あまり効果はなかった
          i64 center;
#ifdef B_64
          long long mask;
#else
          int mask;
#endif
          do {
            center = (left+right)>>1;
            cp = child(DB->U, w, center);
            c1 = B_VAL_MASK&done(bddp(cp<<1));
            if (c1==f1) return bddp(cp<<1);
#ifdef B_64
            mask = ((long long)(f1)-(long long)(c1))>>63;
#else
            mask = (int(f1)-int(c1))>>31;
#endif
            left = ((left^(center+1))&mask)^left;
            --center;
            right = ((center^right)&mask)^center;
          } while (left<=right);
//          miss_count++;
//          i64 center = (left+right)/2;
//          while (left<=right) {
//            cp = child(DB->U, w, center);
//            bddp cone = B_VAL_MASK&done((bddp)cp<<1);
//            if (cone>f1) {
//              left = center+1;
//            }
//            else if (cone<f1) {
//              right = center-1;
//            }
//            else {
//              return (bddp)cp<<1;
//            }
//            center = (left+right)/2;
//          }
        }
      }
    }
  }
  return bddnull;
}

static bddp getbddp(bddvar v, bddp f0, bddp f1)
/* Returns bddnull if not enough memory */
{
  struct B_NodeTable *fp;
  
  /* Check elimination rule */
  if(f0 == f1)
  {
    if(!B_CST(f0)) { fp = B_NP(f0); B_RFC_DEC_NP(fp); }
    return f0;
  }
  
  /* Negative edge constraint */
  if(B_NEG(f0))
  {
    bddp h;
    
    h = getnode(v, B_NOT(f0), B_NOT(f1));
    if(h == bddnull) return bddnull;
    return B_NOT(h);
  }
  return getnode(v, f0, f1);
}

static bddp apply(bddp f, bddp g, unsigned char op, unsigned char skip)
/* Returns bddnull if not enough memory */
{
  struct B_NodeTable *fp, *gp;
  struct B_CacheTable *cachep;
  bddp key, f0, f1, g0, g1, h0, h1, h;
  bddvar v, flev, glev;
  char z; /* flag to check ZBDD node */
  
  /* Check terminal case */
  if(!skip) switch(op)
  {
    case BC_AND:
      /* Check trivial cases */
      if(f == bddfalse || g == bddfalse || f == B_NOT(g))
        return bddfalse;
      if(f == g)
      {
        if(f != bddtrue) { fp = B_NP(f); B_RFC_INC_NP(fp); }
        return f;
      }
      if(f == bddtrue) { fp = B_NP(g); B_RFC_INC_NP(fp); return g; }
      if(g == bddtrue) { fp = B_NP(f); B_RFC_INC_NP(fp); return f; }
      /* Check operand swap */
      if(f < g) { h = f; f = g; g = h; } /* swap (f, g) */
      break;
      
    case BC_XOR:
      /* Check trivial cases */
      if(f == g) return bddfalse;
      if(f == B_NOT(g)) return bddtrue;
      if(f == bddfalse) { fp = B_NP(g); B_RFC_INC_NP(fp); return g; }
      if(g == bddfalse) { fp = B_NP(f); B_RFC_INC_NP(fp); return f; }
      if(f == bddtrue) {fp=B_NP(g); B_RFC_INC_NP(fp); return B_NOT(g);}
      if(g == bddtrue) {fp=B_NP(f); B_RFC_INC_NP(fp); return B_NOT(f);}
      /* Check negation */
      if(B_NEG(f) && B_NEG(g)) { f = B_NOT(f); g = B_NOT(g); }
      else if(B_NEG(f) || B_NEG(g))
      {
        f = B_ABS(f); g = B_ABS(g);
        /* Check operand swap */
        h = (f < g)? apply(g, f, op, 1): apply(f, g, op, 1);
        if(h == bddnull) return bddnull;
        return B_NOT(h);
      }
      /* Check operand swap */
      if(f < g) { h = f; f = g; g = h; } /* swap (f, g) */
      break;
      
    case BC_COFACTOR:
      /* Check trivial cases */
      if(B_CST(f)) return f;
      if(g == bddfalse || f == B_NOT(g)) return bddfalse;
      if(f == g) return bddtrue;
      if(g == bddtrue) { fp = B_NP(f); B_RFC_INC_NP(fp); return f; }
      break;
      
    case BC_UNIV:
      /* Check trivial cases */
      if(B_CST(f)) return f;
      if(B_CST(g)) { fp = B_NP(f); B_RFC_INC_NP(fp); return f; }
      if(B_NEG(g)) g = B_NOT(g);
      break;
      
    case BC_SUPPORT:
      if(B_CST(f)) return bddfalse;
      if(B_NEG(f)) f = B_NOT(f);
      break;
      
    case BC_INTERSEC:
      /* Check trivial cases */
      if(f == bddfalse || g == bddfalse) return bddfalse;
      if(f == bddtrue) return B_NEG(g)? bddtrue: bddfalse;
      if(g == bddtrue) return B_NEG(f)? bddtrue: bddfalse;
      //          added by Lee
      if (f==g) {
        if (!IN_DENSE(f)) {
          fp = B_NP(f);
          B_RFC_INC_NP(fp);
        }
        return f;
      }
      //          end of addition
      //    if(f == g) { fp = B_NP(f); B_RFC_INC_NP(fp); return f; }
      
      //          added by Lee
      if (f == B_NOT(g)) {
        if (!IN_DENSE(f)) {
          fp=B_NP(f);
          B_RFC_INC_NP(fp);
        }
        return B_ABS(f);
      }
      //          end of addition
      
      //    if(f == B_NOT(g)) {fp=B_NP(f); B_RFC_INC_NP(fp); return B_ABS(f); }
      
      /* Check operand swap */
      if(f < g) { h = f; f = g; g = h; } /* swap (f, g) */
      break;
      
    case BC_UNION:
      /* Check trivial cases */
      if(f == bddfalse)
      {
        //        added by Lee
        if (!B_CST(g)) {
          if (!IN_DENSE(g)) {
            fp=B_NP(g);
            B_RFC_INC_NP(fp);
          }
        }
        return g;
        //        end of addition
        
        //      if(!B_CST(g)) {fp=B_NP(g); B_RFC_INC_NP(fp); }
        //      return g;
      }
      if(f == bddtrue)
      {
        //        added by Lee
        if (!B_CST(g)) {
          if (!IN_DENSE(g)) {
            fp=B_NP(g);
            B_RFC_INC_NP(fp);
          }
        }
        return B_NEG(g)? g: B_NOT(g);
        //        end of addition
        
        //      if(!B_CST(g)) {fp=B_NP(g); B_RFC_INC_NP(fp); }
        //      return B_NEG(g)? g: B_NOT(g);
      }
      if(g == bddfalse || f == g)
        //        added by Lee
      {
        if (!IN_DENSE(f)) {
          fp = B_NP(f);
          B_RFC_INC_NP(fp);
        }
        return f;
      }
      //          end of addition
      //      { fp=B_NP(f); B_RFC_INC_NP(fp); return f; }
      
      if(g == bddtrue || f == B_NOT(g))
      {
        //         added by Lee
        if (!IN_DENSE(f)) {
          fp = B_NP(f);
          B_RFC_INC_NP(fp);
        }
        return B_NEG(f)? f: B_NOT(f);
        //        end of addition
        //      fp=B_NP(f); B_RFC_INC_NP(fp);
        //      return B_NEG(f)? f: B_NOT(f);
      }
      /* Check operand swap */
      if(f < g) { h = f; f = g; g = h; } /* swap (f, g) */
      break;
      
    case BC_SUBTRACT:
      /* Check trivial cases */
      if(f == bddfalse || f == g) return bddfalse;
      if(f == bddtrue || f == B_NOT(g))
        return B_NEG(g)? bddfalse: bddtrue;
      //          added by Lee
      if (g==bddfalse) {
        if (!IN_DENSE(f)) {
          fp = B_NP(f);
          B_RFC_INC_NP(fp);
        }
        return f;
      }
      //          end of addition
      
      //    if(g == bddfalse) { fp=B_NP(f); B_RFC_INC_NP(fp); return f; }
      
      //          added by Lee
      if (g==bddtrue) {
        if (!IN_DENSE(f)) {
          fp = B_NP(f);
          B_RFC_INC_NP(fp);
        }
        return B_ABS(f);
      }
      //          end of addition
      
      //    if(g == bddtrue) { fp=B_NP(f); B_RFC_INC_NP(fp); return B_ABS(f); }
      break;
      
    case BC_AT0:
    case BC_AT1:
    case BC_OFFSET:
      /* Check trivial cases */
      if(B_CST(f)) return f;
      /* special cases */
      //          added by Lee
      fp = B_NP(f);
      if (IN_DENSE(f)) {
        flev = dlev(f);
      }
      else {
        flev = Var[B_VAR_NP(fp)].lev;
      }
      //          end of addition
      
      //    fp = B_NP(f); flev = Var[B_VAR_NP(fp)].lev;
      glev = Var[(bddvar)g].lev;
      //          added by Lee
      if (flev<glev) {
        if (!IN_DENSE(f)) {
          B_RFC_INC_NP(fp);
        }
        return f;
      }
      //          end of addition
      //    if(flev < glev) { B_RFC_INC_NP(fp); return f; }
      if(flev == glev)
      {
        if(op != BC_AT1)
        {
          //          added by Lee
          if (IN_DENSE(f)) {
            h = dzero(f);
          }
          else {
            h = B_GET_BDDP(fp->f0);
          }
          //          end of addition
          
          //        h = B_GET_BDDP(fp->f0);
          if(B_NEG(f)^B_NEG(h)) h = B_NOT(h);
        }
        else
        {
          //          for BDD (BC_AT1)
          h = B_GET_BDDP(fp->f1);
          if(B_NEG(f)) h = B_NOT(h);
        }
        //        added by Lee
        if (!(B_CST(h)||IN_DENSE(h))) {
          fp = B_NP(h);
          B_RFC_INC_NP(fp);
        }
        //        end of addition
        //      if(!B_CST(h)) { fp = B_NP(h); B_RFC_INC_NP(fp); }
        return h;
      }
      /* Check negation */
      if(B_NEG(f))
      {
        h = apply(B_NOT(f), g, op, 1);
        if(h == bddnull) return bddnull;
        return B_NOT(h);
      }
      break;
      
    case BC_ONSET:
      /* Check trivial cases */
      if(B_CST(f)) return bddfalse;
      /* special cases */
      //          added by Lee
      fp = B_NP(f);
      if (IN_DENSE(f)) {
        flev = dlev(f);
      }
      else {
        flev = Var[B_VAR_NP(fp)].lev;
      }
      //          end of addition
      //    fp = B_NP(f); flev = Var[B_VAR_NP(fp)].lev;
      glev = Var[(bddvar)g].lev;
      if(flev < glev)  return bddfalse;
      if(flev == glev)
      {
        //        added by Lee
        if (IN_DENSE(f)) {
          h = done(f);
        }
        else {
          h = B_GET_BDDP(fp->f1);
        }
        //        end of addition
        //      h = B_GET_BDDP(fp->f1);
        
        //        added by Lee
        if (!(B_CST(h)||IN_DENSE(h))) {
          fp = B_NP(h);
          B_RFC_INC_NP(fp);
        }
        //        end of addition
        //      if(!B_CST(h)) { fp = B_NP(h); B_RFC_INC_NP(fp); }
        return h;
      }
      /* Check negation */
      if(B_NEG(f)) f = B_NOT(f); //empty set cannot be in onset(f, g)
      break;
      
    case BC_CHANGE:
      /* Check trivial cases */
      if(f == bddfalse) return f;
      if(B_CST(f)) return getzbddp((bddvar)g, bddfalse, f);
      /* special cases */
      //          added by Lee
      fp = B_NP(f);
      if (IN_DENSE(f)) {
        flev = dlev(f);
      }
      else {
        flev = Var[B_VAR_NP(fp)].lev;
      }
      //          end of addition
      //    fp = B_NP(f); flev = Var[B_VAR_NP(fp)].lev;
      glev = Var[(bddvar)g].lev;
      if(flev < glev)
      {
        //        added by Lee
        if (!IN_DENSE(f)) B_RFC_INC_NP(fp);
        //        end of addition
        //      B_RFC_INC_NP(fp);
        h = getzbddp((bddvar)g, bddfalse, f);
        if(h == bddnull) bddfree(f);
        return h;
      }
      if(flev == glev)
      {
        //        added by Lee
        if (IN_DENSE(f)) {
          h0 = done(f);
          h1 = dzero(f);
        }
        else {
          h0 = B_GET_BDDP(fp->f1);
          h1 = B_GET_BDDP(fp->f0);
        }
        //        end of addition
        //      h0 = B_GET_BDDP(fp->f1);
        //      h1 = B_GET_BDDP(fp->f0);
        
        //            added by Lee
        if(B_NEG(f)^B_NEG(h1)) h1 = B_NOT(h1);
        if(!(B_CST(h0)||IN_DENSE(h0))) { fp = B_NP(h0); B_RFC_INC_NP(fp); }
        if(!(B_CST(h1)||IN_DENSE(h1))) { fp = B_NP(h1); B_RFC_INC_NP(fp); }
        h = getzbddp((bddvar)g, h0, h1);
        if(h == bddnull) { bddfree(h0); bddfree(h1); }
        //        end of addition
        
        //      if(B_NEG(f)^B_NEG(h1)) h1 = B_NOT(h1);
        //      if(!B_CST(h0)) { fp = B_NP(h0); B_RFC_INC_NP(fp); }
        //      if(!B_CST(h1)) { fp = B_NP(h1); B_RFC_INC_NP(fp); }
        //      h = getzbddp((bddvar)g, h0, h1);
        //      if(h == bddnull) { bddfree(h0); bddfree(h1); }
        return h;
      }
      break;
      
    case BC_LSHIFT:
    case BC_RSHIFT:
      /* Check trivial cases */
      if(B_CST(f)) return f;
      
      /* Check negation */
      if(B_NEG(f))
      {
        h = apply(B_NOT(f), g, op, 1);
        if(h == bddnull) return bddnull;
        return B_NOT(h);
      }
      break;
      
    case BC_CARD:
      if(B_CST(f)) return (f == bddempty)? 0: 1;
      if(B_NEG(f)) return apply(B_NOT(f), bddempty, op, 1) + 1;
      break;
      
    case BC_LIT:
      if(B_CST(f)) return 0;
      if(B_NEG(f)) f = B_NOT(f);
      break;
      
    case BC_LEN:
      if(B_CST(f)) return 0;
      if(B_NEG(f)) f = B_NOT(f);
      break;
      
    default:
      err("apply: unknown opcode", op);
      break;
  }
  
  /* Non-trivial operations */
  switch(op)
  {
      /* binary operation */
    case BC_AND:
    case BC_XOR:
    case BC_COFACTOR:
    case BC_UNIV:
    case BC_INTERSEC:
    case BC_UNION:
    case BC_SUBTRACT:
      /* Try cache? */
      //          added by Lee
      if((B_CST(f) || ((!IN_DENSE(f))&&B_RFC_ONE_NP(B_NP(f)))) &&
         (B_CST(g) || ((!IN_DENSE(g))&&B_RFC_ONE_NP(B_NP(g))))) key = bddnull;
      //          end of addition
      //    if((B_CST(f) || B_RFC_ONE_NP(B_NP(f))) &&
      //       (B_CST(g) || B_RFC_ONE_NP(B_NP(g)))) key = bddnull;
      else
      {
        /* Checking Cache */
        key = B_CACHEKEY(op, f, g);
        cachep = Cache + key;
        if(cachep->op == op &&
           f == B_GET_BDDP(cachep->f) &&
           g == B_GET_BDDP(cachep->g))
        {
          /* Hit */
          h = B_GET_BDDP(cachep->h);
          //          added by Lee
          if(!B_CST(h) && !IN_DENSE(h) && h != bddnull) { fp = B_NP(h); B_RFC_INC_NP(fp); }
          //          end of addition
          
          //        if(!B_CST(h) && h != bddnull) { fp = B_NP(h); B_RFC_INC_NP(fp); }
          return h;
        }
      }
      /* Get (f0, f1) and (g0, g1)*/
      z = 0;
      fp = B_NP(f);
      //          added by Lee
      if (IN_DENSE(f)) {
        flev = dlev(f);
      }
      else {
        fp = B_NP(f);
        flev = B_CST(f)? 0: Var[B_VAR_NP(fp)].lev;
      }
      //          end of addition
      //    flev = B_CST(f)? 0: Var[B_VAR_NP(fp)].lev;
      
      gp = B_NP(g);
      //          added by Lee
      if (IN_DENSE(g)) {
        glev = dlev(g);
      }
      else {
        glev = B_CST(g)? 0: Var[B_VAR_NP(gp)].lev;
      }
      //          end of addition
      //    glev = B_CST(g)? 0: Var[B_VAR_NP(gp)].lev;
      f0 = f; f1 = f;
      g0 = g; g1 = g;
      
      if(flev <= glev)
      {
        //        added by Lee
        if (IN_DENSE(g)) {
          v = dvar(g);
          if (dzero(g)&1) {
            z = 1;
            if (flev < glev) f1 = bddfalse;
          }
          g0 = dzero(g);
          g1 = done(g);
        }
        else {
          v = B_VAR_NP(gp);
          if(B_Z_NP(gp))
          {
            z = 1;
            if(flev < glev) f1 = bddfalse;
          }
          g0 = B_GET_BDDP(gp->f0);
          g1 = B_GET_BDDP(gp->f1);
        }
        //        end of addition
        
        //      v = B_VAR_NP(gp);
        //      if(B_Z_NP(gp))
        //      {
        //        z = 1;
        //        if(flev < glev) f1 = bddfalse;
        //      }
        //      g0 = B_GET_BDDP(gp->f0);
        //      g1 = B_GET_BDDP(gp->f1);
        if(B_NEG(g)^B_NEG(g0)) g0 = B_NOT(g0);
        if(B_NEG(g) && !z) g1 = B_NOT(g1);
      }
      
      if(flev >= glev)
      {
        //        added by Lee
        if (IN_DENSE(f)) {
          v = dvar(f);
          if (dzero(f)&1) {
            z = 1;
            if (flev > glev) g1 = bddfalse;
          }
          f0 = dzero(f);
          f1 = done(f);
        }
        else {
          v = B_VAR_NP(fp);
          if(B_Z_NP(fp))
          {
            z = 1;
            if(flev > glev) g1 = bddfalse;
          }
          f0 = B_GET_BDDP(fp->f0);
          f1 = B_GET_BDDP(fp->f1);
        }
        //      v = B_VAR_NP(fp);
        //      if(B_Z_NP(fp))
        //      {
        //        z = 1;
        //        if(flev > glev) g1 = bddfalse;
        //      }
        //      f0 = B_GET_BDDP(fp->f0);
        //      f1 = B_GET_BDDP(fp->f1);
        if(B_NEG(f)^B_NEG(f0)) f0 = B_NOT(f0);
        if(B_NEG(f) && !z) f1 = B_NOT(f1);
      }
      break;
      
      /* unary operation */
    case BC_AT0:
    case BC_AT1:
    case BC_LSHIFT:
    case BC_RSHIFT:
    case BC_SUPPORT:
    case BC_OFFSET:
    case BC_ONSET:
    case BC_CHANGE:
      fp = B_NP(f);
      //          added by Lee
      if(!IN_DENSE(f) && B_RFC_ONE_NP(fp)) key = bddnull;
      //          end of addition
      //    if(B_RFC_ONE_NP(fp)) key = bddnull;
      else
      {
        /* Checking Cache */
        key = B_CACHEKEY(op, f, g);
        cachep = Cache + key;
        if(cachep->op == op &&
           f == B_GET_BDDP(cachep->f) &&
           g == B_GET_BDDP(cachep->g))
        {
          /* Hit */
          h = B_GET_BDDP(cachep->h);
          //          added by Lee
          if(!B_CST(h) && !IN_DENSE(h) && h != bddnull) { fp = B_NP(h); B_RFC_INC_NP(fp); }
          //          end of addition
          //        if(!B_CST(h) && h != bddnull) { fp = B_NP(h); B_RFC_INC_NP(fp); }
          return h;
        }
      }
      /* Get (f0, f1)*/
      //          added by Lee
      if (IN_DENSE(f)) {
        v = dvar(f);
        f0 = dzero(f);
        f1 = done(f);
        z = (f0&1)?1:0;
      }
      else {
        v = B_VAR_NP(fp);
        z = B_Z_NP(fp)? 1: 0;
        f0 = B_GET_BDDP(fp->f0);
        f1 = B_GET_BDDP(fp->f1);
      }
      //          end of addition
      //    v = B_VAR_NP(fp);
      //    z = B_Z_NP(fp)? 1: 0;
      //    f0 = B_GET_BDDP(fp->f0);
      //    f1 = B_GET_BDDP(fp->f1);
      if(B_NEG(f)^B_NEG(f0)) f0 = B_NOT(f0);
      if(B_NEG(f) && !z) f1 = B_NOT(f1);
      break;
      
    case BC_CARD:
    case BC_LIT:
    case BC_LEN:
      fp = B_NP(f);
      //          added by Lee
      if (!IN_DENSE(f) && B_RFC_ONE_NP(fp)) key = bddnull;
      //          end of addition
      //    if(B_RFC_ONE_NP(fp)) key = bddnull;
      else
      {
        /* Checking Cache */
        key = B_CACHEKEY(op, f, bddempty);
        cachep = Cache + key;
        if(cachep->op == op &&
           f == B_GET_BDDP(cachep->f) &&
           bddempty == B_GET_BDDP(cachep->g))
        {
          /* Hit */
          return B_GET_BDDP(cachep->h);
        }
      }
      /* Get (f0, f1)*/
      //          added by Lee
      if (IN_DENSE(f)) {
        f0 = dzero(f);
        f1 = done(f);
      }
      else {
        f0 = B_GET_BDDP(fp->f0);
        f1 = B_GET_BDDP(fp->f1);
      }
      //    f0 = B_GET_BDDP(fp->f0);
      //    f1 = B_GET_BDDP(fp->f1);
      if(B_NEG(f)^B_NEG(f0)) f0 = B_NOT(f0);
      break;
      
    default:
      err("apply: unknown opcode", op);
  }
  
  /* Stack overflow limitter */
  BDD_RECUR_INC;
  
  /* Get result node */
  switch(op)
  {
    case BC_AND:
    case BC_XOR:
    case BC_INTERSEC:
    case BC_UNION:
    case BC_SUBTRACT:
      h0 = apply(f0, g0, op, 0);
      if(h0 == bddnull) { h = h0; break; } /* Overflow */
      h1 = apply(f1, g1, op, 0);
      if(h1 == bddnull) { bddfree(h0); h = h1; break; } /* Overflow */
      h = z? getzbddp(v, h0, h1): getbddp(v, h0, h1);
      if(h == bddnull) { bddfree(h0); bddfree(h1); } /* Overflow */
      break;
      
    case BC_COFACTOR:
      if(g0 == bddfalse && g1 != bddfalse)
      {
        h = apply(f1, g1, op, 0);
      }
      else if(g1 == bddfalse && g0 != bddfalse)
      {
        h = apply(f0, g0, op, 0);
      }
      else
      {
        h0 = apply(f0, g0, op, 0);
        if(h0 == bddnull) { h = h0; break; } /* Overflow */
        h1 = apply(f1, g1, op, 0);
        if(h1 == bddnull) { bddfree(h0); h = h1; break; } /* Overflow */
        h = getbddp(v, h0, h1);
        if(h == bddnull) { bddfree(h0); bddfree(h1); } /* Overflow */
      }
      break;
      
    case BC_UNIV:
      if(g0 != g1)
      {
        h0 = apply(f0, g0, op, 0);
        if(h0 == bddnull) { h = h0; break; } /* Overflow */
        h1 = apply(f1, g0, op, 0);
        if(h1 == bddnull) { bddfree(h0); h = h1; break; } /* Overflow */
        h = apply(h0, h1, BC_AND, 0);
        bddfree(h0); bddfree(h1);
      }
      else
      {
        h0 = apply(f0, g0, op, 0);
        if(h0 == bddnull) { h = h0; break; } /* Overflow */
        h1 = apply(f1, g0, op, 0);
        if(h1 == bddnull) { bddfree(h0); h = h1; break; } /* Overflow */
        h = getbddp(v, h0, h1);
        if(h == bddnull) { bddfree(h0); bddfree(h1); } /* Overflow */
      }
      break;
      
    case BC_AT0:
    case BC_AT1:
    case BC_OFFSET:
    case BC_ONSET:
    case BC_CHANGE:
      h0 = apply(f0, g, op, 0);
      if(h0 == bddnull) { h = h0; break; } /* Overflow */
      h1 = apply(f1, g, op, 0);
      if(h1 == bddnull) { bddfree(h0); h = h1; break; } /* Overflow */
      h = z? getzbddp(v, h0, h1): getbddp(v, h0, h1);
      if(h == bddnull) { bddfree(h0); bddfree(h1); } /* Overflow */
      break;
      
    case BC_SUPPORT:
      h0 = apply(f0, bddfalse, op, 0);
      if(h0 == bddnull) { h = h0; break; } /* Overflow */
      h1 = apply(f1, bddfalse, op, 0);
      if(h1 == bddnull) { bddfree(h0); h = h1; break; } /* Overflow */
      h = z? apply(h0, h1, BC_UNION, 0):
      apply(B_NOT(h0), B_NOT(h1), BC_AND, 0);
      bddfree(h0); bddfree(h1);
      if(h == bddnull) break; /* Overflow */
      h0 = h;
      h = z? getzbddp(v, h0, bddtrue):
      getbddp(v, B_NOT(h0), bddtrue);
      if(h == bddnull) bddfree(h0); /* Overflow */
      break;
      
    case BC_LSHIFT:
    case BC_RSHIFT:
      /* Get VarID of new level */
    {
      bddvar flev, newlev;
      
      flev = bddlevofvar(v);
      if(op == BC_LSHIFT)
      {
        newlev = flev + (bddvar)g;
        if(newlev > VarUsed || newlev < flev)
          err("apply: Invald shift", newlev);
      }
      else
      {
        newlev = flev - (bddvar)g;
        if(newlev == 0 || newlev > flev)
          err("apply: Invald shift", newlev);
      }
      v = bddvaroflev(newlev);
    }
      h0 = apply(f0, g, op, 0);
      if(h0 == bddnull) { h = h0; break; } /* Overflow */
      h1 = apply(f1, g, op, 0);
      if(h1 == bddnull) { bddfree(h0); h = h1; break; } /* Overflow */
      h = z? getzbddp(v, h0, h1): getbddp(v, h0, h1);
      if(h == bddnull) { bddfree(h0); bddfree(h1); } /* Overflow */
      break;
      
    case BC_CARD:
      h = apply(f0, bddempty, op, 0)
      + apply(f1, bddempty, op, 0);
      if(h >= bddnull) h = bddnull - 1;
      break;
      
    case BC_LIT:
      h = apply(f0, bddempty, op, 0)
      + apply(f1, bddempty, op, 0);
      if(h >= bddnull) h = bddnull - 1;
      h += apply(f1, bddempty, BC_CARD, 0);
      if(h >= bddnull) h = bddnull - 1;
      break;
      
    case BC_LEN:
      h0 = apply(f0, bddempty, op, 0);
      h1 = apply(f1, bddempty, op, 0) + 1;
      h = (h0 < h1)? h1: h0;
      break;
      
    default:
      err("apply: unknown opcode", op);
      break;
  }
  
  /* Stack overflow limitter */
  BDD_RECUR_DEC;
  
  /* Saving to Cache */
  if(key != bddnull && h != bddnull)
  {
    cachep = Cache + key;
    cachep->op = op;
    B_SET_BDDP(cachep->f, f);
    B_SET_BDDP(cachep->g, g);
    B_SET_BDDP(cachep->h, h);
    if(h == f) switch(op)
    {
      case BC_AT0:
        key = B_CACHEKEY(BC_AT1, f, g);
        cachep = Cache + key;
        cachep->op = BC_AT1;
        B_SET_BDDP(cachep->f, f);
        B_SET_BDDP(cachep->g, g);
        B_SET_BDDP(cachep->h, h);
        break;
      case BC_AT1:
        key = B_CACHEKEY(BC_AT0, f, g);
        cachep = Cache + key;
        cachep->op = BC_AT0;
        B_SET_BDDP(cachep->f, f);
        B_SET_BDDP(cachep->g, g);
        B_SET_BDDP(cachep->h, h);
        break;
      case BC_OFFSET:
        key = B_CACHEKEY(BC_ONSET, f, g);
        cachep = Cache + key;
        cachep->op = BC_ONSET;
        B_SET_BDDP(cachep->f, f);
        B_SET_BDDP(cachep->g, g);
        B_SET_BDDP(cachep->h, bddfalse);
        break;
      default:
        break;
    }
    if(h == bddfalse && op == BC_ONSET)
    {
      key = B_CACHEKEY(BC_OFFSET, f, g);
      cachep = Cache + key;
      cachep->op = BC_OFFSET;
      B_SET_BDDP(cachep->f, f);
      B_SET_BDDP(cachep->g, g);
      B_SET_BDDP(cachep->h, f);
    }
  }
  return h;
}

static void gc1(struct B_NodeTable *np)
{
  /* np is a node ptr to be collected. (refc == 0) */
  bddp key, nx1, f0, f1;
  struct B_VarTable *varp;
  struct B_NodeTable *np1, *np2;
  bddp_32 *p_32;
#ifdef B_64
  bddp_h8 *p_h8;
#endif
  
  /* remove the node from hash list */
  varp = Var + B_VAR_NP(np);
  f0 = B_GET_BDDP(np->f0);
  f1 = B_GET_BDDP(np->f1);
  key = B_HASHKEY(f0, f1, varp->hashSpc);
  B_SET_NXP(p, varp->hash, key);
  nx1 = B_GET_BDDP(*p);
  np1 = Node + nx1;
  
  if(np1 == np) B_CPY_BDDP(*p, np->nx);
  else
  {
    while(np1 != np)
    {
      if(nx1 == bddnull) err("gc1: Fail to find the node to be deleted", np-Node);
      np2 = np1;
      nx1 = B_GET_BDDP(np2->nx);
      np1 = Node + nx1;
    }
    B_CPY_BDDP(np2->nx, np->nx);
  }
  varp->hashUsed--;
  
  /* append the node to avail list */
  B_SET_BDDP(np->nx, Avail);
  Avail = np - Node;
  
  NodeUsed--;
  np->varrfc = 0;
  
  /* Check sub-graphs recursively */
  if(!(B_CST(f0)||IN_DENSE(f0)))
  {
    np1 = B_NP(f0);
    B_RFC_DEC_NP(np1);
    if(B_RFC_ZERO_NP(np1))
    {  BDD_RECUR_INC; gc1(np1); BDD_RECUR_DEC; }
  }
  if(!(B_CST(f1)||IN_DENSE(f1)))
  {
    np1 = B_NP(f1);
    B_RFC_DEC_NP(np1);
    if(B_RFC_ZERO_NP(np1))
    {  BDD_RECUR_INC; gc1(np1); BDD_RECUR_DEC; }
  }
}

static bddp count(bddp f)
{
  bddp nx;
  bddp c, g;
  bddvar flev, glev;
  struct B_NodeTable *fp;
  struct B_NodeTable *gp;
  
  /* Check consistensy
   if(f == bddnull)
   err("count: bddnull found", bddnull);
   */
  
  if(B_CST(f)) return 0; /* Constant */
  fp = B_NP(f);
  
  /* Check visit flag */
  //    added by Lee
  if (IN_DENSE(f)) {
    if (dvisit(f)) return 0;
  }
  else {
    nx = B_GET_BDDP(fp->nx);
    if(nx & B_CST_MASK) return 0;
  }
  //    end of addition
  //  nx = B_GET_BDDP(fp->nx);
  //  if(nx & B_CST_MASK) return 0;
  
  /* Check consistensy
   flev = Var[B_VAR_NP(fp)].lev;
   g = B_GET_BDDP(fp->f0);
   if(!B_CST(g))
   {
   gp = B_NP(g); glev = Var[B_VAR_NP(gp)].lev;
   if(flev <= glev)
   err("count: inconsistensy found at f0", fp-Node);
   }
   g = B_GET_BDDP(fp->f1);
   if(!B_CST(g))
   {
   gp = B_NP(g); glev = Var[B_VAR_NP(gp)].lev;
   if(flev <= glev)
   err("count: inconsistensy found at f1", fp-Node);
   }
   */
  
  BDD_RECUR_INC;
  //    added by Lee
  if (IN_DENSE(f)) {
    c = count(dzero(f)) + count(done(f)) + 1U;
  }
  else {
    c = count(B_GET_BDDP(fp->f0)) + count(B_GET_BDDP(fp->f1)) + 1U ;
  }
  //    end of addition
  //  c = count(B_GET_BDDP(fp->f0)) + count(B_GET_BDDP(fp->f1)) + 1U ;
  BDD_RECUR_DEC;
  
  /* Set visit flag */
  //    added by Lee
  if (IN_DENSE(f)) {
    setdvisit(f, 1);
  }
  else {
    B_SET_BDDP(fp->nx, nx | B_CST_MASK);
  }
  //    end of addition
  //  B_SET_BDDP(fp->nx, nx | B_CST_MASK);
  
  return c;
}

static void dump(bddp f)
{
  bddp nx, f0, f1;
  bddvar v;
  struct B_NodeTable *fp;
  
  if(B_CST(f)) return; /* Constant */
  fp = B_NP(f);
  
  /* Check visit flag */
  //    added by Lee
  nx = B_GET_BDDP(fp->nx);
  if (IN_DENSE(f)) {
    if (dvisit(f)) return;
  }
  else {
    if(nx & B_CST_MASK) return;
  }
  //    end of addition
  //  nx = B_GET_BDDP(fp->nx);
  //  if(nx & B_CST_MASK) return;
  
  /* Set visit flag */
  //    added by Lee
  if (IN_DENSE(f)) {
    setdvisit(f, 1);
  }
  else {
    B_SET_BDDP(fp->nx, nx | B_CST_MASK);
  }
  //    end of addition
  //  B_SET_BDDP(fp->nx, nx | B_CST_MASK);
  
  /* Dump its subgraphs recursively */
  //    added by Lee
  if (IN_DENSE(f)) {
    v = dvar(f);
    f0 = dzero(f);
    f1 = done(f);
  }
  else {
    v = B_VAR_NP(fp);
    f0 = B_GET_BDDP(fp->f0);
    f0 = B_ABS(f0);
    f1 = B_GET_BDDP(fp->f1);
  }
  //    end of addition
  //  v = B_VAR_NP(fp);
  //  f0 = B_GET_BDDP(fp->f0);
  //  f0 = B_ABS(f0);
  //  f1 = B_GET_BDDP(fp->f1);
  BDD_RECUR_INC;
  dump(f0);
  dump(f1);
  BDD_RECUR_DEC;
  
  /* Dump this node */
  printf("N");
  //    added by Lee
  printf(B_BDDP_FD, B_ABS(f)>>1);
  //    end of addition
  //  printf(B_BDDP_FD, B_NDX(f));
  printf(" = [V%d(%d), ", v, Var[v].lev);
  if(B_CST(f0)) printf(B_BDDP_FD, B_VAL(f0));
  //    added by Lee
  else { printf("N"); printf(B_BDDP_FD, B_ABS(f0)>>1); }
  //    end of addition
  //  else { printf("N"); printf(B_BDDP_FD, B_NDX(f0)); }
  printf(", ");
  if(B_NEG(f1)) putchar('~');
  if(B_CST(f1)) printf(B_BDDP_FD, B_ABS(B_VAL(f1)));
  //    added by Lee
  else { printf("N"); printf(B_BDDP_FD, B_ABS(f1)>>1); }
  //    end of addition
  //  else { printf("N"); printf(B_BDDP_FD, B_NDX(f1)); }
  printf("]");
  //    added by Lee
  if (IN_DENSE(f)) {
    if (dzero(f)&1) printf(" #Z");
  }
  else {
    if(B_Z_NP(fp)) printf(" #Z");
  }
  //    end of addition
  //  if(B_Z_NP(fp)) printf(" #Z");
  printf("\n");
}

static void reset(bddp f)
{
  bddp nx;
  struct B_NodeTable *fp;
  
  if(B_CST(f)) return; /* Constant */
  fp = B_NP(f);
  
  /* Check visit flag */
  //    added by Lee
  if (IN_DENSE(f)) {
    if (dvisit(f)) {
      setdvisit(f, 0);
      BDD_RECUR_INC;
      reset(dzero(f));
      reset(done(f));
      BDD_RECUR_DEC;
    }
  }
  else {
    nx = B_GET_BDDP(fp->nx);
    if(nx & B_CST_MASK)
    {
      /* Reset visit flag */
      B_SET_BDDP(fp->nx, nx & ~B_CST_MASK);
      BDD_RECUR_INC;
      reset(B_GET_BDDP(fp->f0));
      reset(B_GET_BDDP(fp->f1));
      BDD_RECUR_DEC;
    }
  }
  //  nx = B_GET_BDDP(fp->nx);
  //  if(nx & B_CST_MASK)
  //  {
  //    /* Reset visit flag */
  //    B_SET_BDDP(fp->nx, nx & ~B_CST_MASK);
  //    BDD_RECUR_INC;
  //    reset(B_GET_BDDP(fp->f0));
  //    reset(B_GET_BDDP(fp->f1));
  //    BDD_RECUR_DEC;
  //  }
}

static bddp getzbddp(bddvar v, bddp f0, bddp f1)
/* Returns bddnull if not enough memory */
{
  struct B_NodeTable *fp;
  
  /* Check elimination rule */
  if(f1 == bddfalse) return f0;
  
  /* Negative edge constraint */
  if(B_NEG(f0))
  {
    bddp h;
    
    h = getnode(v, f0, f1);
    if(h == bddnull) return bddnull;
    return B_NOT(h);
  }
  return getnode(v, B_NOT(f0), f1);
}

static bddp dzero(bddp f)
/* return zero-child node of f in dense */
{
  bddp g = (bddp)DB->M->pred(DB->M, parent(DB->U, f>>1), 1);
  if (g==-1) g = 0;
  //    bddp g = (bddp)DB->M->rank(DB->M, parent(DB->U, f>>1), 1);
  return (g==0)?bddfalse|1:(g<<1)|1; //for ZBDD node
  //    return (g==0)?bddfalse:(g<<1); //for BDD node
}

static bddp done(bddp f)
/* return one-child node of f in dense */
{
  f = (bddp)DB->M->rank(DB->M, f>>1, 1);
  bddp g = (bddp)DB->I->getbits(DB->I, DB->width*(f-1), (int)DB->width);
  return (g==0)?bddfalse:(g==1)?bddtrue:g;
}

static bddp done_r(bddp f)
/* return one-child node of f in dense */
{
  bddp g = (bddp)DB->I->getbits(DB->I, DB->width*(f-1), (int)DB->width);
  return (g==0)?bddfalse:(g==1)?bddtrue:g;
}

static bddvar dvar(bddp f)
{
  return VarID[dlev(f)];
}

static bddvar dlev(bddp f)
{
  return (bddvar)depth(DB->U, f>>1)-1;
}

static bddp dvisit(bddp f)
{
  f = (bddp)DB->M->rank(DB->M, f>>1, 1);
  return (bddp)DB->visit->getbit(DB->visit, f-1);
}

static void setdvisit(bddp f, int flag)
{
  f = (bddp)DB->M->rank(DB->M, f>>1, 1);
  DB->visit->setbit(DB->visit, f-1, flag);
}

static int err(const char msg[], bddp num)
{
  fprintf(stderr,"***** ERROR  %s ( ", msg);
  fprintf(stderr, B_BDDP_FX, num);
  fprintf(stderr," ) *****\n");
  fprintf(stderr," NodeLimit : ");
  fprintf(stderr, B_BDDP_FD, NodeLimit);
  fprintf(stderr,"\t NodeSpc : ");
  fprintf(stderr, B_BDDP_FD, NodeSpc);
  fprintf(stderr,"\t VarSpc : %d",VarSpc);
  fprintf(stderr,"\n CacheSpc : ");
  fprintf(stderr, B_BDDP_FD, CacheSpc);
  fprintf(stderr,"\t NodeUsed : ");
  fprintf(stderr, B_BDDP_FD, NodeUsed);
  fprintf(stderr,"\t VarUsed : %d\n",VarUsed);
  exit(1);
  return 1;
}

static int rfc_inc_ovf(struct B_NodeTable *np)
{
  bddp ix, nx, nx2, key, rfc, oldSpc;
  bddp *p, *p2;
  struct B_RFC_Table *oldRFCT;
  
  /* printf("rfc_inc %d (u:%d)\n", np-Node, RFCT_Used); */
  if(RFCT_Spc == 0)
  {
    /* Create RFC-table */
    RFCT = B_MALLOC(struct B_RFC_Table, B_RFCT_SPC0);
    if(!RFCT)
    {
      err("B_RFC_INC_NP: rfc memory over flow", np-Node);
      return 1;
    }
    for(ix=0; ix<B_RFCT_SPC0; ix++)
    {
      B_SET_BDDP((RFCT+ix)->nx, bddnull);
      B_SET_BDDP((RFCT+ix)->rfc, (bddp)0);
    }
    RFCT_Spc = B_RFCT_SPC0;
  }
  
  nx = np - Node;
  key = nx & (RFCT_Spc-1);
  nx2 = B_GET_BDDP((RFCT+key)->nx);
  while(nx2 != bddnull)
  {
    if(nx == nx2)
    {
      if(np->varrfc < B_RFC_MASK)
      {
        rfc = 0;
        np->varrfc += B_RFC_UNIT;
      }
      else rfc = B_GET_BDDP((RFCT+key)->rfc) + 1;
      B_SET_BDDP((RFCT+key)->rfc, rfc);
      return 0;
    }
    key = (key+1) & (RFCT_Spc-1);
    nx2 = B_GET_BDDP((RFCT+key)->nx);
  }
  
  /* new rfc entry */
  B_SET_BDDP((RFCT+key)->nx, nx);
  B_SET_BDDP((RFCT+key)->rfc, (bddp)0);
  np->varrfc += B_RFC_UNIT;
  RFCT_Used++;
  
  if((RFCT_Used << 1) >= RFCT_Spc)
  {
    /* Enlarge RFC-table */
    oldSpc = RFCT_Spc;
    RFCT_Spc <<= 2;
    
    oldRFCT = RFCT;
    RFCT = B_MALLOC(struct B_RFC_Table, RFCT_Spc);
    if(!RFCT)
    {
      err("B_RFC_INC_NP: rfc memory over flow", np-Node);
      return 1;
    }
    for(ix=0; ix<RFCT_Spc; ix++)
    {
      B_SET_BDDP((RFCT+ix)->nx, bddnull);
      B_SET_BDDP((RFCT+ix)->rfc, (bddp)0);
    }
    for(ix=0; ix<oldSpc; ix++)
    {
      nx = B_GET_BDDP((oldRFCT+ix)->nx);
      if(nx == bddnull) continue;
      key = nx & (RFCT_Spc-1);
      nx2 = B_GET_BDDP((RFCT+key)->nx);
      while(nx2 != bddnull)
      {
        key = (key+1) & (RFCT_Spc-1);
        nx2 = B_GET_BDDP((RFCT+key)->nx);
      }
      B_SET_BDDP((RFCT+key)->nx, nx);
      rfc = B_GET_BDDP((oldRFCT+ix)->rfc);
      B_SET_BDDP((RFCT+key)->rfc, rfc);
    }
    free(oldRFCT);
  }
  
  return 0;
}

static int rfc_dec_ovf(B_NodeTable *np)
{
  bddp nx, key, nx2, rfc;
  
  /* printf("rfc_dec %d (u:%d)\n", np-Node, RFCT_Used); */
  nx = np - Node;
  key = nx & (RFCT_Spc-1);
  nx2 = B_GET_BDDP((RFCT+key)->nx);
  while(nx2 != bddnull)
  {
    if(nx == nx2)
    {
      rfc = B_GET_BDDP((RFCT+key)->rfc);
      if(rfc == 0)
      {
        np->varrfc -= B_RFC_UNIT;
        return 0;
      }
      B_SET_BDDP((RFCT+key)->rfc, rfc-1);
      return 0;
    }
    key = (key+1) & (RFCT_Spc-1);
    nx2 = B_GET_BDDP((RFCT+key)->nx);
  }
  return 0;
}

//rootを含める
static bddp sizeofdense()
{
  if (DB) {
    return (bddp)DB->M->rank(DB->M, DB->M->n, 1) + 1;
  }
  return 1;
}

///////////for compress///////////
typedef struct PrankRange {
  vector<bddp> znodes;
  bddp dzleft;
  bddp dzright;
  bddp from;
  bddp prev;
  bool firstrange;
} PrankRange;
typedef vector<PrankRange *> PrankRanges;
bddp *dzmemo;
bitvector *junctionflags;
bitvector *dzrefflags;
bitvector *zrefflags;
bitvector *order;
bddp *zmemo1;
bddp *zmemo2;
bddp *stmemo;
PrankRanges *PrankList;
//////////////////////////////////
static void countDZJunc(bddp f, bddp parent, int edge);
static void trans(bddp f, bddp parent, int edge);
static bool isJunction(bddp f);
static void setJunction(bddp f);
static bddp calcSubtreeSize(bddp f, bddp root);
static bddp subtreeSize(bddp f);
class CompareRevZeros;
class CompareSameLevels;
class ComparePranks;
static void calcPrerank();
static void calcPrerankFrom(bddp f, bddp *rank);
static void setPrerank(bddp f, bddp rank); //fにrankを割り当てる(ZDDのノードのみに用いる)
static bddp prerank(bddp f);
static void setLastDZ(bddp f, bddp dz, bool isParent);
static bddp lastDZ(bddp f);
static void convertToDense(bddp f, bool hasznode, FBitIO &uio, FBitIO &mio, FBitIO &iio);
static bool canCopy(bddp f);
static bool hasZnode(bddp f);
static void getnewdense(bddp *roots, FBitIO &uio, FBitIO &mio, FBitIO &iio, unsigned long lenU, unsigned long lenM, unsigned long lenI, int opt);

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
  void compress(bddp *roots, int opt) {
    clock_t t0, t1;
    if (DB) {
      i64 n = sizeofdense();
      i64 m = NodeUsed;
//      if (m<Z_LIMIT) {
//        return;
//      }
      float th = (float)m/(float)n;
      if (th<THETA) {
        return;
      }
    }
    bool flag = true;
    t0 = clock();
    for (int i=0; roots[i]!=bddnull; i++) {
      if (!(B_CST(roots[i]))) {
        if (flag) {
          free(Cache);
          Cache = NULL;
          free(RFCT);
          RFCT = NULL;
          bddp m = sizeofdense();
          junctionflags = bitvector_new(m+NodeSpc);
          dzrefflags = bitvector_new(m);
          zrefflags = bitvector_new(NodeSpc);
          flag = false;
        }
        countDZJunc(roots[i], bddnull, 1);
      }
    }
    if (flag) {
      return;
    }
    bitvector_makeindex(junctionflags, 128, SDARRAY_RANK1);
    i64 dzjuncnum = junctionflags->rank(junctionflags, junctionflags->n, 1);
    dzmemo = new bddp[dzjuncnum];
    for (int i=0; i<dzjuncnum; i++) {
      dzmemo[i] = bddnull;
    }
    bitvector_free(zrefflags);
    zrefflags = bitvector_new(NodeSpc);
    for (int i=0; roots[i]!=bddnull; i++) {
      if (!(B_CST(roots[i]))) {
        trans(roots[i], bddnull, 1);
      }
    }
    t1 = clock();
    trans_time += t1 - t0;
    compress_count++;
    //    bitvector_makeindex(junctionflags, 512, SDARRAY_RANK1|SDARRAY_SELECT1);
    //    bitvector_makeindex(dzrefflags, 512, SDARRAY_RANK1|SDARRAY_SELECT1);
    //    bitvector_makeindex(zrefflags, 512, SDARRAY_RANK1/*|SDARRAY_SELECT1*/);
    bitvector_makeindex(junctionflags, 128, SDARRAY_RANK1|SDARRAY_SELECT1); // sada 4/14
    bitvector_makeindex(dzrefflags, 128, SDARRAY_RANK1|SDARRAY_SELECT1);
    bitvector_makeindex(zrefflags, 128, SDARRAY_RANK1/*|SDARRAY_SELECT1*/);
    bddp dzrefnum = (bddp)dzrefflags->rank(dzrefflags, dzrefflags->n, 1);
    bddp zrefnum = (bddp)zrefflags->rank(zrefflags, zrefflags->n, 1);
    order = bitvector_new(dzrefnum+zrefnum);
    
    
    zmemo1 = new bddp[zrefnum];
    zmemo2 = new bddp[zrefnum];
    ///////////////////////////////////////////////////////////
    // Algorithm 1: How to make stmemo
    ///////////////////////////////////////////////////////////
    stmemo = new bddp[dzjuncnum];
    bddp sum = 0;
    for (int i=0; i<dzjuncnum; i++) {
      for (bddp rzero=dzmemo[i]; rzero!=bddnull; rzero=B_GET_BDDP((B_NP(rzero))->nx)) {
        sum += calcSubtreeSize(rzero, rzero);
      }
      stmemo[i] = sum;
    }
    t0 = clock();
    calcPrerank();
    t1 = clock();
    calcPrerank_time += t1 - t0;

    delete [] zmemo2;
    delete [] stmemo;
    //    bitvector_makeindex(order, 512, SDARRAY_SELECT0);
    bitvector_makeindex(order, 128, SDARRAY_SELECT0); // sada 4/14
    bddp width = 1;
    while ((1<<width)-1<order->n-1) {
      width++;
    }
    
    for (int i=0; roots[i]!=bddnull; i++) {
      if (!B_CST(roots[i])) {
        bddp r = prerank(roots[i]);
        roots[i] = (roots[i]&1)?(r<<1)|1:r<<1;
      }
    }
    
    FBitIO uio("./u.densebdd");
    FBitIO mio("./m.densebdd");
    FBitIO iio("./i.densebdd");
    
    uio.start_writing(1);
    mio.start_writing(1);
    iio.start_writing(int(width+1));
    uio.write(1);
    mio.write(0);
    
    t0 = clock();
    convertToDense(bddempty, canCopy(bddempty), uio, mio, iio);
    t1 = clock();
    convertToDense_time += t1 - t0;
    
    uio.write(0);
    mio.write(0);
    
    delete [] zmemo1;
    delete [] dzmemo;
    bitvector_free(junctionflags);
    bitvector_free(dzrefflags);
    bitvector_free(zrefflags);
    bitvector_free(order);
    
    unsigned long lenU = uio.finish_writing();
    unsigned long lenM = mio.finish_writing();
    unsigned long lenI = iio.finish_writing();
    
    uio.start_reading(sizeof(pb)*8);
    mio.start_reading(sizeof(bitvec_t)*8);
    iio.start_reading(int(width+1));
    
    t0 = clock();
    getnewdense(roots, uio, mio, iio, lenU, lenM, lenI, opt);
    t1 = clock();
    getnewdense_time += t1 - t0;
    
    uio.finish_writing();
    mio.finish_writing();
    iio.finish_writing();
  }
  
  void showresult(void) {
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage)==0) {
      cout <<"maxrss " << usage.ru_maxrss << endl;
    }
    cout << "compress count " << compress_count << endl;
//    cout << "z limit " << Z_LIMIT << endl;
    cout << "THETA " << THETA << endl;
    cout << "total nodes " << sizeofdense() + NodeUsed << endl;
    cout << "trans time " << ((double)trans_time/CLOCKS_PER_SEC) << endl;
    cout << "calcPrerank time " << ((double)calcPrerank_time/CLOCKS_PER_SEC) << endl;
    cout << "convertToDense time " << ((double)convertToDense_time/CLOCKS_PER_SEC) << endl;
    cout << "getnewdense time " << ((double)getnewdense_time/CLOCKS_PER_SEC) << endl;
//    cout << "miss rate " << ((double)miss_count/getdnode_count) << endl;
  }
  
  
#ifdef __cplusplus
}
#endif /* __cplusplus */

void densezdd_init(void)
{
  bitvector_make_selecttbl();
  darray_make_selecttbl();
  bp_make_matchtbl();
#ifdef USE_CACHE
    //block_cache_init(50);
#endif
}

static void countDZJunc(bddp f, bddp parent, int edge) {
  if (IN_DENSE(f)||B_CST(f)) {
    if (!(parent==bddnull||edge)) {
      if (isJunction(f)) {
        return;
      }
      else {
        setJunction(f);
        return;
      }
    }
  }
  else {
    bddp fidx = B_NDX(f);
    B_NodeTable *fnode = B_NP(f);
    if (!zrefflags->getbit(zrefflags, fidx)) {
      zrefflags->setbit(zrefflags, fidx, 1);
      countDZJunc(B_GET_BDDP(fnode->f0), f, 0);
      countDZJunc(B_GET_BDDP(fnode->f1), f, 1);
      return;
    }
  }
}

////////////////////////////////////////////////////////////////
// ZDDのノードをたどって使用中フラグを立てる
// junction node を記憶する
////////////////////////////////////////////////////////////////
static void trans(bddp f, bddp parent, int edge) {
  if (IN_DENSE(f)||B_CST(f)) {
    bddp fidx = (IN_DENSE(f))?(bddp)DB->M->rank(DB->M, f>>1, 1):0;
    if ((!edge)&&(!(parent==bddnull||IN_DENSE(parent)))) {
      // f and parent are junction nodes
      i64 idx = (IN_DENSE(f))?junctionflags->rank(junctionflags, fidx, 1) - 1:0;
      bddp eldest_son = dzmemo[idx];
      B_SET_BDDP(B_NP(parent)->nx, eldest_son);
      dzmemo[idx] = parent;
      if (!isJunction(parent)) {
        setJunction(parent);
      }
    }
    if (!dzrefflags->getbit(dzrefflags, fidx)) {
      dzrefflags->setbit(dzrefflags, fidx, 1);
      if (IN_DENSE(f)) {
        trans(dzero(f), f, 0);
        trans(done(f), f, 1);
      }
    }
  }
  else {
    if (!zrefflags->getbit(zrefflags, B_NDX(f))) {
      zrefflags->setbit(zrefflags, B_NDX(f), 1);
      B_SET_BDDP(B_NP(f)->nx, bddnull);
      trans(B_GET_BDDP(B_NP(f)->f0), f, 0);
      trans(B_GET_BDDP(B_NP(f)->f1), f, 1);
      B_SET_BDDP(B_NP(f)->f0, (edge)?bddnull:parent);
    }
    else {
      if (!edge) {
        B_SET_BDDP(B_NP(parent)->nx, B_GET_BDDP(B_NP(f)->f0));
        B_SET_BDDP(B_NP(f)->f0, parent);
      }
    }
  }
}

static bool isJunction(bddp f) {
  if (B_CST(f)) {
    return junctionflags->getbit(junctionflags, 0);
  }
  else if (IN_DENSE(f)) {
    return junctionflags->getbit(junctionflags, DB->M->rank(DB->M, f>>1, 1));
  }
  else {
    return junctionflags->getbit(junctionflags, B_NDX(f)+sizeofdense());
  }
}

static void setJunction(bddp f) {
  if (B_CST(f)) {
    junctionflags->setbit(junctionflags, 0, 1);
  }
  else if (IN_DENSE(f)) {
    junctionflags->setbit(junctionflags, DB->M->rank(DB->M, f>>1, 1), 1);
  }
  else {
    junctionflags->setbit(junctionflags, B_NDX(f)+sizeofdense(), 1);
  }
}

static bddp calcSubtreeSize(bddp f, bddp root) {
  bddp sum = 1;
  for (bddp rzero=B_GET_BDDP(B_NP(f)->f0); rzero!=bddnull; rzero=B_GET_BDDP(B_NP(rzero)->nx)) {
    sum += calcSubtreeSize(rzero, root);
  }
  bddp idx = (bddp)zrefflags->rank(zrefflags, B_NDX(f), 1) - 1;
  zmemo1[idx] = sum;
  zmemo2[idx] = root;
  return sum;
}

//////////////////////////////////////////////////////
// Algorithm 2: How to obtain the subtree size rooted at f
//////////////////////////////////////////////////////
static bddp subtreeSize(bddp f) {
  if (B_CST(f)||IN_DENSE(f)) {
    bddp fidx = (B_CST(f))?0:(bddp)DB->M->rank(DB->M, f>>1, 1);
    i64 fp = (B_CST(f))?0:f>>1;
    i64 fcp = find_close(DB->U, fp);
    bddp lastidx = (bddp)DB->M->rank(DB->M, fcp, 1);
    bddp from;
    if (B_CST(f)) {
      from = 0;
    }
    else {
      i64 idx = junctionflags->rank(junctionflags, fidx-1, 1);
      // idx: fの前にあるjunction nodeの数
      if (idx) {
        from = stmemo[idx-1];
      }
      else {
        from = 0;
      }
    }
    bddp to;
    {
      i64 idx = junctionflags->rank(junctionflags, lastidx, 1);
      if (idx) {
        to = stmemo[idx-1];
      }
      else {
        to = 0;
      }
    }
    return bddp(dzrefflags->rank(dzrefflags, lastidx, 1)-dzrefflags->rank(dzrefflags, fidx, 1)+1+to-from);
  }
  else {
    return zmemo1[zrefflags->rank(zrefflags, B_NDX(f), 1)-1];
  }
}

////////////////////////////////////////////////////////////////
// Algorithm 3: compareLevels
////////////////////////////////////////////////////////////////
class CompareRevZeros {
public:
  bool operator() (const bddp &zero1, const bddp &zero2) {
    return bddlevofvar(bddtop(zero1)) > bddlevofvar(bddtop(zero2));
  }
};

////////////////////////////////////////////////////////////////
// Algorithm 4: compareSameLevels
////////////////////////////////////////////////////////////////
class CompareSameLevels {
public:
  bool operator() (const bddp &zero1, const bddp &zero2) {
    //  trivial case
    if (IN_DENSE(zero1)&&IN_DENSE(zero2)) {
      return zero1<zero2;
    }
      
    //  get 1-child
    bddp one1 = (IN_DENSE(zero1))?done(zero1):B_GET_BDDP(B_NP(zero1)->f1);
    bddp one2 = (IN_DENSE(zero2))?done(zero2):B_GET_BDDP(B_NP(zero2)->f1);
    if (B_CST(one1)&&B_CST(one2)) {
      return one1>one2;
    }
    else if (B_CST(one1)) {
      return 0;
    }
    else if (B_CST(one2)) {
      return 1;
    }
    else if (IN_DENSE(one1)&&IN_DENSE(one2)) {
      return one1>one2;
    }
    else if (IN_DENSE(one1)) {
      bddp last = lastDZ(one2);
      if (B_CST(last)) {
        last = 0;
      }
      return B_ABS(one1)>B_ABS(last);
    }
    else if (IN_DENSE(one2)) {
      bddp last = lastDZ(one1);
      if (B_CST(last)) {
        last = 0;
      }
      return B_ABS(last)>=B_ABS(one2);
    }
    else {
      return ((prerank(one1)<<1)|(one1&1))>((prerank(one2)<<1)|(one2&1));
    }
  }
};

class ComparePranks {
public:
  bool operator() (const bddp &zero1, const bddp &zero2) {
    if (IN_DENSE(zero1)&&IN_DENSE(zero2)) {
      return zero1<zero2;
    }
    return prerank(zero1)<prerank(zero2);
  }
};

///////////////////////////////////////////////////////////////////////
// Algorithm 6: calcPrerank
///////////////////////////////////////////////////////////////////////
static void calcPrerank() {
  bddvar maxLev = bddvarused();
  PrankList = new PrankRanges[maxLev];
  bddp rank = 0;
  calcPrerankFrom(bddempty, &rank);
  for (int i=0; i<maxLev; i++) {
    PrankRanges ranges = PrankList[i];
    for (vector<PrankRange *>::iterator it=ranges.begin(); it!=ranges.end(); it++) {
      vector<bddp> &znodes = (*it)->znodes;
      bddp nextz = B_GET_BDDP(B_NP(znodes.back())->nx);
      sort(znodes.begin(), znodes.end(), CompareSameLevels());
      
      // リストのつなげ直し
      if ((*it)->firstrange) {
        bddp f = (*it)->prev;
        if (IN_DENSE(f)) {
          bddp fidx = (bddp)DB->M->rank(DB->M, f>>1, 1);
          i64 idx = junctionflags->rank(junctionflags, fidx, 1) - 1;
          dzmemo[idx] = znodes.front();
        }
        else if (B_CST(f)) {
          dzmemo[0] = znodes.front();
        }
        else {
          B_SET_BDDP(B_NP(f)->f0, znodes.front());
        }
      }
      else {
        B_SET_BDDP(B_NP((*it)->prev)->nx, znodes.front());
      }
      for (vector<bddp>::iterator zp=znodes.begin(); zp!=znodes.end(); zp++) {
        B_SET_BDDP(B_NP(*zp)->nx, ((zp+1)!=znodes.end())?*(zp+1):nextz);
      }
      
      bddp lastdz = bddnull;
      rank = (*it)->from;
      unsigned int idx = 0;
      bddp dz = (*it)->dzleft;
      while (idx<znodes.size()||dz!=bddnull) {
        bddp f;
        if (idx<znodes.size()&&dz!=bddnull) {
          if (CompareSameLevels()(znodes[idx], dz)) {
            // Z-nodeを進める
            f = znodes[idx++];
            if (lastdz!=bddnull) {
              setLastDZ(f, lastdz, false);
            }
          }
          else {
            // DZ-nodeを進める
            f = lastdz = dz;
            if (dz==(*it)->dzright) {
              dz = bddnull;
            }
            else {
              dz = bddp(DB->M->select(DB->M, (bddp)dzrefflags->succ(dzrefflags, DB->M->rank(DB->M, find_close(DB->U, dz>>1), 1)+1, 1), 1)<<1);
            }
          }
        }
        else if (dz!=bddnull) {
          // DZ-nodeを進める
          f = dz;
          if (dz==(*it)->dzright) {
            dz = bddnull;
          }
          else {
            dz = bddp(DB->M->select(DB->M, (bddp)dzrefflags->succ(dzrefflags, DB->M->rank(DB->M, find_close(DB->U, dz>>1), 1)+1, 1), 1)<<1);
          }
        }
        else {
          // Z-nodeを進める
          f = znodes[idx++];
          if (lastdz!=bddnull) {
            setLastDZ(f, lastdz, false);
          }
        }
        calcPrerankFrom(f, &rank);
      }
      delete *it;
    }
  }
  delete [] PrankList;
}

//////////////////////////////////////////////////////////////
// Algorithm 7: calcPrerankFrom
//////////////////////////////////////////////////////////////
static void calcPrerankFrom(bddp f, bddp *rank) {
  if ((IN_DENSE(f)||(DB&&B_CST(f)))&&(!isJunction(f))) {
    // この時は深さ優先探索を続ける
    i64 cp, fcp;
    // get cp and fcp
    {
      i64 fp = (B_CST(f))?0:f>>1; // fのBP上の位置
      bddp fidx = (B_CST(f))?0:(bddp)DB->M->rank(DB->M, fp, 1);
      fcp = find_close(DB->U, fp); // fの閉じカッコの位置
      cp = DB->M->select(DB->M, (bddp)dzrefflags->succ(dzrefflags, fidx+1, 1), 1); // 現在注目しているDZ側の子ノードのBP上の位置
    }
    
    if (!hasZnode(f)) {
      // fのsubtree内にZDDのノードがないので，
      // subtreeのサイズを足せば良い
      *rank += subtreeSize(f);
    }
    else {
      (*rank)++; // fに(*rank)を割当て
      while (cp<fcp) {
        calcPrerankFrom(bddp(cp<<1), rank);
        cp = DB->M->select(DB->M, (bddp)dzrefflags->succ(dzrefflags, DB->M->rank(DB->M, find_close(DB->U, cp), 1)+1, 1), 1);
      }
    }
  }
  else {
    // この時はfのreverse childrenにZDDのノードが含まれる
    vector<bddp> revzzeros; // reverse zero nodes in the ZDD part
    // start making revzzeros
    {
      bddp firstchild;
      i64 idx = 0;
      if (IN_DENSE(f)) {
        bddp fidx = (bddp)DB->M->rank(DB->M, f>>1, 1);
        idx = junctionflags->rank(junctionflags, fidx, 1) - 1;
        firstchild = dzmemo[idx];
      }
      else if (B_CST(f)) {
        idx = 0;
        firstchild = dzmemo[0];
      }
      else {
        firstchild = B_GET_BDDP((B_NP(f))->f0);
      }
      int size = 0;
      for (bddp revzero=firstchild; revzero!=bddnull; revzero=B_GET_BDDP((B_NP(revzero))->nx)) {
        size++;
      }
      if (size) {
        revzzeros.reserve(size);
        for (bddp revzero=firstchild; revzero!=bddnull; revzero=B_GET_BDDP((B_NP(revzero))->nx)) {
          revzzeros.push_back(revzero);
        }
        sort(revzzeros.begin(), revzzeros.end(), CompareRevZeros());
        if (IN_DENSE(f)||B_CST(f)) {
          dzmemo[idx] = revzzeros.front();
        }
        else {
          B_SET_BDDP((B_NP(f))->f0, revzzeros.front());
        }
        for (vector<bddp>::iterator it=revzzeros.begin(); it!=revzzeros.end(); it++) {
          B_SET_BDDP((B_NP(*it))->nx, (it+1!=revzzeros.end())?*(it+1):bddnull);
        }
      }
    }
    // finish making revzzeros

    if (IN_DENSE(f)||(DB&&B_CST(f))) {
      // f is a DZ-junction node
      // このケースが最も厄介
      // revzzerosをイテレートしながら，DZ-nodeも同時に見ていく
      // DZ-nodeはすでに正しい順に並んでいることに注意
      i64 cp, fcp;
      // get cp and fcp
      {
        i64 fp = (B_CST(f))?0:f>>1; // fのBP上の位置
        bddp fidx = (B_CST(f))?0:(bddp)DB->M->rank(DB->M, fp, 1);
        fcp = find_close(DB->U, fp); // fの閉じカッコの位置
        cp = DB->M->select(DB->M, (bddp)dzrefflags->succ(dzrefflags, fidx+1, 1), 1); // 現在注目しているDZ側の子ノードのBP上の位置
      }
      (*rank)++; // fに(*rank)を割当て
      int from = 0;
      int to = 0;
      bddp lastdz = f;
      bddp last = f;
      while (from<revzzeros.size()) {
        bddvar lev = bddlevofvar(bddtop(revzzeros[from]));
        // DZ-nodeをイテレートしていく
        while (cp<fcp&&bddlevofvar(bddtop(bddp(cp<<1)))>lev) {
          calcPrerankFrom(bddp(cp<<1), rank);
          lastdz = bddp(cp<<1); // lastdzの更新
          cp = DB->M->select(DB->M, (bddp)dzrefflags->succ(dzrefflags, DB->M->rank(DB->M, find_close(DB->U, cp), 1)+1, 1), 1);
        }
        
        bddp stsize = 0;
        do {
          setLastDZ(revzzeros[to], lastdz, f==lastdz);
          stsize += subtreeSize(revzzeros[to]);
        } while ((to+1<revzzeros.size()&&bddlevofvar(bddtop(revzzeros[to+1]))==lev)&&(++to));
        // from~toのlevが同じ
        
        PrankRange *range;
        if (from!=to||(cp<fcp&&bddlevofvar(bddtop(bddp(cp<<1)))==lev)) {
          range = new PrankRange();
          range->znodes.reserve(to-from+1);
          for (int i=from; i<=to; i++) {
            range->znodes.push_back(revzzeros[i]);
          }
          range->dzleft = range->dzright = (cp<fcp&&bddlevofvar(bddtop(bddp(cp<<1)))==lev)?bddp(cp<<1):bddnull;
          range->from = *rank;
          range->prev = last;
          range->firstrange = (f==last);
        }
        else {
          range = nullptr;
        }
        // update last
        last = revzzeros[to];
        
        while (cp<fcp&&bddlevofvar(bddtop(bddp(cp<<1)))==lev) {
          // 同じレベルなので，rangeに追加
          lastdz = range->dzright = bddp(cp<<1);
          stsize += subtreeSize(bddp(cp<<1));
          cp = DB->M->select(DB->M, (bddp)dzrefflags->succ(dzrefflags, DB->M->rank(DB->M, find_close(DB->U, cp), 1)+1, 1), 1);
        }
        if (range) {
          PrankList[lev-1].push_back(range);
          *rank += stsize;
        }
        else {
          // continue depth-first search
          calcPrerankFrom(revzzeros[from], rank);
        }
        from = ++to;
      }
      while (cp<fcp) {
        calcPrerankFrom(bddp(cp<<1), rank);
        cp = DB->M->select(DB->M, (bddp)dzrefflags->succ(dzrefflags, DB->M->rank(DB->M, find_close(DB->U, cp), 1)+1, 1), 1);
      }
    }
    else {
      // f is a Z-node
      if (!B_CST(f)) {
        setPrerank(f, (*rank)); // fに(*rank)を割当て
      }
      (*rank)++;
      int from = 0;
      int to = 0;
      bddp last = f;
      while (from<revzzeros.size()) {
        bddp stsize = 0;
        bddvar lev = bddlevofvar(bddtop(revzzeros[from]));
        do {
          stsize += subtreeSize(revzzeros[to]);
        } while ((to+1<revzzeros.size()&&bddlevofvar(bddtop(revzzeros[to+1]))==lev)&&(++to));
        // from~toのlevが同じ
        if (from!=to) {
          PrankRange *range = new PrankRange();
          range->znodes.reserve(to-from+1);
          for (int i=from; i<=to; i++) {
            range->znodes.push_back(revzzeros[i]);
          }
          range->dzleft = range->dzright = bddnull;
          range->from = *rank;
          range->prev = last;
          range->firstrange = (last==f);
          PrankList[lev-1].push_back(range);
          *rank += stsize;
        }
        else {
          // continue depth-first search
          calcPrerankFrom(revzzeros[from], rank);
        }
        // update last
        last = revzzeros[to];
        from = ++to;
      }
    }
  }
}

static void setPrerank(bddp f, bddp rank) {
  zmemo1[zrefflags->rank(zrefflags, B_NDX(f), 1) - 1] = rank;
  order->setbit(order, rank, 1);
}

static bddp prerank(bddp f) {
  if (B_CST(f)) {
    return 0;
  }
  else if (IN_DENSE(f)) {
    f = (bddp)dzrefflags->rank(dzrefflags, DB->M->rank(DB->M, f>>1, 1), 1);
    return (bddp)order->select(order, f, 0);
  }
  else {
    return zmemo1[zrefflags->rank(zrefflags, B_NDX(f), 1) - 1];
  }
}

//zjuncnodeに対してのみ用いる
static void setLastDZ(bddp f, bddp dz, bool isParent) {
  if (!isParent) {
    i64 dzp = (B_CST(dz))?0:dz>>1;
    i64 dzcp = find_close(DB->U, dzp);
    dz = (bddp)DB->M->pred(DB->M, dzcp, 1)<<1;
  }
  zmemo2[zrefflags->rank(zrefflags, B_NDX(f), 1) - 1] = dz;
}


///////////////////////////////////////////////////////////////////
// Algorithm 5: How to obatin the last visited DZ-node before
//              the Z-node f in preorder.
///////////////////////////////////////////////////////////////////
static bddp lastDZ(bddp f) {
  if (isJunction(f)) {
    return zmemo2[(bddp)zrefflags->rank(zrefflags, B_NDX(f), 1) - 1];
  }
  else {
    bddp root = zmemo2[zrefflags->rank(zrefflags, B_NDX(f), 1)-1];
    return zmemo2[zrefflags->rank(zrefflags, B_NDX(root), 1)-1];
  }
}

static void convertToDense(bddp f, bool copy, FBitIO &uio, FBitIO &mio, FBitIO &iio) {
  if (copy) {
    // この時は，元のzero-edge treeをコピーすれば良い
    i64 fp = (B_CST(f))?0:f>>1;
    i64 fcp = find_close(DB->U, fp);
    bddp last = (bddp)DB->M->pred(DB->M, fcp-1, 1);
    for (i64 i=fp+1; i<fcp; i++) { // 何ビットかまとめて出力できる
      //uio.write(getbit(DB->U->B, i));
      uio.write(darray_getbit(DB->U->da, i));
      mio.write(DB->M->getbit(DB->M, i));
    }
    // one-child arrayを順番に見ているだけだから，毎回rankを計算する必要は無い
    bddp r1 = bddp(DB->M->rank(DB->M, fp, 1)+1);
    bddp r2 = bddp(DB->M->rank(DB->M, last, 1));
    for (bddp r=r1; r<=r2; r++) {
      bddp one = done_r(r);
      bddp d = prerank(one);
      iio.write((one&1)?(d<<1)|1:d<<1);
    }
    return;
  }

  bddvar flev = bddlevofvar(bddtop(f));
  bddvar i = flev;
  
  if ((IN_DENSE(f)||(DB&&B_CST(f)))&&(!isJunction(f))) {
    i64 cp, fcp;
    // get cp and fcp
    {
      i64 fp = (B_CST(f))?0:f>>1; // fのBP上の位置
      bddp fidx = (B_CST(f))?0:(bddp)DB->M->rank(DB->M, fp, 1);
      fcp = find_close(DB->U, fp); // fの閉じカッコの位置
      cp = DB->M->select(DB->M, (bddp)dzrefflags->succ(dzrefflags, fidx+1, 1), 1); // 現在注目しているDZ側の子ノードのBP上の位置
    }
    while (cp<fcp) {
      bddp revzero = bddp(cp<<1);
      cp = DB->M->select(DB->M, (bddp)dzrefflags->succ(dzrefflags, DB->M->rank(DB->M, find_close(DB->U, cp), 1)+1, 1), 1);
      bddvar lev = bddlevofvar(bddtop(revzero));
      if (i+1<lev) {
        // first_childのとき
        while (i+1<lev) {
          uio.write(1);
          mio.write(0);
          i++;
        }
      }
      else if (i+1>lev) {
        while (i+1>lev) {
          uio.write(0);
          mio.write(0);
          i--;
        }
      }
      uio.write(1);
      mio.write(1);
      bddp one = done(revzero);
      bddp d = (one&1)?(prerank(one)<<1)|1:prerank(one)<<1;
      iio.write(d);
      
      convertToDense(revzero, canCopy(revzero), uio, mio, iio);
      
      uio.write(0);
      mio.write(0);
    }
  }
  else {
    // この時はfのreverse childrenにZDDのノードが含まれる
    
    bddp revzzero;
    if (IN_DENSE(f)) {
      bddp fidx = (bddp)DB->M->rank(DB->M, f>>1, 1);
      i64 idx = junctionflags->rank(junctionflags, fidx, 1) - 1;
      revzzero = dzmemo[idx];
    }
    else if (B_CST(f)) {
      revzzero = dzmemo[0];
    }
    else {
      revzzero = B_GET_BDDP((B_NP(f))->f0);
    }
    
    if (IN_DENSE(f)||(DB&&B_CST(f))) {
      i64 cp, fcp;
      // get cp and fcp
      {
        i64 fp = (B_CST(f))?0:f>>1; // fのBP上の位置
        bddp fidx = (B_CST(f))?0:(bddp)DB->M->rank(DB->M, fp, 1);
        fcp = find_close(DB->U, fp); // fの閉じカッコの位置
        cp = DB->M->select(DB->M, (bddp)dzrefflags->succ(dzrefflags, fidx+1, 1), 1); // 現在注目しているDZ側の子ノードのBP上の位置
      }
      while (revzzero!=bddnull||cp<fcp) {
        bddp revzero;
        if (revzzero!=bddnull&&cp<fcp) {
          if (ComparePranks()(revzzero, bddp(cp<<1))) {
            // Z-nodeを進める
            revzero = revzzero;
            revzzero=B_GET_BDDP((B_NP(revzzero))->nx);
          }
          else {
            // DZ-nodeを進める
            revzero = bddp(cp<<1);
            cp = bddp(DB->M->select(DB->M, (bddp)dzrefflags->succ(dzrefflags, DB->M->rank(DB->M, find_close(DB->U, cp), 1)+1, 1), 1));
          }
        }
        else if (cp<fcp) {
          // DZ-nodeを進める
          revzero = bddp(cp<<1);
          cp = bddp(DB->M->select(DB->M, (bddp)dzrefflags->succ(dzrefflags, DB->M->rank(DB->M, find_close(DB->U, cp), 1)+1, 1), 1));
        }
        else {
          // Z-nodeを進める
          revzero = revzzero;
          revzzero=B_GET_BDDP((B_NP(revzzero))->nx);
        }
        bddvar lev = bddlevofvar(bddtop(revzero));
        if (i+1<lev) {
          // first_childのとき
          while (i+1<lev) {
            uio.write(1);
            mio.write(0);
            i++;
          }
        }
        else if (i+1>lev) {
          while (i+1>lev) {
            uio.write(0);
            mio.write(0);
            i--;
          }
        }
        uio.write(1);
        mio.write(1);
        bddp one;
        if (IN_DENSE(revzero)) {
          one = done(revzero);
        }
        else {
          one = B_GET_BDDP(B_NP(revzero)->f1);
        }
        bddp d = prerank(one);
        d = (one&1)?(d<<1)|1:d<<1;
        iio.write(d);
        
        convertToDense(revzero, canCopy(revzero), uio, mio, iio);
        
        uio.write(0);
        mio.write(0);
      }
    }
    else {
      // f is a Z-node
      for (; revzzero!=bddnull; revzzero=B_GET_BDDP((B_NP(revzzero))->nx)) {
        bddvar lev = bddlevofvar(bddtop(revzzero));
        if (i+1<lev) {
          // first_childのとき
          while (i+1<lev) {
            uio.write(1);
            mio.write(0);
            i++;
          }
        }
        else if (i+1>lev) {
          while (i+1>lev) {
            uio.write(0);
            mio.write(0);
            i--;
          }
        }
        uio.write(1);
        mio.write(1);
        bddp one = B_GET_BDDP(B_NP(revzzero)->f1);
        bddp d = prerank(one);
        d = (one&1)?(d<<1)|1:d<<1;
        iio.write(d);
        
        convertToDense(revzzero, canCopy(revzzero), uio, mio, iio);
        
        uio.write(0);
        mio.write(0);
      }
    }
  }
  while (i>flev) {
    uio.write(0);
    mio.write(0);
    i--;
  }
}

static bool canCopy(bddp f) {
  if (hasZnode(f)) {
    return false;
  }
  else {
    bddp fidx, lastidx;
    i64 fp, fcp;
    if (B_CST(f)) {
      fidx = 0;
      fp = 0;
    }
    else {
      fidx = (bddp)DB->M->rank(DB->M, f>>1, 1);
      fp = f>>1;
    }
    fcp = find_close(DB->U, fp);
    lastidx = (bddp)DB->M->rank(DB->M, fcp, 1);
    return (bddp)(dzrefflags->rank(dzrefflags, lastidx, 1)-dzrefflags->rank(dzrefflags, fidx, 1))==lastidx-fidx;
  }
}

static bool hasZnode(bddp f) {
  if (B_CST(f)) {
    if (junctionflags->rank(junctionflags, junctionflags->n, 1)) {
      return true;
    }
    return false;
  }
  else if (IN_DENSE(f)) {
    bddp fidx = (bddp)DB->M->rank(DB->M, f>>1, 1);
    i64 fp = f>>1;
    i64 fcp = find_close(DB->U, fp);
    bddp lastidx = (bddp)DB->M->rank(DB->M, fcp, 1);
    if (junctionflags->rank(junctionflags, fidx-1, 1)!=junctionflags->rank(junctionflags, lastidx, 1)) {
      return true;
    }
    return false;
  }
  else {
    return true;
  }
}

static void getnewdense(bddp *roots, FBitIO &uio, FBitIO &mio, FBitIO &iio, unsigned long lenU, unsigned long lenM, unsigned long lenI, int opt) {
  if (DB) {
    // free old U, M, I
    bp_free(DB->U);
    delete DB->U; // sada 2015/1/16
    bitvector_free(DB->M);
    bitvector_free(DB->I);
    bitvector_free(DB->visit);
  }
  else {
    DB = (struct DenseBDD *)malloc(sizeof(struct DenseBDD));
  }
  // free old BDD.
  bddvar varNum = VarUsed;
  free(Node);
  Node = NULL;
  for(int i=0; i<VarSpc; i++)
  {
    if(Var[i].hash_32) free(Var[i].hash_32);
#ifdef B_64
    if(Var[i].hash_h8) free(Var[i].hash_h8);
#endif
  }
  free(Var);
  Var = NULL;
  free(VarID);
  VarID = NULL;
  
  unsigned long buf;
  
  // create new dummy node vector
  bitvector *M = bitvector_new(lenM);
  unsigned long countM = (lenM+sizeof(bitvec_t)*8-1)/(sizeof(bitvec_t)*8); //危なそう
  for (unsigned long i=0; i<countM; i++) {
    mio.read(&buf);
    M->buf[i] = (bitvec_t)buf;
  }
  //    bitvector_makeindex(M, 512, SDARRAY_RANK1 | SDARRAY_SELECT1);
  if (opt & COMP_M) {
    bitvector_makeindex(M, 128, SDARRAY_RANK1 | SDARRAY_SELECT1 | SDARRAY_SPARSE);
  }
  else {
    bitvector_makeindex(M, 128, SDARRAY_RANK1 | SDARRAY_SELECT1); // sada 4/14
  }
  
//  i64 numOfNodes = M->rank(M, M->n, 1)+1;
//  bddp *rank_memo = new bddp[numOfNodes];
//  int i=0;
//  for (bddp pos=0; pos<M->n; pos=(bddp)M->succ(M, pos+1, 1)) {
//    rank_memo[i++] = pos;
//  }
  //    create new one child array
  bddp last = (bddp)M->pred(M, M->n-1, 1);
  bddp width = 1;
  while ((1<<width)-1<last) {
    width++;
  }
  width++;
  bitvector *I = bitvector_new(lenI*width);
  for (int i=0; i<lenI; i++) {
    iio.read(&buf);
    if (buf>>1) {
      //            rootじゃない時
      buf = (M->select(M, buf>>1, 1)<<1)|(buf&1); //BP上の位置に変換
//      buf = (rank_memo[buf>>1]<<1)|(buf&1);
    }
    I->setbits(I, i*width, (int)width, (bddp)buf);
  }
  for (int i=0; roots[i]!=bddnull; i++) {
    if (!B_CST(roots[i])) {
      roots[i] = (bddp)(M->select(M, roots[i]>>1, 1)<<1)|(roots[i]&1);
//      roots[i] = (rank_memo[roots[i]>>1]<<1)|(roots[i]&1);
    }
  }
//  delete [] rank_memo;
  
  //    create new BP
  unsigned long lenB = (lenU+D-1)/D;
  pb *B = (pb *)malloc(lenB*sizeof(pb));
  B[lenB-1] = 0;
  unsigned long countB = (lenU+sizeof(pb)*8-1)/(sizeof(pb)*8);
  for (unsigned long i=0; i<countB; i++) {
    uio.read(&buf);
    B[i] = (pb)buf;
  }
  bp *U = new bp;
  if (opt & COMP_U) {
    bp_construct(U, lenU, B, OPT_FAST_LCA | OPT_DEGREE | OPT_COMP_RLE);
  }
  else {
    bp_construct(U, lenU, B, OPT_FAST_LCA | OPT_DEGREE);
  }
  
  // get new dense
  DB->U = U;
  DB->M = M;
  DB->I = I;
  DB->width = width;
  DB->visit = bitvector_new(lenI);
  FirstIdx = last + 1;
  
  bddinit(B_NODE_SPC0, B_NODE_MAX);
  for (int i=0; i<varNum; i++) bddnewvar();
}

