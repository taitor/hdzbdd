/*****************************************
 *  BDD Package (SAPPORO-1.55)   - Header *
 *  (C) Shin-ichi MINATO  (Dac. 11, 2012) *
 ******************************************/

#ifndef bddc_h
#define bddc_h

/***************** Internal macro for index *****************/
#define B_VAR_WIDTH 16U  /* Width of variable index */
#define B_VAR_MASK       ((1U << B_VAR_WIDTH) - 1U)

/***************** Internal macro for bddp *****************/
//#define B_64
#ifdef B_64
#  define B_MSB_POS   39ULL
#  define B_LSB_MASK  1ULL
#else
#  define B_MSB_POS   31U
#  define B_LSB_MASK  1U
#endif
#define B_MSB_MASK  (B_LSB_MASK << B_MSB_POS)
#define B_INV_MASK  B_LSB_MASK /* Mask of inverter-flag */
#define B_CST_MASK  B_MSB_MASK /* Mask of constant-flag */
#define B_VAL_MASK  (B_MSB_MASK - 1U) /* Mask of value-field */

/***************** For stack overflow limit *****************/
extern const int BDD_RecurLimit;
extern int BDD_RecurCount;

/***************** External typedef *****************/
typedef unsigned int bddvar;
#ifdef B_64
typedef unsigned long long bddp;
#else
typedef unsigned int bddp;
#endif

/***************** External Macro *****************/
#define bddvarmax B_VAR_MASK /* Max value of variable index */
#define bddnull   B_VAL_MASK /* Special value for null pointer */
#define bddfalse  B_CST_MASK /* bddp of constant false (0) */
#define bddtrue   (bddfalse ^ B_INV_MASK)
/* bddp of constant true (1) */
#define bddempty  bddfalse /* bddp of empty ZBDD (0) */
#define bddsingle bddtrue  /* bddp of single unit ZBDD (1) */
#define bddconst(c) (((c) & B_VAL_MASK) | B_CST_MASK)
/* bddp of a constant valued node */
#define bddvalmax B_VAL_MASK  /* Max constant value */

/***************** Option for DenseZDD *****************/
#define COMP_U 1
#define COMP_M (1<<1)

/***************** External operations *****************/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/***************** Init. and config. ****************/
extern int    bddinit(bddp initsize, bddp limitsize);
extern bddvar bddnewvar();
extern bddvar bddnewvaroflev(bddvar lev);
extern bddvar bddlevofvar(bddvar v);
extern bddvar bddvaroflev(bddvar lev);
extern bddvar bddvarused();

/************** Basic logic operations *************/
extern bddp   bddprime(bddvar v);
extern bddvar bddtop(bddp f);
extern bddp   bddcopy(bddp f);
extern bddp   bddnot(bddp f);
extern bddp   bddand(bddp f, bddp g);
extern bddp   bddor(bddp f, bddp g);
extern bddp   bddxor(bddp f, bddp g);
extern bddp   bddnand(bddp f, bddp g);
extern bddp   bddnor(bddp f, bddp g);
extern bddp   bddxnor(bddp f, bddp g);
extern bddp   bddat0(bddp f, bddvar v);
extern bddp   bddat1(bddp f, bddvar v);

/********** Memory management and observation ***********/
extern void   bddfree(bddp f);
extern bddp   bddused();
extern int    bddgc();
extern bddp   bddsize(bddp f);
extern bddp   bddvsize(bddp *p, int lim);
extern void   bddexport(FILE *strm, bddp *p, int lim);
extern int    bddimport(FILE *strm, bddp *p, int lim);
extern void   bdddump(bddp f);
extern void   bddvdump(bddp *p, int lim);
extern void   bddgraph(bddp f);
extern void   bddgraph0(bddp f);
extern void   bddvgraph(bddp *p, int lim);
extern void   bddvgraph0(bddp *p, int lim);

/************** Advanced logic operations *************/
extern bddp   bddlshift(bddp f, bddvar shift);
extern bddp   bddrshift(bddp f, bddvar shift);
extern bddp   bddsupport(bddp f);
extern bddp   bdduniv(bddp f, bddp g);
extern bddp   bddexist(bddp f, bddp g);
extern bddp   bddcofactor(bddp f, bddp g);
extern bddp   bddimply(bddp f, bddp g);
extern bddp   bddrcache(unsigned char op, bddp f, bddp g);
extern void   bddwcache(unsigned char op, bddp f, bddp g, bddp h);

/************** ZBDD operations *************/
extern bddp   bddoffset(bddp f, bddvar v);
extern bddp   bddonset(bddp f, bddvar v);
extern bddp   bddonset0(bddp f, bddvar v);
extern bddp   bddchange(bddp f, bddvar v);
extern bddp   bddintersec(bddp f, bddp g);
extern bddp   bddunion(bddp f, bddp g);
extern bddp   bddsubtract(bddp f, bddp g);
extern bddp   bddcard(bddp f);
extern bddp   bddlit(bddp f);
extern bddp   bddlen(bddp f);

//added by Lee
extern bddp   bddmult(bddp f, bddp g);
extern bddp   bdddiv(bddp f, bddp p);
extern bddp   bddquot(bddp f, bddp p);
extern int    bddimportz(FILE *strm, bddp *p, int lim);

/************** SeqBDD operations *************/
//extern bddp   bddpush B_ARG((bddp f, bddvar v));

/************** Hybrid DenseBDD operations, added by Lee **************/
extern void   compress(bddp *roots, int opt);

/************** For Debug, added by Lee ****************/
extern void showresult();

#ifdef __cplusplus
}
#endif /* __cplusplus */
extern void densezdd_init(void);
#endif /* bddc_h */
