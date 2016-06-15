#include <stdio.h>
#include <sys/resource.h>
#ifdef USE_DENSE
#include "hdbdd.h"
#else
#include "bddc.h"
#endif
#include "ZBDD.h"
#include "CtoI.h"
#include "darray.h"

void showmaxrss(void)
{
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage)==0) {
        cout << "maxrss " << usage.ru_maxrss << endl;
    }
}

#ifdef USE_DENSE
extern double THETA;
extern int opt;
//extern clock_t C_T;
//extern int bp_opt;
#endif

int main(int argc, char *argv[])
{
  int cache_size = 10;
  clock_t t01 = clock();
#ifdef USE_DENSE
  densezdd_init();
  THETA = 0.1;
  opt = 0;
//  C_T = 0;
//  bp_opt = 1;
  if (argc>3) THETA = atof(argv[3]);
  printf("THETA = %lf\n", THETA);
  if (argc>4) opt = atoi(argv[4]) ? COMP_M | COMP_U : 0;
    printf("COMPRESS: %s\n", opt ? "YES" : "NO");
#endif
  BDD_Init(256, 100000000);

  CtoI c;
  c = CtoI_LcmA(argv[1], NULL, atoi(argv[2]));

  clock_t t02 = clock();
  cout << "Total time:" << (double)(t02-t01)/CLOCKS_PER_SEC << endl;
#ifdef USE_DENSE
//  cout << "Compress time:" << (double)C_T/CLOCKS_PER_SEC << endl;
#endif
  showmaxrss();
  ZBDD z;
  z = c.GetZBDD();
  printf("Size %ld Card %ld\n", z.Size(), z.Card());
#ifdef USE_DENSE
  showresult();
  //hdbdd_info();
#ifdef USE_CACHE
  if (cache_size > 0) {
    printf("cache hit %ld miss %ld cache size %d\n", block_cache_hit, block_cache_miss, cache_size);
    block_cache_free();
  }
#endif
#endif
  return 0;
}

