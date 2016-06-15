//
//  powerset_union.cpp
//  hdzbdd
//
//  Created by Taito Lee on 2015/10/28.
//  Copyright © 2015年 otita. All rights reserved.
//

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include <time.h>
#include <vector>
#include <algorithm>

#ifdef USE_DENSE
#include "hdbdd.h"
#else
#include "bddc.h"
#endif

using namespace std;

#define N_DEFAULT 1000
#define K_DEFAULT 150
#define M_DEFAULT 30

#ifdef USE_DENSE
extern double THETA;
#endif

#ifndef USE_DENSE
#include <sys/resource.h>
  void showresult(void) {
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage)==0) {
      cout << "maxrss " << usage.ru_maxrss << endl;
    }
  }
#endif

int main(int argc, const char * argv[]) {
  srand(1);
  int N, K, M;
  if (argc>1) N = atoi(argv[1]);
  else N = N_DEFAULT;
  if (argc>2) K = atoi(argv[2]);
  else K = K_DEFAULT;
  if (argc>3) M = atoi(argv[3]);
  else M = M_DEFAULT;
#ifdef USE_DENSE
  int opt = 0;
  if (argc>4) THETA = atof(argv[4]);
  printf("THETA=%lf\n", THETA);
  if (argc>5) opt = atoi(argv[5]) ? COMP_M | COMP_U : 0;
  printf("COMPRESS: %s\n", opt ? "YES" : "NO");
#endif
  printf("powerset union (N,K,M)=(%d,%d,%d)\n", N, K, M);
#ifdef USE_DENSE
  densezdd_init();
#endif
  clock_t t01 = clock();
  clock_t c_t = 0;
  bddinit(0, UINT_MAX);
  vector<bddvar> x;
  x.reserve(N);
  srand(3);
  for (int i=0; i<N; i++) {
    x.push_back(bddnewvar());
  }
  bddp roots[2] = {bddempty, bddnull};
  for (int i=0; i<M; i++) {
    random_shuffle(x.begin(), x.end());
    bddp temp = bddsingle;
    for (int j=0; j<K; j++) {
      temp = bddunion(temp, bddchange(temp, x[j]));
    }
    roots[0] = bddunion(roots[0], temp);
#ifdef USE_DENSE
        clock_t t1 = clock();
        compress(roots, opt);
        clock_t t2 = clock();
        c_t += t2 - t1;
#endif
  }
  
  clock_t t02 = clock();
  cout << "Total time:" << (double)(t02-t01)/CLOCKS_PER_SEC << endl;
#ifdef USE_DENSE
  cout << "Compress time:" << (double)c_t/CLOCKS_PER_SEC << endl;
#endif
  showresult();
  cout << bddcard(roots[0]) << endl;
  
  return 0;
}
