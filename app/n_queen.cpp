//
//  n_queen.cpp
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
#elif USE_DENSE_OLD
#include "hdbdd_old.h"
#else
#include "bddc.h"
#endif

using namespace std;

#define N_DEFAULT 13

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
  int N;
  if (argc>1) N = atoi(argv[1]);
  else N = N_DEFAULT;
#ifdef USE_DENSE
  int opt = 0;
  if (argc>2) THETA = atof(argv[2]);
  printf("THETA=%lf\n", THETA);
  if (argc>3) opt = atoi(argv[3]) ? COMP_M | COMP_U : 0;
  printf("COMPRESS: %s\n", opt ? "YES" : "NO");
#endif
  printf("%d Queen\n", N);
#ifdef USE_DENSE
  densezdd_init();
#endif
  
  clock_t t01 = clock();
  clock_t c_t = 0;
  
  bddinit(0, UINT_MAX);
  bddvar x[N][N];
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      x[i][j] = bddnewvar();
    }
  }
  bddp roots[3] = {bddnull, bddnull, bddnull};
  for (int i=0; i<N; i++) {
    roots[1] = roots[0];
    roots[0] = bddempty;
    for (int j=0; j<N; j++) {
      if (i>0) {
        bddp temp = roots[1];
        for (int k=1; k<=i; k++) {
          int left = j-k;
          int center = j;
          int right = j+k;
          int row = i-k;
          
          if (left>=0) {
            temp = bddquot(temp, bddchange(bddsingle, x[row][left]));
          }
          
          temp = bddquot(temp, bddchange(bddsingle, x[row][center]));
          
          if (right<N) {
            temp = bddquot(temp, bddchange(bddsingle, x[row][right]));
          }
        }
        roots[0] = bddunion(roots[0], bddmult(bddchange(bddsingle, x[i][j]), temp));
#ifdef USE_DENSE
        clock_t t1 = clock();
        compress(roots, opt);
        clock_t t2 = clock();
        c_t += t2 - t1;
#endif
      }
      else {
        roots[0] = bddunion(roots[0], bddchange(bddsingle, x[i][j]));
      }
    }
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
