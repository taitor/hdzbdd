//
//  main.cpp
//  HDZBDD
//
//  Created by Taito Lee on 2015/08/18.
//  Copyright (c) 2015年 otita. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <vector>

#include <time.h>
#ifdef USE_DENSE
clock_t c_t = 0;
#endif

#include "SOP.h"
using namespace std;

#ifdef USE_DENSE
extern double THETA;
int opt=0;
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

ZBDD *pla_to_zbdd(istream &istream, int *out)
{
  ZBDD *f;
  string line;
  char b[100], b2[100];
  int ni, no, nl;
  ni = 0; no = 0; nl = 0;
  int i, j;
  int header=1;

  nl = 100000000;

  while (header) {
    getline(istream, line);
    if (line[0] == '.') {
      sscanf(line.c_str(), "%s", b);
      if (strcmp(b, ".i") == 0) {
        sscanf(line.c_str(), "%s %d", b2, &ni);
      } else 
      if (strcmp(b, ".o") == 0) {
        sscanf(line.c_str(), "%s %d", b2, &no);
      } else 
      if (strcmp(b, ".p") == 0) {
        sscanf(line.c_str(), "%s %d", b2, &nl);
//        header = 0;
      }
    } else 
    if (line[0] == '#') {
      continue;
    } else
    if (!isalnum(line[0])) {
      continue;
    } else {
      break;
    }
  }

  printf("number of inputs = %d\n", ni);
  printf("number of outputs = %d\n", no);

  *out = no;

  for (i=0; i<ni*2; i++) {
    bddnewvar();
  }

  f = new ZBDD[no];
#ifdef USE_DENSE
  bddword *roots = new bddword[no+1];
  roots[no] = bddnull;
#endif
  for (i=0; i<no; i++) {
    f[i] = ZBDD();
#ifdef USE_DENSE
    roots[i] = f[i]._zbdd;
#endif
  }

  for (j=0; j<nl; j++) {
    ZBDD cube = 1;
    for (i=0; i<ni; i++) {
      switch (line[i]) {
      case '1':
        cube = cube.Change((i+1)*2);
        break;
      case '0':
        cube = cube.Change((i+1)*2-1);
        break;
      case '-':
        break;
      default:
        printf("??? [%c][%d]\n", line[i], line[i]);
        exit(1);
      }
    }

    for (i=0; i<no; i++) {
      if (line[ni+1+i] == '1') {
        f[i] += cube;
#ifdef USE_DENSE
        roots[i] = f[i]._zbdd;
#endif
      }
    }

#ifdef USE_DENSE
    clock_t t1 = clock();
    compress(roots, opt);
    clock_t t2 = clock();
    c_t += t2 - t1;
    for (i=0; i<no; i++) {
      f[i]._zbdd = roots[i];
    }
#endif
    getline(istream, line);
    if (line[0] == 0 || line[0] == '.') break;
  }

#ifdef USE_DENSE
  delete [] roots;
#endif
  return f;
}


int main(int argc, const char * argv[]) {
  int i, out;
  clock_t t01 = clock();
  BDD_Init(256, 100000000);
#ifdef USE_DENSE
  bddword *roots;
  densezdd_init();
  if (argc>2) THETA = atof(argv[2]);
  printf("THETA = %lf\n", THETA);
  if (argc>3) opt = atoi(argv[3]) ? COMP_M | COMP_U : 0;
  printf("COMPRESS: %s\n", opt ? "YES" : "NO");
#endif
  ifstream input(argv[1]);
  if (input.fail()) {
    std::cerr << "cannot open " << argv[1] << std::endl;
    return -1;
  }
  ZBDD *fz = pla_to_zbdd(input, &out);

#ifdef USE_DENSE
  roots = new bddword[out+2];
  for (i=0; i<out; i++) roots[i] = fz[i]._zbdd;
  roots[out] = roots[out+1] = bddnull;
#endif

  for (i=0; i<out; i++) {
#ifdef USE_DENSE
    SOP f = SOP();
    f._zbdd._zbdd = roots[i];
#else
    SOP f = SOP(fz[i]);
#endif
//    cout << "i Size " << f.Size() << " Cube " << f.Cube() << " Lit " << f.Lit() << endl;
    SOP p = f.Divisor();
    while (p!=f) {
      // pが新たな変数
      int v = SOP_NewVar(); // pを表す変数
      f = SOP(1).And1(v)*(f/p) + (f%p);
      p = f.Divisor();
#ifdef USE_DENSE
      roots[i] = f._zbdd._zbdd;
      roots[out] = p._zbdd._zbdd;
      clock_t t1 = clock();
      compress(roots, opt);
      clock_t t2 = clock();
      c_t += t2 - t1;
      f._zbdd._zbdd = roots[i];
      p._zbdd._zbdd = roots[out];
#endif
//      cout << "o Size " << f.Size() << " Cube " << f.Cube() << " Lit " << f.Lit() << endl;
    }
  }

#ifdef USE_DENSE
  free(roots);
#endif
  clock_t t02 = clock();
  cout << "Total time:" << (double)(t02-t01)/CLOCKS_PER_SEC << endl;
#ifdef USE_DENSE
  cout << "Compress time:" << (double)c_t/CLOCKS_PER_SEC << endl;
#endif
  showresult();
  
  return 0;
}
