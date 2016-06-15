//
//  bitio.h
//  hdzbdd
//
//  Created by Taito Lee on 2015/10/28.
//  Copyright © 2015年 otita. All rights reserved.
//

#ifndef bitio_h
#define bitio_h

#include <stdio.h>

class FBitIO {
private:
  FILE *fp_;
  unsigned long c_;
  int i_;
  int size_;
  unsigned long count_;
  int width_;
  unsigned long mask_;
public:
  FBitIO(const char *fname, const char *mode="w+b");
  virtual ~FBitIO();
  void start_writing(int width);
  void write(unsigned long num);
  unsigned long finish_writing();
  void start_reading(int width);
  int read(unsigned long *buf);
  void finish_reading();
  void rewind();
};


#endif /* bitio_h */
