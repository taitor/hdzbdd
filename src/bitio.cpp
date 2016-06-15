//
//  bitio.cpp
//  hdzbdd
//
//  Created by Taito Lee on 2015/10/28.
//  Copyright © 2015年 otita. All rights reserved.
//
#include <limits.h>
#include "bitio.h"

FBitIO::FBitIO(const char fname[], const char *mode) {
  fp_ = fopen(fname, mode);
}

FBitIO::~FBitIO() {
  fclose(fp_);
}

void FBitIO::start_writing(int width) {
  rewind();
  width_ = width;
  c_ = 0;
  i_ = 0;
  size_ =sizeof(c_)*8/width_;
  mask_ = (1UL<<width_) - 1;
  count_ = 0;
}

void FBitIO::write(unsigned long num) {
  num = num&mask_;
  c_ |= num<<(width_*(size_-(++i_)));
  if (i_==size_) {
    fwrite(&c_, sizeof(c_), 1, fp_);
    count_ += size_;
    c_ = 0;
    i_ = 0;
  }
}

unsigned long FBitIO::finish_writing() {
  //    残りを書き込む
  if (i_) {
    count_ += i_;
    fwrite(&c_, sizeof(c_), 1, fp_);
  }
  return count_;
}

void FBitIO::start_reading(int width) {
  rewind();
  width_ = width;
  c_ = 0;
  size_ = sizeof(c_)*8/width_;
  mask_ = (width_==sizeof(c_)*8) ? ULONG_MAX : (1UL<<width_) - 1;
  i_ = size_;
}

int FBitIO::read(unsigned long *buf) {
  if (i_==size_) {
    size_t n = fread(&c_, sizeof(c_), 1, fp_);
    if (n==0) {
      return -1;
    }
    i_ = 0;
  }
  unsigned long num = (c_ >> (width_ * (size_ - ++i_))) & mask_;
  *buf = num;
  return 0;
}

void FBitIO::finish_reading() {
  
}

void FBitIO::rewind() {
  ::rewind(fp_);
}