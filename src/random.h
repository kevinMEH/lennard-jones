#ifndef RANDOM
#define RANDOM

#include "inttypes.h"

typedef struct XorshiftState {
    uint64_t x[2];
} XorshiftState;

void setXorshiftState(uint64_t first, uint64_t second);

uint64_t xorshift64();
uint64_t xorshift64State(XorshiftState* state);

// Random between 0 and 1 not including 1
double xorshiftDouble();
double xorshiftDoubleState(XorshiftState* state);

double xorshiftNormal(double center, double deviation);
double xorshiftNormalState(XorshiftState* state, double center, double deviation);

#endif