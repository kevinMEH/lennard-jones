#include <math.h>
#include "random.h"

static XorshiftState state;

void setXorshiftState(uint64_t first, uint64_t second) {
    state.x[0] = first;
    state.x[1] = second;
}

uint64_t xorshift64() {
    uint64_t t = state.x[0];
    uint64_t s = state.x[1];
    state.x[0] = s;
    t ^= t << 23;
    t ^= t >> 18;
    t ^= s ^ (s >> 5);
    state.x[1] = t;
    return t + s;
}

uint64_t xorshift64State(XorshiftState* state) {
    uint64_t t = state->x[0];
    uint64_t s = state->x[1];
    state->x[0] = s;
    t ^= t << 23;
    t ^= t >> 18;
    t ^= s ^ (s >> 5);
    state->x[1] = t;
    return t + s;
}

double xorshiftDouble() {
    return (double) xorshift64() / (0x8000000000000000u * 2.0);
}

double xorshiftDoubleState(XorshiftState* state) {
    return (double) xorshift64State(state) / (0x8000000000000000u * 2.0);
}

double xorshiftNormal(double center, double deviation) {
    return center + deviation *
    sqrt(-2 * log(xorshiftDouble())) * cos(2 * M_PI * xorshiftDouble());
}

double xorshiftNormalState(XorshiftState* state, double center, double deviation) {
    return center + deviation *
    sqrt(-2 * log(xorshiftDoubleState(state))) * cos(2 * M_PI * xorshiftDoubleState(state));
}