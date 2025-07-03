#include <simd/simd.h>

int isnanf(float v) {
    return isnan(v);
}
int isnand(double v) {
    return isnan(v);
}
simd_int2 isnanf2(simd_float2 v) {
    return simd_make_int2(isnan(v.x),isnan(v.y));
}
simd_int2 isnand2(simd_double2 v) {
    return simd_make_int2(isnan(v.x),isnan(v.y));
}

int isinff(float v) {
    return isinf(v);
}
int isinfd(double v) {
    return isinf(v);
}
simd_int2 isinff2(simd_float2 v) {
    return simd_make_int2(isinf(v.x),isinf(v.y));
}

simd_int2 isinfd2(simd_double2 v) {
    return simd_make_int2(isinf(v.x),isinf(v.y));
}
