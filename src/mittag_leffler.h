#ifndef MITTAG_LEFFLER_H
#define MITTAG_LEFFLER_H

#include "base.h"

enum { MLQuadK = 30 };
struct MLQuad
{
	F64 A[MLQuadK + 1], B[MLQuadK + 1], C[MLQuadK + 1], D[MLQuadK + 1];
	F64 Alpha, Beta, AtZero;
};

internal void MLQuadInit(F64 Alpha, F64 Beta, MLQuad *Q);
internal F64 MLQuadEval(F64 X, const MLQuad *Q);


#endif // MITTAG_LEFFLER_H

