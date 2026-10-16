#include <math.h>
#include <complex.h>
#include <float.h>


internal void
MLQuadInit(F64 Alpha, F64 Beta, MLQuad *Q)
{
	*Q = (MLQuad){.Alpha = Alpha, .Beta = Beta, .AtZero = 1/tgamma(Beta)};
	if (Alpha == 1 && (Beta == 1 || Beta == 2)) { return; }
	for (I32 K = 0; K <= MLQuadK; ++K)
	{
		complex64 W = 1 + _Complex_I*(K*0.2);
		complex64 A = cexp(W*W)*cpow(W, 2*(Alpha - Beta) + 1)/PI;
		complex64 P = cpow(W, 2*Alpha);
		Q->A[K] = creal(A); Q->B[K] = cimag(A);
		Q->C[K] = creal(P); Q->D[K] = cimag(P);
	}
}

internal F64
MLQuadEval(F64 X, const MLQuad *Q)
{
	if (!isfinite(X) || X < 0) { return NAN; }
	if (!X) { return Q->AtZero; }
	if (Q->Alpha == 1)
	{
		if (Q->Beta == 1) { return exp(-X); }
		if (Q->Beta == 2) { return -expm1(-X)/X; }
	}
	F64 Sum = Q->A[0]/(Q->C[0] + X);
	for (I32 K = 1; K <= MLQuadK; ++K)
	{
		F64 CX = Q->C[K] + X;
		F64 Num = Q->A[K]*CX + Q->B[K]*Q->D[K];
		F64 Den = CX*CX + Q->D[K]*Q->D[K];
		Sum += 2*Num/Den;
	}
	return 0.2*Sum;
}


