#include <math.h>


internal void
MLQuadInit(F64 Alpha, F64 Beta, MLQuad *Q)
{
	*Q = (MLQuad){.Alpha = Alpha, .Beta = Beta, .AtZero = 1/tgamma(Beta)};
	if (Alpha == 1 && (Beta == 1 || Beta == 2)) { return; }
	for (I32 K = 0; K <= MLQuadK; ++K)
	{
		F64 U = K*0.2;
		F64 Theta = atan(U);
		F64 LogR = 0.5*log1p(U*U);
		F64 Power = 2*(Alpha - Beta) + 1;
		F64 Scale = exp(1 - U*U + Power*LogR)/PI;
		F64 Phase = 2*U + Power*Theta;
		Q->A[K] = Scale*cos(Phase);
		Q->B[K] = Scale*sin(Phase);
		F64 Radius = exp(2*Alpha*LogR);
		F64 Angle = 2*Alpha*Theta;
		Q->C[K] = Radius*cos(Angle);
		Q->D[K] = Radius*sin(Angle);
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


