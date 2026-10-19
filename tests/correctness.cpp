#include <stdio.h>
#include "mittag_leffler.h"
#include "mittag_leffler.cpp"

internal I32
Sample(F64 Alpha, F64 Beta, I32 Steps, F64 Low, F64 High)
{
	MLQuad Q;
	MLQuadInit(Alpha, Beta, &Q);
	if (MLQuadEval(0, &Q) != 1/tgamma(Beta) ||
		MLQuadEval(-0.0, &Q) != Q.AtZero ||
		!isnan(MLQuadEval(-1, &Q)) || !isnan(MLQuadEval(NAN, &Q)) ||
		!isnan(MLQuadEval(INFINITY, &Q)) || !isnan(MLQuadEval(-INFINITY, &Q)))
	{
		fprintf(stderr, "Boundary failure: alpha=%g beta=%g\n", Alpha, Beta);
		return 1;
	}
	printf("%.17g,%.17g,0,%.17g\n", Alpha, Beta, MLQuadEval(0, &Q));
	for (I32 Index = 0; Index <= Steps; ++Index)
	{
		F64 X = pow(10.0, Low + (High - Low)*(F64)Index/Steps);
		F64 Y = MLQuadEval(X, &Q);
		if (!isfinite(Y))
		{
			fprintf(stderr, "Nonfinite result: alpha=%g beta=%g x=%g\n", Alpha, Beta, X);
			return 1;
		}
		printf("%.17g,%.17g,%.17g,%.17g\n", Alpha, Beta, X, Y);
	}
	return 0;
}

I32
main()
{
	const F64 Original[][2] = {
		{0.55, 0.55}, {0.55, 1}, {0.55, 1.55},
		{0.75, 0.75}, {0.75, 1}, {0.75, 1.75},
		{0.95, 0.95}, {0.95, 1}, {0.95, 1.95},
		{1, 1}, {1, 2}, {0.5, 1},
	};
	I32 Failures = 0;
	printf("alpha,beta,x,value\n");
	for (const auto &Pair : Original)
	{
		Failures += Sample(Pair[0], Pair[1], 300, -8, 4);
	}
	for (I32 Index = 1; Index <= 23; ++Index)
	{
		const F64 Edges[] = {0.001, 0.01, 0.99, 0.999};
		F64 Alpha = Index <= 19 ? (F64)Index/20 : Edges[Index - 20];
		Failures += Sample(Alpha, Alpha, 1000, -12, 6);
		Failures += Sample(Alpha, 1, 1000, -12, 6);
		Failures += Sample(Alpha, Alpha + 1, 1000, -12, 6);
	}
	fprintf(stderr, "Boundary and finite-value checks: %s\n", Failures ? "FAIL" : "PASS");
	return Failures;
}
