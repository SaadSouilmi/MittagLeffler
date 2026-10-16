#include <complex.h>
#include <stdio.h>

#include "../mittag_leffler.h"
#include "../mittag_leffler.cpp"

I32
main()
{
	const F64 Parameters[][2] = {
		{0.55, 0.55}, {0.55, 1.0}, {0.55, 1.55},
		{0.75, 0.75}, {0.75, 1.0}, {0.75, 1.75},
		{0.95, 0.95}, {0.95, 1.0}, {0.95, 1.95},
		{1.0, 1.0}, {1.0, 2.0}, {0.5, 1.0},
	};
	printf("alpha,beta,x,value\n");
	for (const auto &Pair : Parameters)
	{
		MLQuad Q;
		MLQuadInit(Pair[0], Pair[1], &Q);
		printf("%.17g,%.17g,0,%.17g\n", Pair[0], Pair[1], MLQuadEval(0, &Q));
		for (I32 Index = 0; Index <= 300; ++Index)
		{
			F64 X = pow(10.0, -8.0 + 12.0*(F64)Index/300.0);
			printf("%.17g,%.17g,%.17g,%.17g\n", Pair[0], Pair[1], X, MLQuadEval(X, &Q));
		}
	}
	return 0;
}
