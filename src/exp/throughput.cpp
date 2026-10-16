#include <complex.h>
#include <stdio.h>

#include "../base.h"
#include "../mittag_leffler.h"
#include "../profile.h"

#include "../mittag_leffler.cpp"
#include "../profile.cpp"

I32
main()
{
	const F64 Parameters[][2] = {
		{0.55, 0.55}, {0.55, 1.0}, {0.55, 1.55},
		{0.75, 0.75}, {0.75, 1.0}, {0.75, 1.75},
		{0.95, 0.95}, {0.95, 1.0}, {0.95, 1.95},
		{1.0, 1.0}, {1.0, 2.0},
	};
	const U64 BatchSize = 4096;
	const U64 MinNanos = 250000000;
	F64 Inputs[BatchSize];
	for (U64 Index = 0; Index < BatchSize; ++Index)
	{
		Inputs[Index] = pow(10.0, -6.0 + 8.0*(F64)Index/(F64)(BatchSize - 1));
	}

	printf("MLQuadEval: logarithmically spaced X in [1e-6, 100], initialization excluded\n");
	printf("Minimum 0.25 s per pair, %llu evaluations per batch\n\n",
		(unsigned long long)BatchSize);
	printf("%5s %5s %10s %14s %12s %14s %9s %18s\n",
		"alpha", "beta", "path", "eval/s", "ns/eval", "evaluations", "seconds", "checksum");

	for (const auto &Pair : Parameters)
	{
		MLQuad Q;
		MLQuadInit(Pair[0], Pair[1], &Q);

		F64 Warmup = 0;
		for (U64 Index = 0; Index < BatchSize; ++Index)
		{
			Warmup += MLQuadEval(Inputs[Index], &Q);
		}
		__asm__ volatile("" : : "m"(Warmup) : "memory");

		F64 Checksum = 0;
		U64 Count = 0;
		U64 Elapsed;
		U64 Start = ReadOSTime();
		do
		{
			// Keep evaluations inside each timed batch.
			__asm__ volatile("" : : "m"(Inputs), "m"(Q) : "memory");
			for (U64 Index = 0; Index < BatchSize; ++Index)
			{
				Checksum += MLQuadEval(Inputs[Index], &Q);
			}
			Count += BatchSize;
			Elapsed = ReadOSTime() - Start;
		} while (Elapsed < MinNanos);

		if (!isfinite(Warmup) || !isfinite(Checksum))
		{
			fprintf(stderr, "Nonfinite result for alpha=%g, beta=%g\n", Pair[0], Pair[1]);
			return 1;
		}
		F64 Seconds = (F64)Elapsed*1e-9;
		const char *Path = Pair[0] == 1.0 ? (Pair[1] == 1.0 ? "exp" : "expm1") : "quadrature";
		printf("%5.2f %5.2f %10s %14.3f %12.3f %14llu %9.4f %18.10e\n",
			Pair[0], Pair[1], Path, (F64)Count/Seconds, (F64)Elapsed/(F64)Count,
			(unsigned long long)Count, Seconds, Checksum);
		fflush(stdout);
	}
	return 0;
}
