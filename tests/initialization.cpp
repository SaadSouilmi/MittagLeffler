#include <stdio.h>
#include "mittag_leffler.h"
#include "profile.h"
#include "mittag_leffler.cpp"
#include "profile.cpp"

I32
main()
{
	F64 Parameters[][2] = {
		{0.55, 0.55}, {0.55, 1}, {0.55, 1.55},
		{0.75, 0.75}, {0.75, 1}, {0.75, 1.75},
		{0.95, 0.95}, {0.95, 1}, {0.95, 1.95},
	};
	const U64 Count = 65536;
	printf("alpha,beta,init_s,ns_init,ticks_init,floor_ns_init,checksum\n");
	for (auto &Pair : Parameters)
	{
		MLQuad Q;
		for (U64 Index = 0; Index < 512; ++Index)
		{
			MLQuadInit(Pair[0], Pair[1], &Q);
			__asm__ volatile("" : : "m"(Q), "m"(Pair) : "memory");
		}
		F64 Checksum = 0;
		U64 Start = ReadOSTime();
		U64 Ticks = ReadTSC();
		for (U64 Index = 0; Index < Count; ++Index)
		{
			MLQuadInit(Pair[0], Pair[1], &Q);
			__asm__ volatile("" : : "m"(Q), "m"(Pair) : "memory");
			Checksum += Q.AtZero;
		}
		U64 ElapsedTicks = ReadTSC() - Ticks;
		U64 Elapsed = ReadOSTime() - Start;
		Start = ReadOSTime();
		Ticks = ReadTSC();
		for (U64 Index = 0; Index < Count; ++Index)
		{
			__asm__ volatile("" : : "m"(Q), "m"(Pair) : "memory");
			Checksum += Q.AtZero;
		}
		U64 FloorTicks = ReadTSC() - Ticks;
		U64 Floor = ReadOSTime() - Start;
		if (!isfinite(Checksum) || Elapsed <= Floor || ElapsedTicks <= FloorTicks) { return 1; }
		F64 Nanos = (F64)(Elapsed - Floor)/Count;
		printf("%.17g,%.17g,%.9g,%.9g,%.9g,%.9g,%.17g\n",
			Pair[0], Pair[1], 1e9/Nanos, Nanos, (F64)(ElapsedTicks - FloorTicks)/Count,
			(F64)Floor/Count, Checksum);
		fflush(stdout);
	}
	return 0;
}
