#define _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_DEPRECATE

#include <ctime>
#include <stdlib.h>
#include <stdio.h>

#include "../include/namespace.h"
#include "../include/multi_threads_run.h"

using namespace deardia;

int main(int argc, char* argv[])
{
	time_t startTime = time(NULL);

	srand(((unsigned int)time(NULL)));

	deardia::MultiThreadsRun(argc, argv);

	float elapseTime = (float)(time(NULL) - startTime) / 60.0f;

	printf("Elapse % .3f minutes.\n", elapseTime);

	return 0;
}