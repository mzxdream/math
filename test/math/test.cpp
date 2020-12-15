#include <random>
#include <ctime>
#include <iostream>
#include "test_math.h"

#define OUT_TO_FILE

int main(int argc, char *argv[])
{
#ifdef OUT_TO_FILE
    if (freopen("1.log", "w+", stdout) == nullptr)
    {
        fprintf(stderr, "freopen failed");
        return 0;
    }
#endif
    std::cout.precision(17);
    srand((unsigned)time(NULL));
    TestMath(100000);
#ifdef OUT_TO_FILE
    fclose(stdout);
#endif
    return 0;
}