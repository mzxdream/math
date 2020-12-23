#include <random>
#include <ctime>
#include <iostream>
#include "test_math_util.h"
#include "test_vector3.h"
#include "test_quaternion.h"

//#define OUT_TO_FILE

int main(int argc, char *argv[])
{
#ifdef OUT_TO_FILE
    if (freopen("1.log", "w+", stdout) == nullptr)
    {
        fprintf(stderr, "freopen failed");
        return 0;
    }
#endif
    std::cout << std::numeric_limits<int64_t>::min() << std::endl;
    std::cout << std::abs(std::numeric_limits<int64_t>::min()) << std::endl;
    std::cout << std::boolalpha << std::is_signed_v<bool> << std::endl;
    std::cout.precision(17);
    srand((unsigned)time(NULL));
    TestMathUtil(100000);
    TestVector3(10000);
    TestQuaternion(10000);
#ifdef OUT_TO_FILE
    fclose(stdout);
#endif
    return 0;
}