#include <cmath>
#include <cfenv>
#include <iostream>
#include <limits>
#include <cfloat>


int main()
{
    const double PI2IVAL = 6.28318530717958647693;
    int q1, q2, q3;
    double r1, r2, r3;
    char angle[256];
    std::cout << "Enter the angle (rad):\n";
    fgets(angle, 256, stdin);
    // std::cout.precision(65);
    
    r1 = remquo(atof(angle), PI2IVAL, &q1);
    std::cout << r1 << ',' << q1 << '\n';
    
    std::fesetround(FE_DOWNWARD);
    r2 = remquo(atof(angle), PI2IVAL, &q2);
    std::cout << r2 << ',' << q2 << '\n';
    
    std::fesetround(FE_UPWARD);
    r3 = remquo(atof(angle), PI2IVAL, &q3);
    std::cout << r3 << ',' << q3 << '\n';

    std::cout << std::numeric_limits<double>::min() << std::endl;
    std::cout << std::numeric_limits<double>::denorm_min() << std::endl;

    std::cout << DBL_MIN << std::endl;
    // double eta0, eta = std::numeric_limits<double>::epsilon();
    // while (eta > 0.0) {
    // 	eta0 = eta;
    // 	eta /= 2.0;
    // }
    // std::cout << eta0 << std::endl;

    return 1;
}
