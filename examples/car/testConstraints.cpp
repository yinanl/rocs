#include <iostream>
#include "src/definitions.h"
#include "car.hpp"

int main() {

    rocs::Rn xr= {1.4, 3, -0.5};
    rocs::Rn yr= cst_v5<rocs::Rn>(xr);

    for (int iter=0; iter<yr.size(); ++iter)
	std::cout << yr[iter] << ", ";
    std::cout << std::endl;

    rocs::ivec xi(3);
    xi[0] = rocs::interval(1.3, 1.5);
    xi[1] = rocs::interval(2.9, 3.1);
    xi[2] = rocs::interval(-0.51, -0.49);

    rocs::ivec yi= cst_v5<rocs::ivec>(xi);
    for (int iter=0; iter<yi.getdim(); ++iter)
	std::cout << yi[iter] << ", ";
    std::cout << std::endl;
    
    return 0;
}
