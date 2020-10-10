#include "src/hdf5io.h"
#include "src/interval_vector.h"

int main()
{
    std::vector<double> x{0.1, 0.2, 0.3, 1.1, 1.2, 1.3};
    rocs::h5FileHandler h5f("testdata.h5");
    std::string rv("x");
    std::string rn("ts");
    h5f.write_real_array(x, rv);

    double ts = 0.004;
    h5f.write_real_number(ts, rn);

    rocs::ivec statespace{rocs::interval(-1., 3.),
			  rocs::interval(0.57, 10.4),
			  rocs::interval(-20, -13)};
    h5f.write_state_space(statespace, std::string("X"));
    return 0;
}
