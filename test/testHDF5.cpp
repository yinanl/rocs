#include "src/hdf5io.h"
#include "src/interval_vector.h"

int main()
{
    rocs::h5FileHandler h5f("testdata.h5", H5F_ACC_TRUNC);
    
    /* Test write number and arrays */
    double ts = 0.004;
    h5f.write_number<double>(ts, "ts");
    
    rocs::ivec statespace{rocs::interval(-1., 3.),
			  rocs::interval(0.57, 10.4),
			  rocs::interval(-20, -13)};
    h5f.write_state_space(statespace, std::string("X"));
    
    std::vector<double> xdb{0.1, 0.2, 0.3, 1.1, 1.2, 1.3};
    h5f.write_array<double>(xdb, "Xdb");
    
    std::vector<int> xint{0, 3, 5, -18, 9};
    h5f.write_array<int>(xint, "Xint");
    
    std::vector< std::vector<size_t> > xs{{0, 3, 5, 18, 9},
    					  {3, 100, 77, 186, 1},
    					  {0, 0, 0, 0, 0}};
    size_t xarr2d[]{0, 3, 5, 18, 9, 3, 100, 77, 186, 1, 0, 0, 0, 0, 0};
    size_t dim[2]{3, 5};
    h5f.write_2d_array<size_t>(xarr2d, dim, "Xarr2d");
    h5f.write_2d_array<size_t>(xs, "Xsv");


    /* Test read and write arrays */
    rocs::h5FileHandler reader("data_invset_0.8.h5", H5F_ACC_RDONLY);
    std::vector<double> pavings;
    std::vector<int> tags;
    boost::dynamic_bitset<> ctlr;
    size_t pdims[2], cdims[2];
    reader.read_sptree_controller(pavings, pdims, tags, ctlr, cdims);
    
    /* write the data that are just read */
    h5f.write_2d_array<double>(pavings, pdims, "pavings");
    h5f.write_array<int>(tags, "tags");
    std::vector<unsigned char> cc(ctlr.size());
    h5f.write_2d_array<unsigned char>(cc, cdims, "ctlr");
    
    return 0;
}
