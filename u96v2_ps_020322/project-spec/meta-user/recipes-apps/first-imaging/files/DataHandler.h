#include <iostream>
#include <vector>
#include <arm_neon.h>
#include <string>
#include <fstream>
#include <time.h>

#include "Structures.h"

namespace iQoLi{
	class DataHandler{
	public:
		DataHandler();
		virtual ~DataHandler();
		int load_bin_file(std::string fName, bool debugFlag);
        int split_iq_data(bool debugFlag);
		int convert_neon(bool debugFlag);
		int load_pl_data();
    
    public:
        int16_t* raw_data;
        int16_t* raw_data_I;
        int16_t* raw_data_Q;

		int16x4_t* neon_data_I;
		int16x4_t* neon_data_Q;
	};
}
