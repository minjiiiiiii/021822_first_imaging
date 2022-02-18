#include <iostream>
#include <vector>
#include <arm_neon.h>
#include <string>
#include <fstream>
#include <time.h>
#include <cstring>
#include <cmath>

#include "Structures.h"


namespace iQoLi{
	class ImageProcessor{
	public:
		ImageProcessor();
		virtual ~ImageProcessor();

        int normal_process_iq(bool debugFlag, int16_t* data_i, int16_t* data_q);
        int normal_log_compression(bool debugFlag);

        int process_iq(bool debugFlag, int16x4_t* data_i, int16x4_t* data_q);
        int log_compression(bool debugFlag);
        int scan_conversion(bool debugFlag);

        str_SC generate_Param_SC();
    
    public:
		int16x4_t* neon_data_signed;
        uint16x4_t* neon_data_unsigned;
        uint32x4_t* neon_expanded;
        uint8_t* result_image;

        float* non_neon_data_unsigned;
        uint16_t* non_neon_log_data;

        str_SC sc_params;
	};
}
