#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <arm_neon.h>

#include "DataHandler.h"
#include "ImageProcessor.h"

using namespace std;

int main(int argc, char *argv[])
{
	int result_code = 0;
	bool debug = true;

	if (debug) {
		std::cout << "FIRST IMAGING TEST" << std::endl;
	}
	
	iQoLi::DataHandler* dHand = new iQoLi::DataHandler();
	iQoLi::ImageProcessor* iProc = new iQoLi::ImageProcessor();

	int16_t* display_buffer = (int16_t*)malloc(320 * 240 * sizeof(int16_t));
	std::ifstream file("/usr/bin/data/main.bin", std::ifstream::binary);
	if (true == file.is_open()){
		if (file){
			file.seekg(0, file.end);
			int length = (int)file.tellg();
			file.seekg(0, file.beg);
			file.read((char*)display_buffer, length);
			file.close();
		}
	}

	std::string file_dir = "/usr/bin/data/iqdata_v9_211021.bin";

	result_code = dHand->load_bin_file(file_dir, debug);
	result_code = dHand->split_iq_data(debug);
	result_code = dHand->convert_neon(debug);

	result_code = iProc->normal_process_iq(debug, dHand->raw_data_I, dHand->raw_data_Q);
	result_code = iProc->normal_log_compression(debug);
	result_code = iProc->process_iq(debug, dHand->neon_data_I, dHand->neon_data_Q);
	result_code = iProc->log_compression(debug);
	result_code = iProc->scan_conversion(debug);


    std::ofstream fp;
    fp.open("/usr/bin/data/result_image.bin", ios::out | ios :: binary );
    fp.write((char*)iProc->result_image, sizeof(iProc->result_image));
	fp.close();


	delete dHand;
	delete iProc;
	
	return 0;
}
