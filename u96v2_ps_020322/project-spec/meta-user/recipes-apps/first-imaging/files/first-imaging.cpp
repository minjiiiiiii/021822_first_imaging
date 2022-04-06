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

	result_code = iProc->compare_result(debug);
	result_code = iProc->scan_conversion(debug);

	/*
	// non_neon_data_unsigned[i] -non neon iq data
	std::ofstream non_neon_iq_fp;
	non_neon_iq_fp.open("/usr/bin/data/non_neon_iq_image.bin", ios::out | ios :: binary );
    non_neon_iq_fp.write((char*)iProc->non_neon_data_unsigned, sizeof(*(iProc->non_neon_data_unsigned)) * 512 * 256);
	non_neon_iq_fp.close(); 
	*/


	std::ofstream non_neon_fp;
	non_neon_fp.open("/usr/bin/data/non_neon_log_image.bin", ios::out | ios :: binary );
    non_neon_fp.write((char*)iProc->non_neon_log_data, sizeof(*(iProc->non_neon_log_data)) * 512 * 256);
	non_neon_fp.close(); 
	
	/*
   //non_neon_log_data_xnor
	std::ofstream non_neon_xnor_fp;
	non_neon_xnor_fp.open("/usr/bin/data/non_neon_log_xnor_image.bin", ios::out | ios :: binary );
    non_neon_xnor_fp.write((char*)iProc->non_neon_log_data_xnor, sizeof(*(iProc->non_neon_log_data_xnor)) * 512 * 256);
	non_neon_xnor_fp.close(); 	
	*/


	std::ofstream neon_fp;
	neon_fp.open("/usr/bin/data/neon_log_image.bin", ios::out | ios :: binary );
    neon_fp.write((char*)iProc-> neon_stored_data, sizeof(*(iProc-> neon_stored_data)) *512 * 256);
	neon_fp.close(); 


    std::ofstream fp;
    
	fp.open("/usr/bin/data/result_image.bin", ios::out | ios :: binary );
    fp.write((char*)iProc->result_image, sizeof(*(iProc->result_image)) * 320 * 240);
	fp.close(); 
	
	
	delete dHand;
	delete iProc;
	
	return 0;
}
