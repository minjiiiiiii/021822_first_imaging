#include "DataHandler.h"

using namespace iQoLi;

DataHandler::DataHandler(){
    raw_data_I = (int16_t*)malloc(512 * 256 * sizeof(int16_t)); // TODO : CHANGE THIS PARAMETER RECONFIGURABLE AND ROBUST
    raw_data_Q = (int16_t*)malloc(512 * 256 * sizeof(int16_t));

    neon_data_I = (int16x4_t*)malloc(512 * 256 * sizeof(int16x4_t));
    neon_data_Q = (int16x4_t*)malloc(512 * 256 * sizeof(int16x4_t));
}

DataHandler::~DataHandler(){

}

/*
 * Class : DataHandler
 * Function : load_bin_file()
 * Parameter :
 * Date : 2021-09-03
 * Writer : Hyunwoo Cho
 * Description :
 * 		Load data from .bin format file.
 * 		This function is test-only purpose for PS stand-alone application.
 */

int DataHandler::load_bin_file(std::string fName, bool debugFlag) {
    clock_t start, end;

    if (debugFlag) {
        start = clock();
        std::cout << "------------------------------" << std::endl;
        std::cout << "FILE READ START..." << std::endl;
        std::cout << fName << std::endl;
    }

    std::ifstream file(fName, std::ifstream::binary);
	if (true == file.is_open()){
		if (file){
			file.seekg(0, file.end);
			int length = (int)file.tellg();
			file.seekg(0, file.beg);

			raw_data = (int16_t*)malloc(length / sizeof(int16_t));
			file.read((char*)raw_data, length);
			file.close();
		}
	}
	else {
        if (debugFlag) {
            std::cout << "FILE READ FAIL" << std::endl;
            std::cout << "------------------------------" << std::endl;
        }   
        return -1;
	}

    if (debugFlag) {
        std::cout << "FILE READ SUCCESS" << std::endl;
        end = clock();
        std::cout << "1st DATA I : " << raw_data[0] << std::endl;
        std::cout << "2nd DATA Q : " << raw_data[1] << std::endl;
        std::cout << "3rd DATA I : " << raw_data[2] << std::endl;
        std::cout << "4th DATA Q : " << raw_data[3] << std::endl;
        std::cout << "EXECUTION TIME : " << ((double)(end - start) / CLOCKS_PER_SEC) * 1000 << " ms" << std::endl;
        std::cout << "------------------------------" << std::endl;
    }
	return 0;
}

int DataHandler::split_iq_data(bool debugFlag) {
    clock_t start, end;
    if (debugFlag) {
        start = clock();
        std::cout << "------------------------------" << std::endl;
        std::cout << "DATA SPLIT START..." << std::endl;
    }
    
    for (int i = 1; i < 128+1; i++) { // TODO : CHANGE THIS PARAMETER RECONFIGURABLE AND ROBUST
        for (int j = 1; j < 512+1; j++) {
            if (j < 3) { // TODO : PROCESS HEADER
                raw_data_I[1024*(i-1)+((j-1))] = 0;
                raw_data_Q[1024*(i-1)+((j-1))] = 0;
                raw_data_I[1024*(i-1)+((j-1)+512)] = 0;
                raw_data_Q[1024*(i-1)+((j-1)+512)] = 0;
            }
            else {
                raw_data_I[1024*(i-1)+((j-1))] = raw_data[2048*(i-1) + 4*(j-1)];
                raw_data_Q[1024*(i-1)+((j-1))] = raw_data[2048*(i-1) + 4*(j-1)+1];
                raw_data_I[1024*(i-1)+((j-1)+512)] = raw_data[2048*(i-1) + 4*(j-1)+2];
                raw_data_Q[1024*(i-1)+((j-1)+512)] = raw_data[2048*(i-1) + 4*(j-1)+3];
            }
            
        }
    }

    if (debugFlag) {
        end = clock();
        std::cout << "FIRST DATA I : " << raw_data_I[0] << std::endl;
        std::cout << "FIRST DATA Q : " << raw_data_Q[0] << std::endl;
        std::cout << "100th DATA I : " << raw_data_I[99] << std::endl;
        std::cout << "100th DATA Q : " << raw_data_Q[99] << std::endl;
        std::cout << "512th DATA I : " << raw_data_I[511] << std::endl;
        std::cout << "512th DATA Q : " << raw_data_Q[511] << std::endl;
        std::cout << "513th DATA I : " << raw_data_I[512] << std::endl;
        std::cout << "513th DATA Q : " << raw_data_Q[512] << std::endl;
        std::cout << "DATA SPLIT SUCCESS" << std::endl;
        std::cout << "EXECUTION TIME : " << ((double)(end - start) / CLOCKS_PER_SEC) * 1000 << " ms" << std::endl;
        std::cout << "------------------------------" << std::endl;
    }

    return 0;
}

int DataHandler::convert_neon(bool debugFlag) {
    clock_t start, end;
    if (debugFlag) {
        start = clock();
        std::cout << "------------------------------" << std::endl;
        std::cout << "NEON DATA CONVERSION START..." << std::endl;
    }


    int num_iter = (512*256) >> 2;

    for (int i=0; i < num_iter; i++) {
        neon_data_I[i] = vld1_s16(raw_data_I + 4 * i);
        neon_data_Q[i] = vld1_s16(raw_data_Q + 4 * i);
    }
   
    if (debugFlag) {
        end = clock();
        std::cout << "NEON DATA CONVERSION SUCCESS" << std::endl;
        std::cout << "EXECUTION TIME : " << ((double)(end - start) / CLOCKS_PER_SEC) * 1000 << " ms" << std::endl;
        std::cout << "------------------------------" << std::endl;
    }
    return 0;
}

