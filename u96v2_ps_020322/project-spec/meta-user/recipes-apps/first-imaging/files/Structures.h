#pragma once

#include <iostream>
#include <vector>
#include <arm_neon.h>
#include "Macros.h"

struct str_SC{
	uint16_t NumPix_X;
	uint16_t NumPix_Z;
	uint8_t RoughScaleBit;
	uint8_t InputBitResol;
	uint32_t NSamples;
	uint16_t NScanlines;

	int X_Real[Param_SC_ScreenSize];
	int Z_Real[Param_SC_ScreenSize];
	int X_Rough[Param_SC_ScreenSize];
	int Z_Rough[Param_SC_ScreenSize];
	int LUT[2][Param_SC_NumPix_Z];
	int LUT_Znum ;
	uint8_t dBval = 0;

	int16_t* input_ptr;
	int8_t* output_ptr;
	uint16_t* Uinput_ptr;
	uint8_t* Uoutput_ptr;
};

