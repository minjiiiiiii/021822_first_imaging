#include "ImageProcessor.h"
#include "NeonMath.h"
#include "Macros.h"

using namespace iQoLi;

ImageProcessor::ImageProcessor(){
    neon_data_signed = (int16x4_t*)malloc(512 * 256 * sizeof(int16x4_t));
    neon_data_unsigned = (uint16x4_t*)malloc(512 * 256 * sizeof(uint16x4_t));
    //neon_expanded = (uint32x4_t*)malloc(512 * 256 * sizeof(uint32x4_t));
    neon_iq_data = (uint32x4_t*)malloc(512 * 256 * sizeof(uint32x4_t)); 
    neon_log_data = (uint16x4_t*)malloc(512 * 256 * sizeof(uint16x4_t));   //vectorized log data 


    result_image = (uint8_t*)malloc(320 * 240 * sizeof(uint8_t));
   
    non_neon_data_unsigned = (float*)malloc(512 * 256 * sizeof(float));
    non_neon_data_unsigned_xsqrt = (float*)malloc(512 * 256 * sizeof(float));
    non_neon_log_data = (uint16_t*)malloc(512 * 256 * sizeof(uint16_t));
    non_neon_log_data_xnor = (uint16_t*)malloc(512 * 256 * sizeof(uint16_t)); //normalization x 
    neon_iq_stored_data =  (uint32_t*)malloc(512 * 256 * sizeof(uint32_t));  //neon data-> store
    neon_stored_data =  (uint16_t*)malloc(512 * 256 * sizeof(uint16_t));     //neon data-> store
    err_iq =  (uint32_t*)malloc(512 * 256 * sizeof(uint32_t));               //for error calculation 
    err_log =  (uint16_t*)malloc(512 * 256 * sizeof(uint16_t));              //for error calculation


    sc_params = generate_Param_SC();
}

ImageProcessor::~ImageProcessor(){

}

int ImageProcessor::normal_process_iq(bool debugFlag, int16_t* data_i, int16_t* data_q) {
    clock_t start, end;

    if (debugFlag) {
        start = clock();
        std::cout << "------------------------------" << std::endl;
        std::cout << "NON SIMD PREPROCESSING IQ START..." << std::endl;
        
    }
                                                                                                                          
    int raw_i, raw_q;
   

    for (int i=0; i < 512 * 256; i++) {
        raw_i = data_i[i] * data_i[i];
        raw_q = data_q[i] * data_q[i];
        non_neon_data_unsigned_xsqrt[i]= raw_i + raw_q; //execution time -> delete this
        non_neon_data_unsigned[i] = std::sqrt(raw_i + raw_q);
        
    }
    
    if (debugFlag) {
        std::cout << "NON SIMD PREPROCESSING IQ SUCCESS" << std::endl;
        end = clock();
        std::cout << "EXECUTION TIME : " << ((double)(end - start) / CLOCKS_PER_SEC) * 1000 << " ms" << std::endl;
        std::cout << "------------------------------" << std::endl;
    }
	return 0;
}


int ImageProcessor::normal_log_compression(bool debugFlag) {
    clock_t start, end;

    if (debugFlag) {
        start = clock();
        std::cout << "------------------------------" << std::endl;
        std::cout << "NON SIMD LOG COMPRESSION START..." << std::endl;
    }
   
    int num_iter = (512 * 256);

    float max = 0;

    for(int i=0; i<num_iter; i++)
	{
        if (max < non_neon_data_unsigned[i]) max = non_neon_data_unsigned[i];
	}

    for(int i=0; i<num_iter; i++)
	{
        non_neon_log_data[i] = (uint16_t) (20* (std::log10( (std::abs(non_neon_data_unsigned[i]) / max) + 0.001)) + 60.0);
        //non_neon_log_data_xnor[i] = (uint16_t) (20* (std::log10( (std::abs(non_neon_data_unsigned[i])) + 0.001)) );


	}

    if (debugFlag) {
        std::cout << "NON SIMD LOG COMPRESSION SUCCESS" << std::endl;
        end = clock();
        std::cout << "EXECUTION TIME : " << ((double)(end - start) / CLOCKS_PER_SEC) * 1000 << " ms" << std::endl;
        std::cout << "------------------------------" << std::endl;
    }

    return 0;
}

int ImageProcessor::process_iq(bool debugFlag, int16x4_t* data_i, int16x4_t* data_q) {
    clock_t start, end;

    if (debugFlag) {
        start = clock();
        std::cout << "------------------------------" << std::endl;
        std::cout << "PREPROCESSING IQ START..." << std::endl;
    }

    int num_iter = (512 * 256) >>2;

    uint16x4_t u16I, u16Q;
	uint32x4_t u32I, u32Q;
   
    
    for(int i=0; i<num_iter; i++)
	{
		u16I = vreinterpret_u16_s16(vabs_s16(data_i[i]));
        u16Q = vreinterpret_u16_s16(vabs_s16(data_q[i]));
        u32I = vmull_u16(u16I, u16I);
        u32I = vmlal_u16(u32I, u16Q, u16Q);
        neon_iq_data[i]=u32I;

    }
   

    if (debugFlag) {
        std::cout << "PREPROCESSING IQ SUCCESS" << std::endl;
        end = clock();
        std::cout << "EXECUTION TIME : " << ((double)(end - start) / CLOCKS_PER_SEC) * 1000 << " ms" << std::endl;
        std::cout << "------------------------------" << std::endl;
    }
	return 0;
}


int ImageProcessor::log_compression(bool debugFlag) {
    clock_t start, end;
    float scale = float(Param_Log_ScaleupFactor);
    if (debugFlag) {
        start = clock();
        std::cout << "------------------------------" << std::endl;
        std::cout << "LOG COMPRESSION START..." << std::endl;
        std::cout << "LOG SCALE UP FACTOR : " <<  Param_Log_ScaleupFactor * Param_Log_Nat2Ten << std::endl;
        std::cout << "scale : " << scale << std::endl; //scale : 157.534
    }
  
    uint32x4_t tmp_uint32;
	float32x4_t tmp_float32; 
    
    /*
    float32x4_t v_test =vmulq_n_f32((log_ps(vdupq_n_f32((float32_t)12))),float32_t(0.43429)); // log_ps-> ln
    float32_t* v_testP; 
    v_testP =  (float32_t*)malloc(sizeof(float32_t));
    */

	uint32x4_t v_one = vdupq_n_u32((uint32_t)1);

	int num_iter = (512 * 256) >> 2;

	for(int i=0; i<num_iter; i++)
	{
		tmp_uint32 = vaddq_u32(v_one, neon_iq_data[i]);
        tmp_float32 =log_ps(vcvtq_f32_u32(tmp_uint32));  
        //tmp_float32 = 10*log_ps(vcvtq_f32_u32(tmp_uint32)); 
        //tmp_float32 = vmulq_n_f32(tmp_float32, scale);  //-> max error of log data :64451
        //tmp_float32 = vmulq_n_f32(tmp_float32, float32_t(0.43429));
        tmp_float32 = vmulq_n_f32(tmp_float32, float32_t(4.3429));
        tmp_uint32 = vcvtq_u32_f32(tmp_float32);
        neon_log_data[i] = vmovn_u32(tmp_uint32);
       
	}


    if (debugFlag) {
        
        std::cout << "SIMD LOG COMPRESSION SUCCESS" << std::endl;
        end = clock();
        std::cout << "EXECUTION TIME : " << ((double)(end - start) / CLOCKS_PER_SEC) * 1000 << " ms" << std::endl;
        std::cout << "------------------------------" << std::endl;
        
       /*
        vst1q_f32(v_testP,v_test); 
         for (int i=0; i<4; i=i+1) {
            
            std::cout << "log test" << v_testP[i] <<std::endl;
        }
        */
        

    }
	return 0;
}

int ImageProcessor::compare_result(bool debugFlag) {
    
    if(debugFlag) {
        
        
        //neon_log_data[i] = vmovn_u32(tmp_uint32);
        int num_iter = (512 * 256); 
        int err; 
        for (int i=0; i< num_iter; i=i+4) {
            vst1q_u32(&neon_iq_stored_data[i], neon_iq_data[i/4]);
        }

        for (int i=0; i< num_iter; i=i+4) {
            vst1_u16(&neon_stored_data[i], neon_log_data[i/4]);
        }
        std::cout << "------------------------------" << std::endl;
        std::cout << "Comparing NEON/NON NEON IQ DATA" << std::endl; 
        std::cout << "------------------------------" << std::endl;
        
        std::cout << "1st NON NEON IQ DATA: " << non_neon_data_unsigned_xsqrt[0] <<std::endl;
        std::cout << "2nd NON NEON IQ DATA: " << non_neon_data_unsigned_xsqrt[1] <<std::endl;
        std::cout << "3rd NON NEON IQ DATA: " << non_neon_data_unsigned_xsqrt[2] <<std::endl;
        std::cout << "4th NON NEON IQ DATA: " << non_neon_data_unsigned_xsqrt[3] <<std::endl; 

        std::cout << "------------------------------" << std::endl;
        std::cout << "1st NEON IQ DATA: " << neon_iq_stored_data[0] <<std::endl;
        std::cout << "2nd NEON IQ DATA: " << neon_iq_stored_data[1] <<std::endl;
        std::cout << "3rd NEON IQ DATA: " << neon_iq_stored_data[2] <<std::endl;
        std::cout << "4th NEON IQ DATA: " << neon_iq_stored_data[3] <<std::endl; 

        std::cout << "------------------------------" << std::endl;
        std::cout << "101st NON NEON IQ DATA: " << non_neon_data_unsigned_xsqrt[100] <<std::endl;
        std::cout << "102nd NON NEON IQ DATA: " << non_neon_data_unsigned_xsqrt[101] <<std::endl;
        std::cout << "103rd NON NEON IQ DATA: " << non_neon_data_unsigned_xsqrt[102] <<std::endl;
        std::cout << "104th NON NEON IQ DATA: " << non_neon_data_unsigned_xsqrt[103] <<std::endl; 
        std::cout << "------------------------------" << std::endl;
        std::cout << "101st NEON IQ DATA: " << neon_iq_stored_data[100] <<std::endl;
        std::cout << "102nd NEON IQ DATA: " << neon_iq_stored_data[101] <<std::endl;
        std::cout << "103rd NEON IQ DATA: " << neon_iq_stored_data[102] <<std::endl;
        std::cout << "104th NEON IQ DATA: " << neon_iq_stored_data[103] <<std::endl; 

        std::cout << "------------------------------" << std::endl;
        std::cout << "Comparing NEON/NON NEON LOG COMPRESSION DATA" << std::endl; 
        std::cout << "------------------------------" << std::endl;
        std::cout << "1st NON NEON LOG COMPRESSION DATA : " <<  non_neon_log_data[0] <<std::endl; 
        std::cout << "2nd NON NEON LOG COMPRESSION DATA : " <<  non_neon_log_data[1] <<std::endl;
        std::cout << "3rd NON NEON LOG COMPRESSION DATA : " <<  non_neon_log_data[2] <<std::endl;
        std::cout << "4th NON NEON LOG COMPRESSION DATA : " <<  non_neon_log_data[3] <<std::endl;
       
        std::cout << "------------------------------" << std::endl;
        std::cout << "1st NEON LOG COMPRESSION DATA : " << neon_stored_data[0] <<std::endl; 
        std::cout << "2nd NEON LOG COMPRESSION DATA : " <<  neon_stored_data[1] <<std::endl;
        std::cout << "3rd NEON LOG COMPRESSION DATA : " <<  neon_stored_data[2] <<std::endl;
        std::cout << "4th NEON LOG COMPRESSION DATA : " <<  neon_stored_data[3] <<std::endl;

        std::cout << "------------------------------" << std::endl;
        std::cout << "101th NON NEON LOG COMPRESSION DATA : " << non_neon_log_data[100] <<std::endl;
        std::cout << "102nd NON NEON LOG COMPRESSION DATA : " <<  non_neon_log_data[101] <<std::endl;
        std::cout << "103rd NON NEON LOG COMPRESSION DATA : " <<  non_neon_log_data[102] <<std::endl;
        std::cout << "104th NON NEON LOG COMPRESSION DATA : " <<  non_neon_log_data[103] <<std::endl;

        std::cout << "------------------------------" << std::endl;
        std::cout << "101th NEON LOG COMPRESSION DATA : " << neon_stored_data[100] <<std::endl;
        std::cout << "102nd NEON LOG COMPRESSION DATA : " << neon_stored_data[101] <<std::endl;
        std::cout << "103rd NEON LOG COMPRESSION DATA : " << neon_stored_data[102] <<std::endl;
        std::cout << "104th NEON LOG COMPRESSION DATA : " << neon_stored_data[103] <<std::endl;
     
       uint32_t max32 = 0; 
       uint16_t max16 = 0;
      
       for (int i=0; i< num_iter; i=i+4) {
           err_iq[i] = ((uint32_t)non_neon_data_unsigned_xsqrt[i]) - neon_iq_stored_data[i];
           err_log[i] =  non_neon_log_data[i] -  neon_stored_data[i];
                if (max32 < err_iq[i]) max32 = err_iq[i];
                if (max16 < err_log[i]) max16 = err_log[i];
               
        }
        std::cout << "max error of iq data :" << max32 <<std::endl;
        std::cout << "max error of log data :" << max16 <<std::endl;
        

    }
    return 0; 


}


int ImageProcessor::scan_conversion(bool debugFlag) {
    clock_t start, end;

    if (debugFlag) {
        start = clock();
        std::cout << "------------------------------" << std::endl;
        std::cout << "SCAN CONVERSION START..." << std::endl;
    }

    memset(result_image, 0, 320 * 240);

    int* X_real;
    int* Z_real;
    int* X_rough;
    int* Z_rough;
    int addr_add;
    int num_iter = (512 * 256) >> 2;


    uint16_t* input = non_neon_log_data;
    // uint16_t* input = (uint16_t*)malloc(512 * 256 * sizeof(uint16_t));

    // for (int i=0; i < num_iter; i++) {
    //     vst1_u16((input + 4 * i), *(neon_data_unsigned + i));
    // }

    X_real = sc_params.X_Real;
    Z_real = sc_params.Z_Real;
    X_rough = sc_params.X_Rough;
    Z_rough = sc_params.Z_Rough;

    uint rough_scale_factor = (1<<Param_SC_RoughScaleBit);
    int doubleroubit = (Param_SC_RoughScaleBit << 1);
    int GrayScaleBit = Param_SC_InputBitResol - 8;

    uint temp;
    uint MM = 0;

	int realX, realZ, roughX, roughZ;

    for (int pix_z = 0; pix_z < sc_params.LUT_Znum; pix_z++){
		for (int pix_x = sc_params.LUT[0][pix_z]+1; pix_x < sc_params.LUT[1][pix_z] ; pix_x ++ ){

			addr_add = Param_SC_NumPix_Z* (pix_x-1) + pix_z;

			realX =  *(X_real + addr_add);
			realZ =  *(Z_real + addr_add);
			roughX = *(X_rough + addr_add);
			roughZ = *(Z_rough + addr_add);

			temp = ((*(input+(realX-1)*512 + realZ  ))
					+(*(input+(realX  )*512 + realZ  ))
					+(*(input+(realX-1)*512 + realZ+1))
					+(*(input+(realX  )*512 + realZ+1)));
			// temp = (temp>>(doubleroubit+GrayScaleBit+1));

			MM = MAX(temp, MM);
            *(result_image+addr_add) = (uint8_t)(temp);
			//*(result_image+addr_add) = (uint8_t)((temp & 0x000000FF));
		}
	}
	float factor = (float)255/(float)MM;

	for (int pix_z = 1 ; pix_z <sc_params.LUT_Znum; pix_z ++){
		for (int pix_x = sc_params.LUT[0][pix_z] ; pix_x < sc_params.LUT[1][pix_z]; pix_x++){
			addr_add = Param_SC_NumPix_Z* (pix_x-1) + pix_z;
			*(result_image+addr_add) = *(result_image+addr_add) * factor;
		}
	}

    if (debugFlag) {
        std::cout << "SCAN CONVERSION SUCCESS" << std::endl;
        end = clock();
        std::cout << "EXECUTION TIME : " << ((double)(end - start) / CLOCKS_PER_SEC) * 1000 << " ms" << std::endl;
        std::cout << "------------------------------" << std::endl;
    }

    return 0;
}


str_SC ImageProcessor::generate_Param_SC() {
    str_SC param_sc;

    int num_samples = 512;
    int num_scalines = 256;
    int rough_scale_factor = (1<<Param_SC_RoughScaleBit);

    float center_x = (float)(Param_SC_NumPix_X-1)/2;
    float center_z = 0;
    double sample_space = ((double)Param_SC_View_Depth / (num_samples - 1));
	double pix_len = (double) (SQRT2 * Param_SC_View_Depth / Param_SC_NumPix_X);
    float z_num = (float)Param_SC_View_Depth/pix_len;

    param_sc.LUT_Znum = z_num;

	int leftcounter, validcounter;

	double coord_x, coord_z;
	double theta, distance;
	double idx_sc, idx_samp;

	int real_sc, real_samp, rough_sc, rough_samp;

	int pix_z, pix_x;
	for (pix_z = 0 ; pix_z < 240 ; pix_z ++){
		leftcounter = 0;
		validcounter = 0;
		for (pix_x = 0 ; pix_x < 320 ; pix_x ++){
			coord_x = pix_x - center_x;
			coord_z = pix_z - center_z;

			theta = atan2(coord_x, coord_z) + Param_SC_View_Angle/2;
			distance = double(sqrt(coord_x * coord_x + coord_z * coord_z) * pix_len);

			idx_sc = (double)(num_scalines - 1) * theta / Param_SC_View_Angle;
			idx_samp = double(distance / sample_space);

			real_sc = int(idx_sc);
			real_samp = int(idx_samp);
			rough_sc = int((idx_sc - real_sc)*rough_scale_factor);
			rough_samp = int((idx_samp - real_samp)*rough_scale_factor);

			param_sc.X_Real[pix_x*Param_SC_NumPix_Z + pix_z] = real_sc;
			param_sc.Z_Real[pix_x*Param_SC_NumPix_Z + pix_z] = real_samp;
			param_sc.X_Rough[pix_x*Param_SC_NumPix_Z + pix_z] = rough_sc;
			param_sc.Z_Rough[pix_x*Param_SC_NumPix_Z + pix_z] = rough_samp;

			if (((real_sc < 1) && (real_samp < num_samples)) || ((1 <= real_sc) && (real_samp > num_samples) && (validcounter == 0)))
				leftcounter += 1;
			if ((1 <= real_sc) && (real_sc < num_scalines) && (real_samp < num_samples))
				validcounter += 1;
		}
		param_sc.LUT[0][pix_z] = leftcounter + 1;
		param_sc.LUT[1][pix_z] = leftcounter + validcounter -1;
	}
	return param_sc;
}