/*
 * Macro.h
 *
 *  Created on: 2020. 3. 26.
 *      Author: Hyeyun
 */

#ifndef SRC_MACRO_H_
#define SRC_MACRO_H_


#define Byte_16		2
#define Byte_32 	4
#define Byte_64 	8
#define pi			3.141592653589793
#define SQRT2 		1.414213562373095
#define MAX(x, y) 	(x>y)? x:y
#define MIN(x, y) 	(x<y)? x:y
#define MAXEQ(x, y)	(x>=y)? x:y
#define MINEQ(x, y)	(x<=y)? x:y
///////////////////////////////////////////////////////////////////////////////////////////////
/* -------------------------------------------------------------------------- */
/*                         System Physical Parameters                         */
/* -------------------------------------------------------------------------- */
#define Param_Pre_ELEMENT_PTICH  0.3e-3L
#define Param_Pre_CENTER_FREQ    2500000
#define Param_Pre_CHANNEL_NUM    32
#define Param_Pre_F_NUMBER       1
#define Param_Pre_TX_FOCAL_DEPTH 7e-3L

/* -------------------------------------------------------------------------- */
/*                               Pre Processing                               */
/* -------------------------------------------------------------------------- */
// Initial ADC Settings
#define Param_ADC_Spi_Clock_Div 40
#define Param_ADC_Dout_Cross    1
#define Param_ADC_DTGC_Mode     0b01  // RAMP
//#define Param_ADC_DTGC_Mode          0b10  // Fixed Gain
#define Param_ADC_DTGC_Profile       0b00
#define Param_ADC_DTGC_Pos_Gain      3
#define Param_ADC_DTGC_Pos_Step_Time 1
#define Param_ADC_DTGC_Neg_Gain      0xFF
#define Param_ADC_DTGC_Neg_Step_Time 1
#define Param_ADC_DTGC_Manual_Gain   200
#define Param_ADC_DTGC_Start_Time    0
#define Param_ADC_DTGC_Start_Gain    0
#define Param_ADC_DTGC_Stop_Time     319
#define Param_ADC_DTGC_Stop_Gain     319
#define Param_ADC_Tx_Trig            16
#define Param_ADC_High_Power         0
#define Param_ADC_LNA_Prog           3
#define Param_ADC_LNA_Dis            0
#define Param_ADC_LPF_Prog           1
#define Param_ADC_DIG_EN             0
#define Param_ADC_Res_Sel            0
#define Param_ADC_Flip_Atten         0
#define Param_ADC_Disable_Atten      0

#define Param_ADC_Fix_Atten 0

#define Param_ADC_Fix_Atten_Enable_Profile0 0
#define Param_ADC_Fix_Atten_Value_Profile0  0
#define Param_ADC_Fix_Atten_Enable_Profile1 0
#define Param_ADC_Fix_Atten_Value_Profile1  0
#define Param_ADC_Fix_Atten_Enable_Profile2 0
#define Param_ADC_Fix_Atten_Value_Profile2  0
#define Param_ADC_Fix_Atten_Enable_Profile3 0
#define Param_ADC_Fix_Atten_Value_Profile3  0

// Initial Pulser Settings
#define Param_Pulser_Freq         80000000
#define Param_Pulser_PLL_EN       1
#define Param_Pulser_Wave_ProFile 1
#define Param_Pulser_Invert       0
#define Param_Pulser_Discharge    1
#define Param_Pulser_Delay_Only   1
#define Param_Pulser_REP   		  6 // REP * half cycle

// Initial Bladder RTC Settings
#define Param_RTC_Time_Interval 128
#define Param_RTC_Delay_Count   128
#define Param_RTC_LNA_Count     64
#define Param_RTC_TGC_Start     64
#define Param_RTC_Decimation    1

// Initial RxBF Settings
#define Param_RxBF_Aperture_Open     1
#define Param_RxBF_Apodization_Coeff 32
#define Param_RxBF_Apodization_Open  1

/* -------------------------------------------------------------------------- */
/*                               Mid Processing                               */
/* -------------------------------------------------------------------------- */

//BF
//#define BF_InputBit				15
//#define Param_Init_NScanlines	128
//#define Param_BF_SamplingFreq 	40000000 //10MHz
//#define Param_BF_CenterFreq	  	2500000 // 2.5 MHz
//#define Param_BF_NScanlines	  	128
//#define Param_BF_NSamples   	10240
//#define Param_BF_C 			  	1540000 //[mm/s]
//#define Param_BF_ViewDepth		200		//[mm]


// BF
#define Param_BF_SamplingFreq 40000000  // 40MHz
#define Param_BF_CenterFreq   2500000   // 2.5 MHz
#define Param_BF_NScanlines   128
#define Param_BF_NSamples     10240
#define Param_BF_C            1540
#define Param_BF_ViewDepth		200		//[mm]


//DCR
#define Param_DCR_intBit	 	14	// Resolution (Q1.14) total 15 bits
#define Param_DCR_FilterLength 	65
#define DCR_HalfFL 				(Param_DCR_FilterLength>>1)
#define Param_DCR_CutoffFreqOrg 2000000 // for 0.05
//#define Param_DCR_CutoffFreqOrg 500000 // 500kHz

//PD
#define Param_PD_intBit			12
#define Param_PD_FilterLength 	33
#define Param_PD_DecimFactor	3

//RFI
#define Param_RFI_intBit 		12
#define Param_RFI_InterpFactor 	1
#define Param_RFI_PPFLength 	5
#define Param_RFI_AddNumSc		2
#define Param_RFI_NewNScanlines	(Param_RFI_InterpFactor*Param_BF_NScanlines)

//QBP
#define Param_QBP_intBit 		 12	// Resolution (Q1.12) total 13 bits
#define Param_QBP_cosintBit 	 11 // Resolution (Q1.11) total 12 bits
#define Param_QBP_NumberofZone 	 1
#define Param_QBP_FilterLength 	 33
#define QDM_HalfFL 				 (Param_QBP_FilterLength>>1)

#define Param_QBP_FreqStart 	 2500000 //
#define Param_QBP_FreqEnd 		 2500000 //
#define Param_QBP_BandwidthStart 2500000 //
#define Param_QBP_BandwidthEnd 	 2500000 //
#define Param_QBP_DepthStart 	 0
#define Param_QBP_NumElemPeriod (int)((float)Param_BF_SamplingFreq/(float)Param_BF_CenterFreq)

//QDM
#define Param_QDM_CosIter		(4+(Param_BF_NSamples+DCR_HalfFL)*Param_RFI_NewNScanlines-DCR_HalfFL)

//////////////////////////////////////////////////////////////////////////////////////
//LOG
#define Param_Log_OutputBit			15 // Unsigned -2021.06.07
#define Param_Log_scaleupBit		6 // for fitting in 12 bit output of LOG
#define Param_Log_scaleup 			(int)(1<<((int)Param_Log_scaleupBit)) //2^(log_scaleupbit -1)
#define Param_Log_MaxdB				(10*20.8)//log(2^30-1)//22 //10*loge(2^32-1) --> log(2^32-1)
#define Param_Log_ScaleupFactor		((float)((1<<Param_Log_OutputBit)-1) / Param_Log_MaxdB)
#define Param_Log_Ten2Nat			2.3026
#define Param_Log_Nat2Ten			0.4343
//#define Param_Log_OutputBit			13 // Unsigned

//ABHF
#define Param_ABHF0_Offset_double 	0.3
#define Param_ABHF0_ScaleupFactor	Param_Log_ScaleupFactor
#define Param_ABHF0_Ten2Nat			Param_Log_Ten2Nat

//Gain
#define Param_Gain_Unit 			0.5
#define Param_Gain_Default 			45// 25
#define Param_Gain_User 			45

//DCC
#define Param_DCC_DRHigh			50
#define Param_DCC_DRCompr 			45

//#define Param_DCC_ScaleBit 			12
#define Param_DCC_ScaleBit 			15  // unsigned 2021/06/08
#define Param_DCC_Dynamic_StartIdx 	8000
#define Param_DCC_Dynamic_EndIdx 	10000
//#define Param_DCC_Dynamic_StartIdx 	2200
//#define Param_DCC_Dynamic_EndIdx 	2500

//LBHF
#define Param_LBHF_Offset 			0.3
#define Param_LBHF_ScaleupFactor	Param_Log_ScaleupFactor
#define Param_LBHF_Ten2Nat			Param_ABHF0_Ten2Nat

//AEEF
#define Param_AEEF_intBit	 		 13 // filter coefficient scale-up bit
//#define Param_AEEF_intBit	 		 15 // 2021/06/08
#define Param_AEEF_HPFLength 		 17
#define Param_AEEF_CutFreqNorm		 0.05
#define Param_AEEF_EdgeFilterFilling 4	// 4: Default
#define Param_AEEF_AEEFRatio 		 0.2
#define Param_AEEF_RatioTotStepBit	 5 // 2^5 = 32
#define Param_AEEF_RatioSelect		 7  // 7/32 = 0.218 close to 0.2(AEEFRatio)
#define Param_AEEF_AEEFRatioBit 	 2
#define Param_AEEF_EdgeFilOrgLeng	 3

//LSF
//#define Param_LSF_intBit 	   		 12
#define Param_LSF_intBit 	   		 15  // 2021/06/08
#define Param_LSF_Opt 		  		 1
#define Param_LSF_FilterLength 		 3


//REJ
#define Param_REJ_RejectThreshold 	2
#define Param_REJ_ScaleBit 		  	Param_DCC_ScaleBit
//#define Param_REJ_OutputBitResol  	12
#define Param_REJ_OutputBitResol  	15 // 2021.06.07
#define Param_REJ_BitShift			(16-Param_REJ_OutputBitResol)

//SC
#define Param_SC_NumPix_X 			320
#define Param_SC_NumPix_Z 	 		240
#define Param_SC_NumPix_ZH 	 		240
#define Param_SC_ScreenSize 	 	Param_SC_NumPix_X * Param_SC_NumPix_Z
#define Param_SC_RoughScaleBit 	 	8
#define Param_SC_InputBitResol	 	Param_REJ_OutputBitResol
#define Param_SC_View_Depth 	 	Param_BF_ViewDepth //[mm]
#define Param_SC_View_Angle 		1.570796326794897
#define Frame_Num 			   		10


///////////////////////////////////////////////////////////////////////////////////////////////
// Array Size Setting

#define Total_ArrSize	((Param_BF_NSamples + DCR_HalfFL)*Param_Init_NScanlines + DCR_HalfFL)

#define InputArraySize 	((Param_BF_NSamples + DCR_HalfFL)*Param_BF_NScanlines + DCR_HalfFL)
#define MaxArraySize 	((Param_BF_NSamples + Param_DCR_FilterLength)*Param_RFI_NewNScanlines)
#define H_MaxArraySize	(MaxArraySize>>1)

#define MaxArray_Aeef 	(((Param_BF_NSamples + Param_AEEF_HPFLength+1)*Param_RFI_NewNScanlines)>>2)
#define SC_ScreenSize 	(Param_SC_NumPix_X * Param_SC_NumPix_Z)
#define B_EArraySize	(Param_BF_NSamples*Param_RFI_NewNScanlines)//back end array size
#define B_EArraySizeD	(Param_BF_NSamples*Param_RFI_NewNScanlines*2)//back end array size
#define QuartBE_ArrSize	(B_EArraySize>>2)
///////////////////////////////////////////////////////////////////////////////////////////////
// LOG value setting

#define c_inv_mant_mask 	~0x7f800000u
#define c_cephes_SQRTHF 	0.707106781186547524
#define c_cephes_log_p0 	7.0376836292E-2
#define c_cephes_log_p1 	-1.1514610310E-1
#define c_cephes_log_p2 	1.1676998740E-1
#define c_cephes_log_p3 	-1.2420140846E-1
#define c_cephes_log_p4 	+1.4249322787E-1
#define c_cephes_log_p5 	-1.6668057665E-1
#define c_cephes_log_p6 	+2.0000714765E-1
#define c_cephes_log_p7 	-2.4999993993E-1
#define c_cephes_log_p8 	+3.3333331174E-1
#define c_cephes_log_q1 	-2.12194440e-4
#define c_cephes_log_q2 	0.693359375
///////////////////////////////////////////////////////////////////////////////////////////////
// address in DDR3 ////////////////////////////////////////////////////////////////////////////

// restore data
#define addr_data_BF	0x10000000
//#define addr_data_BF	0x10100000
//#define addr_data_SC 	0x16C00000

#define BlockMemSize 	0x00800000
#define BlockMemSizeH	(0x00800000>>1)

#define addr_Vec_CPU0_ping	addr_data_BF 		+ BlockMemSize
#define addr_Vec_CPU0_pong	addr_Vec_CPU0_ping 	+ BlockMemSize
#define addr_Vec_CPU0_pang	addr_Vec_CPU0_pong 	+ BlockMemSize

#define addr_Vec_CPU1_ping	addr_Vec_CPU0_ping 	+ (BlockMemSize>>1)
#define addr_Vec_CPU1_pong	addr_Vec_CPU0_pong 	+ (BlockMemSize>>1)
#define addr_Vec_CPU1_pang	addr_Vec_CPU0_pang 	+ (BlockMemSize>>1)



#endif /* SRC_MACRO_H_ */
