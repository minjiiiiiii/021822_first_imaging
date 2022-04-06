#include "NeonMath.h"


float32x4_t log_ps(float32x4_t x) {
	float32x4_t one = vdupq_n_f32(1);
	x = vmaxq_f32(x, vdupq_n_f32(0)); /* force flush to zero on denormal values */
	uint32x4_t invalid_mask = vcleq_f32(x, vdupq_n_f32(0));
	int32x4_t ux = vreinterpretq_s32_f32(x);
	int32x4_t emm0 = vshrq_n_s32(ux, 23);
	ux = vandq_s32(ux, vdupq_n_s32(c_inv_mant_mask));
	ux = vorrq_s32(ux, vreinterpretq_s32_f32(vdupq_n_f32(0.5f)));
	x = vreinterpretq_f32_s32(ux);

	emm0 = vsubq_s32(emm0, vdupq_n_s32(0x7f));
	float32x4_t e = vcvtq_f32_s32(emm0);

	e = vaddq_f32(e, one);

  /* part2:
	 if( x < SQRTHF ) {
	   e -= 1;
	   x = x + x - 1.0;
	 } else { x = x - 1.0; }
  */
	uint32x4_t mask = vcltq_f32(x, vdupq_n_f32(c_cephes_SQRTHF));
	float32x4_t tmp = vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(x), mask));
  	x = vsubq_f32(x, one);
  	e = vsubq_f32(e, vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(one), mask)));
  	x = vaddq_f32(x, tmp);

  	float32x4_t z = vmulq_f32(x,x); // x^2

  	float32x4_t y = vdupq_n_f32(c_cephes_log_p0);   //samples  
  	y = vmlaq_f32(vdupq_n_f32(c_cephes_log_p1), y, x); // y = p1 + y*x
  	y = vmlaq_f32(vdupq_n_f32(c_cephes_log_p2), y, x); // y = p2 + (p1+y*x)*x
  	y = vmlaq_f32(vdupq_n_f32(c_cephes_log_p3), y, x); // y = p3 + y*x
	/*
	y = vmlaq_f32(vdupq_n_f32(c_cephes_log_p4), y, x); // y = p4 + y*x
  	y = vmlaq_f32(vdupq_n_f32(c_cephes_log_p5), y, x); // y = p5+ y*x
	y = vmlaq_f32(vdupq_n_f32(c_cephes_log_p6), y, x); // y = p6 + y*x
  	y = vmlaq_f32(vdupq_n_f32(c_cephes_log_p7), y, x); // y = p7 + y*x
  	y = vmlaq_f32(vdupq_n_f32(c_cephes_log_p8), y, x); // y = p8 + y*x
  	*/
	y = vmulq_f32(y, x); // y = x * y
  	y = vmulq_f32(y, z); // y = y * x^2

  	tmp = vmulq_f32(e, vdupq_n_f32(c_cephes_log_q1));
  	y = vaddq_f32(y, tmp); //y = y + tmp

  	tmp = vmulq_f32(z, vdupq_n_f32(0.5f)); // tmp = z * 0.5 = 0.5 * x^2
  	y = vsubq_f32(y, tmp); // y = y - tmp

  	tmp = vmulq_f32(e, vdupq_n_f32(c_cephes_log_q2));
  	x = vaddq_f32(x, y);
  	x = vaddq_f32(x, tmp);
  	x = vreinterpretq_f32_u32(vorrq_u32(vreinterpretq_u32_f32(x), invalid_mask)); // negative arg will be NAN

  	return x;
}

const float __log10f_rng =  0.3010299957f;

const float __log10f_lut[8] = {
	-0.99697286229624, 		//p0
	-1.07301643912502, 		//p4
	-2.46980061535534, 		//p2
	-0.07176870463131, 		//p6
	2.247870219989470, 		//p1
	0.366547581117400, 		//p5
	1.991005185100089, 		//p3
	0.006135635201050,		//p7
};

float log10f_c(float x)
{
	float a, b, c, d, xx;
	int m;
	
	union {
		float   f;
		int 	i;
	} r;
	
	//extract exponent
	r.f = x;
	m = (r.i >> 23);
	m = m - 127;
	r.i = r.i - (m << 23);
		
	//Taylor Polynomial (Estrins)
	xx = r.f * r.f;
	a = (__log10f_lut[4] * r.f) + (__log10f_lut[0]);
	b = (__log10f_lut[6] * r.f) + (__log10f_lut[2]);
	c = (__log10f_lut[5] * r.f) + (__log10f_lut[1]);
	d = (__log10f_lut[7] * r.f) + (__log10f_lut[3]);
	a = a + b * xx;
	c = c + d * xx;
	xx = xx * xx;
	r.f = a + c * xx;

	//add exponent
	r.f = r.f + ((float) m) * __log10f_rng;

	return r.f;
}


float log10f_neon_hfp(float x)
{
#ifdef __MATH_NEON
	asm volatile (
	
	"vdup.f32		d0, d0[0]				\n\t"	//d0 = {x,x};
	
	//extract exponent
	"vmov.i32		d2, #127				\n\t"	//d2 = 127;
	"vshr.u32		d6, d0, #23				\n\t"	//d6 = d0 >> 23;
	"vsub.i32		d6, d6, d2				\n\t"	//d6 = d6 - d2;
	"vshl.u32		d1, d6, #23				\n\t"	//d1 = d6 << 23;
	"vsub.i32		d0, d0, d1				\n\t"	//d0 = d0 + d1;

	//polynomial:
	"vmul.f32 		d1, d0, d0				\n\t"	//d1 = d0*d0 = {x^2, x^2}	
	"vld1.32 		{d2, d3, d4, d5}, [%1]	\n\t"	//q1 = {p0, p4, p2, p6}, q2 = {p1, p5, p3, p7} ;
	"vmla.f32 		q1, q2, d0[0]			\n\t"	//q1 = q1 + q2 * d0[0]		
	"vmla.f32 		d2, d3, d1[0]			\n\t"	//d2 = d2 + d3 * d1[0]		
	"vmul.f32 		d1, d1, d1				\n\t"	//d1 = d1 * d1 = {x^4, x^4}	
	"vmla.f32 		d2, d1, d2[1]			\n\t"	//d2 = d2 + d1 * d2[1]		

	//add exponent 	
	"vdup.32 		d7, %0					\n\t"	//d7 = {rng, rng}
	"vcvt.f32.s32 	d6, d6					\n\t"	//d6 = (float) d6
	"vmla.f32 		d2, d6, d7				\n\t"	//d2 = d2 + d6 * d7		

	"vmov.f32 		s0, s4					\n\t"	//s0 = s4

	:: "r"(__log10f_rng), "r"(__log10f_lut) 
    : "d0", "d1", "q1", "q2", "d6", "d7"
	);
#endif
}


float log10f_neon_sfp(float x)
{
#ifdef __MATH_NEON
	asm volatile ("vmov.f32 s0, r0 		\n\t");
	log10f_neon_hfp(x);
	asm volatile ("vmov.f32 r0, s0 		\n\t");
#else
	return log10f_c(x);
#endif
};