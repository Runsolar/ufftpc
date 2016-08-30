/*
ufftpc v.1.0
Uno Fast Fourier Transform for any PC projects where need a real time spectrum analisys for speech recognition. 
This opensource program for engineering and scientific purposes.
It uses float type pointers in complex computation with SSE4 instruction.
License: GNU GPL.
Maked by Danijar Wolf, 2016.
*/

#include<stdio.h>
#include<math.h>

#define L 256

#define K L/4
#define LOG2(x) log(x)/log(2)
#define PI 3.14159265359

float  Cosinewave[L];

struct SSEREGS4
{
	float a[4];
	float b[4];
};

struct SSEREGS4 sse_mul4(struct SSEREGS4 X) {
	_asm 
	{
		movups xmm0, X.a    ;  // load a in xmm0
		movups xmm1, X.b    ;  // load b in xmm1
		mulps xmm0, xmm1    ;  // Mul: xmm0 = xmm0 * xmm1
		;  // xmm10 = xmm10 * xmm00
		;  // xmm11 = xmm11 * xmm01
		;  // xmm12 = xmm12 * xmm02
		;  // xmm13 = xmm13 * xmm03					  
		movups X.a, xmm0    ;  // load results from xmm0 int a
	};

	return X;
}

//A test vector
int Vector[256] = {
	2,1,0,0,2,5,3,0,
	-2,0,2,3,4,7,6,
	-5,-8,0,5,1,-6,-7,
	-3,11,24,18,-2,-23,-35,
	-37,-16,14,29,27,15,-8,
	-28,-30,-22,-8,4,7,0,
	2,12,18,11,-2,-12,-10,
	12,30,24,13,21,-5,-46,
	-22,33,40,-1,-27,-35,-10,
	60,91,31,-45,-82,-105,-86,
	-1,71,67,36,2,-35,-47,
	-39,-33,-30,-6,11,5,13,
	42,55,27,-14,-44,-21,46,
	84,70,41,-21,-98,-71,40,
	106,58,-15,-56,-49,35,128,
	99,-31,-125,-154,-130,-32,85,
	98,34,-13,-46,-57,-39,-22,
	-35,-39,-20,-9,9,45,67,
	45,12,-18,-39,-5,66,110,
	102,39,-92,-136,-28,98,107,
	24,-32,-46,-1,93,117,6,
	-112,-155,-137,-53,66,110,50,
	-19,-54,-62,-41,-16,-14,-27,
	-21,-5,5,17,33,32,12,
	1,-11,-1,28,61,75,54,
	-30,-99,-57,37,77,40,-7,
	-34,-19,39,76,34,-45,-94,
	-99,-55,18,63,45,0,-29,
	-36,-24,-10,-10,-17,-18,-9,
	1,14,21,15,1,-5,-6,
	-3,12,25,32,30,1,-36,
	-33,3,28,21,5,-6,-7,
	7,21,14,-8,-26,-30,-20,
	0,17,16,1,-10,-11,-6,
	-2,-1,-4,-6,-4,-1,2,
	6,6,2,-2,-3,-2,1,
	6,9,9};

class Complex
{
public:
	float re;
	float im;   
public:
	Complex(const float r=0, const float i=0) : re(r), im(i) {};
	float real() { return re; };
	float imag() { return im; };          
	~Complex() {}

	//MATH
	Complex operator + (const Complex &c)
	{
		return Complex(re + c.re, im + c.im);
	}

	Complex operator - (const Complex &c)
	{
		return Complex(re - c.re, im - c.im);
	}

	Complex operator * (const Complex &c)
	{
		SSEREGS4 X; 
		X.a[0]=re;
		X.b[0]=c.re;
		X.a[1]=im;
		X.b[1]=c.im;
		X.a[2]=re;
		X.b[2]=c.im;
		X.a[3]=im;
		X.b[3]=c.re;
		X=sse_mul4(X);
		float r = X.a[0] - X.a[1];
		float i = X.a[2] + X.a[3];     
		return Complex(r, i);  
	} 
};

static Complex w(int z, int p)
{ 
	int m = L / p;
	int i = m * z;

	float c = Cosinewave[i];
	float s = Cosinewave[i + K];
	return Complex(c, s);
}

//Init cosine wave data
static void InitCosinewave()
{
	for(int i=0; i<L; i++) 
	{
		Cosinewave[i] = (float)cos(-2*PI*i/L);
	}
}

static void ufftpc(float *Re, float *Im, int N)
{
	if (N == 2)
	{
		float tmp=Re[0];
		Re[0] = Re[0] + Re[1];
		Re[1] = tmp - Re[1];
	}
	else
	{
		int LEN = N >> 1;

		float *Re_even;
		Re_even = Re;
		float *Re_odd;
		Re_odd = &Re[LEN];

		float *Im_even;
		Im_even = Im;
		float *Im_odd;
		Im_odd = &Im[LEN];

		float tmpRe = 0;
		float tmpIm = 0;
		for (int i = 1; i < LEN; i++)
		{
			tmpRe = Re[i];
			tmpIm = Im[i];

			for(int j = i; j < N - 1; j++)
			{
				Re[j] = Re[j + 1];
				Im[j] = Im[j + 1];
			}
			Re[N - 1] = tmpRe;
			Im[N - 1] = tmpIm;
		}

		ufftpc(Re_even, Im_even, LEN);
		ufftpc(Re_odd, Im_odd, LEN);

		for (int i = 0; i < LEN; i++)
		{ 
			Complex X_even(Re_even[i], Im_even[i]);        
			Complex X_odd(Re_odd[i], Im_odd[i]);
			Complex X_tmp=0;

			X_tmp = X_even;

			X_odd=w(i, N) * X_odd;

			X_even = X_even + X_odd;
			X_odd = X_tmp - X_odd;        

			Re_even[i] = X_even.real();
			Re_odd[i] =  X_odd.real();   

			Im_even[i] = X_even.imag();
			Im_odd[i] =  X_odd.imag();     
		}   
	}    
}

void main()
{
	//Example
	int N=L;
	float Re[L];
	float Im[L];

	for(int i=0; i<N; i++)
	{
		Re[i]=Vector[i];
		Im[i]=0;  
	}

	InitCosinewave();
	ufftpc(Re, Im, N);

	for(int i=0; i<N; i++)
	{
		float R=(float)Re[i]*Re[i];
		float I=(float)Im[i]*Im[i];
		printf("%f\n", sqrt(R+I));
	}
}
