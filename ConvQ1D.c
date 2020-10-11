/*
	ConvQ1D.c
		Direct BackConvolution
		Assumption: TSF is symmetric
*/

#include	"TOF_ParaBeam2D.h"
void 	ConvQ1D(
	float 		*f,			// N + N_TSF - 1
	float		*y,				
	int			N, 
	float		*tsf,
	int			N_TSF)
{
	//	Variables
		int		i, j;

	//	Initilization
		for(j=0; j<N; j++)
			f[j] = (float) 0.0;

	//	Convolution
		for(i=0; i<N; i++)
			for(j=-N_TSF/2; j<=N_TSF/2; j++)
				f[i] += y[i+j+N_TSF/2] * tsf[j+N_TSF/2];
}
