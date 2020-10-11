/*
	ConvP1D.c
		Direct Convolution
		Assumption: TSF is symmetric
*/

#include	"TOF_ParaBeam2D.h"
void 	ConvP1D(
	float 		*y,			// N + N_TSF - 1
	float		*f,				
	int			N, 
	float		*tsf,
	int			N_TSF)
{
	//	Variables
		int		i, j;

	//	Convolution
		for(i=0; i<N+N_TSF-1; i++)
			y[i] = (float) 0.0;

		for(i=0; i<N; i++)
			for(j=-N_TSF/2; j<=N_TSF/2; j++)
				y[i+j+N_TSF/2] += f[i] * tsf[j+N_TSF/2];
}
