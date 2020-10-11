/*
	ConvP1D.c
*/

#include	"TOF_ParaBeam2D.h"
void 	ConvP1D(
	float 		*y,			// N +N_TSF - 1
	float		*f,				
	int			N, 
	float		*tsf,
	int			N_TSF)
{
	//	Variables
		int		i, j;
		float	sum;

	//	Convolution
		for(i=0; i<N+N_TSF-1; i++){
			sum = (float) 0.0;
			for(j=0; j<N_TSF; j++)
				if(i+j-N_TSF/2 >= 0 && i+j-N_TSF/2 < N)
					y[i] += f[i+j-N_TSF/2] * tsf[j];
		}

	//	Free
		free(tmp_f);
}
