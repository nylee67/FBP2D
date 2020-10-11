/*
	TOFSinogram_Phantom2D.c
*/

#include	"TOFFBP2D.h"

void 	Rotation2D(
	float		**g,				// Rotated One:		N x N
	float		**f,				// Original:		N x N
	int			N,
	double		theta,
	double		dd);

void TOFSinogram_Phantom2D(
	float	***TOF_Y,
	float	**f, 
	int		Na,
	int		N,
	double	dd,
	int		N_TOFPSF)
{
	int		ia, iu, it, k;
	int		N_Wing = (N_TOFPSF-1)/2;

	double 	theta;

	float	**rot_f;

	rot_f = (float **) calloc(N, sizeof(float *));
	for(iu=0; iu<N; iu++)		
		rot_f[iu] = (float *) calloc(N, sizeof(float));

    for(ia=0; ia<Na; ia++) {   		
		theta = (ia+0.5) * MY_PI/Na;	
		Rotation2D(rot_f, f, N, theta, dd);

		for(it=0; it<N; it++){		
			for(iu=0; iu<N_Wing; iu++){
				TOF_Y[ia][iu][it] = (float) 0.0;
				rot_f[iu][it]	  = (float) 0.0;
			}
			for(iu=N_Wing; iu<N-N_Wing; iu++){
				TOF_Y[ia][iu][it] = (float) 0.0;
				for(k=0; k<N_TOFPSF; k++)
					TOF_Y[ia][iu][it] += rot_f[iu-N_Wing+k][it]*dd;
				TOF_Y[ia][iu][it] = TOF_Y[ia][iu][it]/N_TOFPSF;
			}
			for(iu=N-N_Wing; iu<N; iu++){
				TOF_Y[ia][iu][it] = (float) 0.0;
				rot_f[iu][it]	  = (float) 0.0;
			}
		}		
	}// loop for ia

	for(iu=0; iu<N; iu++)
		free(rot_f[iu]);
	free(rot_f);
}
