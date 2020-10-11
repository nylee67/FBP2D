/*
	ParaBeamQ2D.c
*/

#include	"TOF_ParaBeam2D.h"
void 	Rotate2D(float **rot_f, float **f, int N, double theta);

void 	ParaBeamQ2D(
	float	**f,
	float	**Y, 
	int		Na,
	int 	N,
	double 	dd)
{
	int	    ia, iu, i;
	double 	theta; 

	float 	**rot_f, **tmp_f;
	
	rot_f = (float **) calloc(N, sizeof(float *));
	for(i=0; i<N; i++)		
		rot_f[i] = (float *) calloc(N, sizeof(float));
		
	tmp_f = (float **) calloc(N, sizeof(float *));
	for(i=0; i<N; i++)		
		tmp_f[i] = (float *) calloc(N, sizeof(float));
		
	for(iu=0; iu<N; iu++)
		for(i=0; i<N; i++)
			f[iu][i] = (float) 0.0;
		
	for(ia=0; ia<Na; ia++){
		for(iu=0; iu<N; iu++)
			for(i=0; i<N; i++)
				rot_f[i][iu] = Y[ia][iu]; 

		theta = ia*MY_PI/Na;
		
		Rotate2D(tmp_f, rot_f, N, theta);

		for(iu=0; iu<N; iu++)
			for(i=0; i<N; i++)
				f[iu][i] += tmp_f[iu][i];
	}
	
	for(iu=0; iu<N; iu++)
		for(i=0; i<N; i++)
			f[iu][i] = f[iu][i]/Na;

	for(i=0; i<N; i++)
		free(rot_f[i]);
	free(rot_f);
}
