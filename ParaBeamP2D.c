/*
	ParaBeamP2D.c
		Y(a,u)		Na x Nu projection
			a = ia * PI/Na			
			u = (iu - Nu/2.0 + 0.5) * du: 		equispaced detectors 
*/

#include	"TOF_ParaBeam2D.h"
void 	Rotate2D(float **rot_f, float **f, int N, double theta);

void 	ParaBeamP2D(
	float	**Y,
	float	**f, 
	int		Na,
	int 	N,
	double 	dd)
{
	int	    ia, iu, i;
	double 	theta; 

	float 	**rot_f;
	
	rot_f = (float **) calloc(N, sizeof(float *));
	for(i=0; i<N; i++)		
		rot_f[i] = (float *) calloc(N, sizeof(float));
		
	for(ia=0; ia<Na; ia++){
		theta = ia*MY_PI/Na;
		
		Rotate2D(rot_f, f, N, -theta);
		
		for(iu=0; iu<N; iu++){
			Y[ia][iu] = (float) 0.0;
			for(i=0; i<N; i++)
				Y[ia][iu] += rot_f[i][iu]*dd; 
		}
	}
	
	for(i=0; i<N; i++)
		free(rot_f[i]);
	free(rot_f);
}
