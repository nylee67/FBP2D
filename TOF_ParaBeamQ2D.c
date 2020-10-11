/*
	TOF_ParaBeamQ2D.c
*/

#include	"TOF_ParaBeam2D.h"
void 	Rotate2D(float **rot_f, float **f, int N, double theta);
void	ConvQ1D(float *f, float *y, int N, float *tsf, int N_TSF); 

void 	TOF_ParaBeamQ2D(
	float	**f, 
	float	***p,
	int		Na,
	int 	N,
	double 	dd,
	int 	N_TSF,
	double 	Sigma_TSF)
{
	int	    ia, iu, i;
	double 	theta, x; 
	float 	**rot_f, **tmp_f;
	float   *lor, *tsf, sum;
	
	rot_f = (float **) calloc(N, sizeof(float *));
	for(i=0; i<N; i++)		
		rot_f[i] = (float *) calloc(N, sizeof(float));
	
	tmp_f = (float **) calloc(N, sizeof(float *));
	for(i=0; i<N; i++)		
		tmp_f[i] = (float *) calloc(N, sizeof(float));
	
	lor = (float *) calloc(N, sizeof(float));
	tsf = (float *) calloc(N_TSF, sizeof(float));
	
	tsf[N_TSF/2] = (float)1.0;	
	for(i=1; i<N_TSF/2; i++){
		x = i*dd;
		tsf[N_TSF/2 + i] = (float) (exp(-x*x/(2*Sigma_TSF*Sigma_TSF)));
		tsf[N_TSF/2 - i] = (float) (exp(-x*x/(2*Sigma_TSF*Sigma_TSF)));
	}

	sum = 0.0;
	for(i=0; i<N_TSF; i++)
		sum += tsf[i];

	for(i=0; i<N_TSF; i++)
		tsf[i] = tsf[i]/sum;

	for(iu=0; iu<N; iu++)
		for(i=0; i<N; i++)
			f[iu][i] = (float) 0.0;
		
	for(ia=0; ia<Na; ia++){
		for(iu=0; iu<N; iu++){
			ConvQ1D(lor, p[ia][iu], N, tsf, N_TSF); 
			
			for(i=0; i<N; i++)
				rot_f[i][iu] = lor[i]; 
		}
		
		theta = ia*MY_PI/Na;
		
		Rotate2D(tmp_f, rot_f, N, theta);
		
		for(iu=0; iu<N; iu++)
			for(i=0; i<N; i++)
				f[iu][i] += tmp_f[iu][i];
	}
	
	for(i=0; i<N; i++)
		free(rot_f[i]);
	free(rot_f);
	
	for(i=0; i<N; i++)
		free(tmp_f[i]);
	free(tmp_f);
	
	free(lor);
	free(tsf);
}
