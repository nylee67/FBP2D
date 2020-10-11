/*
	DeConv1D_LR.c
*/

#include	"TOF_ParaBeam2D.h"

void 	PConv1D(float *g, float *f, float *tsf, int N, int N_TSF);
void 	QConv1D(float *f, float *g, float *tsf, int N, int N_TSF);

void	dfft(double *, double *, int, int);

void 	DeConv1D_LR(
	float		*f,				// Original -> Convolved One
	float		*y,
	int			N, 
	float		*tsf,
	int			N_TSF)
{
	//	Variables
		int		Nt = N + N_TSF - 1;
		int		N_Iter = 20;
		int 	i, j;
		float	*old_g, *new_g, *tmp_s, *rate_g, *all_1;

		old_g = (float *) calloc(N, sizeof(float));
		new_g = (float *) calloc(N, sizeof(float));
		rate_g = (float *) calloc(N, sizeof(float));
		
		tmp_s = (float *) calloc(Nt, sizeof(float));
		all_1 = (float *) calloc(Nt, sizeof(float));
		
		for(i=0; i<Nt; i++)
			all_1[i] = (float) 1.0;
		
		QConv1D(rate_g, all_1, tsf, N, N_TSF);
		
		for(i=0; i<N; i++)
			old_g[i] = (float) 1.0; 

		for(i=0; i<N_Iter; i++){
			PConv1D(tmp_s, old_g, tsf, N, N_TSF);
		
			for(j=0; j<Nt; j++)
				if(tmp_s[j] > 0.0)
					tmp_s[j] = y[j]/tmp_s[j];
				else
					tmp_s[j] = (float) 0.0;
			
			QConv1D(new_g, tmp_s, tsf, N, N_TSF);
		
			for(j=0; j<N; j++)
				if(rate_g[j] != 0.0)
					new_g[j] = old_g[j]*new_g[j]/rate_g[j];
				else
					new_g[j] = (float) 0.0;
			
			for(j=0; j<N; j++)
				old_g[j] = new_g[j];
		}

		for(i=0; i<N; i++)
			f[i] = new_g[i];
		
	//	Free
		free(tmp_s);
		free(old_g);
		free(new_g);
}

void 	PConv1D(float *g, float *f, float *tsf, int N, int N_TSF)
{
	int 	i, j;
	int		Nt = N + N_TSF - 1;
	
	for(j=0; j<Nt; j++)
		g[j] = (float) 0.0;
	
	for(i=0; i<N; i++)
		for(j=-N_TSF/2; j<=N_TSF/2; j++)
			g[i+j+N_TSF/2] += f[i]*tsf[j+N_TSF/2]; 
}

void 	QConv1D(float *f, float *g, float *tsf, int N, int N_TSF)
{
	int 	i, j;
	
	for(j=0; j<N; j++)
		f[j] = (float) 0.0;
	
	for(i=0; i<N; i++)
		for(j=-N_TSF/2; j<=N_TSF/2; j++)
			f[i] += g[i+j+N_TSF/2] * tsf[j+N_TSF/2];
}
