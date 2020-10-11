/*
	TOFFBP2D.c
*/

#include	"TOF_ParaBeam2D.h"

void	RampFilter1D(double *real_f, double *imag_f, int BigNu, double CutOff);
void 	TOF_ParaBeamQ2D(float **f, float ***y, int Na, int N, double dd, int N_TSF, double Sigma_TSF);

void	TOFFBP2D(
	float		**f,
	float		***y,
	int			Na,
	int			N,
	double		dd,
	int 		N_TSF,
	double 		Sigma_TSF,
	double		CutOff)
{
	//
	//	Variables
	//
		int		i, j, k, BigNu;
		int		Nt = N + N_TSF - 1;
		float	***tmp_y;
		double	*realF, *imagF, x1, x2;

		BigNu = 1;
		while(BigNu < N)
			BigNu = 2*BigNu;

		tmp_y = (float ***) calloc(Na, sizeof(float **));
        for(i=0; i<Na; i++){
			tmp_y[i] = (float **) calloc(N, sizeof(float *));
			for(j=0; j<N; j++)
				tmp_y[i][j] = (float *) calloc(Nt, sizeof(float));
		}
		
        realF = (double *) calloc(BigNu, sizeof(double));
        imagF = (double *) calloc(BigNu, sizeof(double));

	//	RampFilter1D
		for(i=0; i<Na; i++){
			for(k=0; k<Nt; k++){
				for(j=0; j<BigNu; j++){
					realF[j] = 0.0;
					imagF[j] = 0.0;
				}
				
				for(j=0; j<N; j++)
					realF[j] = (double) y[i][j][k];
				
				RampFilter1D(realF, imagF, BigNu, CutOff); 

				for(j=0; j<N; j++)
					tmp_y[i][j][k] = (float) realF[j];
			}
		}

	//	Backprojection
		TOF_ParaBeamQ2D(f, tmp_y, Na, N, dd, N_TSF, Sigma_TSF); 

	//	adjust reconstructed image function inside circle
		if(1){
			for(i=0;i<N;i++) {
				x2 = (double) (N/2.0-i-0.5) / (double) N;
				for(j=0;j<N;j++) {
					x1 = (double) (j-N/2.0+0.5) / (double) N;
					if(f[i][j] < 0.0)
						f[i][j] = 0.0;
					if(x1*x1 + x2*x2 > 0.25)
						f[i][j] = 0.0;
				}
			}
		}

	//	Free
        for(i=0; i<Na; i++){
			for(j=0; j<N; j++)
				free(tmp_y[i][j]);
			free(tmp_y[i]);
		}
		free(tmp_y);
		
        free(realF); 	free(imagF);	
}
