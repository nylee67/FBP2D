/*
	FBP2D.c
*/

#include	"TOF_ParaBeam2D.h"

void	RampFilter1D(double *real_f, double *imag_f, int BigNu, double CutOff);
void 	ParaBeamQ2D(float **f, float **Y, int Na, int N, double	dd);

void	FBP2D(
	float		**f,
	float		**Y,
	int			Na,
	int			N,
	double		dd,
	double		CutOff)
{
	//
	//	Variables
	//
		int		i, j, BigNu;
		float	**tmp_Y;
		double	**realF, **imagF, x1, x2;

		BigNu = 1;
		while(BigNu < N)
			BigNu = 2*BigNu;

		tmp_Y = (float **) calloc(Na, sizeof(float *));
        for(i=0; i<Na; i++)
			tmp_Y[i] = (float *) calloc(N, sizeof(float));

        realF = (double **) calloc(Na, sizeof(double *));
        imagF = (double **) calloc(Na, sizeof(double *));
        for(i=0; i<Na; i++){
            realF[i] = (double *) calloc(BigNu, sizeof(double));
            imagF[i] = (double *) calloc(BigNu, sizeof(double));
        }

	//	RampFilter1D
		for(i=0; i<Na; i++){
			for(j=0; j<N; j++)
				realF[i][j] = (double) Y[i][j];

			RampFilter1D(realF[i], imagF[i], BigNu, CutOff); 

			for(j=0; j<N; j++)
				tmp_Y[i][j] = (float) realF[i][j];
		}

	//	Backprojection
		ParaBeamQ2D(f, tmp_Y, Na, N, dd); 

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
            free(realF[i]); 	free(imagF[i]);		free(tmp_Y[i]);
        }
        free(realF); 	free(imagF);	free(tmp_Y);
}
