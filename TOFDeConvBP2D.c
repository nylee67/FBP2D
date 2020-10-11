/*
	TOFDeConvBP2D.c
*/

#include	"TOF_ParaBeam2D.h"

void	int2str(char *, int, int);
void	pgm_write_scaled_float2d(char *, float **, int, int);

void	dfft(double *, double *, int, int);

void	DeConv1D_LR(float *f, float *y, int N, float *tsf, int N_TSF);
void	DeConv1D_Wiener(float *f, float *y, double *rTSF, double *iTSF, int N, int N_TSF, int BigN, double Sigma_Wiener);

void 	Rotate2D(float **rot_f, float **f, int N, double theta);

void	TOFDeConvBP2D(
	float		**f,
	float		***y,
	int			Na,
	int			N,
	double		dd,
	int 		N_TSF,
	double 		Sigma_TSF,
	double 		Sigma_Wiener)
{
	//
	//	Variables
	//
		int		Nt = N + N_TSF - 1;
		int		i, j, k, BigNu;
		float	**rot_f, **tmp_f;
		float 	*lor_f, *lor, *tsf, x1, x2, temp, sum;
		double 	x, theta;
		double	*rTSF, *iTSF;
		char 	stemp[100], *temp_fn;
		
		BigNu = 1;
		while(BigNu < Nt)
			BigNu = 2*BigNu;

		tmp_f = (float **) calloc(N, sizeof(float *));
		for(i=0; i<N; i++)		
			tmp_f[i] = (float *) calloc(N, sizeof(float));
	
		rot_f = (float **) calloc(N, sizeof(float *));
		for(i=0; i<N; i++)		
			rot_f[i] = (float *) calloc(N, sizeof(float));
	
		lor_f = (float *) calloc(N, sizeof(float));
		lor = (float *) calloc(Nt, sizeof(float));
		tsf = (float *) calloc(N_TSF, sizeof(float));
	
		temp = (float) (1.0/(sqrt(2*MY_PI) * Sigma_TSF));	// Assume that N_TSF is odd

		tsf[N_TSF/2] = temp;	
		for(i=1; i<N_TSF/2; i++){
			x = i*dd;
			tsf[N_TSF/2 + i] = (float) (temp*exp(-x*x/(2*Sigma_TSF*Sigma_TSF)));
			tsf[N_TSF/2 - i] = (float) (temp*exp(-x*x/(2*Sigma_TSF*Sigma_TSF)));
		}

		sum = 0.0;
		for(i=0; i<N_TSF; i++)
			sum += tsf[i];
	
		for(i=0; i<N_TSF; i++)
			tsf[i] = tsf[i]/sum;
		
		rTSF = (double *) calloc(BigNu, sizeof(double));
		iTSF = (double *) calloc(BigNu, sizeof(double));

		rTSF[0] = tsf[N_TSF/2];
		for(i=1; i<=N_TSF/2; i++){		
			rTSF[i] = tsf[i+N_TSF/2];
			rTSF[BigNu-i] = tsf[-i+N_TSF/2];
		}

		dfft(rTSF, iTSF, BigNu, -1);

	//	DeConv1D
		for(i=0; i<Na; i++){
			for(j=0; j<N; j++){
				for(k=0; k<Nt; k++)
					lor[k] = y[i][j][k];
				
				//DeConv1D_LR(lor_f, lor, N, tsf, N_TSF); 
				DeConv1D_Wiener(lor_f, lor, rTSF, iTSF, N, N_TSF, BigNu, Sigma_Wiener); 

				for(k=0; k<N; k++)
					rot_f[j][k] = lor_f[k];
			}
			
			if(1){
				temp_fn = (char *) calloc(200, sizeof(char)); 	
				strcat(temp_fn, "../../3_Output/Test"); 
				strcat(temp_fn, "_DeConved_TOF_p_sharp_");
				int2str(stemp, 3, i);
				strcat(temp_fn, stemp);
				strcat(temp_fn, ".pgm");
				pgm_write_scaled_float2d(temp_fn, rot_f, N, N);
				free(temp_fn);
			}
			
			for(j=0; j<N; j++)
				for(k=0; k<N; k++)
					tmp_f[j][k] = rot_f[k][j];

			if(0){
				temp_fn = (char *) calloc(200, sizeof(char)); 	
				strcat(temp_fn, "../../3_Output/Test"); 
				strcat(temp_fn, "_DeConved_TOF_p_sharp_switched_");
				int2str(stemp, 3, i);
				strcat(temp_fn, stemp);
				strcat(temp_fn, ".pgm");
				pgm_write_scaled_float2d(temp_fn, tmp_f, N, N);
				free(temp_fn);
			}
			
			theta = i*MY_PI/Na;
		
			Rotate2D(rot_f, tmp_f, N, theta);
		
			if(1){
				temp_fn = (char *) calloc(200, sizeof(char)); 	
				strcat(temp_fn, "../../3_Output/Test"); 
				strcat(temp_fn, "_DeConved_TOF_p_sharp_rotated_");
				int2str(stemp, 3, i);
				strcat(temp_fn, stemp);
				strcat(temp_fn, ".pgm");
				pgm_write_scaled_float2d(temp_fn, rot_f, N, N);
				free(temp_fn);
			}
			
			for(j=0; j<N; j++)
				for(k=0; k<N; k++)
					f[j][k] += rot_f[j][k];
		}
		
		for(j=0; j<N; j++)
			for(k=0; k<N; k++)
				f[j][k] = f[j][k]/Na;
		
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
		for(i=0; i<N; i++)
			free(rot_f[i]);
		free(rot_f);
	
		for(i=0; i<N; i++)
			free(tmp_f[i]);
		free(tmp_f);
	
        free(lor); 
		free(tsf);
}
