/*
	TOFEM2D.c
*/
#include	"TOF_ParaBeam2D.h"

void		int2str(char *, int, int);
void		pgm16_write_scaled_float2d(char *, float **, int, int);

void		TOF_ParaBeamP2D(float ***y, float **f, int Na, int N, double dd, int N_TSF, double Sigma_TSF); 
void		TOF_ParaBeamQ2D(float **f, float ***y, int Na, int N, double dd, int N_TSF, double Sigma_TSF); 

void		TOFEM2D(
	float 			**f, 
	float 			***y, 
	int 			Na, 
	int 			N, 
	double 			dd, 
	int 			N_TSF, 
	double 			Sigma_TSF,
	char 			*PGM_Header)
{
	//
	//	Declaration of Input Parameters
	//
		float		**old_f, **new_f, **rate_f, ***OS_y, ***all_1_y;
		char		*temp_fn, stemp[50];
		
		int			i, j, k, n;
		int			Nt = N + N_TSF - 1;
		int 		N_Iter = 20; 
	
	//
	// 	Memory Allocation
	//
		old_f = (float **) calloc(N, sizeof(float *));
		for(i=0; i<N; i++)		
			old_f[i] = (float *) calloc(N, sizeof(float));
	
			
		new_f = (float **) calloc(N, sizeof(float *));
		for(i=0; i<N; i++)		
			new_f[i] = (float *) calloc(N, sizeof(float));
	
		rate_f = (float **) calloc(N, sizeof(float *));
		for(i=0; i<N; i++)		
			rate_f[i] = (float *) calloc(N, sizeof(float));
	
		OS_y = (float ***) calloc(Na, sizeof(float **));
		for(i=0; i<Na; i++){		
			OS_y[i] = (float **) calloc(N, sizeof(float *));
			for(j=0; j<N; j++)		
				OS_y[i][j] = (float *) calloc(Nt, sizeof(float));
		}
		
		all_1_y = (float ***) calloc(Na, sizeof(float **));
		for(i=0; i<Na; i++){		
			all_1_y[i] = (float **) calloc(N, sizeof(float *));
			for(j=0; j<N; j++)		
				all_1_y[i][j] = (float *) calloc(Nt, sizeof(float));
		}
		
		for(i=0; i<Na; i++)		
			for(j=0; j<N; j++)		
				for(k=0; k<Nt; k++)		
					all_1_y[i][j][k] = (float) 1.0;
					
	//
	//	Ready for Iteration[
	//
		TOF_ParaBeamQ2D(rate_f, all_1_y, Na, N, dd, N_TSF, Sigma_TSF); 
		
		for(i=0; i<N; i++)
			for(j=0; j<N; j++)
				if((i-N/2)*(i-N/2) + (j-N/2)*(j-N/2) > N*N/4)
					rate_f[i][j] = (float) 0.0;

		if(1){
			temp_fn = (char *) calloc(200, sizeof(char)); 	
			strcat(temp_fn, PGM_Header); 
			strcat(temp_fn, "_EPD.pgm");
			pgm16_write_scaled_float2d(temp_fn, rate_f, N, N);
			free(temp_fn);
		}
		
	//
	//	EM Iteration
	//
		for(i=0; i<N; i++)		
			for(j=0; j<N; j++)
				old_f[i][j] = (float) 1.0;

		for(n=0; n<N_Iter; n++){

			printf("TOFEM at %4d-th Iteration\n", n);
			
			TOF_ParaBeamP2D(OS_y, old_f, Na, N, dd, N_TSF, Sigma_TSF); 

			for(i=0; i<Na; i++)		
				for(j=0; j<N; j++)		
					for(k=0; k<Nt; k++)
						if(OS_y[i][j][k] > EPSILON)
							OS_y[i][j][k] = y[i][j][k] / OS_y[i][j][k];
						else
							OS_y[i][j][k] = (float) 0.0;
						
			TOF_ParaBeamQ2D(new_f, OS_y, Na, N, dd, N_TSF, Sigma_TSF); 

			for(i=0; i<N; i++)
				for(j=0; j<N; j++)
					if(rate_f[i][j] > EPSILON && (i-N/2)*(i-N/2) + (j-N/2)*(j-N/2) < N*N/4)
						new_f[i][j] = old_f[i][j] * new_f[i][j] / rate_f[i][j];
					else
						new_f[i][j] = 0.0;

			for(i=0; i<N; i++)
				for(j=0; j<N; j++)
					old_f[i][j] = new_f[i][j];

			if(1){
				temp_fn = (char *) calloc(200, sizeof(char)); 	
				strcat(temp_fn, PGM_Header); 
				strcat(temp_fn, "_TOFEM-6_");
				int2str(stemp, 3, n);	
				strcat(temp_fn, stemp); 	
				strcat(temp_fn, ".pgm");
				pgm16_write_scaled_float2d(temp_fn, new_f, N, N);
				free(temp_fn);
			}
		}// n-loop
		
		for(i=0; i<N; i++)		
			for(j=0; j<N; j++)
				f[i][j] = new_f[i][j];
			
	// Free
		for(i=0; i<N; i++){
			free(old_f[i]); free(new_f[i]); free(rate_f[i]);
		}
		free(old_f); free(new_f); free(rate_f);
		
		for(i=0; i<Na; i++){
			for(j=0; j<N; j++){
				free(OS_y[i][j]); free(all_1_y[i][j]);
			}			
			free(OS_y[i]); free(all_1_y[i]);
		}
		free(OS_y); free(all_1_y);
}
