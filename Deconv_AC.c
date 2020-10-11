/*
	Deconv_AC.c
*/
#include	"TOF_ParaBeam2D.h"

void	int2str(char *, int, int);
void	pgm_write_scaled_float2d(char *, float **, int, int);

void	dfft(double *, double *, int, int);

void	DeConv1D_LR(float *f, float *y, int N, float *tsf, int N_TSF);
void	DeConv1D_Wiener(float *f, float *y, double *rTSF, double *iTSF, int N, int N_TSF, int BigN, double Sigma_Wiener);

void 	Rotate2D(float **rot_f, float **f, int N, double theta);

void	Deconv_AC(
	int 	ia_0,
	int 	iu_0,
	int 	it_0,
	float	**a,
	float 	***m,
	int 	Na,
	int 	N,
	double	dd,
	int 	N_TSF,
	double 	Sigma_TSF,
	double 	Sigma_Wiener)
{
	int		Nu = N;
	int 	Nt = N + N_TSF - 1;
	int 	ia, iu, i, j, k, BigNu;
	
	float 	***m_sharp, ***m_sharp_deattenuated, ***m_sharp_rotated, **tmp_f, **ref_f, *tsf, temp, *vec_x, *vec_y, **est_M, **est_P, sum_xy, sum_x2, sum_mean, sum_each;
	int		**idx_A, N_vec;
	double 	*rTSF, *iTSF, sum;
	
	double	d_theta = MY_PI/Na;

	double  ref_value, ref_center;
	
	double	theta;
	double  x;
	
	char 	stemp[100], *temp_fn;

	float	eps = (float) 0.001;
	float	atten_eps = (float) 0.03;

	//
	//	Memory Allocation
	//
		m_sharp = (float ***) calloc(Na, sizeof(float **));
		for(ia=0; ia<Na; ia++){
			m_sharp[ia] = (float **) calloc(Nu, sizeof(float *));
			for(iu=0; iu<Nu; iu++)
				m_sharp[ia][iu] = (float *) calloc(N, sizeof(float));
		}
			
		m_sharp_deattenuated = (float ***) calloc(Na, sizeof(float **));
		for(ia=0; ia<Na; ia++){
			m_sharp_deattenuated[ia] = (float **) calloc(Nu, sizeof(float *));
			for(iu=0; iu<Nu; iu++)
				m_sharp_deattenuated[ia][iu] = (float *) calloc(N, sizeof(float));
		}
			
		m_sharp_rotated = (float ***) calloc(Na, sizeof(float **));
		for(ia=0; ia<Na; ia++){
			m_sharp_rotated[ia] = (float **) calloc(Nu, sizeof(float *));
			for(iu=0; iu<Nu; iu++)
				m_sharp_rotated[ia][iu] = (float *) calloc(N, sizeof(float));
		}
			
		est_M = (float **) calloc(Na, sizeof(float *));
		for(i=0; i<N; i++)
			est_M[i] = (float *) calloc(Nu, sizeof(float));

		est_P = (float **) calloc(Na, sizeof(float *));
		for(i=0; i<N; i++)
			est_P[i] = (float *) calloc(Nu, sizeof(float));

		tmp_f = (float **) calloc(N, sizeof(float *));
		for(i=0; i<N; i++)
			tmp_f[i] = (float *) calloc(N, sizeof(float));

		ref_f = (float **) calloc(N, sizeof(float *));
		for(i=0; i<N; i++)
			ref_f[i] = (float *) calloc(N, sizeof(float));

		idx_A = (int **) calloc(Na, sizeof(int *));
		for(i=0; i<Na; i++)
			idx_A[i] = (int *) calloc(Nu, sizeof(int));

		vec_x = (float *) calloc(N, sizeof(float));
		vec_y = (float *) calloc(N, sizeof(float));

	//
	//	Deconvolution: m --> m_sharp
	//
		BigNu = 1;
		while(BigNu < Nt)
			BigNu = 2*BigNu;
		//BigNu = 4*BigNu;

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
		
		for(ia=0; ia<Na; ia++){
			for(iu=0; iu<Nu; iu++){
				//DeConv1D_LR(m_sharp[ia][iu], m[ia][iu], N, tsf, N_TSF); 
				DeConv1D_Wiener(m_sharp[ia][iu], m[ia][iu], rTSF, iTSF, N, N_TSF, BigNu, Sigma_Wiener); 

				for(i=0; i<N; i++){
					if(m_sharp[ia][iu][i] < 0.0)
						m_sharp[ia][iu][i] = (float) 0.0;
					m_sharp_deattenuated[ia][iu][i] = m_sharp[ia][iu][i];
				}			
			}

			if(1){
				temp_fn = (char *) calloc(200, sizeof(char)); 	
				strcat(temp_fn, "../../3_Output/Test"); 
				strcat(temp_fn, "_DeConved_TOF_m_sharp_");
				int2str(stemp, 3, ia);
				strcat(temp_fn, stemp);
				strcat(temp_fn, ".pgm");
				pgm_write_scaled_float2d(temp_fn, m_sharp[ia], N, N);
				free(temp_fn);
			}
		}

		ref_value 	= m_sharp[ia_0][iu_0][it_0];
		ref_center 	= m_sharp[ia_0-Na/2][N/2][N-iu_0-1];
		if(1){
			printf("m_sharp at (MY_PI/2, u_0, 0) = %20.8f\n", ref_value);
			printf("m_sharp at (0, 0, u_0)       = %20.8f\n", ref_center);
		}
		
		free(tsf);
		free(rTSF);
		free(iTSF);
	
		if(ref_center <= 0.0){
			puts("ref_center is not positive\n");
			exit(1);
		}

	//
	//	Deconv_AC - Initialization: Compute m_sharp at (0,0,0) and (pi/2,0,0) 	
	//
		if(0){
			for(ia=0; ia<Na; ia++)
				printf("Center Value at %4d = %20.12f\n", ia, m_sharp[ia][N/2][N/2]); 
			exit(1);
		}

		// Initial	
		if(ia_0 != Na/2){
			puts("ia_0 must be Na/2\n");
			exit(1);
		}

		ref_center = m_sharp_deattenuated[0][N/2][N/2]/ref_center * ref_value;
		if(1){
			printf("ref_center = %20.12f\n", ref_center);
		}

		for(ia=0; ia<Na/2; ia++){

			// vertical center point
			idx_A[ia][N/2] = 1;
			a[ia][N/2] = m_sharp[ia][N/2][N/2]/ref_center; 

			if(0)
				printf("a[ia][N/2] = %20.12f\n",a[ia][N/2]);
 
			// vertical center line
			for(k=0; k<N; k++){
				vec_x[k] = m_sharp_deattenuated[ia][N/2][k];
				m_sharp_deattenuated[ia][N/2][k] = m_sharp_deattenuated[ia][N/2][k] / a[ia][N/2]; 
			}

			if(1){
				temp_fn = (char *) calloc(200, sizeof(char)); 	
				strcat(temp_fn, "../../3_Output/Test"); 
				strcat(temp_fn, "_DeAttenuated_m_sharp_deattenuated_v0_");
				int2str(stemp, 3, ia);
				strcat(temp_fn, stemp);
				strcat(temp_fn, ".pgm");
				pgm_write_scaled_float2d(temp_fn, m_sharp_deattenuated[ia], N, N);
				free(temp_fn);
			}

			// horizontal lineS
			for(k=0; k<N; k++)
				if(m_sharp_deattenuated[ia][N/2][k] > eps && idx_A[ia+Na/2][N-k-1] == 0){
					idx_A[ia+Na/2][N-k-1] = 1;
					if(m_sharp_deattenuated[ia+Na/2][N-k-1][N/2] > eps)
						a[ia+Na/2][N-k-1] = m_sharp_deattenuated[ia+Na/2][N-k-1][N/2] / m_sharp_deattenuated[ia][N/2][k];
					else
						idx_A[ia+Na/2][N-k-1] = 0;
						
					if(a[ia+Na/2][N-k-1] < atten_eps)
						idx_A[ia+Na/2][N-k-1] = 0;

					for(j=0; j<N; j++)
						if(idx_A[ia+Na/2][N-k-1])
							m_sharp_deattenuated[ia+Na/2][N-k-1][j] = m_sharp_deattenuated[ia+Na/2][N-k-1][j] / a[ia+Na/2][N-k-1];
				}
			
			if(1){
				temp_fn = (char *) calloc(200, sizeof(char)); 	
				strcat(temp_fn, "../../3_Output/Test"); 
				strcat(temp_fn, "_DeAttenuated_m_sharp_deattenuated_h1_");
				int2str(stemp, 3, ia+Na/2);
				strcat(temp_fn, stemp);
				strcat(temp_fn, ".pgm");
				pgm_write_scaled_float2d(temp_fn, m_sharp_deattenuated[ia+Na/2], N, N);
				free(temp_fn);
			}

			if(0)
				for(j=0; j<N; j++)
					printf("h1 a[ia+Na/2][j] = %20.12f at %4d\n", a[ia+Na/2][j], j);

			// multi-vertical lineS
			for(j=0; j<N; j++)
				m_sharp_deattenuated[ia][N/2][j] = vec_x[j];

			for(k=0; k<N; k++){

				N_vec = 0;
				for(j=0; j<N; j++)
					if(m_sharp_deattenuated[ia+Na/2][j][k] > eps && m_sharp_deattenuated[ia][k][N-j-1] > eps){
						vec_y[N_vec] = m_sharp_deattenuated[ia][k][N-j-1];
						vec_x[N_vec] = m_sharp_deattenuated[ia+Na/2][j][k];
						N_vec += 1;
					}

				if(N_vec > 0){
					idx_A[ia][k] = 1;
					sum_x2 = (float) 0.0;
					sum_xy = (float) 0.0;
					for(j=0; j<N_vec; j++){
						sum_xy += vec_x[j]*vec_y[j];
						sum_x2 += vec_x[j]*vec_x[j];
					}
					a[ia][k] = sum_xy/sum_x2;

					//printf("v1 a[ia][k] = %4d %20.12f\n", k, a[ia][k]);

					for(j=0; j<N; j++)
						m_sharp_deattenuated[ia][k][j] = m_sharp_deattenuated[ia][k][j] / a[ia][k];	
				}
			}

			if(1){
				temp_fn = (char *) calloc(200, sizeof(char)); 	
				strcat(temp_fn, "../../3_Output/Test"); 
				strcat(temp_fn, "_DeAttenuated_m_sharp_deattenuated_v1_");
				int2str(stemp, 3, ia);
				strcat(temp_fn, stemp);
				strcat(temp_fn, ".pgm");
				pgm_write_scaled_float2d(temp_fn, m_sharp_deattenuated[ia], N, N);
				free(temp_fn);
			}

			// mutil-multi-horizontal lineSS
			for(k=0; k<N; k++){

				N_vec = 0;
				for(j=0; j<N; j++)
					if(m_sharp_deattenuated[ia+Na/2][N-k-1][j] > eps && m_sharp_deattenuated[ia][j][k] > eps){
						vec_x[N_vec] = m_sharp_deattenuated[ia][j][k];
						vec_y[N_vec] = m_sharp_deattenuated[ia+Na/2][N-k-1][j];
						N_vec += 1;
					}

				if(N_vec > 0){
					idx_A[ia+Na/2][N-k-1] = 1;
					sum_x2 = (float) 0.0;
					sum_xy = (float) 0.0;
					for(j=0; j<N_vec; j++){
						sum_xy += vec_x[j]*vec_y[j];
						sum_x2 += vec_x[j]*vec_x[j];
					}
					a[ia+Na/2][N-k-1] = sum_xy/sum_x2;

					//printf("h2 a[ia+Na/2][N-k-1] = %4d %20.12f\n", k, a[ia+Na/2][N-k-1]);

					for(j=0; j<N; j++)
						m_sharp_deattenuated[ia+Na/2][N-k-1][j] = m_sharp_deattenuated[ia+Na/2][N-k-1][j] / a[ia+Na/2][N-k-1];	
				}
			}

			if(1){
				temp_fn = (char *) calloc(200, sizeof(char)); 	
				strcat(temp_fn, "../../3_Output/Test"); 
				strcat(temp_fn, "_DeAttenuated_m_sharp_deattenuated_h2_");
				int2str(stemp, 3, ia+Na/2);
				strcat(temp_fn, stemp);
				strcat(temp_fn, ".pgm");
				pgm_write_scaled_float2d(temp_fn, m_sharp_deattenuated[ia+Na/2], N, N);
				free(temp_fn);
			}

			// multi-vertical lineS
			for(k=0; k<N; k++){

				N_vec = 0;
				for(j=0; j<N; j++)
					if(m_sharp_deattenuated[ia+Na/2][j][k] > eps && m_sharp_deattenuated[ia][k][N-j-1] > eps){
						vec_y[N_vec] = m_sharp_deattenuated[ia][k][N-j-1];
						vec_x[N_vec] = m_sharp_deattenuated[ia+Na/2][j][k];
						N_vec += 1;
					}

				if(N_vec > 0){
					idx_A[ia][k] = 1;
					sum_x2 = (float) 0.0;
					sum_xy = (float) 0.0;
					for(j=0; j<N_vec; j++){
						sum_xy += vec_x[j]*vec_y[j];
						sum_x2 += vec_x[j]*vec_x[j];
					}
					a[ia][k] = sum_xy/sum_x2;

					//printf("v1 a[ia][k] = %4d %20.12f\n", k, a[ia][k]);

					for(j=0; j<N; j++)
						m_sharp_deattenuated[ia][k][j] = m_sharp_deattenuated[ia][k][j] / a[ia][k];	
				}
			}

			if(1){
				temp_fn = (char *) calloc(200, sizeof(char)); 	
				strcat(temp_fn, "../../3_Output/Test"); 
				strcat(temp_fn, "_DeAttenuated_m_sharp_deattenuated_v2_");
				int2str(stemp, 3, ia);
				strcat(temp_fn, stemp);
				strcat(temp_fn, ".pgm");
				pgm_write_scaled_float2d(temp_fn, m_sharp_deattenuated[ia], N, N);
				free(temp_fn);
			}

			// mutil-multi-horizontal lineSS
			for(k=0; k<N; k++){

				N_vec = 0;
				for(j=0; j<N; j++)
					if(m_sharp_deattenuated[ia+Na/2][N-k-1][j] > eps && m_sharp_deattenuated[ia][j][k] > eps){
						vec_x[N_vec] = m_sharp_deattenuated[ia][j][k];
						vec_y[N_vec] = m_sharp_deattenuated[ia+Na/2][N-k-1][j];
						N_vec += 1;
					}

				if(N_vec > 0){
					idx_A[ia+Na/2][N-k-1] = 1;
					sum_x2 = (float) 0.0;
					sum_xy = (float) 0.0;
					for(j=0; j<N_vec; j++){
						sum_xy += vec_x[j]*vec_y[j];
						sum_x2 += vec_x[j]*vec_x[j];
					}
					a[ia+Na/2][N-k-1] = sum_xy/sum_x2;

					//printf("h2 a[ia+Na/2][N-k-1] = %4d %20.12f\n", k, a[ia+Na/2][N-k-1]);

					for(j=0; j<N; j++)
						m_sharp_deattenuated[ia+Na/2][N-k-1][j] = m_sharp_deattenuated[ia+Na/2][N-k-1][j] / a[ia+Na/2][N-k-1];	
				}
			}

			if(1){
				temp_fn = (char *) calloc(200, sizeof(char)); 	
				strcat(temp_fn, "../../3_Output/Test"); 
				strcat(temp_fn, "_DeAttenuated_m_sharp_deattenuated_h3_");
				int2str(stemp, 3, ia+Na/2);
				strcat(temp_fn, stemp);
				strcat(temp_fn, ".pgm");
				pgm_write_scaled_float2d(temp_fn, m_sharp_deattenuated[ia+Na/2], N, N);
				free(temp_fn);
			}

			// multi-vertical lineS
			for(k=0; k<N; k++){

				N_vec = 0;
				for(j=0; j<N; j++)
					if(m_sharp_deattenuated[ia+Na/2][j][k] > eps && m_sharp_deattenuated[ia][k][N-j-1] > eps){
						vec_y[N_vec] = m_sharp_deattenuated[ia][k][N-j-1];
						vec_x[N_vec] = m_sharp_deattenuated[ia+Na/2][j][k];
						N_vec += 1;
					}

				if(N_vec > 0){
					idx_A[ia][k] = 1;
					sum_x2 = (float) 0.0;
					sum_xy = (float) 0.0;
					for(j=0; j<N_vec; j++){
						sum_xy += vec_x[j]*vec_y[j];
						sum_x2 += vec_x[j]*vec_x[j];
					}
					a[ia][k] = sum_xy/sum_x2;

					//printf("v1 a[ia][k] = %4d %20.12f\n", k, a[ia][k]);

					for(j=0; j<N; j++)
						m_sharp_deattenuated[ia][k][j] = m_sharp_deattenuated[ia][k][j] / a[ia][k];	
				}
			}

			if(1){
				temp_fn = (char *) calloc(200, sizeof(char)); 	
				strcat(temp_fn, "../../3_Output/Test"); 
				strcat(temp_fn, "_DeAttenuated_m_sharp_deattenuated_v3_");
				int2str(stemp, 3, ia);
				strcat(temp_fn, stemp);
				strcat(temp_fn, ".pgm");
				pgm_write_scaled_float2d(temp_fn, m_sharp_deattenuated[ia], N, N);
				free(temp_fn);
			}

			ref_center = (m_sharp_deattenuated[ia][N/2][N/2] + m_sharp_deattenuated[ia+Na/2][N/2][N/2])/2.0; 
		}

		for(ia=0; ia<Na; ia++){

			theta = ia*d_theta;

			for(j=0; j<N; j++)
				for(k=0; k<N; k++)
					tmp_f[j][k] = m_sharp[ia][k][j];
	
			Rotate2D(m_sharp_rotated[ia], tmp_f, N, theta);

			for(iu=0; iu<Nu; iu++)
				for(k=0; k<N; k++)
					ref_f[iu][k] += m_sharp_rotated[ia][iu][k];
		}

		if(1){
			temp_fn = (char *) calloc(200, sizeof(char)); 	
			strcat(temp_fn, "../../3_Output/Test"); 
			strcat(temp_fn, "_ref_f_SumUp_UnDeAttenuated.pgm");
			pgm_write_scaled_float2d(temp_fn, ref_f, N, N);
			free(temp_fn);
		}

		for(iu=0; iu<Nu; iu++)
			for(k=0; k<N; k++)
				ref_f[iu][k] = (float) 0.0;

		for(ia=0; ia<Na; ia++){

			theta = ia*d_theta;

			for(j=0; j<N; j++)
				for(k=0; k<N; k++)
					tmp_f[j][k] = m_sharp_deattenuated[ia][k][j];
	
			Rotate2D(m_sharp_rotated[ia], tmp_f, N, theta);

			if(1){
				temp_fn = (char *) calloc(200, sizeof(char)); 	
				strcat(temp_fn, "../../3_Output/Test"); 
				strcat(temp_fn, "_m_sharp_rotated_");
				int2str(stemp, 3, ia);
				strcat(temp_fn, stemp);
				strcat(temp_fn, ".pgm");
				pgm_write_scaled_float2d(temp_fn, m_sharp_rotated[ia], N, N);
				free(temp_fn);
			}

			for(iu=0; iu<Nu; iu++)
				for(k=0; k<N; k++){
					est_M[ia][iu] += m_sharp[ia][iu][k];
					est_P[ia][iu] += m_sharp_deattenuated[ia][iu][k];
				}

			for(iu=0; iu<Nu; iu++)
				for(k=0; k<N; k++)
					ref_f[iu][k] += m_sharp_rotated[ia][iu][k];
		}

		if(1){
			temp_fn = (char *) calloc(200, sizeof(char)); 	
			strcat(temp_fn, "../../3_Output/Test"); 
			strcat(temp_fn, "_ref_f_SumUp.pgm");
			pgm_write_scaled_float2d(temp_fn, ref_f, N, N);
			free(temp_fn);
		}

		for(ia=0; ia<Na; ia++)
			for(iu=0; iu<Nu; iu++){
				if(est_P[ia][iu] > 0.0)
					a[ia][iu] = est_M[ia][iu]/est_P[ia][iu];
				else
					a[ia][iu] = (float) 1.0;

				if(a[ia][iu] >= 1.0)
					a[ia][iu] = (float) 1.0;

				if(a[ia][iu] <= 0.01)
					a[ia][iu] = (float) 0.01;

				a[ia][iu] = - log(a[ia][iu]);
			}

		sum_mean = (float) 0.0;
		for(ia=0; ia<Na; ia++)
			for(iu=0; iu<Nu; iu++)
				sum_mean += a[ia][iu];

		sum_mean = sum_mean/Na;
	
		for(ia=0; ia<Na; ia++){

			sum_each = (float) 0.0;
			for(iu=0; iu<Nu; iu++)
				sum_each += a[ia][iu];

			for(iu=0; iu<Nu; iu++)
				a[ia][iu] = a[ia][iu]/sum_each * sum_mean;
		}

		for(ia=0; ia<Na; ia++){
			for(iu=0; iu<Nu; iu++)
				free(m_sharp[ia][iu]);
			free(m_sharp[ia]);
		}
		free(m_sharp);

		for(ia=0; ia<Na; ia++){
			for(iu=0; iu<Nu; iu++)
				free(m_sharp_rotated[ia][iu]);
			free(m_sharp_rotated[ia]);
		}
		free(m_sharp_rotated);

		for(i=0; i<Na; i++)
			free(est_M[i]);
		free(est_M);

		for(i=0; i<Na; i++)
			free(est_P[i]);
		free(est_P);

		for(i=0; i<N; i++)
			free(tmp_f[i]);
		free(tmp_f);

		for(i=0; i<N; i++)
			free(ref_f[i]);
		free(ref_f);

		for(ia=0; ia<Na; ia++)
			free(idx_A[ia]);
		free(idx_A);

		free(vec_x);
		free(vec_y);
}
