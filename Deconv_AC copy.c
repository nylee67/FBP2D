/*
	Deconv_AC.c
*/
#include	"TOF_ParaBeam2D.h"

void	int2str(char *, int, int);
void	pgm_write_scaled_float2d(char *, float **, int, int);

void	dfft(double *, double *, int, int);

void	DeConv1D_LR(float *f, float *y, int N, float *tsf, int N_TSF);
void	DeConv1D_Wiener(float *f, float *y, double *rTSF, double *iTSF, int N, int N_TSF, int BigN, double Sigma_Wiener);

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
	int 	ia, iu, it, i, BigNu;
	
	float 	***m_sharp, *ref_f, *tsf, temp;
	int		**idx_A;
	double 	*rTSF, *iTSF, sum;
	
	double	du = dd;
	double	d_theta = MY_PI/Na;

	double 	theta_0 = ia_0*d_theta;
	double 	u_0 = (iu_0-Nu/2)*du;
	double 	t_0 = (it_0-N/2)*dd;
	
	double 	x_0 = u_0*cos(theta_0) - t_0*sin(theta_0);
	double  y_0 = u_0*sin(theta_0) + t_0*cos(theta_0);
	
	double  ref_value;
	
	double	theta, u;
	double  t, x, y;
	
	double 	theta_ast, u_ast;
	int 	ia_ast, iu_ast;
	
	double 	t_0ast, t_ast;
	int 	it_0ast, it_ast;
	
	double 	max;
	char 	stemp[100], *temp_fn;

	float	tmp_a, tmp_b, tmp_c;

	//
	//	Memory Allocation
	//
		m_sharp = (float ***) calloc(Na, sizeof(float **));
		for(ia=0; ia<Na; ia++){
			m_sharp[ia] = (float **) calloc(Nu, sizeof(float *));
			for(iu=0; iu<Nu; iu++)
				m_sharp[ia][iu] = (float *) calloc(Nu, sizeof(float));
		}
			
		ref_f = (float **) calloc(N, sizeof(float *));
		for(i=0; i<N; i++)
			ref_f[i] = (float *) calloc(N, sizeof(float));

		idx_A = (int **) calloc(Na, sizeof(int *));
		for(i=0; i<Na; i++)
			idx_A[i] = (int *) calloc(Nu, sizeof(int));

	//
	//	Deconvolution: m --> m_sharp
	//
		BigNu = 1;
		while(BigNu < Nt)
			BigNu = 2*BigNu;

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

				for(i=0; i<N; i++)
					if(m_sharp[ia][iu][i] < 0.0)
						m_sharp[ia][iu][i] = (float) 0.0;			
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

		ref_value = m_sharp[ia_0][iu_0][it_0];
		if(1)
			printf("m_sharp at ref. point = %20.8f\n", ref_value); 
		
		free(tsf);
		free(rTSF);
		free(iTSF);
	
	//
	//	
	//
		for(ia=0; ia<Na; ia++){
			theta = ia*d_theta;

			for(iu=0; iu<Nu; iu++){
				u = (-Nu/2.0 + iu + 0.5)*du;

				max = 0.0;
				for(i=0; i<N; i++)
					if(m_sharp[ia][iu][i] > max){
						max = m_sharp[ia][iu][i];
						it = i;
					}
				
				if(max > 0.0){// if max. activity time is found

					t = (N/2-it)*dd;
					x = u*cos(theta) - t*sin(theta);
					y = u*sin(theta) + t*cos(theta);
 
					// find ell_(theta_ast, u_ast) connects (x_0,y_0) and (x,y)
					if(x != x_0)
						theta_ast = MY_PI/2 + atan((y-y_0)/(x-x_0));
					else
						theta_ast = 0.0;
					
					u_ast = x*cos(theta_ast) + y*sin(theta_ast);

					t_ast = -x*sin(theta_ast) + y*cos(theta_ast);
					t_0ast = -x_0*sin(theta_ast) + y_0*cos(theta_ast);

					ia_ast = (int) (theta_ast/d_theta);
					if(ia_ast >= Na || ia_ast < 0)
						ia_ast = 0;
					iu_ast = (int) (u_ast/dd + N/2);
					it_ast = (int) (N/2 - t_ast/dd);
					it_0ast = (int) (N/2 - t_0ast/dd);

					tmp_a = m_sharp[ia][iu][it];

					if(0 < it_0ast && it_0ast < N-1 && 0 < iu_ast && iu_ast < N-1){ 
						tmp_b  = m_sharp[ia_ast][iu_ast-1][it_0ast-1] + m_sharp[ia_ast][iu_ast-1][it_0ast] + m_sharp[ia_ast][iu_ast-1][it_0ast+1];
						tmp_b += m_sharp[ia_ast][iu_ast  ][it_0ast-1] + m_sharp[ia_ast][iu_ast  ][it_0ast] + m_sharp[ia_ast][iu_ast  ][it_0ast+1];
						tmp_b += m_sharp[ia_ast][iu_ast+1][it_0ast-1] + m_sharp[ia_ast][iu_ast+1][it_0ast] + m_sharp[ia_ast][iu_ast+1][it_0ast+1];
						tmp_b = tmp_b/9;
					}
					else
						tmp_b = m_sharp[ia_ast][iu_ast][it_0ast];
						//tmp_b = ref_value;

					if(0 < it_ast && it_ast < N-1 && 0 < iu_ast && iu_ast < N-1){ 
						tmp_c  = m_sharp[ia_ast][iu_ast-1][it_ast-1] + m_sharp[ia_ast][iu_ast-1][it_ast] + m_sharp[ia_ast][iu_ast-1][it_ast+1];
						tmp_c += m_sharp[ia_ast][iu_ast  ][it_ast-1] + m_sharp[ia_ast][iu_ast  ][it_ast] + m_sharp[ia_ast][iu_ast  ][it_ast+1];
						tmp_c += m_sharp[ia_ast][iu_ast+1][it_ast-1] + m_sharp[ia_ast][iu_ast+1][it_ast] + m_sharp[ia_ast][iu_ast+1][it_ast+1];
						tmp_c = tmp_c/9;
					}
					else
						tmp_c = m_sharp[ia_ast][iu_ast][it_ast];
						//tmp_c = ref_value;

					if(tmp_c > 0.0)
						a[ia][iu] = tmp_a/ref_value * tmp_b/tmp_c;
					else
						a[ia][iu] = (float) 1.0;

					if(a[ia][iu] >= 1.0)
						a[ia][iu] = (float) 1.0;

					if(a[ia][iu] <= 0.05)
						a[ia][iu] = (float) 0.05;

				}
				else
					a[ia][iu] = (float) 1.0;

				if(0)
					printf("Here %12.4f %4d %4d %4d %4d %4d %4d %4d\n", a[ia][iu], ia, iu, ia_ast, iu_ast, it_ast, it_0ast, it);
			}
		}
		a[ia_0][iu_0] = (float) 1.0;
		
		for(ia=0; ia<Na; ia++){
			for(iu=0; iu<Nu; iu++)
				free(m_sharp[ia][iu]);
			free(m_sharp[ia]);
		}
		free(m_sharp);

		for(i=0; i<N; i++)
			free(ref_f[i]);
		free(ref_f);

		for(ia=0; ia<Na; ia++)
			free(idx_A[ia]);
		free(idx_A);

		for(ia=0; ia<Na; ia++)
			for(iu=0; iu<Nu; iu++)
				a[ia][iu] = - log(a[ia][iu]);
}
