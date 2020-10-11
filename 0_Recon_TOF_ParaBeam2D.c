/*
	0_Recon_TOF_ParaBeam2D.c
*/

#include	"TOF_ParaBeam2D.h"

float	gasdev(long *);

void	int2str(char *, int, int);
void	pgm16_write_fixed_float2d(char *, float **, int, int, float, float);

void	Read_Image2D(char *FileName_Image, char *PGMHeader, float **f, int N_Pixel, int Flag_PGM);
void	ParaBeamP2D(float **R_mu, float **mu, int Na, int N, double dd); 

void	TOF_ParaBeamP2D(float ***y, float **f, int Na, int N, double dd, int N_TSF, double Sigma_TSF); 

void 	Diff_AC(float **R_mu, float **a, float ***m, int Na, int N, double dd, int N_TSF, double Sigma_TSF);
void 	Deconv_AC(int ia_0, int iu_0, int it_0, float **a, float ***m, int Na, int N, double dd, int N_TSF, double Sigma_TSF, double Sigma_Wiener);

void	FBP2D(float **f, float **Y, int Na, int N, double dd, double CutOff);
void	TOFFBP2D(float **f, float ***y, int Na, int N, double dd, int N_TSF, double Sigma_TSF, double CutOff);
void	TOFDeConvBP2D(float **f, float ***y, int Na, int N, double dd, int N_TSF, double Sigma_TSF, double Sigma_Wiener);
void	TOFEM2D(float **f, float ***y, int Na, int N, double dd, int N_TSF, double Sigma_TSF, char *PGM_Header);

int		main()
{
	//
	//	Declaration of Input Parameters
	//
		float 		**orig_f, **f, **mu, **R_mu, ***p, ***m, **Y, **a, ***da_m;
		int			i, j, k; 
		char		stemp[100], *temp_fn;
		float		max;

		long		seed = (long) 777, *idum;

	//
	//	Parameter Setup
	//
		char		FileName_Phantom[100] 		= "../../2_Phantom/ColdRod_001x201.pht3d";
		//char		FileName_Phantom[100] 		= "../../2_Phantom/HotRod_001x201.pht3d";
		char		FileName_AttenMap[100]		= "../../2_Phantom/AttenMap_001x201.pht3d";

		char		PGM_Header[100]				= "../../3_Output/ColdRod_8_FBP_2";

		int			N				= 201;
		int 		Na				= 8;
		double 		dd 				= 2.0;
		double 		mean;
		clock_t		start, end;

		double 		AttenRate 		= 0.0;
		
		int			Choice_TSF;					// 0: 300ps, 		1: 600ps, 			2: 3ns
		int			Choice_Noise;				// 0: 0%,			1: 1%,				2: 5%
		int			Choice_Attenuation;			// 0: NoDeAtten		1: Diff,			2: Deconv
		int			Choice_Method;				// 0: FBP,			1: TOFFBP,			2: TOFDeconvBP,		3: TOFEM

		double		Sigma_TSF, Sigma_Noise_InPercent;
		int			N_TSF, Nt;
		
		int 		Flag_EPPS		= 0;		// Add External Point Positron Source

		double 		CutOff 			= 0.75;
		double 		Sigma_Wiener 	= 0.0001;
 
		float		max_orig		= (float) 1.0;
		float		max_Rm			= (float) 1.6275;
		//float		max_FBP			= (float) 75.0;
		//float		max_TOFFBP		= (float) 375.0;
		//float		max_TOFDeconvBP	= (float) 1.0;
		//float		max_TOFEM		= (float) 1.0;
		float		max_P			= (float) 42.0;
		float		max_M			= (float) 42.0;
		float		max_p			= (float) 0.85;
		float		max_m			= (float) 0.85;

		idum = (long *) calloc(1, sizeof(long));
		*idum = seed;

		int			iy_0		= 15;
		int			ix_0 		= N/2;
		
		int 		ia_0		= Na/2;
		int 		iu_0 		= N-iy_0-1;
		int 		it_0		= N/2;
		
	//
	//	Choice
	//
		Choice_TSF 			= 2;
		Choice_Noise 		= 0;
		Choice_Attenuation	= 0;
		Choice_Method		= 0;
		
	//
	//	Choice-dependent Parameters
	//
		if(Choice_TSF == 0)
			Sigma_TSF = 38.22;
		if(Choice_TSF == 1)
			Sigma_TSF = 76.44;
		if(Choice_TSF == 2)
			Sigma_TSF = 382.2;
		
		N_TSF = 2*ceil(2.5*Sigma_TSF/2.0) + 1;
		printf("N_TSF = %8d\n", N_TSF);
		Nt = N + N_TSF - 1;
					
		if(Choice_Noise == 0)
			Sigma_Noise_InPercent = (double) 0.0;
		if(Choice_Noise == 1)
			Sigma_Noise_InPercent = (double) 0.01;
		if(Choice_Noise == 5)
			Sigma_Noise_InPercent = (double) 0.05;

	//
	//	Memory Allocation
	//
		orig_f = (float **) calloc(N, sizeof(float *));
		for(i=0; i<N; i++)		
			orig_f[i] = (float *) calloc(N, sizeof(float));
	
		f = (float **) calloc(N, sizeof(float *));
		for(i=0; i<N; i++)		
			f[i] = (float *) calloc(N, sizeof(float));
	
		mu = (float **) calloc(N, sizeof(float *));
		for(i=0; i<N; i++)		
			mu[i] = (float *) calloc(N, sizeof(float));
	
		R_mu = (float **) calloc(Na, sizeof(float *));
		for(i=0; i<Na; i++)		
			R_mu[i] = (float *) calloc(N, sizeof(float));
	
		p = (float ***) calloc(Na, sizeof(float **));
		for(i=0; i<Na; i++){		
			p[i] = (float **) calloc(N, sizeof(float *));
			for(j=0; j<N; j++)		
				p[i][j] = (float *) calloc(Nt, sizeof(float));
		}
	
		m = (float ***) calloc(Na, sizeof(float **));
		for(i=0; i<Na; i++){		
			m[i] = (float **) calloc(N, sizeof(float *));
			for(j=0; j<N; j++)		
				m[i][j] = (float *) calloc(Nt, sizeof(float));
		}
	
		a = (float **) calloc(Na, sizeof(float *));
		for(i=0; i<Na; i++)		
			a[i] = (float *) calloc(N, sizeof(float));
	
		da_m = (float ***) calloc(Na, sizeof(float **));
		for(i=0; i<Na; i++){		
			da_m[i] = (float **) calloc(N, sizeof(float *));
			for(j=0; j<N; j++)		
				da_m[i][j] = (float *) calloc(Nt, sizeof(float));
		}
	
		Y = (float **) calloc(Na, sizeof(float *));
		for(i=0; i<Na; i++)		
			Y[i] = (float *) calloc(N, sizeof(float));
	
	//
	//	Read Phantom
	//
		strcpy(stemp, PGM_Header);
		strcat(stemp, "_Phantom");
		Read_Image2D(FileName_Phantom, stemp, orig_f, N, 1);

		if(Flag_EPPS){// Add External Positron Point Source
			orig_f[iy_0][ix_0] = (float) 1.0;
		
			max = (float) 0.0;
			for(i=0; i<N; i++)
				for(j=0; j<N; j++)
					if(max < orig_f[i][j])
						max = orig_f[i][j];

			strcpy(stemp, PGM_Header);
			strcat(stemp, "_PhantomEPPS.pgm");
			pgm16_write_fixed_float2d(stemp, orig_f, N, N, (float) 0.0, max_orig);

			printf("Max of PhantomEPPS %12.4f\n", max_orig);
		}
		
	//
	//	Read AttenMap
	//
		strcpy(stemp, PGM_Header);
		strcat(stemp, "_AttenMap");
		Read_Image2D(FileName_AttenMap, stemp, mu, N, 1);

		if(0){
			for(i=0; i<N; i++)
				for(j=0; j<N; j++){
					orig_f[i][j] = (float) 0.0;
					if((N/2-i-40)*(N/2-i-40) + (-N/2+j-40)*(-N/2+j-40) < N*N/4*0.01)
						orig_f[i][j] = (float) 1.0;
				}

			orig_f[iy_0][ix_0] = (float) 1.0;

			strcpy(stemp, PGM_Header);
			strcat(stemp, "_Phantom.pgm");
			pgm16_write_fixed_float2d(stemp, orig_f, N, N, (float) 0.0, (float) 1.0);
		}

		if(0){
			for(i=0; i<N; i++)
				for(j=0; j<N; j++){
					mu[i][j] = (float) 0.0;
					if((N/2-i)*(N/2-i) + (-N/2+j)*(-N/2+j) < N*N/4*0.49)
						mu[i][j] = (float) 1.0;
				}

			strcpy(stemp, PGM_Header);
			strcat(stemp, "_AttenMap.pgm");
			pgm16_write_fixed_float2d(stemp, mu, N, N, (float) 0.0, (float) 1.0);
		}

	//
	//	Sinogram_Atten, AttenRate
	//
		ParaBeamP2D(R_mu, mu, Na, N, dd); 

		for(i=0; i<Na; i++)
			for(j=0; j<N; j++)
				R_mu[i][j] = AttenRate * R_mu[i][j];

		if(1){
			max = (float) 0.0;
			for(i=0; i<Na; i++)
				for(j=0; j<N; j++)
					if(max < R_mu[i][j])
						max = R_mu[i][j];

			max_Rm = max;
			strcpy(stemp, PGM_Header);
			strcat(stemp, "_AttenRate.pgm");
			pgm16_write_fixed_float2d(stemp, R_mu, Na, N, (float) 0.0, max_Rm);

			printf("Max of Sinogram_Atten %12.4f\n", max_Rm);
		}

		if(0){
			max = (float) 1.0;
			for(i=0; i<Na; i++)
				for(j=0; j<N; j++)
					if(max < exp(R_mu[i][j]))
						max = exp(R_mu[i][j]);

			printf("Max of Attenuation Effect %12.4f\n", max);
		}

	//
	//	Observation
	//
		TOF_ParaBeamP2D(p, orig_f, Na, N, dd, N_TSF, Sigma_TSF);

		if(1){
			max = (float) 0.0;
			for(i=0; i<Na; i++)
				for(j=0; j<N; j++)
					for(k=0; k<Nt; k++)
						if(max < p[i][j][k])
							max = p[i][j][k];

			printf("Max of Unattenuated p %12.4f\n", max);
			max_p = max;

			for(i=0; i<Na; i++){
				temp_fn = (char *) calloc(200, sizeof(char)); 	
				strcat(temp_fn, PGM_Header); 
				strcat(temp_fn, "_Unattenuated_TOF_p_");
				int2str(stemp, 3, i);
				strcat(temp_fn, stemp);
				strcat(temp_fn, ".pgm");
       			pgm16_write_fixed_float2d(temp_fn, p[i], N, Nt, (float) 0.0, max_p);
				free(temp_fn);
			}
		}

		if(1){
			for(i=0; i<Na; i++)
				for(j=0; j<N; j++){
					Y[i][j] = (float) 0.0;
					for(k=0; k<Nt; k++)
						Y[i][j] += p[i][j][k];
				}

			max = (float) 0.0;
			for(i=0; i<Na; i++)
				for(j=0; j<N; j++)
					if(max < Y[i][j])
						max = Y[i][j];

			max_P = max;
			printf("Max of UnAttenuated Y %12.4f\n", max_P);

			temp_fn = (char *) calloc(200, sizeof(char)); 	
			strcat(temp_fn, PGM_Header); 
			strcat(temp_fn, "_UnAttenuated_P.pgm");
       		pgm16_write_fixed_float2d(temp_fn, Y, Na, N, (float) 0.0, max_P);
			free(temp_fn);
		}

		for(i=0; i<Na; i++)
			for(j=0; j<N; j++)
				for(k=0; k<Nt; k++)
					m[i][j][k] = exp(-R_mu[i][j]) * p[i][j][k];

	//
	//	Adding Noise
	//
		if(1){
			mean = 0.0;
			for(i=0; i<Na; i++)
				for(j=0; j<N; j++)
					for(k=0; k<Nt; k++)
						mean += m[i][j][k];

			mean = mean/(Na*N*N);
		
			for(i=0; i<Na; i++)
				for(j=0; j<N; j++)
					for(k=0; k<Nt; k++) {
						m[i][j][k] = m[i][j][k] + Sigma_Noise_InPercent * mean * gasdev(idum);
						if(m[i][j][k] < 0.0)
							m[i][j][k] = (float) 0.0;
					}
		}

		if(1){
			max = (float) 0.0;
			for(i=0; i<Na; i++)
				for(j=0; j<N; j++)
					for(k=0; k<Nt; k++)
						if(max < m[i][j][k])
							max = m[i][j][k];

			printf("Max of Attenuated m %12.4f\n", max);

			for(i=0; i<Na; i++){
				temp_fn = (char *) calloc(200, sizeof(char)); 	
				strcat(temp_fn, PGM_Header); 
				strcat(temp_fn, "_Attenuated_TOF_m_");
				int2str(stemp, 3, i);
				strcat(temp_fn, stemp);
				strcat(temp_fn, ".pgm");
       			pgm16_write_fixed_float2d(temp_fn, m[i], N, Nt, (float) 0.0, max_p);
				free(temp_fn);
			}
		}

		if(1){
			for(i=0; i<Na; i++)
				for(j=0; j<N; j++){
					Y[i][j] = (float) 0.0;
					for(k=0; k<Nt; k++)
						Y[i][j] += m[i][j][k];
				}

			max = (float) 0.0;
			for(i=0; i<Na; i++)
				for(j=0; j<N; j++)
					if(max < Y[i][j])
						max = Y[i][j];

			max_M = max;
			printf("Max of Attenuated M %12.4f\n", max_M);

			temp_fn = (char *) calloc(200, sizeof(char)); 	
			strcat(temp_fn, PGM_Header); 
			strcat(temp_fn, "_Attenuated_M.pgm");
       		pgm16_write_fixed_float2d(temp_fn, Y, Na, N, (float) 0.0, max_P);
			free(temp_fn);
		}

	//
	//	Attenuation Correction: m --> a
	//
		if(Choice_Attenuation == 0){
			puts("No Attenuation Correction\n");

			for(i=0; i<Na; i++)
				for(j=0; j<N; j++)
					a[i][j] = (float) 0.0;
		}
		
		if(Choice_Attenuation == 1){
			puts("Diff-based AC\n");
			
			Diff_AC(R_mu, a, m, Na, N, dd, N_TSF, Sigma_TSF);
		}
		
		if(Choice_Attenuation == 2){
			puts("Deconv-based AC\n");

			Deconv_AC(ia_0, iu_0, it_0, a, m, Na, N, dd, N_TSF, Sigma_TSF, Sigma_Wiener);
		}

		if(0){
			max = (float) 0.0;
			for(i=0; i<Na; i++)
				for(j=0; j<N; j++)
					if(max < a[i][j])
						max = a[i][j];

			printf("Max of Atten Sino a %12.4f\n", max);

			temp_fn = (char *) calloc(200, sizeof(char)); 	
			strcat(temp_fn, PGM_Header); 
			if(Choice_Attenuation == 0)
				strcat(temp_fn, "_AttenSino_0_a.pgm");
			if(Choice_Attenuation == 1)
				strcat(temp_fn, "_AttenSino_1_a.pgm");
			if(Choice_Attenuation == 2)
				strcat(temp_fn, "_AttenSino_2_a.pgm");
       		pgm16_write_fixed_float2d(temp_fn, a, Na, N, (float) 0.0, max);
			free(temp_fn);
		}

		for(i=0; i<Na; i++)
			for(j=0; j<N; j++)
				for(k=0; k<Nt; k++)
					da_m[i][j][k] = exp(a[i][j])*m[i][j][k];
				
		if(1){
			max = (float) 0.0;
			for(i=0; i<Na; i++)
				for(j=0; j<N; j++)
					for(k=0; k<Nt; k++)
						if(max < m[i][j][k])
							max = da_m[i][j][k];

			printf("Max of DeAttenuated m %12.4f\n", max);

			for(i=0; i<Na; i++){
				temp_fn = (char *) calloc(200, sizeof(char)); 	
				strcat(temp_fn, PGM_Header); 
				strcat(temp_fn, "_DeAttenuated_TOF_m_");
				int2str(stemp, 3, i);
				strcat(temp_fn, stemp);
				strcat(temp_fn, ".pgm");
       			pgm16_write_fixed_float2d(temp_fn, da_m[i], N, Nt, (float) 0.0, max_m);
				free(temp_fn);
			}
		}

	//
	//	Reconstruction
	//
		if(Choice_Method == 0){// FBP on Y
			start = clock();
			for(i=0; i<Na; i++)
				for(j=0; j<N; j++){
					Y[i][j] = (float) 0.0;
					for(k=0; k<Nt; k++)
						Y[i][j] += da_m[i][j][k];
				}
			
			FBP2D(f, Y, Na, N, dd, CutOff);
			end = clock();

			if(1){
				max = (float) 0.0;
				for(i=0; i<N; i++)
					for(j=0; j<N; j++)
						if(max < f[i][j])
							max = f[i][j];

				printf("Max of FBP f %12.4f and time = %20.12f\n", max, (end-start)/(double) CLOCKS_PER_SEC);

				temp_fn = (char *) calloc(200, sizeof(char)); 	
				strcat(temp_fn, PGM_Header); 
				if(Choice_Attenuation == 0)
					strcat(temp_fn, "_FBP_0.pgm");
				if(Choice_Attenuation == 1)
					strcat(temp_fn, "_FBP_1.pgm");
				if(Choice_Attenuation == 2)
					strcat(temp_fn, "_FBP_2.pgm");
       			pgm16_write_fixed_float2d(temp_fn, f, N, N, (float) 0.0, max);
				free(temp_fn);
			}
		}

		if(Choice_Method == 1){// TOFFBP on da_m
			start = clock();
			TOFFBP2D(f, da_m, Na, N, dd, N_TSF, Sigma_TSF, CutOff);
			end = clock();

			if(1){
				max = (float) 0.0;
				for(i=0; i<N; i++)
					for(j=0; j<N; j++)
						if(max < f[i][j])
							max = f[i][j];

				printf("Max of TOFFBP f %12.4f and time = %20.12f\n", max, (end-start)/(double) CLOCKS_PER_SEC);

				temp_fn = (char *) calloc(200, sizeof(char)); 	
				strcat(temp_fn, PGM_Header); 
				if(Choice_Attenuation == 0)
					strcat(temp_fn, "_TOFFBP_0.pgm");
				if(Choice_Attenuation == 1)
					strcat(temp_fn, "_TOFFBP_1.pgm");
				if(Choice_Attenuation == 2)
					strcat(temp_fn, "_TOFFBP_2.pgm");
       			pgm16_write_fixed_float2d(temp_fn, f, N, N, (float) 0.0, max);
				free(temp_fn);
			}
		}
	
		if(Choice_Method == 2){// TOFDeConvBP on da_m
			start = clock();
			TOFDeConvBP2D(f, da_m, Na, N, dd, N_TSF, Sigma_TSF, Sigma_Wiener);
			end = clock();

			if(1){
				max = (float) 0.0;
				for(i=0; i<N; i++)
					for(j=0; j<N; j++)
						if(max < f[i][j])
							max = f[i][j];

				printf("Max of TOFDeConvBP f %12.4f and time = %20.12f\n", max, (end-start)/(double) CLOCKS_PER_SEC);

				temp_fn = (char *) calloc(200, sizeof(char)); 	
				strcat(temp_fn, PGM_Header); 
				if(Choice_Attenuation == 0)
					strcat(temp_fn, "_TOFDeconvBP_0.pgm");
				if(Choice_Attenuation == 1)
					strcat(temp_fn, "_TOFDeconvBP_1.pgm");
				if(Choice_Attenuation == 2)
					strcat(temp_fn, "_TOFDeconvBP_2.pgm");
       			pgm16_write_fixed_float2d(temp_fn, f, N, N, (float) 0.0, max);
				free(temp_fn);
			}
		}
		
		if(Choice_Method == 3){// TOFEM on da_m
			start = clock();
			TOFEM2D(f, da_m, Na, N, dd, N_TSF, Sigma_TSF, PGM_Header);
			end = clock();

			if(1){
				max = (float) 0.0;
				for(i=0; i<N; i++)
					for(j=0; j<N; j++)
						if(max < f[i][j])
							max = f[i][j];

				printf("Max of TOFEM f %12.4f and time = %20.12f\n", max, (end-start)/(double) CLOCKS_PER_SEC);

				temp_fn = (char *) calloc(200, sizeof(char)); 	
				strcat(temp_fn, PGM_Header); 
				if(Choice_Attenuation == 0)
					strcat(temp_fn, "_TOFEM_0.pgm");
				if(Choice_Attenuation == 1)
					strcat(temp_fn, "_TOFEM_1.pgm");
				if(Choice_Attenuation == 2)
					strcat(temp_fn, "_TOFEM_2.pgm");
       			pgm16_write_fixed_float2d(temp_fn, f, N, N, (float) 0.0, max);
				free(temp_fn);
			}
		}
		
	return 0;
}
