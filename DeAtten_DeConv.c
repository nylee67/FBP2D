/*
	EstimateAtten.c
*/
#include	"../1_Header/TOFPET.h"

void	BackSinogram_Rotation2D(double ***, double ***, Parameter);
void	dfft(double *, double *, int, int);

void 	pgm_write_scaled_double2d(char	*, double **, int, int);

void	Smooth(double *, double *, double *, int, int);
void	Median3(double *, int);

void	EstimateAtten(
	double		***g,
	double		****Y,
	Parameter	par)
{
	int			ip, ia, iu, it, BigNu, Nsupp, Nsupp_min, u_istep;
	double		sum, da, du, dt, t, u, Dm, temp, *line, *smooth_line, *u_psf, ****Smooth_Y, ***Du_Y, ***Da_Y, ***Dt_m, ***Da_m, ***Du_m, ***DDtu_m, Ja, Huu, Ju, Hua, Haa, *Real_F, *Imag_F, *weight;
	double		cutoff = 1.0;
	
	//
	//	Memory Allocation
	//
		BigNu = 1;
		while(BigNu < par.N_Sinogram_Bin)
			BigNu = 2*BigNu;
		
		if(cutoff < 100.0)
			Nsupp = (int) (cutoff*BigNu/2);		// cutoff: 0.0 - 1.0
		else
			Nsupp = BigNu/2;

		Nsupp_min = BigNu/2;
		if(Nsupp_min > Nsupp)
			Nsupp_min = Nsupp;

        Real_F = (double *) calloc(BigNu, sizeof(double));
        Imag_F = (double *) calloc(BigNu, sizeof(double));
		
        weight = (double *) calloc(BigNu, sizeof(double));
		
		if(0){
			for(iu=1; iu<Nsupp_min; iu++) {
				weight[iu] = (0.5 + 0.5 * cos(iu * MY_PI / Nsupp));
				weight[BigNu-iu] = (0.5 + 0.5 * cos(iu * MY_PI / Nsupp));
			}
			
			for(iu=Nsupp_min; iu<BigNu/2; iu++){
				weight[iu] = 0.0;
				weight[BigNu-iu] = 0.0;
			}
			weight[0] = 1.0; weight[BigNu/2] = 0.0;
		}
		else
			for(iu=0; iu<BigNu; iu++)
				weight[iu] = 1.0;
			
		Smooth_Y = (double ****) calloc(par.N_Sinogram_Plate, sizeof(double ***));
		for(ip=0; ip<par.N_Sinogram_Plate; ip++){
			Smooth_Y[ip] = (double ***) calloc(par.N_Sinogram_Angle, sizeof(double **));
			for(ia=0; ia<par.N_Sinogram_Angle; ia++){
				Smooth_Y[ip][ia] = (double **) calloc(par.N_Sinogram_Bin, sizeof(double *));
				for(iu=0; iu<par.N_Sinogram_Bin; iu++)
					Smooth_Y[ip][ia][iu] = (double *) calloc(par.N_TOF_TimeSlot, sizeof(double));
			}
		}
		
		Du_Y = (double ***) calloc(par.N_Sinogram_Plate, sizeof(double **));
		for(ip=0; ip<par.N_Sinogram_Plate; ip++){
			Du_Y[ip] = (double **) calloc(par.N_Sinogram_Angle, sizeof(double *));
			for(ia=0; ia<par.N_Sinogram_Angle; ia++)
				Du_Y[ip][ia] = (double *) calloc(par.N_Sinogram_Bin, sizeof(double));	
		}

		Da_Y = (double ***) calloc(par.N_Sinogram_Plate, sizeof(double **));
		for(ip=0; ip<par.N_Sinogram_Plate; ip++){
			Da_Y[ip] = (double **) calloc(par.N_Sinogram_Angle, sizeof(double *));
			for(ia=0; ia<par.N_Sinogram_Angle; ia++)
				Da_Y[ip][ia] = (double *) calloc(par.N_Sinogram_Bin, sizeof(double));	
		}
		
		Dt_m = (double ***) calloc(par.N_Sinogram_Angle, sizeof(double **));
		for(ia=0; ia<par.N_Sinogram_Angle; ia++){
			Dt_m[ia] = (double **) calloc(par.N_Sinogram_Bin, sizeof(double *));
			for(iu=0; iu<par.N_Sinogram_Bin; iu++)
				Dt_m[ia][iu] = (double *) calloc(par.N_TOF_TimeSlot, sizeof(double));	
		}
	
		Da_m = (double ***) calloc(par.N_Sinogram_Angle, sizeof(double **));
		for(ia=0; ia<par.N_Sinogram_Angle; ia++){
			Da_m[ia] = (double **) calloc(par.N_Sinogram_Bin, sizeof(double *));
			for(iu=0; iu<par.N_Sinogram_Bin; iu++)
				Da_m[ia][iu] = (double *) calloc(par.N_TOF_TimeSlot, sizeof(double));	
		}
	
		Du_m = (double ***) calloc(par.N_Sinogram_Angle, sizeof(double **));
		for(ia=0; ia<par.N_Sinogram_Angle; ia++){
			Du_m[ia] = (double **) calloc(par.N_Sinogram_Bin, sizeof(double *));
			for(iu=0; iu<par.N_Sinogram_Bin; iu++)
				Du_m[ia][iu] = (double *) calloc(par.N_TOF_TimeSlot, sizeof(double));	
		}
	
		DDtu_m = (double ***) calloc(par.N_Sinogram_Angle, sizeof(double **));
		for(ia=0; ia<par.N_Sinogram_Angle; ia++){
			DDtu_m[ia] = (double **) calloc(par.N_Sinogram_Bin, sizeof(double *));
			for(iu=0; iu<par.N_Sinogram_Bin; iu++)
				DDtu_m[ia][iu] = (double *) calloc(par.N_TOF_TimeSlot, sizeof(double));	
		}
	
	//
	//	Compute ***Du_Y from ****Y
	//
		da = MY_PI/par.N_Sinogram_Angle;
		du = 2*par.FOV_Radius/par.N_Sinogram_Bin;
		dt = par.TOF_time_step;
		dt = du;
		
		pgm_write_scaled_double2d("Before Smoothing Y[0][10][*][*].pgm", Y[0][10], par.N_Sinogram_Bin, par.N_TOF_TimeSlot);
		
		// Smoothing Y[][][][]
		u_istep = 2;
		u_psf = (double *)calloc(2*u_istep+1, sizeof(double));
		
		sum = 0.0;
		for(iu=-u_istep; iu<=u_istep; iu++){
			u_psf[iu+u_istep] = exp(-iu*iu/(2*0.3*u_istep*0.3*u_istep)); 
			sum += u_psf[iu+u_istep];
		}
		for(iu=-u_istep; iu<=u_istep; iu++)
			u_psf[iu+u_istep] = u_psf[iu+u_istep]/sum; 
				
		line = (double *)calloc(par.N_Sinogram_Bin+2*u_istep, sizeof(double));
		smooth_line = (double *)calloc(par.N_Sinogram_Bin, sizeof(double));
		
		for(ip=0; ip<par.N_Sinogram_Plate; ip++)
			for(ia=0; ia<par.N_Sinogram_Angle; ia++)
				for(it=0; it<par.N_TOF_TimeSlot; it++){
					for(iu=0; iu<u_istep; iu++){
						line[iu] = 0.0;
						line[par.N_Sinogram_Bin + 2*u_istep - 1 - iu] = 0.0;
					}
					for(iu=0; iu<par.N_Sinogram_Bin; iu++){
						line[iu+u_istep] = Y[ip][ia][iu][it];
						smooth_line[iu] = 0.0;
					}
					
					Smooth(smooth_line, line, u_psf, u_istep, par.N_Sinogram_Bin); 
					
					for(iu=0; iu<par.N_Sinogram_Bin; iu++)
						Smooth_Y[ip][ia][iu][it] = smooth_line[iu];
				}
		free(line);
		free(smooth_line);
		
		/*
		a_istep = 1;
		a_psf = (double *)calloc(2*a_istep+1, sizeof(double));
		
		sum = 0.0;
		for(iu=-a_istep; iu<=a_istep; iu++){
			a_psf[iu+a_istep] = exp(-iu*iu/(2*0.3*a_istep*0.3*a_istep)); 
			sum += a_psf[iu+a_istep];
		}
		for(iu=-a_istep; iu<=a_istep; iu++)
			a_psf[iu+a_istep] = a_psf[iu+u_istep]/sum; 
		
		line = (double *)calloc(par.N_Sinogram_Angle+2*a_istep, sizeof(double));
		smooth_line = (double *)calloc(par.N_Sinogram_Angle, sizeof(double));
		for(ip=0; ip<par.N_Sinogram_Plate; ip++)
			for(iu=0; iu<par.N_Sinogram_Bin; iu++)
				for(it=0; it<par.N_TOF_TimeSlot; it++){
					for(ia=0; ia<a_istep; ia++){
						line[ia] = Y[ip][par.N_Sinogram_Angle-a_istep+ia][par.N_Sinogram_Bin-1-iu][par.N_Sinogram_Bin-1-it];
						line[par.N_Sinogram_Angle+2*a_istep-1-ia] = Y[ip][a_istep-ia][par.N_Sinogram_Bin-1-iu][par.N_Sinogram_Bin-1-it];
					}
					for(ia=0; ia<par.N_Sinogram_Angle; ia++){
						line[ia+a_istep] = Y[ip][ia][iu][it];
						smooth_line[ia] = 0.0;
					}
					
					Smooth(smooth_line, line, a_psf, a_istep, par.N_Sinogram_Angle); 
					
					for(ia=0; ia<par.N_Sinogram_Angle; ia++)
						Smooth_Y[ip][ia][iu][it] = smooth_line[ia];
				}
		free(line);
		free(smooth_line);
		*/
		
		pgm_write_scaled_double2d("After Smoothing Y[0][10][*][*].pgm", Smooth_Y[0][10], par.N_Sinogram_Bin, par.N_TOF_TimeSlot);
		
		for(ip=0; ip<par.N_Sinogram_Plate; ip++){
			
			// Derivative
			for(ia=0; ia<par.N_Sinogram_Angle; ia++)
				for(iu=0; iu<par.N_Sinogram_Bin; iu++){
					Dt_m[ia][iu][0] = 0.0;
					for(it=1; it<par.N_TOF_TimeSlot-1; it++)
						Dt_m[ia][iu][it] = (Smooth_Y[ip][ia][iu][it+1] - Smooth_Y[ip][ia][iu][it-1])/(2*dt);
					Dt_m[ia][iu][par.N_TOF_TimeSlot-1] = 0.0;
				}
			pgm_write_scaled_double2d("Dt_m[0][*][*].pgm", Dt_m[0], par.N_Sinogram_Bin, par.N_TOF_TimeSlot);
				
			for(ia=0; ia<par.N_Sinogram_Angle; ia++)
				for(it=0; it<par.N_TOF_TimeSlot; it++){
					Du_m[ia][0][it] = 0.0;
					for(iu=1; iu<par.N_Sinogram_Bin-1; iu++)
						Du_m[ia][iu][it] = (Smooth_Y[ip][ia][iu+1][it] - Smooth_Y[ip][ia][iu-1][it])/(2*du);
					Du_m[ia][par.N_Sinogram_Bin-1][it] = 0.0;
				}
			pgm_write_scaled_double2d("Du_m[0][*][*].pgm", Du_m[0], par.N_Sinogram_Bin, par.N_TOF_TimeSlot);

			for(iu=0; iu<par.N_Sinogram_Bin; iu++)
				for(it=0; it<par.N_TOF_TimeSlot; it++){
					Da_m[0][iu][it] = (Smooth_Y[ip][1][iu][it] - Smooth_Y[ip][par.N_Sinogram_Angle-1][par.N_Sinogram_Bin-1-iu][par.N_TOF_TimeSlot-1-it])/(2*da);
					for(ia=1; ia<par.N_Sinogram_Angle-1; ia++)
						Da_m[ia][iu][it] = (Smooth_Y[ip][ia+1][iu][it] - Smooth_Y[ip][ia-1][iu][it])/(2*da);
					Da_m[par.N_Sinogram_Angle-1][iu][it] = (Smooth_Y[ip][0][par.N_Sinogram_Bin-1-iu][par.N_TOF_TimeSlot-1-it] - Smooth_Y[ip][par.N_Sinogram_Angle-2][iu][it])/(2*da);
				}
			pgm_write_scaled_double2d("Da_m[0][*][*].pgm", Da_m[0], par.N_Sinogram_Bin, par.N_TOF_TimeSlot);
				
			for(ia=0; ia<par.N_Sinogram_Angle; ia++){
				for(it=1; it<par.N_TOF_TimeSlot-1; it++)
					for(iu=1; iu<par.N_Sinogram_Bin-1; iu++)
						DDtu_m[ia][iu][it] = (Smooth_Y[ip][ia][iu+1][it+1] - Smooth_Y[ip][ia][iu-1][it+1] - Smooth_Y[ip][ia][iu+1][it-1] + Smooth_Y[ip][ia][iu-1][it-1])/(4*du*dt);
					
				for(iu=0; iu<par.N_Sinogram_Bin; iu++){
					DDtu_m[ia][iu][0] = 0.0;
					DDtu_m[ia][iu][par.N_TOF_TimeSlot-1] = 0.0;
				}
				
				for(it=0; it<par.N_TOF_TimeSlot; it++){
					DDtu_m[ia][0][it] = 0.0;
					DDtu_m[ia][par.N_Sinogram_Bin-1][it] = 0.0;
				}
			}
			pgm_write_scaled_double2d("DDtu_m[0][*][*].pgm", DDtu_m[0], par.N_Sinogram_Bin, par.N_TOF_TimeSlot);
			
			for(ia=0; ia<par.N_Sinogram_Angle; ia++)
				for(iu=0; iu<par.N_Sinogram_Bin; iu++){
					
					u = - par.FOV_Radius + (iu+0.5)*du;
					
					Ja = 0.0;
					Ju = 0.0;
					Huu = 0.0;
					Hua = 0.0;
					Haa = 0.0;
					for(it=0; it<par.N_TOF_TimeSlot; it++){
						t = (it-par.N_TOF_TimeSlot/2)*dt;
						
						temp = Smooth_Y[ip][ia][iu][it]*t + par.TOF_Sigma*par.TOF_Sigma*Dt_m[ia][iu][it];
						Dm = t*Du_m[ia][iu][it] + Da_m[ia][iu][it] - u*Dt_m[ia][iu][it] + par.TOF_Sigma*par.TOF_Sigma*DDtu_m[ia][iu][it];
						
						Ja += Dm*Smooth_Y[ip][ia][iu][it]*dt;
						Ju += Dm*temp*dt;
						Huu += temp*temp*dt;
						Hua += Smooth_Y[ip][ia][iu][it]*temp*dt;
						Haa += Smooth_Y[ip][ia][iu][it]*Smooth_Y[ip][ia][iu][it]*dt;
					}
					
					if(Huu*Haa - Hua*Hua != 0.0){
						Du_Y[ip][ia][iu] = - (Ju*Haa - Ja*Hua)/(Huu*Haa - Hua*Hua);
						Da_Y[ip][ia][iu] = - (Ja*Huu - Ju*Hua)/(Huu*Haa - Hua*Hua);
					}
					else{
						Du_Y[ip][ia][iu] = 0.0;
						Da_Y[ip][ia][iu] = 0.0;
					}
				}
		}
		pgm_write_scaled_double2d("Before_Du_Y.pgm", Du_Y[0], par.N_Sinogram_Angle, par.N_Sinogram_Bin);
		pgm_write_scaled_double2d("Before_Da_Y.pgm", Da_Y[0], par.N_Sinogram_Angle, par.N_Sinogram_Bin);
		for(iu=0; iu<par.N_Sinogram_Bin; iu++)
			printf("Before Du_Y %10d %40.20f%40.20f\n", iu, Du_Y[0][20][iu], Du_Y[0][21][iu]);
		for(ia=0; ia<par.N_Sinogram_Angle; ia++)
			printf("Before Da_Y %10d %60.20f\n", ia, Da_Y[0][ia][90]);
		
		// Median3
		for(ip=0; ip<par.N_Sinogram_Plate; ip++)
			for(ia=0; ia<par.N_Sinogram_Angle; ia++){
				Median3(Du_Y[ip][ia], par.N_Sinogram_Bin);	
				Median3(Du_Y[ip][ia], par.N_Sinogram_Bin);	
				Median3(Du_Y[ip][ia], par.N_Sinogram_Bin);	
			}
			
		pgm_write_scaled_double2d("After_Du_Y.pgm", Du_Y[0], par.N_Sinogram_Angle, par.N_Sinogram_Bin);
		pgm_write_scaled_double2d("After_Da_Y.pgm", Da_Y[0], par.N_Sinogram_Angle, par.N_Sinogram_Bin);
		for(iu=0; iu<par.N_Sinogram_Bin; iu++)
			printf("After Du_Y %10d %40.20f%40.20f\n", iu, Du_Y[0][20][iu], Du_Y[0][21][iu]);
		for(ia=0; ia<par.N_Sinogram_Angle; ia++)
			printf("After Da_Y %10d %60.20f\n", ia, Da_Y[0][ia][90]);
		
	//
	//	Compute ***g from ***Du_Y
	//
		for(ip=0; ip<par.N_Sinogram_Plate; ip++)
			for(ia=0; ia<par.N_Sinogram_Angle; ia++){
				for(iu=0; iu<BigNu; iu++){
					Real_F[iu] = 0.0;
					Imag_F[iu] = 0.0;
				}
				
				for(iu=0; iu<par.N_Sinogram_Bin; iu++)
					Real_F[iu] = Du_Y[ip][ia][iu];	
				
				dfft(Real_F, Imag_F, BigNu, -1);
				
				for(iu=1; iu<BigNu/2; iu++){
					temp = - Real_F[iu]/(2*MY_PI)*weight[iu];
					Real_F[iu] = Imag_F[iu]/(2*MY_PI)*weight[iu];
					Imag_F[iu] = temp;
					
					temp = Real_F[BigNu-iu]/(2*MY_PI)*weight[BigNu-iu];
					Real_F[BigNu-iu] = -Imag_F[BigNu-iu]/(2*MY_PI)*weight[BigNu-iu];
					Imag_F[BigNu-iu] = temp;
				}
				
				temp = - Real_F[0]/(2*MY_PI)*weight[0];
				Real_F[0] = Imag_F[0]/(2*MY_PI)*weight[0];
				Imag_F[0] = temp;
				
				temp = -Real_F[BigNu/2]/(2*MY_PI)*weight[BigNu/2];
				Real_F[BigNu/2] = Imag_F[BigNu/2]/(2*MY_PI)*weight[BigNu/2];
				Imag_F[BigNu/2] = temp;
					
				dfft(Real_F, Imag_F, BigNu, 1);
				
				for(iu=0; iu<par.N_Sinogram_Bin; iu++)
					Du_Y[ip][ia][iu] = Real_F[iu] * par.N_Sinogram_Bin;	
			}
		
		BackSinogram_Rotation2D(Du_Y, g, par);
		
		pgm_write_scaled_double2d("Hilbert_Du_Y.pgm", Du_Y[0], par.N_Sinogram_Angle, par.N_Sinogram_Bin);
	//
	//	Freeing Memory
	// ***Du_Y, ***Dt_m, ***Da_m, ***Du_m, ***DDtu_
		for(ip=0; ip<par.N_Sinogram_Plate; ip++){
			for(ia=0; ia<par.N_Sinogram_Angle; ia++)
				free(Du_Y[ip][ia]);
			free(Du_Y[ip]);
		}
		free(Du_Y);
			
		for(ia=0; ia<par.N_Sinogram_Angle; ia++){
			for(iu=0; iu<par.N_Sinogram_Bin; iu++){
				free(Dt_m[ia][iu]);
				free(Da_m[ia][iu]);
				free(Du_m[ia][iu]);
				free(DDtu_m[ia][iu]);
			}
			free(Dt_m[ia]);
			free(Da_m[ia]);
			free(Du_m[ia]);
			free(DDtu_m[ia]);
		}
		free(Dt_m);
		free(Da_m);
		free(Du_m);
		free(DDtu_m);
		
		free(Real_F);	free(Imag_F);	free(weight);	
		
		for(ip=0; ip<par.N_Sinogram_Plate; ip++){
			for(ia=0; ia<par.N_Sinogram_Angle; ia++){
				for(iu=0; iu<par.N_Sinogram_Bin; iu++)
					free(Smooth_Y[ip][ia][iu]);
				free(Smooth_Y[ip][ia]);
			}
			free(Smooth_Y[ip]);
		}
		free(Smooth_Y);
}

void	Smooth(double *Y, double *line_Y, double *PSF, int N_Half, int N)
{
		int			i,j;
		
		for(i=0; i<N; i++)
			for(j=-N_Half; j<=N_Half; j++)
				Y[i] += line_Y[i+j+N_Half]*PSF[j+N_Half];
}

void	Median3(double *Y, int N)
{
		int			i;
		double		sum = 0.0, var = 0.0, *median_Y;
		
		median_Y = (double *)calloc(N, sizeof(double));
		
		for(i=1; i<N/2; i++)
			if(Y[i-1] == 0.0 && Y[i] != 0.0){
				Y[i] = 0.0;
				break;
			}
			
		for(i=1; i<N/2; i++)
			if(Y[N-i-1] != 0.0 && Y[N-i] == 0.0){
				Y[N-i-1] = 0.0;
				break;
			}
			
		median_Y[0] = 0.0;
		for(i=1; i<N-1; i++){
			median_Y[i] = Y[i];
			sum += Y[i];
		}
		median_Y[N-1] = 0.0;
		sum = sum/N;
		
		for(i=1; i<N-1; i++)
			var += fabs(Y[i]); 
		var = var/N;
		
		for(i=1; i<N-1; i++)
			if(Y[i] < -10*var || Y[i] > 10*var)
				median_Y[i] = 0.0;
		for(i=0; i<N; i++)
			Y[i] = median_Y[i];
		free(median_Y);
}