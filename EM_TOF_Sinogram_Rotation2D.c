/*
	EM_TOF_Sinogram_Rotation2D.c
*/
#include	"../1_Header/TOFPET.h"

void		int2str(char *, int, int);
//void		pgm_write_scaled_double2d(char *, double **, int, int);
void		pgm_write_scaled_double3d(char *, double ***, int, int, int, int);
//void		fscanf_double3d(char *, double ***, int, int, int);
//void		fprintf_double3d(char *, double ***, int, int, int);
double		rms_scaled_double3d(double ***, double ***, int, int, int);

void		TOF_Sinogram_Rotation2D(double ****, double ***, Parameter);
void		TOF_BackSinogram_Rotation2D(double ****, double ***, Parameter);

void		EM_TOF_Sinogram_Rotation2D(
	double			****TOF_Sinogram_Y,
	double			***orig_f,
	Parameter		par)
{
	//
	//	Declaration of Input Parameters
	//
		double		***old_f, ***new_f, ***rate_f, ****Y, ****all_1_Y;

		char		*temp_fn, stemp[50];
		
		int			ip, ia, iu, it, i,j,k,n;
	
		double		error;
	
	//
	// 	Memory Allocation
	//
		old_f = (double ***) calloc(par.N_Phantom_Slice, sizeof(double **));
		for(i=0; i<par.N_Phantom_Slice; i++){		
			old_f[i] = (double **) calloc(par.N_Phantom_Pixel, sizeof(double *));
			for(j=0; j<par.N_Phantom_Pixel; j++)
				old_f[i][j] = (double *) calloc(par.N_Phantom_Pixel, sizeof(double));
		}				
		
		new_f = (double ***) calloc(par.N_Phantom_Slice, sizeof(double **));
		for(i=0; i<par.N_Phantom_Slice; i++){		
			new_f[i] = (double **) calloc(par.N_Phantom_Pixel, sizeof(double *));
			for(j=0; j<par.N_Phantom_Pixel; j++)
				new_f[i][j] = (double *) calloc(par.N_Phantom_Pixel, sizeof(double));
		}				
		
		rate_f = (double ***) calloc(par.N_Phantom_Slice, sizeof(double **));
		for(i=0; i<par.N_Phantom_Slice; i++){		
			rate_f[i] = (double **) calloc(par.N_Phantom_Pixel, sizeof(double *));
			for(j=0; j<par.N_Phantom_Pixel; j++)
				rate_f[i][j] = (double *) calloc(par.N_Phantom_Pixel, sizeof(double));
		}				
		
		Y = (double ****) calloc(par.N_Sinogram_Plate, sizeof(double ***));
		for(ip=0; ip<par.N_Sinogram_Plate; ip++){
			Y[ip] = (double ***) calloc(par.N_Sinogram_Angle, sizeof(double **));
			for(ia=0; ia<par.N_Sinogram_Angle; ia++){
				Y[ip][ia] = (double **) calloc(par.N_Sinogram_Bin, sizeof(double *));
				for(iu=0; iu<par.N_Sinogram_Bin; iu++)
					Y[ip][ia][iu] = (double *) calloc(par.N_TOF_TimeSlot, sizeof(double));
			}
		}
		
		all_1_Y = (double ****) calloc(par.N_Sinogram_Plate, sizeof(double ***));
		for(ip=0; ip<par.N_Sinogram_Plate; ip++){
			all_1_Y[ip] = (double ***) calloc(par.N_Sinogram_Angle, sizeof(double **));
			for(ia=0; ia<par.N_Sinogram_Angle; ia++){
				all_1_Y[ip][ia] = (double **) calloc(par.N_Sinogram_Bin, sizeof(double *));
				for(iu=0; iu<par.N_Sinogram_Bin; iu++)
					all_1_Y[ip][ia][iu] = (double *) calloc(par.N_TOF_TimeSlot, sizeof(double));
			}
		}
		
		for(ip=0; ip<par.N_Sinogram_Plate; ip++)
			for(ia=0; ia<par.N_Sinogram_Angle; ia++)
				for(iu=0; iu<par.N_Sinogram_Bin; iu++)
					for(it=0; it<par.N_TOF_TimeSlot; it++)
						all_1_Y[ip][ia][iu][it] = 1.0;
					
	//
	//	Ready for Iteration[
	//
		TOF_BackSinogram_Rotation2D(all_1_Y, rate_f, par); 
		
		if(par.Flag_PGM){
			for(i=0; i<par.N_Phantom_Slice; i++)
				for(j=0; j<par.N_Phantom_Pixel; j++)
					for(k=0; k<par.N_Phantom_Pixel; k++)
						if((j-par.N_Phantom_Pixel/2)*(j-par.N_Phantom_Pixel/2) + (k-par.N_Phantom_Pixel/2)*(k-par.N_Phantom_Pixel/2) >= par.N_Phantom_Pixel*par.N_Phantom_Pixel/4)
							rate_f[i][j][k] = 0.0;
							
			temp_fn = (char *) calloc(200, sizeof(char)); 	
			strcat(temp_fn, par.Output_Header); 	
			strcat(temp_fn, "_TOF_EPD.pgm");
			pgm_write_scaled_double3d(temp_fn, rate_f, par.N_Phantom_SliceRow, par.N_Phantom_SliceCol, par.N_Phantom_Pixel, par.N_Phantom_Pixel); 				
			free(temp_fn);
		}
		
	//
	//	EM Iteration
	//
		for(i=0; i<par.N_Phantom_Slice; i++)
			for(j=0; j<par.N_Phantom_Pixel; j++)
				for(k=0; k<par.N_Phantom_Pixel; k++)
					old_f[i][j][k] = 1.0;

		for(n=0; n<par.N_Iteration; n++){

			TOF_Sinogram_Rotation2D(Y, old_f, par);

			for(ip=0; ip<par.N_Sinogram_Plate; ip++)
				for(ia=0; ia<par.N_Sinogram_Angle; ia++)
					for(iu=0; iu<par.N_Sinogram_Bin; iu++)
						for(it=0; it<par.N_TOF_TimeSlot; it++)
							if(Y[ip][ia][iu][it] > EPSILON)
								Y[ip][ia][iu][it] = TOF_Sinogram_Y[ip][ia][iu][it] / Y[ip][ia][iu][it];
							else
								Y[ip][ia][iu][it] = 0.0;
										
			TOF_BackSinogram_Rotation2D(Y, new_f, par);

			for(i=0; i<par.N_Phantom_Slice; i++)
				for(j=0; j<par.N_Phantom_Pixel; j++)
					for(k=0; k<par.N_Phantom_Pixel; k++)
						if(rate_f[i][j][k] > EPSILON && (j-par.N_Phantom_Pixel/2)*(j-par.N_Phantom_Pixel/2) + (k-par.N_Phantom_Pixel/2)*(k-par.N_Phantom_Pixel/2) < par.N_Phantom_Pixel*par.N_Phantom_Pixel/4)
							new_f[i][j][k] = old_f[i][j][k] * new_f[i][j][k] / rate_f[i][j][k];
						else
							new_f[i][j][k] = 0.0;

			for(i=0; i<par.N_Phantom_Slice; i++)
				for(j=0; j<par.N_Phantom_Pixel; j++)
					for(k=0; k<par.N_Phantom_Pixel; k++)
						old_f[i][j][k] = new_f[i][j][k];

			// Result at Current State
			//if(par.Flag_Phantom){
				error = rms_scaled_double3d(orig_f, new_f, par.N_Phantom_Slice, par.N_Phantom_Pixel, par.N_Phantom_Pixel);
				printf("Error by EM Iteration at %4d-th = %20.12f\n", n, error);
			//}
			//else
				//printf("EM Iteration at %4d-th\n", n);

			temp_fn = (char *) calloc(200, sizeof(char)); 	
			strcat(temp_fn, par.Output_Header); 	
			strcat(temp_fn, "_EM_");
			int2str(stemp, 3, n);	
			strcat(temp_fn, stemp); 	
			strcat(temp_fn, "_");
			int2str(stemp, 6, (int) (error * 10000));	
			strcat(temp_fn, stemp); 	
			strcat(temp_fn, ".pgm");
			pgm_write_scaled_double3d(temp_fn, new_f, par.N_Phantom_SliceRow, par.N_Phantom_SliceCol, par.N_Phantom_Pixel, par.N_Phantom_Pixel); 				
			free(temp_fn);
		}// n-loop
}
