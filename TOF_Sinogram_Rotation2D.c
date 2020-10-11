/*
	TOF_Sinogram_Rotation2D.c
*/
#include	"../1_Header/TOFPET.h"

void	Expansion1D(double *, double *, int, int, double);
void 	Rotation2D(double **, double **, int, double, double); 
void	pgm_write_scaled_double2d(char *, double **, int, int);

void	ValidConv_TOF(double *, double *, double *, int, int);

void	TOF_Sinogram_Rotation2D(
	double		****Y,
	double		***f,
	Parameter	par)
{
	int			ip, ia, iu, it, i;
	double		**mid_f, **slice_f, **tmp_f, *column_f, *line_f, *time_f, theta, du;
	
	du = 2*par.FOV_Radius/par.N_Sinogram_Bin;
	
	mid_f = (double **) calloc(par.N_Phantom_Pixel, sizeof(double *));
	for(i=0; i<par.N_Phantom_Pixel; i++)
		mid_f[i] = (double *) calloc(par.N_Sinogram_Bin, sizeof(double));
	
	slice_f = (double **) calloc(par.N_Sinogram_Bin, sizeof(double *));
	for(iu=0; iu<par.N_Sinogram_Bin; iu++)
		slice_f[iu] = (double *) calloc(par.N_Sinogram_Bin, sizeof(double));
	
	tmp_f = (double **) calloc(par.N_Sinogram_Bin, sizeof(double *));
	for(iu=0; iu<par.N_Sinogram_Bin; iu++)
		tmp_f[iu] = (double *) calloc(par.N_Sinogram_Bin, sizeof(double));
	
	column_f = (double *) calloc(par.N_Phantom_Pixel, sizeof(double));
	
	line_f = (double *) calloc(par.N_Sinogram_Bin, sizeof(double));
	
	time_f = (double *) calloc(par.N_TOF_TimeSlot + 2*par.N_TOF_PSFHalfSupp, sizeof(double));
	
	if(par.N_Sinogram_Plate != par.N_Phantom_Slice){
		puts("par.N_Sinogram_Plate != par.N_Phantom_Slice\n");
		exit(1);
	}
	
	for(ip=0; ip<par.N_Sinogram_Plate; ip++){
		
		for(i=0; i<par.N_Phantom_Pixel; i++)
			Expansion1D(mid_f[i], f[ip][i], par.N_Sinogram_Bin, par.N_Phantom_Pixel, par.FOV_Radius); 	
		
		for(iu=0; iu<par.N_Sinogram_Bin; iu++){
			for(i=0; i<par.N_Phantom_Pixel; i++)
				column_f[i] = mid_f[i][iu];
			
			Expansion1D(line_f, column_f, par.N_Sinogram_Bin, par.N_Phantom_Pixel, par.FOV_Radius); 	
			
			for(i=0; i<par.N_Sinogram_Bin; i++)
				slice_f[i][iu] = line_f[i];
		}
		
		for(ia=0; ia<par.N_Sinogram_Angle; ia++){
			theta = (ia+0.5)*MY_PI/par.N_Sinogram_Angle;
			Rotation2D(tmp_f, slice_f, par.N_Sinogram_Bin, -theta, par.FOV_Radius);
			
			for(iu=0; iu<par.N_Sinogram_Bin; iu++){
				Expansion1D(Y[ip][ia][iu], tmp_f[iu], par.N_TOF_TimeSlot, par.N_Sinogram_Bin, par.FOV_Radius);
				
				for(it=0; it<par.N_TOF_TimeSlot; it++)
					time_f[it+par.N_TOF_PSFHalfSupp] = Y[ip][ia][iu][it] * du;
				
				ValidConv_TOF(Y[ip][ia][iu], time_f, par.TOF_PSF, par.N_TOF_PSFHalfSupp, par.N_TOF_TimeSlot);
			}
		}
	}
	
	for(iu=0; iu<par.N_Sinogram_Bin; iu++){
		free(tmp_f[iu]); 	free(slice_f[iu]);
	}
	free(tmp_f); 	free(slice_f);
	
	for(i=0; i<par.N_Phantom_Pixel; i++)
		free(mid_f[i]);
	free(mid_f);
	
	free(line_f);	free(column_f);		free(time_f);
}

void	ValidConv_TOF(double *Y, double *line_Y, double *PSF, int N_Half, int N)
{
		int			i,j;
		
		for(i=0; i<N; i++)
			for(j=-N_Half; j<=N_Half; j++)
				Y[i] += line_Y[i+j+N_Half]*PSF[j+N_Half];
}