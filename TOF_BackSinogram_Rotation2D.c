/*
	TOF_BackSinogram_Rotation2D.c
*/
#include	"../1_Header/TOFPET.h"

void	FullConv_TOF(double *, double *, double *, int, int);

void	Expansion1D(double *, double *, int, int, double);
void 	Rotation2D(double **, double **, int, double, double); 
void	pgm_write_scaled_double2d(char *, double **, int, int);

void	TOF_BackSinogram_Rotation2D(
	double		****Y,
	double		***f,
	Parameter	par)
{
	int			ip, ia, iu, it, i;
	double		**mid_f, **slice_f, **tmp_slice_f, **tmp_f, *column_f, *line_f, *time_f, theta, du, da;
	
	da = MY_PI/par.N_Sinogram_Angle;
	du = 2*par.FOV_Radius/par.N_Sinogram_Bin;
	
	mid_f = (double **) calloc(par.N_Phantom_Pixel, sizeof(double *));
	for(i=0; i<par.N_Phantom_Pixel; i++)
		mid_f[i] = (double *) calloc(par.N_Sinogram_Bin, sizeof(double));
	
	slice_f = (double **) calloc(par.N_Sinogram_Bin, sizeof(double *));
	for(iu=0; iu<par.N_Sinogram_Bin; iu++)
		slice_f[iu] = (double *) calloc(par.N_Sinogram_Bin, sizeof(double));

	tmp_slice_f = (double **) calloc(par.N_Sinogram_Bin, sizeof(double *));
	for(iu=0; iu<par.N_Sinogram_Bin; iu++)
		tmp_slice_f[iu] = (double *) calloc(par.N_Sinogram_Bin, sizeof(double));
	
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
		
		for(iu=0; iu<par.N_Sinogram_Bin; iu++)
			for(i=0; i<par.N_Sinogram_Bin; i++)
				slice_f[iu][i] = 0.0;
			
		for(ia=0; ia<par.N_Sinogram_Angle; ia++){
			
			//if(ia == 20){
			//	pgm_write_scaled_double2d("before_Y[0][20].pgm", Y[ip][ia], par.N_Sinogram_Bin, par.N_TOF_TimeSlot);
			//}
			
			for(iu=0; iu<par.N_Sinogram_Bin; iu++){
				
				for(it=0; it<par.N_TOF_TimeSlot + 2*par.N_TOF_PSFHalfSupp; it++)
					time_f[it] = 0.0;
			
				FullConv_TOF(Y[ip][ia][iu], time_f, par.TOF_PSF, par.N_TOF_PSFHalfSupp, par.N_TOF_TimeSlot);
				
				for(it=0; it<par.N_TOF_TimeSlot; it++)
					Y[ip][ia][iu][it] = time_f[it+par.N_TOF_PSFHalfSupp];
				
				Expansion1D(tmp_f[iu], Y[ip][ia][iu], par.N_Sinogram_Bin, par.N_TOF_TimeSlot, par.FOV_Radius);
			}
			
			theta = (ia+0.5)*MY_PI/par.N_Sinogram_Angle;
			Rotation2D(tmp_slice_f, tmp_f, par.N_Sinogram_Bin, theta, par.FOV_Radius);
			
			//if(ia == 20){
			//	pgm_write_scaled_double2d("after_Y[0][20].pgm", Y[ip][ia], par.N_Sinogram_Bin, par.N_TOF_TimeSlot);
			//	pgm_write_scaled_double2d("tmp_f.pgm", tmp_f, par.N_Sinogram_Bin, par.N_Sinogram_Bin);
			//	exit(1);
			//}
			for(iu=0; iu<par.N_Sinogram_Bin; iu++)
				for(i=0; i<par.N_Sinogram_Bin; i++)
					slice_f[iu][i] += tmp_slice_f[iu][i];
		}
		
		for(iu=0; iu<par.N_Sinogram_Bin; iu++)
			for(i=0; i<par.N_Sinogram_Bin; i++)
				slice_f[iu][i] = slice_f[iu][i] * da;
		
		for(iu=0; iu<par.N_Sinogram_Bin; iu++){
			for(i=0; i<par.N_Sinogram_Bin; i++)
				line_f[i] = slice_f[i][iu];
			
			Expansion1D(column_f, line_f, par.N_Phantom_Pixel, par.N_Sinogram_Bin, par.FOV_Radius); 	
			
			for(i=0; i<par.N_Phantom_Pixel; i++)
				mid_f[i][iu] = column_f[i];
		}
		
		for(i=0; i<par.N_Phantom_Pixel; i++)
			Expansion1D(f[ip][i], mid_f[i], par.N_Phantom_Pixel, par.N_Sinogram_Bin, par.FOV_Radius); 	
	}
	
	for(iu=0; iu<par.N_Sinogram_Bin; iu++){
		free(tmp_f[iu]); 	free(slice_f[iu]);		free(tmp_slice_f[iu]);
	}
	free(tmp_f); 	free(slice_f);		free(tmp_slice_f);
	
	for(i=0; i<par.N_Phantom_Pixel; i++)
		free(mid_f[i]);
	free(mid_f);
	
	free(line_f);	free(column_f);		free(time_f);
}

void	FullConv_TOF(double *Y, double *line_Y, double *PSF, int N_Half, int N)
{
		int			i,j;
		
		for(i=0; i<N; i++)
			for(j=-N_Half; j<=N_Half; j++)
				line_Y[i+j+N_Half] += Y[i] * PSF[j+N_Half];
}