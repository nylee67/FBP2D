/*
	BackSinogram_Rotation2D.c
*/
#include	"../1_Header/TOFPET.h"

void	Expansion1D(double *, double *, int, int, double);
void 	Rotation2D(double **, double **, int, double, double); 
void	pgm_write_scaled_double2d(char *, double **, int, int);

void	BackSinogram_Rotation2D(
	double		***Y,
	double		***f,
	Parameter	par)
{
	int			ip, ia, iu, i;
	double		**mid_f, **slice_f, **tmp_slice_f, **tmp_f, *column_f, *line_f, theta, du, da;
	
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
	
	if(par.N_Sinogram_Plate != par.N_Phantom_Slice){
		puts("par.N_Sinogram_Plate != par.N_Phantom_Slice\n");
		exit(1);
	}
	
	for(ip=0; ip<par.N_Sinogram_Plate; ip++){
		
		for(iu=0; iu<par.N_Sinogram_Bin; iu++)
			for(i=0; i<par.N_Sinogram_Bin; i++)
				slice_f[iu][i] = 0.0; 
				
		for(ia=0; ia<par.N_Sinogram_Angle; ia++){
			
			for(iu=0; iu<par.N_Sinogram_Bin; iu++)
				for(i=0; i<par.N_Sinogram_Bin; i++)
					tmp_f[iu][i] = Y[ip][ia][iu]; 
			
			theta = (ia+0.5)*MY_PI/par.N_Sinogram_Angle;
			Rotation2D(tmp_slice_f, tmp_f, par.N_Sinogram_Bin, theta, par.FOV_Radius);
			
			for(iu=0; iu<par.N_Sinogram_Bin; iu++)
				for(i=0; i<par.N_Sinogram_Bin; i++)
					slice_f[iu][i] += tmp_slice_f[iu][i];
		}
		
		for(iu=0; iu<par.N_Sinogram_Bin; iu++)
			for(i=0; i<par.N_Sinogram_Bin; i++)
				slice_f[iu][i] = slice_f[iu][i] * da;
		
		//pgm_write_scaled_double2d("slice_f.pgm", slice_f, par.N_Sinogram_Bin, par.N_Sinogram_Bin);
		//for(i=0; i<par.N_Sinogram_Bin; i++)
		//	printf("mid_f[ip][100] %30.20f\n", mid_f[100][i]);
		//exit(1);
		
		for(iu=0; iu<par.N_Sinogram_Bin; iu++){
			for(i=0; i<par.N_Sinogram_Bin; i++)
				line_f[i] = slice_f[i][iu];
			
			Expansion1D(column_f, line_f, par.N_Phantom_Pixel, par.N_Sinogram_Bin, par.FOV_Radius); 	
			
			for(i=0; i<par.N_Phantom_Pixel; i++)
				mid_f[i][iu] = column_f[i];
		}
		//pgm_write_scaled_double2d("orig_f.pgm", f[0], par.N_Phantom_Pixel, par.N_Phantom_Pixel);
		//for(i=0; i<par.N_Phantom_Pixel; i++)
		//	printf("f[ip][100] %30.20f\n", f[ip][100][i]);
		
		for(i=0; i<par.N_Phantom_Pixel; i++)
			Expansion1D(f[ip][i], mid_f[i], par.N_Phantom_Pixel, par.N_Sinogram_Bin, par.FOV_Radius); 	
		//pgm_write_scaled_double2d("mid_f.pgm", mid_f, par.N_Phantom_Pixel, par.N_Sinogram_Bin);
		//for(i=0; i<par.N_Sinogram_Bin; i++)
		//	printf("mid_f[ip][100] %30.20f\n", mid_f[100][i]);
		
		//pgm_write_scaled_double2d("slice_f.pgm", slice_f, par.N_Sinogram_Bin, par.N_Sinogram_Bin);
		//for(i=0; i<par.N_Sinogram_Bin; i++)
		//	printf("slice_f[ip][100] %30.20f\n", slice_f[100][i]);
	}
	
	for(iu=0; iu<par.N_Sinogram_Bin; iu++){
		free(tmp_f[iu]); 	free(slice_f[iu]);		free(tmp_slice_f[iu]);
	}
	free(tmp_f); 	free(slice_f);		free(tmp_slice_f);
	
	for(i=0; i<par.N_Phantom_Pixel; i++)
		free(mid_f[i]);
	free(mid_f);
	
	free(line_f);	free(column_f);
}