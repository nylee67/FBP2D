/*
	Read_Image2D.c 
*/

#include	"TOF_ParaBeam2D.h"

void	fscanf_float2d(char *, float **, int, int);
void	int2str(char *, int, int);
void	pgm16_write_fixed_float2d(char *, float **, int, int, float, float);

void	Read_Image2D(char *FileName_Image, char *PGMHeader, float **f, int N_Pixel, int Flag_PGM)
{
		int 	i,j;
		char	*temp_fn;
		float	max;
		
		fscanf_float2d(FileName_Image, f, N_Pixel, N_Pixel);

		max = (float) 0.0;
		for(i=0; i<N_Pixel; i++)
			for(j=0; j<N_Pixel; j++)
				if(max < f[i][j])
					max = f[i][j];
					
		for(i=0; i<N_Pixel; i++)
			for(j=0; j<N_Pixel; j++)
				f[i][j] = f[i][j]/max;	

		if(Flag_PGM){
			temp_fn = (char *) calloc(200, sizeof(char)); 	
			strcat(temp_fn, PGMHeader); 
			strcat(temp_fn, ".pgm");
        	pgm16_write_fixed_float2d(temp_fn, f, N_Pixel, N_Pixel, (float) 0.0, (float) 1.0);
			free(temp_fn);
		}
}
