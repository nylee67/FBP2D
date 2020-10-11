/*
	TOF_Sinogram_Read.c
*/

#include	"../1_Header/TOFPET.h"

void 	pgm_write_fixed_double2d(char	*, double **, int, int, double, double);
void	int2str(char *,	int, int);

void	TOF_Sinogram_Read(char *FileName, double ****Y, double ***sum_Y, Parameter par)
{
	//
	//	Declaration of Input Parameters
	//
		FILE				*fp;
		int					ip, ia, iu, it;
		double				max;
		char				*temp_fn, stemp[100]; 
		
		if((fp = fopen(FileName, "r")) != NULL){
			for(ip=0; ip<par.N_Sinogram_Plate; ip++)
				for(ia=0; ia<par.N_Sinogram_Angle; ia++)
					for(iu=0; iu<par.N_Sinogram_Bin; iu++)
						for(it=0; it<par.N_TOF_TimeSlot; it++)
							fscanf(fp, "%lf", &Y[ip][ia][iu][it]);
			fclose(fp);
		}
		else{
			printf("File Opeening Error %s\n", FileName);
			exit(1);
		}
		
		for(ip=0; ip<par.N_Sinogram_Plate; ip++)
			for(ia=0; ia<par.N_Sinogram_Angle; ia++)
				for(iu=0; iu<par.N_Sinogram_Bin; iu++){
					sum_Y[ip][ia][iu] = 0.0;
					for(it=0; it<par.N_TOF_TimeSlot; it++)
						sum_Y[ip][ia][iu] += Y[ip][ia][iu][it];
				}
				
		if(par.Flag_PGM){ // sum_Y
			max = 0.0;
			for(ip=0; ip<par.N_Sinogram_Plate; ip++)
				for(ia=0; ia<par.N_Sinogram_Angle; ia++)
					for(iu=0; iu<par.N_Sinogram_Bin; iu++)
						if(max < sum_Y[ip][ia][iu])
							max = sum_Y[ip][ia][iu];
				
			for(ip=0; ip<par.N_Sinogram_Plate; ip++){
				temp_fn = (char *) calloc(200, sizeof(char)); 	
				strcat(temp_fn, par.Output_Header); 
				strcat(temp_fn, "_Sinogram_");
				int2str(stemp, 3, ip);
				strcat(temp_fn, stemp);
				strcat(temp_fn, "_fixed.pgm");
				pgm_write_fixed_double2d(temp_fn, sum_Y[ip], par.N_Sinogram_Angle, par.N_Sinogram_Bin, 0.0, max); 
				free(temp_fn);
			}
		}
		
		if(par.Flag_PGM * par.Flag_TOF_PGM){ // Y
			max = 0.0;
			for(ip=0; ip<par.N_Sinogram_Plate; ip++)
				for(ia=0; ia<par.N_Sinogram_Angle; ia++)
					for(iu=0; iu<par.N_Sinogram_Bin; iu++)
						for(it=0; it<par.N_TOF_TimeSlot; it++)
							if(max < Y[ip][ia][iu][it])
								max = Y[ip][ia][iu][it];
				
			for(ip=0; ip<par.N_Sinogram_Plate; ip++)
				for(ia=0; ia<par.N_Sinogram_Angle; ia++){
					temp_fn = (char *) calloc(200, sizeof(char)); 	
					strcat(temp_fn, par.Output_Header); 
					strcat(temp_fn, "_TOF_Sinogram_");
					int2str(stemp, 3, ip);
					strcat(temp_fn, stemp);
					strcat(temp_fn, "_");
					strcat(temp_fn, "_");
					int2str(stemp, 3, ia);
					strcat(temp_fn, stemp);
					strcat(temp_fn, "_");
					strcat(temp_fn, "_fixed.pgm");
					pgm_write_fixed_double2d(temp_fn, Y[ip][ia], par.N_Sinogram_Bin, par.N_TOF_TimeSlot, 0.0, max); 
					free(temp_fn);
				}
		}
}