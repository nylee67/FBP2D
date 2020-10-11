/*
	Rotate2D.c
*/

#include	"TOF_ParaBeam2D.h"

void 	Rotate2D(
	float		**g,				// Rotated One:		N x N
	float		**f,				// Original:		N x N
	int			N,
	double		theta)
{
	// Variables
		int				i, j;
		double			ch, sh;

		int				i_0, i_1, j_0, j_1;
		double			x_r, y_r, x_o, y_o, c_0, c_1, r_0, r_1;
		float			temp_00, temp_01, temp_10, temp_11;

	// Rotation
		ch = cos(theta);		sh = sin(theta);
	
		for(i=0; i<N; i++){			y_o = (N/2.0 - i - 0.5);
			for(j=0; j<N; j++){		x_o = (j - N/2.0 + 0.5);
				x_r = x_o*ch + y_o*sh; 	
				y_r =  - x_o*sh + y_o*ch;	

				if(x_r*x_r + y_r*y_r < N*N/4){

					i_0 = (int) (N/2.0 - y_r - 0.5); 	i_1 = i_0 + 1;
					r_0 = (N/2.0-i_0-0.5) - y_r; 		r_1 = 1.0 - r_0;

					j_0 = (int) (N/2.0 + x_r - 0.5); 	j_1 = j_0 + 1;
					c_0 = x_r - (-N/2.0+j_0+0.5); 		c_1 = 1.0 - c_0;
			
					if(0 <= i_0 && i_1 < N && 0 <= j_0 && j_1 < N){	
						temp_00 = (float) (r_1*c_1);
						temp_01 = (float) (r_1*c_0);
						temp_10 = (float) (r_0*c_1);
						temp_11 = (float) (r_0*c_0);

						g[i][j]  = f[i_0][j_0]*temp_00;
						g[i][j] += f[i_0][j_1]*temp_01;
						g[i][j] += f[i_1][j_0]*temp_10;
						g[i][j] += f[i_1][j_1]*temp_11;
					}
					else
						g[i][j] = (float) 0.0;	
				}
				else
					g[i][j] = (float) 0.0;
			} //j-loop
		} //i-loop
}
