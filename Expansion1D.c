/*
	Expansion1D.c
*/

#include	"../1_Header/TOFPET.h"
void 	Expansion1D(
	double		*g,				// Expanded One:	Ng
	double		*f,				// Original:		Nf
	int			Ng,
	int			Nf,
	double		FOV_R)
{
	//
	//	Variables
	//
		int				i, i_0, i_1;
		double			dg, df, x, c_0, c_1;

	//
	// 	Rotation
	//
		dg = 2.0*FOV_R/Ng;		df = 2.0*FOV_R/Nf; 		
	
		for(i=0; i<Ng; i++){
			x = (-Ng/2.0 + i + 0.5) * dg;	
			i_0 = (int) (x/df + Nf/2.0 - 0.5);		i_1 = i_0 + 1;
			c_0 = x - (-Nf/2.0+i_0+0.5)*df;			c_1 = (-Nf/2.0+i_1+0.5)*df - x;
			
			if(0 <= i_0 && i_1 < Nf)
				g[i] = f[i_0]*c_1/df + f[i_1]*c_0/df;	
			else
				g[i] = 0.0;
		}
}