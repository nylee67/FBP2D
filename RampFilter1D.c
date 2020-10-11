/*
	RampFilter1D.c
*/

#include	"TOF_ParaBeam2D.h"

void	dfft(double *, double *, int, int);

void	RampFilter1D(
	double		*realF,
	double		*imagF,
	int			BigNu,
	double		CutOff)
{
	//
	//	Declaration of Variables
	//
		int		j, Nsupp, Nsupp_min;
		double	temp;

	//
	//	Parameters Setup
	//
		if(CutOff < 100.0)
			Nsupp = (int) (CutOff*BigNu/2);		
		else
			Nsupp = BigNu/2;		

		Nsupp_min = BigNu/2;
		if(Nsupp_min > Nsupp)
			Nsupp_min = Nsupp;

		dfft(realF, imagF, BigNu, -1);

		realF[0] = imagF[0] = 0.0;			// multiplication by | w |
		for(j=1; j<Nsupp_min; j++) {
			if(CutOff < 100.0)
				temp = (0.5 + 0.5 * cos(j * MY_PI / Nsupp))*j;
			else
				temp = 1.0*j;
			realF[j] *= temp;
			imagF[j] *= temp;
			realF[BigNu-j] *= temp;
			imagF[BigNu-j] *= temp;
		}
		for(j=Nsupp; j<BigNu/2; j++){
			realF[j] = 0.0;
			imagF[j] = 0.0;
			realF[BigNu-j] = 0.0;
			imagF[BigNu-j] = 0.0;
		}
		realF[BigNu/2] = imagF[BigNu/2] = 0.0;

		dfft(realF, imagF, BigNu, 1);
}
