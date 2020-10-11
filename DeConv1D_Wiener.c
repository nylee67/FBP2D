/*
	DeConv1D_Wiener.c
*/

#include	"TOF_ParaBeam2D.h"

void	dfft(double *, double *, int, int);

void 	DeConv1D_Wiener(
	float		*f,
	float		*y,
	double		*rTSF,
	double		*iTSF,
	int			N,
	int			N_TSF,
	int			BigN,
	double 		Sigma_Wiener)
{
	//	Variables
		int 	i;
		double	*real_f, *imag_f, *real_g, *imag_g, temp;
		int		Nt = N + N_TSF - 1;

		real_f = (double *) calloc(BigN, sizeof(double));
		imag_f = (double *) calloc(BigN, sizeof(double));
		real_g = (double *) calloc(BigN, sizeof(double));
		imag_g = (double *) calloc(BigN, sizeof(double));

		for(i=0; i<Nt; i++)
			real_f[i] = y[i];

		dfft(real_f, imag_f, BigN, -1);

		for(i=0; i<BigN; i++){
			temp = rTSF[i]*rTSF[i] + iTSF[i]*iTSF[i] + Sigma_Wiener*Sigma_Wiener;
			real_g[i] = (rTSF[i]*real_f[i] + iTSF[i]*imag_f[i])/temp;
			imag_g[i] = (rTSF[i]*imag_f[i] - iTSF[i]*real_f[i])/temp;
		}
		
		dfft(real_g, imag_g, BigN, 1);

		for(i=0; i<N; i++)
			f[i] = real_g[i+N_TSF/2];

		free(real_f);
		free(imag_f);
		free(real_g);
		free(imag_g);
}
