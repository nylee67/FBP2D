/*
	TOF_ParaBeam2D.h

	Data
		mu[iy][ix]			(attenuation Map)
		R_mu[ia][iu]		(Radon transform of mu)
		estR_mu[ia][iu]		(Estimated Radon transform of mu)
		A[ia][iu]			(exp(-R_mu))
		orig_f[iy][ix]		(phantom)
		f[iy][ix]			(reconstructed image)
		y[ia][iu][it]		(TOFSinogram)
		Y[ia][iu]			(Sinogram, Y[ia][iu] = sum_it y[ia][iu][it]) 
		m[ia][iu][it]		(attenuated TOFSinogram = A .* Y)
		M[ia][iu]			(attenuated Sinogram, M[ia][iu] = sum_it m[ia][iu][it])
		da_m[ia][iu][it]	(deattenuated TOFSinogram)
		da_M[ia][iu]		(deattenuated Sinogram, da_M = A^{-1} .* M)
		dc_m[ia][iu][[it]	(deconvolved attenuated TOFSinogram)
		rf[ia][iu][[it]		(deconvolved deattenuated TOFSinogram)
		
	0_FBP_TOF_ParaBeam2D.c
		Read 		orig_f		Read_Image2D(FileName_Phantom, PGM_Header, orig_f, N, Flag_PGM);
		Read 		mu			Read_Image2D(FileName_AttenMap, PGM_Header, mu, N, Flag_PGM);
		Compute		R_mu		ParaBeamP2D(R_mu, mu, Na, N, dd);
									Rotate2D(tmp_mu, mu, N, theta);
		Compute		y		 	TOF_ParaBeamP2D(y, orig_f, Na, N, dd, N_TSF, Sigma_TSF);
									Rotate2D(y[ia], orig_f, N, theta);
									Conv1D(y[ia][iu], N, tsf, N_TSF);
		Compute		m			m = exp(-R_mu) .* y + noise;
	
		if Choice_Method = FBP of TOF_FBP
		
			Compute		estR_mu		DeAtten_NumDiff(estR_mu, m, Na, N, Sigma_TSF, PGM_Header, Flag_PGM);
			Compute		da_m		da_m = exp(estR_mu) *. m;

			// FBP
			Compute		da_M		da_M = sum(da_m); 
			Compute		f 			FBP2D(f, da_M, N, Na, CutOff, PGM_Header, Flag_PGM)
										RampFilter1D(rampfiltered_da_M[ia], da_M[ia], BigNu, CutOff); 
										ParaBeamQ2D(f, rampfiltered_da_M, N, Na)
			
			// TOF_FBP 
			Compute		f 			TOFFBP2D(f, da_m, N, Na, CutOff, PGM_Header, Flag_PGM)
										RampFilter1D(rampfiltered_dam[ia][.][it], da_m[ia][.][it], CutOff, BigNu); 
										TOF_ParaBeamQ2D(f, rampfiltered_da_m, N, Na, N_TSF, Sigma_TSF);
											Conv1D(ramfiltered_da_m[ia][iu], N, tsf, N_TSF);
											Rotate2D(tmp_f, rampfiltered_da_m[ia], N, theta);
											f += tmp_f;
		
		if Choice_Method = TOF_DeConv						
			Compute		rf			DeAtten_DeConv(rf, m, Na, N, iu_0, it_0, rho);
										DeConv1D(dc_m[ia][iu], m[ia][iu], N_TSF, Sigma_TSF); 
			Compute		f 			Roate2D(tmp_f, rf, N, theta);
									f += tmp_f;
*/
#ifndef		_TOF_ParaBeam2D_H_
#define		_TOF_ParaBeam2D_H_

#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<math.h>
#include	<time.h>

/*
typedef struct {
	char	FileName_Phantom[100];		
	char	FileName_AttenMap[100];		
			
	char	Output_Header[100];
	char	PGM_Header[100];
	
	int		Na;					// = N_TOFSinogramAngle
	int		N;					// = N_PhantomPixel = N_TOFSinogramBin = N_TOFSinogramTime
	double	dd;					// = d_PhantomPixel = d_TOFSinogramBin = d_TOFSinogramTime

	int		N_TSF;
	double	Sigma_TSF;

	double	AttenRate;

	double	FOV_Radius;
	double	CutOff;

	int		Flag_PGM;
} Parameter;

typedef struct {
	char	flag;
	double	x0;
	double	y0;
	double	z0;
	double	x1;
	double	y1;
	double	z1;
	int		ip;
	int		ia;
	int		iu;
} SinoLine;
*/

#define		MY_PI		3.14159265358979
#define		MY_C		300000000000
#define		N_Ignored	20
#define		EPSILON		0.00000000000001

#endif // _TOF_ParaBeam2D_H_
