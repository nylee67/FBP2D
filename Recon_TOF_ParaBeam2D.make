gcc -Wall -c Read_Image2D.c
gcc -Wall -c Rotate2D.c
gcc -Wall -c ParaBeamP2D.c
gcc -Wall -c ParaBeamQ2D.c
gcc -Wall -c ConvP1D.c
gcc -Wall -c ConvQ1D.c
gcc -Wall -c DeConv1D_LR.c
gcc -Wall -c DeConv1D_Wiener.c
gcc -Wall -c TOF_ParaBeamP2D.c
gcc -Wall -c TOF_ParaBeamQ2D.c
gcc -Wall -c RampFilter1D.c
gcc -Wall -c Diff_AC.c
gcc -Wall -c CGM_AttenDI2D.c
gcc -Wall -c Deconv_AC.c
gcc -Wall -c FBP2D.c
gcc -Wall -c TOFFBP2D.c
gcc -Wall -c TOFDeConvBP2D.c
gcc -Wall -c TOFEM2D.c
ar ruv Recon_TOF_ParaBeam2D.a \
	Read_Image2D.o \
	Rotate2D.o \
	ParaBeamP2D.o \
	ParaBeamQ2D.o \
	ConvP1D.o \
	ConvQ1D.o \
	DeConv1D_LR.o \
	DeConv1D_Wiener.o \
	TOF_ParaBeamP2D.o \
	TOF_ParaBeamQ2D.o \
	RampFilter1D.o \
	Diff_AC.o \
	CGM_AttenDI2D.o \
	Deconv_AC.o \
	FBP2D.o \
	TOFFBP2D.o \
	TOFDeConvBP2D.o \
	TOFEM2D.o
gcc -Wall -c 0_Recon_TOF_ParaBeam2D.c
gcc 0_Recon_TOF_ParaBeam2D.o Recon_TOF_ParaBeam2D.a ../../1_Lib/Image_Misc.a ../../1_Lib/Fourier.a -Wall -o Recon_TOF_ParaBeam2D.out -lm
rm *.o
rm Recon_TOF_ParaBeam2D.a
