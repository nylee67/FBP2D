/*
	CGM_AttenDI2D.c
*/
#include	"TOF_ParaBeam2D.h"

void	int2str(char *, int, int);
void	pgm16_write_scaled_float2d(char *, float **, int, int);

void	Expand(float **F, float **f, int Na, int Nu);
void	BackExpand(float **f, float **F, int Na, int Nu);

void	Projection(float **Fa, float **Fu, float **f, int Na, int Nu, float da, float h);					// Included in this file
void	Backprojection(float **f, float **Fa, float **Fu, int Na, int Nu, float da, float h);				// Included in this file

void	Q_transform(float **g, float **f, int Na, int Nu, float da, float h, float lambda);					// Included in this file

void	G_Finder(float **g, float **f, float **b, int Na, int Nu, float da, float h, float lambda);			// Included in this file
void	Linear_Combination(float **f, float **g, float **h, float alpha, float beta, int Na, int Nu); 		// Included in this file
float	Q_Inner_Product(float **f, float **g, int Na, int Nui, float da, float h, float lambda);			// Included in this file
float	Inner_Product(float **f, float **g, int Na, int Nu); 												// Included in this file

void	CGM_AttenDI2D(
	float		**F,
	float		**Fa,
	float		**Fu,
	int			Na,
	int			Nu,
	float 		da, 
	float 		h)
{
	//
	// 	Declaration of Variables	
	//
		int			i, j, k, n, m;
		float		**old_f, **new_f, **old_d, **new_d, **old_g, **new_g, **vec_b, ftemp, alpha, beta;

		int			IterNum = 1000;
		float		lambda = (float) 0.1;
		char		*temp_fn;

	//
	// 	Memory Allocation for CGM
	//
		old_f = (float **) calloc(Na, sizeof(float *));
		for(j=0; j<Na; j++)
			old_f[j] = (float *) calloc(Nu, sizeof(float));

		new_f = (float **) calloc(Na, sizeof(float *));
		for(j=0; j<Na; j++)
			new_f[j] = (float *) calloc(Nu, sizeof(float));

		old_d = (float **) calloc(Na, sizeof(float *));
		for(j=0; j<Na; j++)
			old_d[j] = (float *) calloc(Nu, sizeof(float));

		new_d = (float **) calloc(Na, sizeof(float *));
		for(j=0; j<Na; j++)
			new_d[j] = (float *) calloc(Nu, sizeof(float));

		old_g = (float **) calloc(Na, sizeof(float *));
		for(j=0; j<Na; j++)
			old_g[j] = (float *) calloc(Nu, sizeof(float));

		new_g = (float **) calloc(Na, sizeof(float *));
		for(j=0; j<Na; j++)
			new_g[j] = (float *) calloc(Nu, sizeof(float));

		vec_b = (float **) calloc(Na, sizeof(float *));
		for(j=0; j<Na; j++)
			vec_b[j] = (float *) calloc(Nu, sizeof(float));

	// Initilaization
		Backprojection(vec_b, Fa, Fu, Na, Nu, da, h); 						// b = P^t Y

		if(1){
			temp_fn = (char *) calloc(200, sizeof(char)); 	
			strcat(temp_fn, "../../3_Output/Test_"); 
			strcat(temp_fn, "_vec_b.pgm");
       		pgm16_write_scaled_float2d(temp_fn, vec_b, Na, Nu);
			free(temp_fn);
		}

		for(j=0; j<Na; j++)
			for(k=0; k<Nu; k++)
				old_f[j][k] = (float) 0.0;

		G_Finder(old_g, old_f, vec_b, Na, Nu, da, h, lambda);

		for(j=0; j<Na; j++)
			for(k=0; k<Nu; k++)
				old_d[j][k] = - old_g[j][k];

	// CGM - Interation
	for(i=0; i<IterNum; i++){

		ftemp = Q_Inner_Product(old_d, old_d, Na, Nu, da, h, lambda);
		if(ftemp != 0.0)
			alpha = - Inner_Product(old_g, old_d, Na, Nu)/ftemp;
		else{
			puts("CGM sub-iteration is terminated\n");
			break;			// Stop CGM loop
		}

		Linear_Combination(new_f, old_f, old_d, 1.0, alpha, Na, Nu);
		G_Finder(new_g, new_f, vec_b, Na, Nu, da, h, lambda);
		beta = Q_Inner_Product(new_g, old_d, Na, Nu, da, h, lambda)/ftemp;
 		Linear_Combination(new_d, new_g, old_d, -1.0, beta, Na, Nu);

		for(n=0; n<Na; n++)
			for(m=0; m<Nu; m++){
				old_f[n][m] = new_f[n][m];
				old_g[n][m] = new_g[n][m];
				old_d[n][m] = new_d[n][m];
			}

	}// end of i-loop

	for(n=0; n<Na; n++)
		for(m=0; m<Nu; m++)
			F[n][m] = new_f[n][m];

	for(i=0; i<Na; i++){
		free(old_f[i]); free(new_f[i]); 
		free(old_d[i]); free(new_d[i]); 
		free(old_g[i]); free(new_g[i]); 
		free(vec_b[i]); 
	}
	free(old_f); free(new_f); free(old_d); free(new_d); free(old_g); free(new_g); free(vec_b); 
}

void	Expand(float **F, float **f, int Na, int Nu)
{
	int		j,k;
	
	for(j=1; j<Na+1; j++){
		
		F[j][0] = (float) 0.0;

		for(k=1; k<Nu+1; k++)
			F[j][k] = f[j-1][k-1];

		F[j][Na+1] = (float) 0.0;
	}

	for(k=0; k<Nu+2; k++)
		F[0][k] = F[Na][Nu+2-k-1];

	for(k=0; k<Nu+2; k++)
		F[Na+1][k] = F[1][Nu+2-k-1];
}

void	BackExpand(float **f, float **F, int Na, int Nu)
{
	int		i, j;
	
	for(i=0; i<Na; i++)
		for(j=0; j<Nu; j++)
			f[i][j] = F[i+1][j+1];

	for(j=0; j<Nu; j++)
		f[Na-1][j] += F[0][Nu-j];

	for(j=0; j<Nu; j++)
		f[0][j] += F[Na+1][Nu-j];
}

void	Projection(float **Fa, float **Fu, float **f, int Na, int Nu, float da, float h)										
{
	int		i,j;
	float	 **tmp_F;

	tmp_F = (float **) calloc(Na+2, sizeof(float *));
	for(i=0; i<Na+2; i++)
		tmp_F[i] = (float *) calloc(Nu+2, sizeof(float));

	Expand(tmp_F, f, Na, Nu);

	for(i=0; i<Na; i++)
		for(j=0; j<Nu; j++)
			Fa[i][j] = 1/(2*da) * (tmp_F[i+2][j+1] - tmp_F[i][j+1]);

	for(i=0; i<Na; i++)
		for(j=0; j<Nu; j++)
			Fu[i][j] = 1/(2*h) * (tmp_F[i+1][j+2] - tmp_F[i+1][j]);
  
	for(i=0; i<Na+2; i++)
		free(tmp_F[i]);
	free(tmp_F);
}

void	Backprojection(float **f, float **Fa, float **Fu, int Na, int Nu, float da, float h)									
{
	int		i,j;
	float	 **tmp_F;

	tmp_F = (float **) calloc(Na+2, sizeof(float *));
	for(i=0; i<Na+2; i++)
		tmp_F[i] = (float *) calloc(Nu+2, sizeof(float));

	for(i=2; i<Na; i++)
		for(j=2; j<Nu; j++)
			tmp_F[i][j] = 1/(2*da) * (Fa[i-2][j-1] - Fa[i][j-1]) + 1/(2*h) * (Fu[i-1][j-2] - Fu[i-1][j]);

	for(j=1; j<Nu+1; j++){
		tmp_F[0][j] = -1/(2*da) * Fa[0][j-1];
		tmp_F[Na+1][j] = 1/(2*da) * Fa[Na-1][j-1];
	}

	for(j=2; j<Nu; j++){
		tmp_F[1][j] = -1/(2*da) * Fa[1][j-1] + 1/(2*h) * Fu[0][j-2] - 1/(2*h) * Fu[0][j];
		tmp_F[Na][j] = 1/(2*da) * Fa[Na-2][j-1] + 1/(2*h) * Fu[Na-1][j-2] - 1/(2*h) * Fu[Na-1][j];
	}

	for(i=1; i<Na+1; i++){
		tmp_F[i][0] = -1/(2*h) * Fu[i-1][0];
		tmp_F[i][Na+1] = 1/(2*h) * Fu[i-1][Nu-1];
	}

	for(i=2; i<Na; i++){
		tmp_F[i][1] = 1/(2*da) * Fa[i-2][0] - 1/(2*da) * Fa[i][0] - 1/(2*h) * Fu[i-1][1];
		tmp_F[i][Nu] = 1/(2*da) * Fa[i-2][Nu-1] - 1/(2*da) * Fa[i][Nu-1] + 1/(2*h) * Fu[i-1][Nu-2];
	}

	BackExpand(f, tmp_F, Na, Nu);

	for(i=0; i<Na+2; i++)
		free(tmp_F[i]);
	free(tmp_F);
}

void	Q_transform(float **g, float **f, int Na, int Nu, float da, float h, float lambda)
{
	int		i, j;
	float	**tmp_Fa, **tmp_Fu;

	tmp_Fa = (float **) calloc(Na, sizeof(float *));
	for(i=0; i<Na; i++)
		tmp_Fa[i] = (float *) calloc(Nu, sizeof(float));

	tmp_Fu = (float **) calloc(Na, sizeof(float *));
	for(i=0; i<Na; i++)
		tmp_Fu[i] = (float *) calloc(Nu, sizeof(float));

	Projection(tmp_Fa, tmp_Fu, f, Na, Nu, da, h);
	Backprojection(g, tmp_Fa, tmp_Fu, Na, Nu, da, h);
	
	for(i=0; i<Na; i++)
		for(j=0; j<Nu; j++)
			g[i][j] += lambda*f[i][j];

	for(i=0; i<Na; i++)
		free(tmp_Fa[i]);
	free(tmp_Fa);

	for(i=0; i<Na; i++)
		free(tmp_Fu[i]);
	free(tmp_Fu);
}			

void	G_Finder(float **g, float **f, float **b, int Na, int Nu, float da, float h, float lambda)
{
	int		j, k;

	Q_transform(g, f, Na, Nu, da, h, lambda);

	for(j=0; j<Na; j++)
		for(k=0; k<Nu; k++)
			g[j][k] = g[j][k] - b[j][k];
}

void	Linear_Combination(float **f, float **g, float **h, float alpha, float beta, int Na, int Nu)
{
	int		j, k;

	for(j=0; j<Na; j++)
		for(k=0; k<Nu; k++)
			f[j][k] = alpha*g[j][k] + beta*h[j][k];
}

float	Q_Inner_Product(float **f, float **g, int Na, int Nu, float da, float h, float lambda)
{
	int		j, k; 
	float	**tmp_f;
	float	sum;

	tmp_f = (float **) calloc(Na, sizeof(float *));
	for(j=0; j<Na; j++)
		tmp_f[j] = (float *) calloc(Nu, sizeof(float));

	Q_transform(tmp_f, f, Na, Nu, da, h, lambda);

	sum = (float) 0.0;
	for(j=0; j<Na; j++)
		for(k=0; k<Nu; k++)
			sum += g[j][k]*tmp_f[j][k];

	for(j=0; j<Na; j++)
		free(tmp_f[j]);
	free(tmp_f);

	return	(sum);
}

float	Inner_Product(float **f, float **g, int Na, int Nu)
{
	int		j, k;
	float	sum;

	sum = (float) 0.0;
	for(j=0; j<Na; j++)
		for(k=0; k<Nu; k++)
			sum += g[j][k]*f[j][k];

	return	(sum);
}
