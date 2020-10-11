/*
	Diff_AC.c
*/
#include	"TOF_ParaBeam2D.h"

void	CGM_AttenDI2D(
	float	**a,
	float	**a_a,
	float	**a_u,
	int		Na,
	int		Nu,
	float	da,
	float	h);

void	Diff_AC(
	float	**R_mu,
	float	**a,
	float 	***m,
	int 	Na,
	int 	N,
	double	dd,
	int 	N_TSF,
	double 	Sigma_TSF)
{
	int		ia, iu, it, Nu = N, Nt = N;
	float	***m_u, ***m_a, ***m_t, ***m_ut, ***Dm;
	float	**H_uu, **H_ua, **H_aa, **J_u, **J_a, **a_a, **a_u;
	float	h = dd, dt = dd, du = dt, t, u, d_theta = MY_PI/Na, temp, sum_mean, sum_each;	

	//
	//	Memory Allocation
	//
		m_u = (float ***) calloc(Na, sizeof(float **));
		for(ia=0; ia<Na; ia++){
			m_u[ia] = (float **) calloc(Nu, sizeof(float *));
			for(iu=0; iu<Nu; iu++)
				m_u[ia][iu] = (float *) calloc(Nt, sizeof(float));
		}

		m_a = (float ***) calloc(Na, sizeof(float **));
		for(ia=0; ia<Na; ia++){
			m_a[ia] = (float **) calloc(Nu, sizeof(float *));
			for(iu=0; iu<Nu; iu++)
				m_a[ia][iu] = (float *) calloc(Nt, sizeof(float));
		}

		m_t = (float ***) calloc(Na, sizeof(float **));
		for(ia=0; ia<Na; ia++){
			m_t[ia] = (float **) calloc(Nu, sizeof(float *));
			for(iu=0; iu<Nu; iu++)
				m_t[ia][iu] = (float *) calloc(Nt, sizeof(float));
		}

		m_ut = (float ***) calloc(Na, sizeof(float **));
		for(ia=0; ia<Na; ia++){
			m_ut[ia] = (float **) calloc(Nu, sizeof(float *));
			for(iu=0; iu<Nu; iu++)
				m_ut[ia][iu] = (float *) calloc(Nt, sizeof(float));
		}

		Dm = (float ***) calloc(Na, sizeof(float **));
		for(ia=0; ia<Na; ia++){
			Dm[ia] = (float **) calloc(Nu, sizeof(float *));
			for(iu=0; iu<Nu; iu++)
				Dm[ia][iu] = (float *) calloc(Nt, sizeof(float));
		}

		a_a = (float **) calloc(Na, sizeof(float *));
		for(ia=0; ia<Na; ia++)
			a_a[ia] = (float *) calloc(Nu, sizeof(float));

		a_u = (float **) calloc(Na, sizeof(float *));
		for(ia=0; ia<Na; ia++)
			a_u[ia] = (float *) calloc(Nu, sizeof(float));

		J_u = (float **) calloc(Na, sizeof(float *));
		for(ia=0; ia<Na; ia++)
			J_u[ia] = (float *) calloc(Nu, sizeof(float));

		J_a = (float **) calloc(Na, sizeof(float *));
		for(ia=0; ia<Na; ia++)
			J_a[ia] = (float *) calloc(Nu, sizeof(float));

		H_uu = (float **) calloc(Na, sizeof(float *));
		for(ia=0; ia<Na; ia++)
			H_uu[ia] = (float *) calloc(Nu, sizeof(float));

		H_ua = (float **) calloc(Na, sizeof(float *));
		for(ia=0; ia<Na; ia++)
			H_ua[ia] = (float *) calloc(Nu, sizeof(float));

		H_aa = (float **) calloc(Na, sizeof(float *));
		for(ia=0; ia<Na; ia++)
			H_aa[ia] = (float *) calloc(Nu, sizeof(float));

	//
	//	Computing Deriavtives
	//
		for(ia=0; ia<Na; ia++)
			for(it=0; it<Nt; it++){
				m_u[ia][0][it] = (float) 0.0;
				m_u[ia][Nu-1][it] = (float) 0.0;
				for(iu=1; iu<Nu-1; iu++)
					m_u[ia][iu][it] = (m[ia][iu+1][it] - m[ia][iu-1][it])/(2*h);
			}

		for(ia=1; ia<Na-1; ia++)
			for(iu=0; iu<Nu; iu++)
				for(it=0; it<Nt; it++)
					m_a[ia][iu][it] = (m[ia+1][iu][it] - m[ia-1][iu][it])/(2*d_theta);

		for(iu=0; iu<Nu; iu++)
			for(it=0; it<Nt; it++)
				m_a[0][iu][it] = (m[1][iu][it] - m[Na-1][Nu-iu-1][Nt-it-1])/(2*d_theta);

		for(iu=0; iu<Nu; iu++)
			for(it=0; it<Nt; it++)
				m_a[Na-1][iu][it] = (m[0][Nu-iu-1][Nt-it-1] - m[Na-2][iu][it])/(2*d_theta);

		for(ia=0; ia<Na; ia++)
			for(iu=0; iu<Nu; iu++){
				m_t[ia][iu][0] = (float) 0.0;
				m_t[ia][iu][Nt-1] = (float) 0.0;
				for(it=1; it<Nt-1; it++)
					m_t[ia][iu][it] = (m[ia][iu][it+1] - m[ia][iu][it-1])/(2*h);
			}

		for(ia=0; ia<Na; ia++)
			for(iu=0; iu<Nu; iu++){
				m_ut[ia][iu][0] = (float) 0.0;
				m_ut[ia][iu][Nt-1] = (float) 0.0;
				for(it=1; it<Nt-1; it++)
					m_ut[ia][iu][it] = (m_u[ia][iu][it+1] - m_u[ia][iu][it-1])/(2*h);
			}

		for(ia=0; ia<Na; ia++)
			for(iu=1; iu<Nu-1; iu++){
				u = (iu - Nu/2.0 + 0.5)*du;
				for(it=0; it<Nt; it++){
					t = (it - Nt/2.0 + 0.5)*dt;
					Dm[ia][iu][it] = t*m_u[ia][iu][it] + m_a[ia][iu][it] - u*m_t[ia][iu][it] + Sigma_TSF*Sigma_TSF*m_ut[ia][iu][it];
//printf("Dm[%4d][%4d][%4d] = %12.6f\n", ia, iu, it, Dm[ia][iu][it]);
				}
			}

	//
	//	Computing Integration
	//
		for(ia=0; ia<Na; ia++)
			for(iu=0; iu<Nu; iu++)
				for(it=0; it<Nt; it++){
					t = (it - Nt/2.0 + 0.5)*dt;

					J_u[ia][iu] += Dm[ia][iu][it] * (m[ia][iu][it]*t + Sigma_TSF*Sigma_TSF*m_t[ia][iu][it]) * dt; 
					J_a[ia][iu] += Dm[ia][iu][it] * m[ia][iu][it] * dt; 
					
					H_uu[ia][iu] += (m[ia][iu][it]*t + Sigma_TSF*Sigma_TSF * m_t[ia][iu][it]) * (m[ia][iu][it]*t + Sigma_TSF*Sigma_TSF * m_t[ia][iu][it]) * dt; 
					H_ua[ia][iu] += m[ia][iu][it] * (m[ia][iu][it]*t + Sigma_TSF*Sigma_TSF*m_t[ia][iu][it]) * dt; 
					H_aa[ia][iu] += m[ia][iu][it] * m[ia][iu][it] * dt; 
				}

	//
	//	Estimating Attenuation
	//
		for(ia=0; ia<Na; ia++)
			for(iu=0; iu<Nu; iu++){

				temp = H_uu[ia][iu]*H_aa[ia][iu] - H_ua[ia][iu]*H_ua[ia][iu];
				//printf("temp %4d %4d %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n", ia, iu, temp, J_u[ia][iu], J_a[ia][iu], H_aa[ia][iu], H_ua[ia][iu], H_uu[ia][iu]);
				if(temp != 0.0){ 
					a_a[ia][iu] = (-J_u[ia][iu]*H_aa[ia][iu] + J_a[ia][iu]*H_ua[ia][iu])/temp;
					a_u[ia][iu] = (-J_a[ia][iu]*H_uu[ia][iu] + J_u[ia][iu]*H_ua[ia][iu])/temp;
				}
				else{
					a_a[ia][iu] = (float) 0.0;
					a_u[ia][iu] = (float) 0.0;
				}
			}


//		for(ia=0; ia<Na; ia++)
//			for(iu=0; iu<Nu; iu++)
//				printf("Yhhh %4d %4d %20.12f %20.12f\n", ia, iu, a_a[ia][iu], a_u[ia][iu]);
	//
	//	CGM2D
	//		
if(1){
		for(ia=0; ia<Na; ia++)
			for(iu=1; iu<Nu-1; iu++)
				a_u[ia][iu] = (R_mu[ia][iu+1] - R_mu[ia][iu-1])/(2*h);

		for(ia=0; ia<Na; ia++){
			a_u[ia][0] = (float) 0.0;
			a_u[ia][Nu-1] = (float) 0.0;
		}

		for(ia=1; ia<Na-1; ia++)
			for(iu=0; iu<Nu; iu++)
				a_a[ia][iu] = (R_mu[ia+1][iu] - R_mu[ia-1][iu])/(2*d_theta);

		for(iu=0; iu<Nu; iu++)
			a_a[0][iu] = (R_mu[1][iu] - R_mu[Na-1][Nu-iu-1])/(2*d_theta);

		for(iu=0; iu<Nu; iu++)
			a_a[Na-1][iu] = (R_mu[0][Nu-iu-1] - R_mu[Na-2][iu])/(2*d_theta);
}

		CGM_AttenDI2D(a,a_a,a_u,Na,Nu,d_theta,h);

	//	
	//	Consistence Condition in a-direction
	//
		temp = (float) 0.0; 
		for(ia=0; ia<Na; ia++){
			for(iu=0; iu<5; iu++)
				temp += a[ia][iu];
			for(iu=Nu-5; iu<Nu; iu++)
				temp += a[ia][iu];
		}
		temp = temp/(Na*10);

		printf("min %20.12f\n", temp);

		for(ia=0; ia<Na; ia++)
			for(iu=0; iu<Nu; iu++)
				a[ia][iu] -= -5.5;

		sum_mean = (float) 0.0;
		for(ia=0; ia<Na; ia++)
			for(iu=0; iu<Nu; iu++)
				sum_mean += a[ia][iu];

		sum_mean = sum_mean/Na;

		for(ia=0; ia<Na; ia++){
			sum_each = 0.0;
			for(iu=0; iu<Nu; iu++)
				sum_each += a[ia][iu];
printf("taho %20.12f\n", sum_each);
			for(iu=0; iu<Nu; iu++)
				a[ia][iu] = a[ia][iu]/sum_each * sum_mean;
		}
}
