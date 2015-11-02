/* T		= temperature
 * L		= constant edge length of lattice
 * s[]		= lattice spin configuration (with periodic boudary conditions)
 * e_av 	= <E>
 * esq_av	= <E^2>
 * m_av 	= <M>
 * m_q_av 	= <M^2>
 * H_i		= Hamiltonian along each direction (i = a, b, c)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <cstring>

using namespace std;

#define N (L*L)
#define MCS 10000
#define Transient 5000
#define Jzz 1
#define Jpm 0.915618
#define Tmax 1.65
#define Tmin 1.40
#define Tstep 0.01
#define Pi 3.1415926535897932384626
#define sqr(x) ((x)*(x))
#define pre(x) ((x+L-1)%L)
#define nex(x) ((x+1)%L)
#define ran() ((double)(rand()%32768)/32768.0)

FILE *str_res, *chi_res;
typedef std::complex<double> cdouble;
const cdouble I (0, 1);

int L = 6;

const double frac = 1.0/(N*MCS);
cdouble e_av, esq_av, m_av, m_q_av;
cdouble E;

cdouble Jzp, Jpp;
cdouble H[3][3][3];
cdouble Tc[200][100];

struct spin{
	cdouble x, y, z;
};

struct lat{
	double x, y;
};

struct spin s[60][60];
cdouble S_corr_stripe_av, S_corr_chiral_av;
struct spin s_chiral_phase[3], s_stripe_phase[3];

struct lat r(int i, int j){
	struct lat dl;
	dl.x = (double)i + (double)j/2;
	dl.y = (double)j*sqrt(3)/2;
	return dl;
};

lat q_chiral[3], q_stripe[3];

cdouble dot(int index, struct spin s1, struct spin s2)
{
	cdouble temp = H[index][0][0]*s1.x*s2.x + H[index][0][1]*s1.x*s2.y + H[index][0][2]*s1.x*s2.z + H[index][1][0]*s1.y*s2.x + H[index][1][1]*s1.y*s2.y + H[index][1][2]*s1.y*s2.z + H[index][2][0]*s1.z*s2.x + H[index][2][1]*s1.z*s2.y + H[index][2][2]*s1.z*s2.z;
	return temp;
}

double dp(struct lat l1, struct lat l2)
{
	double temp = l1.x*l2.x + l1.y*l2.y;
	return temp;
}

cdouble dp(struct spin s1, struct spin s2)
{
	cdouble temp = s1.x*s2.x + s1.y*s2.y + s1.z*s2.z;
	return temp;
}

struct spin ad(struct spin s1, struct spin s2)
{
	struct spin ds;
	ds.x = s1.x+s2.x, ds.y = s1.y+s2.y, ds.z = s1.z+s2.z;
	return ds;
}

struct spin sb(struct spin s1, struct spin s2)
{
	struct spin ds;
	ds.x = s1.x-s2.x, ds.y = s1.y-s2.y, ds.z = s1.z-s2.z;
	return ds;
}

struct spin sp(struct spin s, cdouble l)
{
	struct spin ds;
	ds.x = s.x*l, ds.y = s.y*l, ds.z = s.z*l;
	return ds;
}

struct spin cj(struct spin s)
{
	struct spin ds;
	ds.x = conj(s.x), ds.y = conj(s.y), ds.z = conj(s.z);
	return ds;
}

void gen(struct spin * s, double T)
{
	double zeta_1, zeta_2, zeta_q;
	do
	{
		zeta_1 = 1-2*ran();
		zeta_2 = 1-2*ran();
		zeta_q = sqr(zeta_1)+sqr(zeta_2);
		//printf("%12.8f%12.8f%12.8f\n", zeta_1, zeta_2, zeta_q);

	} while (zeta_q > 1);// || (sqr(ss->x - s->x) + sqr(ss->y - s->y) > T));

	s->x = 2*zeta_1*sqrt(1-zeta_q);
	s->y = 2*zeta_2*sqrt(1-zeta_q);
	s->z = 1-2*zeta_q;
}

cdouble energy()
{
	cdouble E = 0;
	int i, j;
	//demo();
	struct spin nb[6], ds;
	for (i = 0; i < L; i++) for (j = 0; j < L; j++)
	{
		nb[0] = s[nex(i)][j], 		nb[1] = s[i][nex(j)];
		nb[2] = s[pre(i)][nex(j)], 	nb[3] = s[pre(i)][j];
		nb[4] = s[i][pre(j)], 		nb[5] = s[nex(i)][pre(j)];

		ds = s[i][j];
		E += dot(0,ds,nb[0]) + dot(1,ds,nb[2]) + dot(2,ds,nb[4]);
	}
	return E;
}

void initialize()
{
	int i, j, k;
	for (i = 0; i < L; i++)	for (j = 0; j < L; j++) gen(&s[i][j], 10);
	e_av = 0;
	m_av = 0;
	esq_av = 0;
	m_q_av = 0;
	S_corr_stripe_av = 0;
	S_corr_chiral_av = 0;
	E = energy();

	memset(s_stripe_phase, 0, sizeof(s_stripe_phase));
	memset(s_chiral_phase, 0, sizeof(s_chiral_phase));
	for (i = 0; i < L; i++) for (j = 0; j < L; j++) for (k = 0; k < 3; k++)
	{
		s_chiral_phase[k] = ad(s_chiral_phase[k], sp(s[i][j], exp(I * dp(r(i, j), q_chiral[k]))));
		s_stripe_phase[k] = ad(s_stripe_phase[k], sp(s[i][j], exp(I * dp(r(i, j), q_stripe[k]))));
	}
}

void sweep(double T)
{
	//printf("SWEEP()\n");
	int i, j, k, count = 0;
	struct spin ss, ds, nb[6];
	cdouble dE;

	for (k = 0; k < N; k++)
	{
		i = rand()%L, j = rand()%L;

		gen(&ss, T); 	ds = sb(ss, s[i][j]);

		nb[0] = s[nex(i)][j], 		nb[1] = s[i][nex(j)];
		nb[2] = s[pre(i)][nex(j)], 	nb[3] = s[pre(i)][j];
		nb[4] = s[i][pre(j)], 		nb[5] = s[nex(i)][pre(j)];

		dE = dot(0,ds,nb[0]) + dot(0,ds,nb[3]) + dot(1,ds,nb[2]) + dot(1,ds,nb[5]) + dot(2,ds,nb[4]) + dot(2,ds,nb[1]);

		if (real(dE) < 1e-8 || ran() < real(exp(-dE/T)))
		{
			//printf("flipped!\n");
			s[i][j] = ss;
			E = E+dE;
			for (int iter = 0; iter < 3; iter++)
			{
				s_chiral_phase[iter] = ad(s_chiral_phase[iter], sp(ds, exp(I * dp(r(i, j), q_chiral[iter]))));
				s_stripe_phase[iter] = ad(s_stripe_phase[iter], sp(ds, exp(I * dp(r(i, j), q_stripe[iter]))));
			}
		}
	}
	e_av += E;
	esq_av += sqr(E);

	for (k = 0; k < 3; k++)
	{
		S_corr_stripe_av += dp(cj(s_stripe_phase[k]), s_stripe_phase[k])-(double)N;
		S_corr_chiral_av += dp(cj(s_chiral_phase[k]), s_chiral_phase[k])-(double)N;
	}
	//printf("swept!\n");
}

double proc()
{
	//printf("PROC()\n");
	double T, C = 0, Ctmp, Tfin;
	int i;
	initialize();

	for (i = 0; i < Transient; i++) sweep(T);

	for (T = Tmax; T >= Tmin; T -= Tstep)
	{
		e_av = 0;
		esq_av = 0;
		S_corr_stripe_av = 0;
		S_corr_chiral_av = 0;
		for (i = 0; i < MCS; i++)
		{
			sweep(T);
			//if (i % 1000 == 0) printf("swept 1000\n");
		}
		e_av *= frac, esq_av *= frac, S_corr_chiral_av *= frac, S_corr_stripe_av *= frac;
		fprintf(str_res, "%12.8f%12.8f\n", T, real(S_corr_stripe_av)/N);
		fprintf(chi_res, "%12.8f%12.8f\n", T, real(S_corr_chiral_av)/N);
	}

	return Tfin;
}

int main()
{
	// Order parameter in the First Brillouin zone
	q_chiral[0].x = 4*Pi/3;  q_chiral[0].y = 0;
	q_chiral[1].x = -2*Pi/3; q_chiral[1].y = +2*Pi/sqrt(3);
	q_chiral[2].x = -2*Pi/3; q_chiral[2].y = -2*Pi/sqrt(3);

	q_stripe[0].x = Pi;		 q_stripe[0].y = Pi/sqrt(3);
	q_stripe[1].x = -Pi;	 q_stripe[1].y = Pi/sqrt(3);
	q_stripe[2].x = 0;		 q_stripe[2].y = -2*Pi/sqrt(3);

	str_res = fopen("str_res.dat", "w");
	chi_res = fopen("chi_res.dat", "w");
    int i, j;
    
    for (j = 0; j <= 10; j++)
        for (i = 0; i <= 20; i++)
        {
            Jzp = (double) j * 0.1;
            Jpp = (double) i * 0.1 - 1;
            fprintf(str_res, "\n%12.8f%12.8f\n\n", real(Jpp), real(Jzp));
            fprintf(chi_res, "\n%12.8f%12.8f\n\n", real(Jpp), real(Jzp));

            H[0][0][0] = 2.0*Jpm+2.0*Jpp; H[0][0][1] = 0; H[0][0][2] = 0;
            H[0][1][0] = 0; H[0][1][1] = 2.0*Jpm-2.0*Jpp; H[0][1][2] = Jzp;
            H[0][2][0] = 0; H[0][2][1] = Jzp; H[0][2][2] = Jzz;

            H[1][0][0] = 2.0*Jpm-Jpp; H[1][0][1] = -sqrt(3)*Jpp; H[1][0][2] = -sqrt(3)/2.0*Jzp;
            H[1][1][0] = -sqrt(3)*Jpp; H[1][1][1] = 2.0*Jpm+Jpp; H[1][1][2] = -Jzp/2.0;
            H[1][2][0] = -sqrt(3)/2.0*Jzp; H[1][2][1] = -Jzp/2.0; H[1][2][2] = Jzz;

            H[2][0][0] = 2.0*Jpm-Jpp; H[2][0][1] = +sqrt(3)*Jpp; H[2][0][2] = +sqrt(3)/2.0*Jzp;
            H[2][1][0] = +sqrt(3)*Jpp; H[2][1][1] = 2.0*Jpm+Jpp; H[2][1][2] = -Jzp/2.0;
            H[2][2][0] = +sqrt(3)/2.0*Jzp; H[2][2][1] = -Jzp/2.0; H[2][2][2] = Jzz;

            Tc[i][j] = proc();
            //printf("{%.8f, %.8f, %.8f}, ", Jpp, Jzp, Tc[i][j]);

            //printf("\n\n");
        }

	return 0;
}
