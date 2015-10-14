/* T		= temperature
 * L		= constant edge length of lattice
 * s[]		= lattice spin configuration (with periodic boudary conditions)
 * e_av 	= <E>
 * e_q_av	= <E^2>
 * m_av 	= <M>
 * m_q_av 	= <M^2>
 * H_i		= Hamiltonian along each direction (i = a, b, c)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define L 6
#define N (L*L)
#define MCS 10000
#define Transient 1000
#define Jzz 1
#define Jpm 0.915618
#define Jzp 0
#define Jpp 0

#define sqr(x) ((x)*(x))
#define pre(x) ((x+L-1)%L)
#define nex(x) ((x+1)%L)
#define ran() ((double)rand()/32767.0)

const double H[3][3][3] = {{{0,0,0},{0,0,0},{0,0,1}},{{0,0,0},{0,0,0},{0,0,1}},{{0,0,0},{0,0,0},{0,0,1}}};
const double frac = 1.0/(N*MCS);
double e_av, e_q_av, m_av, m_q_av;
double E;

struct spin{
	double x, y, z;
};

void display(struct spin s)
{
	printf("%12.8f\n", s.z);
}

struct spin s[L][L];

struct spin diff(struct spin s1, struct spin s2)
{
	struct spin ds;
	ds.x = s1.x-s2.x, ds.y = s1.y-s2.y, ds.z = s1.z-s2.z;
	return ds;
}

double dot(int index, struct spin s1, struct spin s2)
{
	double temp = H[index][0][0]*s1.x*s2.x + H[index][0][1]*s1.x*s2.y + H[index][0][2]*s1.x*s2.z + H[index][1][0]*s1.y*s2.x + H[index][1][1]*s1.y*s2.y + H[index][1][2]*s1.y*s2.z + H[index][2][0]*s1.z*s2.x + H[index][2][1]*s1.z*s2.y + H[index][2][2]*s1.z*s2.z;
	return temp;
}

void gen(struct spin * s)
{
	s->z = -s->z;
}

void demo()
{
	int i, j, k;
	for (i = 0; i < L; i++) for (j = 0; j < L; j++)
	{
		display(s[i][j]);
	}
	/*
	putchar('\n');
	for (k = 0; k < 3; k++)
	{
		for (i = 0; i < 3; i++)
		{
			for (j = 0; j < 3; j++) printf("%12.8f", H[k][i][j]);
			putchar('\n');
		}
		putchar('\n');
	}
	 */
}

double energy()
{
	double E = 0;
	int inow, jnow;
	//demo();
	struct spin nb[6], ds;
	for (inow = 0; inow < L; inow++) for (jnow = 0; jnow < L; jnow++)
	{
		nb[0] = s[nex(inow)][jnow], 		nb[1] = s[inow][nex(jnow)];
		nb[2] = s[pre(inow)][jnow], 		nb[3] = s[inow][pre(jnow)];

		ds = s[inow][jnow];
		E += dot(0,ds,nb[0]) + dot(1,ds,nb[1]);// + dot(2,ds,nb[4]); //+ dot(0,ds,nb[3]) + dot(1,ds,nb[5]) + dot(2,ds,nb[1]);
		//printf("(%d %d): (%d %d) (%d %d) (%d %d)\n", inow, jnow, nex(inow), jnow, pre(inow), nex(jnow), inow, pre(jnow));
	}
	return E;
}

void initialize()
{
	int i, j, inow, jnow;
	for (i = 0; i < L; i++)	for (j = 0; j < L; j++)
	{
		s[i][j].x = 0, s[i][j].y = 0, s[i][j].z = 1-2*(rand()%2);
	}
	e_av = 0;
	m_av = 0;
	e_q_av = 0;
	m_q_av = 0;
	E = energy();
}

void sweep(double T)
{
	int i, j, inow, jnow, count = 0;
	struct spin snow, ds, nb[6];
	double dE;

	//printf("%20.8f\n", E);

	for (i = 0; i < L; i++) for (j = 0; j < L; j++)
	{
		inow = (int)(ran()*L), jnow = (int)(ran()*L), snow = s[inow][jnow];
		gen(&snow); 	ds = diff(snow, s[inow][jnow]);

		nb[0] = s[nex(inow)][jnow], 		nb[1] = s[inow][nex(jnow)];
		nb[2] = s[pre(inow)][jnow], 		nb[3] = s[inow][pre(jnow)];

		dE = dot(0,ds,nb[0]) + dot(0,ds,nb[3]) + dot(1,ds,nb[2]) + dot(2,ds,nb[1]);

		//display(ds);
		//printf("%20.8lf\n\n", dE);


		//printf("%d %d    ", inow, jnow);
		//display(s[inow][jnow]); display(snow); puts("\n");

		if (dE < 0 || dE > 0 && ran() < exp(-dE/T) )
		{
			//printf("%20.8f", exp(-dE/T));
			//display(s[inow][jnow]);
			//display(snow);
			s[inow][jnow] = snow;
			E = E+dE;
			//printf("%60.8f\n", dE);
		}
		e_av += E*frac;
		e_q_av += sqr(E)*frac;
	}
	//printf("%20.8llf\n", E);
	//printf("%20.8llf\n", energy());
}

int main()
{
	freopen("Heat_Capacity.dat", "w", stdout);
	//srand(9);
	double T;
	int i, j;
	initialize();
	//puts("initialized!");
	for (T = 3; T > 2; T -= 0.01)
	{
		//T = 1;
		for (i = 0; i < Transient; i++) sweep(T);

		e_av = 0;
		e_q_av = 0;
		for (i = 0; i < MCS; i++) sweep(T);
		/*
			calculate new E, E^2, M, M^2 and such
		*/
		//printf("%40.8f%40.8f\n%40.8f%40.8f\n%40.8f\n\n", T, E, e_av, e_q_av, e_q_av-sqr(e_av));
		//printf("%40.8f\n", energy());
		printf("{%.8f,%.8f}, ", T, (e_q_av-sqr(e_av))/N);
	}

	return 0;
}
