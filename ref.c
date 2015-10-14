#include <math.h>
#include <stdlib.h>
//random number generator from "Numerical Recipes in C" (Ran1.c)
#include <"random.h">
//file to output data into
//structure for a 2d lattice with coordinates x and y
struct lat_type
{
  int x;
  int y;
};
const int size=2;
//lattice size
const int lsize=size-1;
//array size for lattice
const int n=size*size;
//number of spin points on lattice      
float T=5.0;
//starting point for temperature
const float minT=0.5;
//minimum temperature
float change=0.1;
//size of steps for temperature loop
int lat[size+1][size+1];
//2d lattice for spins
long unsigned int mcs=10000;
//number of Monte Carlo steps
int transient=1000;
//number of transient steps 
double norm=(1.0/float(mcs*n));
//normalization for averaging
long int seed=436675;
   
//seed for random number generator      
//function for random initialization of lattice
void initialize(int lat[size+1][size+1])
{
    for (int y=size; y>=1; y--)
    {
for(int x=1;x<=size;x++)
{
    if(ran1(&seed)>=0.5)
lat[x][y]=1;
    else
lat[x][y]=-1;
}
    }
}
//output of lattice configuration to the screen
void output(int lat[size+1][size+1])
{
    for (int y=size; y>=1; y--)
    {
      for(int x=1;x<=size;x++)
{
  if(lat[x][y]<0)
    puts(" - ");
  else
    puts(" + ");
}
      putchar('\n');
    }
}
//function for choosing random position on lattice
void choose_random_pos_lat(lat_type &pos)
{
  pos.x=(int)ceil(ran1(&seed)*(size));
  pos.y=(int)ceil(ran1(&seed)*(size));
  if(pos.x>size||pos.y>size)
    {
      puts("error in array size allocation for random position on lattice!");
      exit;
    }
}
//function for calculating energy at a particular position on lattice
int energy_pos(lat_type &pos)
{
    //periodic boundary conditions
    int up,down,left,right,e;
    if(pos.y==size)
up=1;
    else
up=pos.y+1;
    if(pos.y==1)
down=size;
    else
down=pos.y-1;
    if(pos.x==1)
left=size;
    else
left=pos.x-1;
    if(pos.x==size)
right=1;
    else
right=pos.x+1;
    //energy for specific position
    e = -1*lat[pos.x][pos.y]
*(lat[left][pos.y]+lat[right][pos.y]+lat[pos.x][up]+lat[pos.x][down]);
    return e;
}
//function for testing the validity of flipping a spin at a selected position
bool test_flip(lat_type pos, int &de)
{
    de = -2*energy_pos(pos);
//change in energy for specific spin
    if(de<0)
return true;
//flip due to lower energy
else if(ran1(&seed)<exp(-de/T))
return true;
//flip due to heat bath
    else 
return false;
//no flip
}
//flip spin at given position 
void flip(lat_type pos)
{
  lat[pos.x][pos.y]=-lat[pos.x][pos.y];
}

//function for disregarding transient results
void transient_results()
{
    lat_type pos;
    int de=0;
    for(int a=1;a<=transient;a++)
    {
for(int b=1;b<=n;b++)
{
    choose_random_pos_lat(pos);
    if(test_flip(pos,de))
    {
flip(pos);
    }
}
    }
}
//function for calculating total magnetization of lattice
int total_magnetization()
{
  int m=0;
  for(int y=size;y>=1;y--)
    {
      for(int x=1;x<=size;x++)
{
  m+=lat[x][y];  
}
    }
  return m;
}  
//function for calculating total energy of lattice
int total_energy()
{
  lat_type pos;
  int e=0;
  for(int y=size;y>=1;y--)
    {
      pos.y=y;
      for(int x=1;x<=size;x++)
{
  pos.x=x;
  e+=energy_pos(pos);
}
    }
  return e;
}
//main program
void main()
{
    //declaring variables to be used in calculating the observables
    double E=0,Esq=0,Esq_avg=0,E_avg=0,etot=0,etotsq=0;
    double M=0,Msq=0,Msq_avg=0,M_avg=0,mtot=0,mtotsq=0;
    double Mabs=0,Mabs_avg=0,Mq_avg=0,mabstot=0,mqtot=0;
    int de=0;
    lat_type pos;
    
    //initialize lattice to random configuration
    initialize(lat);
    
    //Temperature loop
    for(;T>=minT;T=T-change)
    {
//transient function 
transient_results();
//observables adopt equilibrated lattice configurations values 
M=total_magnetization();
Mabs=abs(total_magnetization());
E=total_energy();
//initialize summation variables at each temperature step
etot=0;
etotsq=0;
mtot=0;
mtotsq=0;
mabstot=0;
mqtot=0;
        
        //Monte Carlo loop
       for(int a=1;a<=mcs;a++)
{
            //Metropolis loop
    for(int b=1;b<=n;b++)
    {
choose_random_pos_lat(pos);
if(test_flip(pos,de))
{
    flip(pos);
                    
                    //adjust observables
    E+=2*de;
    M+=2*lat[pos.x][pos.y];
    Mabs+=abs(lat[pos.x][pos.y]);
}
    }
           
    //keep summation of observables
    etot+=E/2.0;//so as not to count the energy for each spin twice
    etotsq+=E/2.0*E/2.0;
    mtot+=M;
    mtotsq+=M*M;
    mqtot+=M*M*M*M;
    mabstot+=(sqrt(M*M));
}
     
//average observables
E_avg=etot*norm;
Esq_avg=etotsq*norm;
M_avg=mtot*norm;
Msq_avg=mtotsq*norm;
Mabs_avg=mabstot*norm; 
Mq_avg=mqtot*norm;
        //output data to file 

}
