/***************************************************************
 * calculate static structure factor S(q)
 * from Fourier transform of radial distribution function g(r)
 * output results in "sq.dat"
 * and directly from coordinate file "movie.xyz"
 * Kai Zhang, Columbia University, 2015
 ****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXConfig 10000
#define MAXEqui 10000
#define MAXNumberOfParticles 10000
#define MAXNumberOfBins 10000
#define MAX_NumberOfWavevectors 100000 // for S(k)

#define SQR(x) ((x)*(x))
#define CUB(x) ((x)*(x)*(x))

int main(void)
{
  int NumberOfConfig,NumberOfEqui,NumberOfParticles;
  int SampleNumber;
  int snapshot,line;
  double L;
  double rho;
  double dr,dq;
  double g[MAXNumberOfBins];
  double ri[MAXNumberOfBins];
  double Sq[MAXNumberOfBins];
  double qi[MAXNumberOfBins];
  int NumberOfrBins;
  int NumberOfqBins;
  double qmax;
  int i,j;
  double dx,dy,dz;
  double rij;
  char AtomID;
  FILE *fpxyz;
  FILE *fpoutput;
  FILE *fpgr;
  double x[MAXNumberOfParticles],y[MAXNumberOfParticles],z[MAXNumberOfParticles];
  
  printf("**************************** 3D S(q) *************************************\n");
  printf("Number of particles in the system ?\n");
  scanf("%d",&NumberOfParticles);
  printf("Simulation box size ?\n");
  scanf("%lf",&L);
  printf("Number of r bins?\n");
  scanf("%d",&NumberOfrBins);
  printf("Number of q bins?\n");
  scanf("%d",&NumberOfqBins);
  printf("qmax ?\n");
  scanf("%lf",&qmax);
  printf("\n");

  rho = NumberOfParticles/CUB(L);
  printf("density %lf\n",rho);  
  dr = 0.5*L/NumberOfrBins;
  dq = (qmax-0.0)/NumberOfqBins;
  fpgr = fopen("gr.dat","r"); 
  for(line=0;line<NumberOfrBins;line++)
  {
   fscanf(fpgr,"r = %lf g(r) = %lf\n",&ri[line],&g[line]);
  }

  for(i=0;i<NumberOfqBins;i++)
  {
   qi[i] = (i+0.5)*dq;
   Sq[i] = 0.0;
   for(j=0;j<NumberOfrBins;j++)
   {
    Sq[i] += ri[j]*(g[j]-1.0)*sin(qi[i]*ri[j])/qi[i];
   }
   Sq[i] *= dr*4.0*M_PI*rho;
   Sq[i] += 1.0;
  }
  
  fpoutput = fopen("sqfourier.dat","w");
  for(i=0;i<NumberOfqBins;i++)
  {
   fprintf(fpoutput,"q = %lf\tS(q) = %lf\n",qi[i],Sq[i]);
  }
  fclose(fpoutput);
 
  /****************** calculate S(k) directly ************************/ 

  printf("Number of Configurations in the .xyz file ?\n");
  scanf("%d",&NumberOfConfig);
  printf("Number of Equilibrations ?\n");
  scanf("%d",&NumberOfEqui);
 
  
  /****************** make wave vectors ************************/ 
  
  int nkbound;
  int dkmultiplier;
  dkmultiplier = 1;
  printf("Number of wavevectors in each dimension ?\n");
  scanf("%d",&nkbound);
  printf("increment dk multiplier ?\n");
  scanf("%d",&dkmultiplier);

  double dkx,dky,dkz;
  double Lx,Ly,Lz;
  Lx=L;Ly=L;Lz=L;
  dkx = 2.0*M_PI/Lx*dkmultiplier;dky = 2.0*M_PI/Ly*dkmultiplier;dkz = 2.0*M_PI/Lz*dkmultiplier;

  int nk,nx,ny,nz;
  typedef struct
  {
        double x;
        double y;
        double z;
  } VECTOR;
  VECTOR Wavevector[MAX_NumberOfWavevectors];
  int NumberOfWavevectors;

  nk = 0;
  for(nx=0;nx<=nkbound;nx++)
  for(ny=0;ny<=nkbound;ny++)
  for(nz=0;nz<=nkbound;nz++)
  {
   Wavevector[nk].x = nx*dkx;
   Wavevector[nk].y = ny*dky;
   Wavevector[nk].z = nz*dkz;
   nk++;
  }
  NumberOfWavevectors = nk;

  printf("number of k's = %d\n",NumberOfWavevectors);
  /*******************************************/ 

  /******************** bubble sort k's ***********************/ 
 int itop,i1;
 double rxi,ryi,rzi;
 int change;
 itop = NumberOfWavevectors-1;
 do
 {
 change = false;
 for(i=0;i<itop;i++)
 {
  i1 = i+1;
  if(Wavevector[i].x*Wavevector[i].x + Wavevector[i].y*Wavevector[i].y + Wavevector[i].z*Wavevector[i].z > Wavevector[i1].x*Wavevector[i1].x+Wavevector[i1].y*Wavevector[i1].y+Wavevector[i1].z*Wavevector[i1].z)
  {
   rxi =  Wavevector[i].x;
   ryi =  Wavevector[i].y;
   rzi =  Wavevector[i].z;

   Wavevector[i].x = Wavevector[i1].x;
   Wavevector[i].y = Wavevector[i1].y;
   Wavevector[i].z = Wavevector[i1].z;

   Wavevector[i1].x = rxi;
   Wavevector[i1].y = ryi;
   Wavevector[i1].z = rzi;
   change = true;
  }
 }//end loop i
 itop--;
 }
 while(change==true && itop > 0);  

  double kmax;
  int nkmax;
  nkmax = NumberOfWavevectors-1;
  kmax = sqrt(Wavevector[nkmax].x*Wavevector[nkmax].x+Wavevector[nkmax].y*Wavevector[nkmax].y+Wavevector[nkmax].z*Wavevector[nkmax].z);
  printf("maximum k = %lf\n",kmax);

  /*******************************************/ 

  double Sk[MAX_NumberOfWavevectors],Sk_cos[MAX_NumberOfWavevectors],Sk_sin[MAX_NumberOfWavevectors];
  VECTOR k_scale;

  SampleNumber = 0;
  for(nk=0;nk<NumberOfWavevectors;nk++)
   Sk[nk] = 0.0;

  SampleNumber = 0;
  fpxyz = fopen("movie.xyz","r"); 
  for(snapshot=0;snapshot<NumberOfConfig;snapshot++)
  {
   printf("snapshot %d\n",snapshot);
   //for(line=0;line<=NumberOfParticles;line++) // N+1 lines
   for(line=0;line<=NumberOfParticles+1;line++) // N+2 lines
   {
     if(line == 0)
     fscanf(fpxyz,"%*d\n");
     else if(line == 1)
     fscanf(fpxyz,"%*s %*s %*d\n");
     else
     //fscanf(fpxyz,"%c %lf %lf %lf\n",&AtomID,&x[line-1],&y[line-1],&z[line-1]);
     //fscanf(fpxyz,"%c %lf %lf %lf\n",&AtomID,&x[line-2],&y[line-2],&z[line-2]);
     fscanf(fpxyz,"%*s %lf %lf %lf\n",&x[line-2],&y[line-2],&z[line-2]);
   }
   if(snapshot >= NumberOfEqui)
   {
    SampleNumber++;
    for(nk=0;nk<NumberOfWavevectors;nk++)
    {
     k_scale.x = Wavevector[nk].x*L/L;k_scale.y = Wavevector[nk].y*L/L; k_scale.z = Wavevector[nk].z*L/L;
     Sk_cos[nk] = 0.0;Sk_sin[nk] = 0.0;
     for(i=0;i<NumberOfParticles;i++) 
     {
      Sk_cos[nk] += cos(k_scale.x*x[i]+k_scale.y*y[i]+k_scale.z*z[i]);
      Sk_sin[nk] += sin(k_scale.x*x[i]+k_scale.y*y[i]+k_scale.z*z[i]);
     }// end loop over particle i
     Sk_cos[nk] = Sk_cos[nk]*Sk_cos[nk];
     Sk_sin[nk] = Sk_sin[nk]*Sk_sin[nk];
     Sk[nk] += Sk_cos[nk]+Sk_sin[nk];
    } // end loop over nk
   } // end if after equilibration
  }//end loop over snapshot
  fclose(fpxyz);

  printf("sampling number %d\n",SampleNumber);  
  for(nk=0;nk<NumberOfWavevectors;nk++)
    Sk[nk] /= SampleNumber*(NumberOfParticles);

  int kcount;
  double knorm,skmean;
  double kminus,kplus;
  fpoutput = fopen("sk.dat","w");
  kcount=0; knorm=0.0; skmean=0.0;
  for(nk=0;nk<NumberOfWavevectors;nk++)
  {
  kcount++;
  knorm += sqrt(Wavevector[nk].x*Wavevector[nk].x+Wavevector[nk].y*Wavevector[nk].y+Wavevector[nk].z*Wavevector[nk].z);
  skmean += Sk[nk];
  kplus = (Wavevector[nk+1].x*Wavevector[nk+1].x+Wavevector[nk+1].y*Wavevector[nk+1].y+Wavevector[nk+1].z*Wavevector[nk+1].z);
  kminus = (Wavevector[nk].x*Wavevector[nk].x+Wavevector[nk].y*Wavevector[nk].y+Wavevector[nk].z*Wavevector[nk].z);
  if( fabs(kplus-kminus) > 0.001 || nk+1 == NumberOfWavevectors )
  {
   knorm /= kcount;  skmean /= kcount;
   fprintf(fpoutput,"k = %lf\tS(k) = %lf\n",knorm,skmean);
   kcount=0; knorm=0.0; skmean=0.0;
  }
 }
  fclose(fpoutput);




 printf("*******************************************************\n");

 return 0;
}

