/****************************************************
 * calculate radial distribution function g(r)
 * of ensemble configurations in "movie.xyz"
 * output results in "gr.dat"
 * Kai Zhang, Duke University, 2011
 *****************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXConfig 10000
#define MAXEqui 10000
#define MAXNumberOfParticles 10000
#define MAXNumberOfBins 10000

#define SQR(x) ((x)*(x))
#define CUB(x) ((x)*(x)*(x))

int main(void)
{
  int NumberOfConfig,NumberOfEqui,NumberOfParticles;
  int SampleNumber;
  int snapshot,line;
  double L;
  double rho;
  double dr;
  double g[MAXNumberOfBins];
  int NumberOfBins;
  int i,j;
  double dx,dy,dz;
  double rij;
  char AtomID;
  FILE *fpxyz;
  FILE *fpoutput;
  double x[MAXNumberOfParticles],y[MAXNumberOfParticles],z[MAXNumberOfParticles];
  
  printf("**************************** 3D g(r) *************************************\n");
  printf("Number of Configurations in the .xyz file ?\n");
  scanf("%d",&NumberOfConfig);
  printf("Number of Equilibrations ?\n");
  scanf("%d",&NumberOfEqui);
  printf("Number of particles in the system ?\n");
  scanf("%d",&NumberOfParticles);
  printf("Simulation box size ?\n");
  scanf("%lf",&L);
  printf("Number of bins of the histogram ?\n");
  scanf("%d",&NumberOfBins);
  printf("\n");

  rho = NumberOfParticles/CUB(L);
  
  dr = 0.5*L/NumberOfBins;
  SampleNumber = 0;
  for(i=0;i<NumberOfBins;i++)
   g[i] = 0;

  SampleNumber = 0;
  fpxyz = fopen("movie.xyz","r"); 
  for(snapshot=0;snapshot<NumberOfConfig;snapshot++)
  {
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
    for(i=0;i<NumberOfParticles-1;i++)
     for(j=i+1;j<NumberOfParticles;j++)
     {
       dx = x[i] - x[j];
       dx = dx - L*round(dx/L);

       dy = y[i] - y[j];
       dy = dy - L*round(dy/L);

       dz = z[i] - z[j];
       dz = dz - L*round(dz/L);
        
       rij = sqrt(SQR(dx)+SQR(dy)+SQR(dz));

       if(rij < 0.5*L)
         g[(int)(rij/dr)] += 2.0;
     }
    }
  }
  fclose(fpxyz);


  fpoutput = fopen("gr.dat","w");
  for(i=0;i<NumberOfBins;i++)
  {
   g[i] /= SampleNumber;
   g[i] /= 4.0*M_PI/3.0*(CUB(i+1)-CUB(i))*CUB(dr)*rho;
   fprintf(fpoutput,"r = %lf\tg(r) = %lf\n",(i+0.5)*dr,g[i]/NumberOfParticles);
  }
  fclose(fpoutput);
 printf("density %lf\n",rho);  
 printf("sampling number %d\n",SampleNumber);  

 printf("*******************************************************\n");

 return 0;
}
