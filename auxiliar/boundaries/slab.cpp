#include <iostream>
#include <math.h>

int main(int argc, char **argv)
{
  if (argc==1)
  {
    printf("Please, start this program with:\n\t %s Ni  > slab.csv \n to create a (Ni,Ni,Ni) slab",argv[0]);
    printf(", or with:\n\t%s Nx Ny Nz  > slab.csv \n to create a (Nx,Ny,Nz) slab\n",argv[0]);
    exit(1);
  }
  
  const int Nx=atoi(argv[1]);
  const int Ny=argc>3?atoi(argv[2]):Nx;
  const int Nz=argc>3?atoi(argv[3]):Nx;
  int pt, i, j, k;
	FILE *output =stdout;
  
//~   Print the header
  fprintf(output,"x,y,z,nx,ny,nz,pt\n");
  
  for( int k= 0; k< Nz; k++)
  {
    for( int j= 0; j< Ny; j++)
    {
      for( int i= 0; i< Nx; i++)
	    {
        if(k==0)
        {
//~           Bottom Surface
          pt=2;
          v[0]= 0.0;
          v[1]= 0.0;
          v[2]=-1.0;
          fprintf(output,"%d,%d,%d,%g,%g,%g,%d\n",i,j,k,v[0],v[1],v[2],pt);
        }else if(k==Nz-1)
        {
//~           Top Surface
          pt=3;
          v[0]= 0.0;
          v[1]= 0.0;
          v[2]= 1.0;
          fprintf(output,"%d,%d,%d,%g,%g,%g,%d\n",i,j,k,v[0],v[1],v[2],pt);
        }
//~         uncomment the line bellow to create a complete .csv file
//~         else { fprintf(output,"%d,%d,%d,%g,%g,%g,%d\n",i,j,k,0,0,0,1); }
      }
    }
  }
  
  
  
	return 0;
}

