#include <iostream>
#include <math.h>

int main(int argc, char **argv)
{
  if (argv<1)
  {
    printf("Please, start this program with %s Nx",argv[0]);
    printf("Please, start this program with %s Nx",argv[0]);
    printf("Please, start this program with %s Nx",argv[0]);
    exit(1);
  }
  
  const int Nx=atoi(argv[1]);
  const int Ny=atoi(argv[1]);
  const int Nz=atoi(argv[1]);
  const double HNx=Nx/2.0;
  const double HNy=Nx/2.0;
  const double HNz=Nx/2.0;
  double delta_x, delta_y, delta_z;
  double v[3], rr;
  int pt, i, j, k;
	FILE *output =fopen("sphere.csv","w");
  fprintf(output,"x,y,z,nx,ny,nz,pt\n");
  for( int k= 0; k< Nz; k++)
  {
    for( int j= 0; j< Ny; j++)
    {
      for( int i= 0; i< Nx; i++)
	    {	          
        delta_x=(i-HNx);
        delta_y=(j-HNy);
        delta_z=(k-HNz);
        rr=sqrt( delta_x*delta_x+delta_y*delta_y+delta_z*delta_z);
        if                  (rr < 49.6*Nx) pt=1;//Setting bulk points
        else if(rr >=49.6*Nx && rr <=50.95*Nx) pt=2;//Setting surface points
        else                            pt=0;//Setting outside points

        v[0]=delta_x/rr;
        v[1]=delta_y/rr;
        v[2]=delta_z/rr;
        if(rr==0){v[0]=v[1]=v[2]=0;}
        
        fprintf(output,"%d,%d,%d,%g,%g,%g,%d\n",i,j,k,v[0],v[1],v[2],pt);
        
      }
    }
  }
  
	return 0;
}

