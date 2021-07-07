#include <iostream>
#include <math.h>

//~ Determine if a the (i,j,k) point is inside of a (a,b,c) sized elpisoid centered in (Cx,Cy,Cz)
bool elipsoid(int i, int j, int k, double Cx, double Cy, double Cz, double a, double b, double c){//test if the point (i,j,k) is inside of the elipsoid of dimension (a,b,c)
  double rx=(i-Cx)/a;
  double ry=(j-Cy)/b;
  double rz=(k-Cz)/c;
  if (rx*rx+ry*ry+rz*rz<=1) return true;
  else return false;
  
}

int main(int argc, char **argv)
{
  if (argc==1)
  {
    printf("Please, start this program with:\n\t %s D > sphere.csv \n to create a sphere of diameter D",argv[0]);
    printf(", or with:\n\t%s Nx Ny Nz > sphere.csv \n to create an elipsoid\n",argv[0]);
    exit(1);
  }
  
  const int Nx=atoi(argv[1]);
  const int Ny=argc>3?atoi(argv[2]):Nx;
  const int Nz=argc>3?atoi(argv[3]):Nx;
  const double Cx=(Nx-1)/2.0;          //get the center x point of the Lattice
  const double Cy=(Ny-1)/2.0;          //get the center y point of the Lattice
  const double Cz=(Nz-1)/2.0;          //get the center z point of the Lattice
  double rx, ry, rz;
  double v[3], rr;
  int pt, i, j, k;
    FILE *output =stdout;
  fprintf(output,"x,y,z,nx,ny,nz,pt\n");//create the .csv headers
  for( int k= 0; k< Nz; k++)
  {
    for( int j= 0; j< Ny; j++)
    {
      for( int i= 0; i< Nx; i++)
        {   
        //evaluate the distance between the point and the center in all coordinates
        rx=(i-Cx);
        ry=(j-Cy);
        rz=(k-Cz);
        if     (elipsoid(i,j,k,Cx,Cy,Cz,Cx-1,Cy-1,Cz-1)) pt=1;      //Setting a bulk point
        else if(elipsoid(i,j,k,Cx,Cy,Cz,Cx+1,Cy+1,Cz+1)) pt=2;      //Setting a surface point
        else                            pt=0;                       //Setting a empty point
        //Evaluate the normals
        v[0]=2*rx/(Cx*Cx);
        v[1]=2*ry/(Cy*Cy);
        v[2]=2*rz/(Cz*Cz);
        rr=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); //evaluate the normalize factor
        if(rr==0)rr = 1;
        //the first line produces a cleaner .csv file and the second a complete one
        if(pt!=1)fprintf(output,"%d,%d,%d,%g,%g,%g,%d\n",i,j,k,v[0]/rr,v[1]/rr,v[2]/rr,pt);
//~         fprintf(output,"%d,%d,%d,%g,%g,%g,%d\n",i,j,k,v[0]/rr,v[1]/rr,v[2]/rr,pt);
        
      }
    }
  }
  
    return 0;
}

