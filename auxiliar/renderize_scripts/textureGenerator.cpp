#include "textureGenerator.h"
#include <iostream>
#include <string.h>
#include <math.h>

void info(char **argv);
void readfile(FILE *input, double *n, int *pt, int *N);
void molp(double* lightInt,double* n,int* pt,int* N);
void printData(double *lightInt, int *N, char fname[]);
void gera_plot (char arquivo[]);

int main(int argc, char **argv){
  if (argc<5 || strcasecmp(argv[1],"-h")==0) info(argv);
  FILE *inputFile = fopen(argv[1],"r");
  if( inputFile == 0) {perror(argv[1]); exit(2);}
  int N[4]={atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),0};
  N[3]=N[0]*N[1]*N[2];
  double *n=(double*) calloc( 3*N[3], sizeof(double));
  double *lightInt=(double*) calloc( N[0]*N[1], sizeof(double));
  int   *pt=(int*)  calloc( N[3], sizeof(int));
  readfile(inputFile, n, pt, N);
  molp(lightInt, n, pt, N);
  printData(lightInt,N,argv[1]);
  
  if(system("which gnuplot>/dev/null")){
    fprintf(stderr,"Gnuplot not found in the system, exiting without renderize the texture!!\n");
    exit(3);
  }
  gera_plot(argv[1]);
  char command[1000];
  sprintf(command,"gnuplot %s.plot && rm %s.plot",argv[1],argv[1]);
  system(command);
  
  
  
  
  return 0;
}

void info(char **argv){

  printf("This code will generate a texture in the .png format using the name pattern given by a .csv input file.\n");
  printf("Execute as:\n %s Input_File Nx Ny Nz\n", argv[0]);
  exit(1);
  
}

void readfile(FILE *input, double *n, int *pt, int *N){
  int x, y, z, ii, tempt;
  double nx, ny, nz;
  
  fscanf(input,"x,y,z,nx,ny,nz,lx,ly,lz,S,P,pt\n")==0;
  while (fscanf(input,"%d,%d,%d,%lf,%lf,%lf,%*f,%*f,%*f,%*f,%*f,%d",&x,&y,&z,&nx,&ny,&nz,&tempt)==7){
    if (x>N[0]||y>N[1]||z>N[2]){
      fprintf(stderr,"The point (%d,%d,%d) is out of the bounds!!\n",x,y,z); exit(3);
    }
    ii=x+N[0]*(y+N[1]*z);
    n[ii+0*N[3]]=nx;
    n[ii+1*N[3]]=ny;
    n[ii+2*N[3]]=nz;
    pt[ii]=tempt;
  }
  FILE *teste=fopen("test.csv","w");
  for (int k = 0; k < N[2]; k++)
  {
    for (int j = 0; j < N[1]; j++)
    {
      for (int i = 0; i < N[0]; i++)
      {
        ii=i+N[0]*(j+N[1]*k);
        fprintf(teste,"%d,%d,%d,%0.3f,%0.3f,%0.3f,%d\n",i,j,k,n[ii+0*N[3]],n[ii+1*N[3]],n[ii+2*N[3]],pt[ii]);
      }
      
    }
    
  }
  
  
}
void molp(double* lightInt,double* n,int* pt,int* N){
  double polIn [16], polOut [16], m_LC[16], S[4], S_in[4]={1,0,0,0}, sTemp[4], mTemp[16], mTot[16];
  double theta, phi, nef, del;
  int i, j, k, ii, m1, m2, m3;
  for (i=0; i<16; i++){
    polIn[i] = 0;
    polOut[i]= 0;
    m_LC[i]  = 0;    
  }
  m_LC[0] =1;
  double n0=1.50;
  double ne=1.72;
  double h=19e-6/N[2];
  double lambda=545e-9;
  
  polIn[0] =polIn[5] =0.5;polIn[3] =polIn[12] =0.5;  //linear polarizer oriented in the x-directions 
  polOut[0]=polOut[5]=0.5;polOut[3]=polOut[12]=-0.5; //linear polarizer oriented in the y-directions
    
  for (i=0;i<N[0];i++){
    for (j=0;j<N[1];j++){
      S[0]=1*polIn[0];
      S[1]=1*polIn[4];
      S[2]=1*polIn[8];
      S[3]=1*polIn[12];
      for (int m1 = 0; m1 < 16; m1++)
      {
        mTot[m1]=0;
      }
      mTot[0]=mTot[5]=mTot[10]=mTot[15]=1;
      for (k=0;k<N[2];k++) 
      {
        if (pt[i+N[0]*(j+N[1]*k)]!=1) continue;   
        ii=i+N[0]*(j+N[1]*k);
        //~ del=pjt=0.0;
        phi=acos(n[2*N[3]+ii]);
        nef=1/sqrt(pow(sin(phi)/ne,2) + pow(cos(phi)/n0,2));
        del=2.0*M_PI*(h/lambda)*(1.0*nef-1.0*n0);
       
        if (fabs(n[ii]) > 0.00001) theta=atan2(n[N[3]+ii],n[ii]); 
        else (theta=M_PI/2.0);
        m_LC[5] = cos(2*theta)*cos(2*theta)+sin(2*theta)*sin(2*theta)*cos(del);
        m_LC[6] = sin(2*theta)*cos(2*theta)*(1-cos(del));
        m_LC[7] = sin(2*theta)*sin(del);
        m_LC[9] = sin(2*theta)*cos(2*theta)*(1-cos(del));
        m_LC[10]= sin(2*theta)*sin(2*theta)+cos(2*theta)*cos(2*theta)*cos(del);
        m_LC[11]=-cos(2*theta)*sin(del);
        m_LC[13]=-sin(2*theta)*sin(del);
        m_LC[14]= cos(2*theta)*sin(del);
        m_LC[15]= cos(del);
        
        for (m1=0;m1<4;m1++)
        {
          for (m2=0;m2<4;m2++){
            mTemp[m1+4*m2]=0;
            for (m3=0;m3<4;m3++) mTemp[m1+4*m2]+=mTot[m1+4*m3]*m_LC[m3+4*m2];
          }
        }
        for (m1=0;m1<16;m1++) mTot[m1]=mTemp[m1];
        
      }
      
      for (m1=0;m1<4;m1++)
      {
        for (m2=0;m2<4;m2++){
          mTemp[m1+4*m2]=0;
          for (m3=0;m3<4;m3++) mTemp[m1+4*m2]+=mTot[m3+4*m2]*polIn[m1+4*m3];
        }
      }
      for (m1=0;m1<16;m1++) mTot[m1]=mTemp[m1];
      
      for (m1=0;m1<4;m1++)
      {
        for (m2=0;m2<4;m2++){
          mTemp[m1+4*m2]=0;
          for (m3=0;m3<4;m3++) mTemp[m1+4*m2]+=mTot[m1+4*m3]*polOut[m3+4*m2];
        }
      }
      for (m1=0;m1<16;m1++) mTot[m1]=mTemp[m1];
      
      for (m1=0;m1<4;m1++)
      {
        S[m1]=0;
        for (m2=0;m2<4;m2++){
          S[m1]+=mTot[m2+4*m1]*S_in[m2];
        }
        
      }
      
      lightInt[i+N[0]*j]=S[0];
    }
  }
}

void printData(double *lightInt, int *N, char fname[]){
  char filename[500];
  sprintf(filename,"%s_lightint.dat",fname); 
  FILE *output=fopen(filename,"w");
  int i, j;
  for (i=0; i<N[1]; i++){
    for (j = 0; j < N[0]; j++)
    {
      fprintf(output,"%g ",lightInt[j+N[0]*i]);
    }
    fprintf(output,"\n");
  }
  fclose(output);
}

void gera_plot (char inputFile[]){
    FILE *arq;
    char nome[100];
    sprintf(nome, "%s.plot",inputFile);
    arq = fopen (nome, "w");
    fprintf(arq,"reset\nset terminal png size 720,720\nset output \"%s.png\"\n\
set lmargin at screen 0.0\nset bmargin at screen 0.0\nset rmargin at screen 1.0\n\
set tmargin at screen 1.0\nset pm3d map\nset palette gray\nset pm3d corners2color c1\n\
set pm3d interpolate 5,5\nunset ytics\nunset xtics\nunset colorbox\nset key off\n\
splot '%s_lightint.dat' matrix w pm3d \nreset",inputFile,inputFile);
    fclose(arq);
    
}
