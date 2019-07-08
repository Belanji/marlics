#ifndef  LC_CONTAINER_
#define  LC_CONTAINER__

class container
{

 public:

  container(struct Simulation_Parameters * );
  ~container(void) {};
  void write_state(double  , const double   *,const int *);
  
 protected:
  
  int  file_number;
  int  director_field_counter;
  int  Nx, Ny, Nz;


};


#endif
