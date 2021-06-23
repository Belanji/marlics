#ifndef  LC_CONTAINER_
#define  LC_CONTAINER_
#include <string>

class container
{

 public:

  container(struct Simulation_Parameters * );
  ~container(void) {};
  void write_state(double  , const double   *,const int *);
  
 protected:
  
  int  file_number;
  char output_folder[200];
  std::string output_fname;
  int  director_field_counter;
  int  Nx, Ny, Nz;


};


#endif
