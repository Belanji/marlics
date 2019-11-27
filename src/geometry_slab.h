#ifndef HOMOGENEOUS_SLAB_
#define HOMOGENEOUS_SLAB_


class slab : public GEOMETRY
{
  
  
 public:


  virtual void fill_ki(double *, const double * ) const ;
    
  slab(const struct Simulation_Parameters *);
  ~slab(void) {};


  
   private:

  virtual int * fill_point_type( void ) const ;  
  

  //void check_bulk_limits(int i, int j, int k) const ;
  //void check_surface_limits(int i, int j, int k) const ;
  
};

#endif
