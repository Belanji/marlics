#ifndef GEOMETRY_BULK_
#define GEOMETRY_BULK_


class Geometry_Bulk : public GEOMETRY
{
  
  
 public:


  virtual void fill_ki(double *, const double * ) const ;
  virtual void compute_forces(double *, const double * ) const ;
  virtual void Energy_calc(double * k_i, const double * Qij) const ;
    
  Geometry_Bulk(const struct Simulation_Parameters *);
  ~Geometry_Bulk(void) {};


  
   private:

  virtual int * fill_point_type( void ) const ;  
  

  //void check_bulk_limits(int i, int j, int k) const ;
  //void check_surface_limits(int i, int j, int k) const ;
  
};

#endif
