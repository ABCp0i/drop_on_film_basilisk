/** Normal impact of a drop on a thin film 
 */
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phaseDOD.h"
#include "tension.h"
//#include "adapt2.h"
#include "tag.h"
//include tag for collecting droplet statistics

#define LARGE 1e36
#define vol_cut 3.0e-5

double max_level = 9;
double L = 4.;
double t_out = 0.002;       
double t_end = 0.990;    

/** dimensionless properties, normalized by scaling variables rhol, D, sigma
 */
double rhog=0.001227;
double mul =0.002487;
double mug =4.6e-5;
double u0  =2.345406;          //falling velocity
double e   =0.013453;            //film thickness
double h   = 0.05;          //initial gap between drop and film

double femax = 0.001;
double uemax = 0.001;

double gravity = 9.8;     //gravity acceleration (m/s2)
double gx;  //dimensionless gravity in x direction  
double spherewidth=1.5;
double sphereheight=.75;


scalar ink[];

u.n[left]  = dirichlet(0);
u.t[left]  = dirichlet(0);
p[left]    = neumann(0);

u.n[right] = neumann(0);
p[right] = dirichlet(0);
pf[right] = dirichlet(0);

int main(int argc, char * argv[])
{

  if (argc > 1)
    max_level = atoi (argv[1]);
  if (argc > 2)
    L         = atof (argv[2]);
  if (argc > 3)
    u0        = atof (argv[3]);
  if (argc > 4)
    t_out     = atof (argv[4]);
  if (argc > 5)
    t_end     = atof (argv[5]);
  if (argc > 6)
    e         = atof (argv[6]);
  if (argc > 7)
    h         = atof (argv[7]);
  if (argc > 8)
    rhog      = atof (argv[8]);
  if (argc > 9)
    mul       = atof (argv[9]);
  if (argc > 10)
    mug       = atof (argv[10]);
  if (argc > 11)
    femax     = atof (argv[11]);
  if (argc > 12)
    uemax     = atof (argv[12]);
  if (argc >13)
    spherewidth=atof(argv[13]);
  if (argc >14)
    sphereheight=atof(argv[14]);
    

  size (L);
  init_grid (256);

  /**
  The liquid phase is water, rho_l=1000 kg/m3, mu_l=1e-3 Ps s; 
  the gas phase is air, rho_g = 1.2 kg/m3, mu_g = 1.7e-5 Pa s; 
  surface tension is sigma=0.0688 N/m; 
  The dimensionless parameters can then be computed based on 
  rho_l, sigma, and D. */

  rho1 = 1., rho2 = rhog; 
  mu1 = mul, mu2 = mug;
  f1.sigma = 1.5;
  f2.sigma = 0.7;
  TOLERANCE = 0.01;  
  //gx  = gravity*1000.*sq(R0)/0.0688;

  run();
}

/**
The initial drop is spherical. */

event init (t = 0)
{
  if (!restore (file = "dump")) {

    double x0 = e+0.5+h; 

    refine ( ( x> 0.8*e && x<1.2*e && level < max_level) || ( sq(y)+sq(x-x0)<sq(0.6) && sq(y)+sq(x-x0)>sq(0.4) && level < max_level ) );
    fraction (f1, -sq(y)/spherewidth-sq(x-x0)/sphereheight+sq(0.5));
    fraction (f2, x<=e);

    /** Initial velocity
    */
    foreach() {
      u.x[] = -u0*f1[];
    }

  }
}


/**

static scalar * interfaces1 = NULL;
event vof (i++) {
 // boundary ({f, ink, uf});
 
 // ink.inverse = false;
 // ink.gradient = minmod2;
 // f.tracers = {ink};

 // vof_advection({f}, i);

 // boundary({f, ink});

  //We set the list of interfaces to NULL so that the default vof() event does nothing (otherwise we would transport f twice).
  interfaces1 = interfaces, interfaces = NULL;
}
*/

//We set the list of interfaces back to its default value.
#if 0
event tracer_advection (i++) {
  interfaces = interfaces1;
}
#endif



event logfile (i+=10)
{
  scalar posy[],posx[],posx_y0_f1[],posx_y0_f2[];
  position (f1,posx,{1,0});
  position (f1,posy,{0,1});
  position (f1,posx_y0_f1,{1,0});
  position (f2,posx_y0_f2,{1,0});


  double area=0.,vol=0.,keg=0.,ked=0.,ud=0.,xd=0.,kef=0.;

  foreach(reduction(+:area) reduction(+:vol)
          reduction(+:keg)  reduction(+:ked)
          reduction(+:ud) reduction(+:xd) reduction(+:kef)) {

    if ( y > Delta ) {
      posx_y0_f1[] = nodata;
      posx_y0_f2[] = nodata;
    }
    

    double dv_axi = pow(Delta, 2.)*2.*pi*y;

    /** kinetic energy */
    double ke = 0.;
    foreach_dimension() {
      ke  += sq(u.x[]);
    }
    ked += dv_axi*    f1[] *ke;
    keg += dv_axi*(1.-f1[]-f2[])*ke;
    kef += dv_axi*f2[]*ke;
    
    /** statistics for drop before impact */
    if ( x> e*2. ) {
      /** Volume */
      vol += dv_axi*f1[];

      /** mean velocity*/
      ud += dv_axi*f1[]*u.x[];

      /** centroid */
      xd += dv_axi*f1[]*x;
    }
  }
  keg /= 2.;
  ked /= 2.;
  ud  /= vol;
  xd  /= vol;

  stats sx = statsf (posx);
  stats sy = statsf (posy);
  stats sx1_y0 = statsf (posx_y0_f1);
  stats sx2_y0 = statsf (posx_y0_f2);

  if ( i == 0 )
    fprintf(ferr,
      "#1: t; 2: dt; 3: xc; 4, uc; 5:ymax_drop; 6:xmax_drop; 7: xmin_drop; 8: xmax_axis_drop; 9:xmin_axis_drop; 10: xmax_axis_film; 11: KE_gas; 12: KE_drop; 13:KE_film; 14: vol_drop; 15: area_drop; 16: n_grid; 17: cput; 18: speed; 19: u_xmax; 20: u_ymax\n");

  if ( i >10 )
    fprintf (ferr, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %ld %g %g %g %g \n",
      t, dt, xd, ud,
      sy.max, sx.max, sx.min, sx1_y0.max, sx1_y0.min, sx2_y0.max, 
      keg, ked, kef, vol, interface_area(f1),
      grid->tn, perf.t, perf.speed, statsf(u.x).max,statsf(u.y).max);

}

#if 0
/** output snapshot of simulation */
event interface (t += t_out; t<= t_end) {

  char name[80];
  sprintf (name, "infc-%05.3f.dat", t);
  FILE * fp1 = fopen (name, "w");
  output_facets (f,fp1);

}
#endif 
  
#if 1	
event snapshot (t += t_out; t<=t_end ) {
  char name[80];
#if 0
  sprintf (name, "snapshot-%05.3f.gfs", t);
  output_gfs (file = name, t=t, list = {f1,u,p,ink});
#endif
  sprintf (name, "dump-%05.3f", t);
  dump (file = name); // so that we can restart

#if 1
  sprintf (name, "infcDrop-%05.3f.dat", t);
  FILE * fp1 = fopen (name, "w");
  output_facets (f1,fp1);
  
  sprintf (name, "infcFilm-%05.3f.dat", t);
  fp1 = fopen (name, "w");
  output_facets (f2,fp1);
#endif 
}
#endif

/**
Adapt mesh based on the volume fraction. */
event adapt (i=10; i++) {
  scalar vliq[];
  foreach() {
    vliq[]=u.y[]*f1[];
  }

  double uemax_liq;
  uemax_liq=uemax*0.01;

  
  //adapt_wavelet ({f2,f1}, (double[]){femax,femax}, minlevel = 5, maxlevel = max_level);
  adapt_wavelet ({f2,f1,u}, (double[]){femax,femax,uemax,uemax}, minlevel = 5, maxlevel = max_level);
  //adapt_wavelet ({f2,f1,u,vliq}, (double[]){femax,femax,uemax,uemax,uemax_liq}, minlevel = 5, maxlevel = max_level);
  //adapt_wavelet2((scalar *){f,u.x,u.y},(double[]){femax,uemax,uemax},(int []){max_level,max_level-3,max_level-3},3);
}

/*
event remove_drops(i+=10)
{
  scalar m[];
  foreach()
    m[]=f1[]+f2[]>0.;
  int n =tag(m);
  double v[n];
  int remove_flag[n];
  for(int j=0;j<n;j++){
    v[j]=0.;
    remove_flag[j]=0.;
  }
  foreach_leaf()
    if(m[]>0){
      int j = m[] - 1;
      v[j]+=dv()*(f1[]+f2[]);
    }
#if _MPI
  MPI_Allreduce(MPI_IN_PLACE,v,n,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif

  for (int j =0;j<n;j++)
    if (v[j]<vol_cut)
      remove_flag[j]=1;
  foreach()
    if(m[]>0){
      int j =m[]-1;
      if (remove_flag[j] == 1)
        f1[]=0.;
      	f2[]=0.;
    }
}

**/


