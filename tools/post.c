/** aerobreakup of a drop 
 */
#include "navier-stokes/centered.h"
#include "two-phaseDOD.h"
#include "tension.h"
//#include "lambda2.h"
#include "view.h"
#include "tag.h"

#define LARGE 1e36

int i_start = 0;
int i_end   = 4998;
int i_gap   = 2;
int fov1 =10;
double tx1=-.25;
double ty1=-.25;
double vmax=3.25;
double pmax=0.5e6;
double l2max=1e-3;
int image_width=1800;
int image_height=1800;

scalar vm[],vml1[],vml2[],vmlboth[];

int main(int argc, char * argv[])
{
  if (argc > 1)
    i_start = atoi (argv[1]);
  if (argc > 2)
    i_end   = atoi (argv[2]);
  if (argc > 3)
    i_gap   = atoi (argv[3]);
  if (argc > 4)
    fov1    = atoi (argv[4]);
  if (argc > 5)
    tx1     = atof (argv[5]);
  if (argc > 6)
    ty1     = atof (argv[6]);
  if (argc > 7)
    vmax    = atof (argv[7]);
  if (argc > 8)
    pmax    = atof (argv[8]);
  if (argc > 9)
    image_width  = atoi (argv[9]);
  if (argc > 10)
    image_height = atoi (argv[10]);
  if (argc > 11)
    l2max   = atof (argv[11]);

  char name[80];

  int i=0; 
  int k=0;
	  //this while loop has been added so that we can proccess dump files with more than 1000 sequence 
	  //like dump-1.243, the integer part is controlled by iterator k
while(i<1000){
  for ( i = i_start; i <= i_end; i+=i_gap ) {
    sprintf (name, "dump-%d.%03d",i/1000,i%1000);
    restore(file = name);

    fprintf (ferr, "Opening file i=%d \n",i);

//    view (fov = fov1, tx=tx1, ty=ty1, quat = {0,0,0,1}, width = image_width, height = image_height);
    view (fov = fov1, tx=tx1, ty=ty1, width = image_width, height = image_height);

// colorless sketch
#if 0
    draw_vof("f1");
    draw_vof("f2");
   // squares("f", linear = true, max = 1., min = 0.);
    box();
    //cells();
    sprintf (name, "vof_%04d.ppm",i);
    save(file=name);
#endif  

#if 1
    foreach(){
      vm[] = sqrt(u.x[]*u.x[]+u.y[]*u.y[]);
      //vml1[] = f1[]*sqrt(u.x[]*u.x[]+u.y[]*u.y[]);
      //vml2[]=f2[]*sqrt(u.x[]*u.x[]+u.y[]*u.y[]);
      vmlboth[]=((f2[]+f1[])>0)*sqrt(u.x[]*u.x[]+u.y[]*u.y[]);
    }
    //boundary ({vm,vml});
    
    clear();
    draw_vof ("f2", color = "vml2", min = 0, max = 1.5);
    draw_vof ("f1", color = "vml1", min =0, max=1.5);

    //squares("vml1", n = {0, 0, 1}, linear = true, max = vmax/5., min = 0.);
    squares("vmlboth",n= {0,0,1},linear=true,max=vmax,min=0);
    sprintf (name, "vof-velLiq_%04d.ppm",i);
    save(file=name);//vml  IS velocity inside liquid

    clear();
    draw_vof ("f1");
    draw_vof("f2");
    squares("vm", n = {0, 0, 1}, linear = true, max = vmax, min = 0);
    sprintf (name, "vof-velmag_%04d.ppm",i);
    save(file=name);// velocity everywhere
/**
    clear();
    draw_vof ("f1");
    draw_vof("f2");
    squares("p", n = {0, 0, 1}, linear = true, max = pmax, min = -pmax);
    sprintf (name, "vof-pre_%04d.ppm",i);
    save(file=name);
   */ 
#endif  
    
#if 0

    view (fov = fov1, tx=tx1, ty=ty1, quat = {0,-0.234,0,0.971}, width = image_width, height = image_height);
    clear();
    draw_vof ("f", color = "vm", min = 0, max = vmax/5.);
    sprintf (name, "vof-vel-backview_%04d.ppm",i);
    save(file=name);
   
    view (fov = fov1, tx=tx1, ty=ty1, quat = {0,0.234,0,0.971}, width = image_width, height = image_height);
    clear();
    draw_vof ("f", color = "vm", min = 0, max = vmax/5.);
    sprintf (name, "vof-vel-frontview_%04d.ppm",i);
    save(file=name);


    view (fov = fov1, tx=tx1, ty=ty1, quat = {0.4,0.55,0.35,0.7}, width = image_width, height = image_height);
    clear();
    draw_vof ("f", color = "vm", min = 0, max = vmax/5.);
    sprintf (name, "vof-vel-front_%04d.ppm",i);
    save(file=name);


    view (fov = fov1, tx=tx1, ty=ty1, quat = {0,-0.234,0,0.971}, width = image_width, height = image_height);
    clear();
    draw_vof ("f", color = "vm", min = 0, max = vmax/5.);
    scalar l2[];
    lambda2 (u, l2);
    isosurface ("l2", -l2max);
    sprintf (name, "vof-l2-backview_%04d.ppm",i);
    save(file=name);
    
    view (fov = fov1, tx=tx1, ty=ty1, quat = {0,0.234,0,0.971}, width = image_width, height = image_height);
    clear();
    draw_vof ("f", color = "vm", min = 0, max = vmax/5.);
    isosurface ("l2", -l2max);
    sprintf (name, "vof-l2-frontview_%04d.ppm",i);
    save(file=name);
 
#endif

#if _GfsView
    double time = i*0.001;
    sprintf (name, "snapshot-%05.3f.gfs", time);
    output_gfs (file = name, t=time, list = {f,u,p});
#endif

#if _dropstats
    scalar m[];
foreach()
      m[] = f1[] > 1e-6;
    int n3 = tag (m);

    double v3[n3];
    coord b3[n3];
    double ux[n3],uy[n3],uz[n3],area[n3];
    for (int j = 0; j < n3; j++)
      v3[j] = b3[j].x = b3[j].y = b3[j].z = ux[j] = uy[j] = uz[j] = area[j] = 0.;

    fprintf(ferr,"We are at point 1\n");
    foreach_leaf()
      if (m[] > 0) {
	      fprintf(ferr,"we are at point 2\n");
         int j = m[] - 1;
         v3[j] += dv()*f[];
         coord p = {x,y,z};
         foreach_dimension()
           b3[j].x += dv()*f[]*p.x;
         ux[j] += dv()*f[]*u.x[];
         uy[j] += dv()*f[]*u.y[];
         uz[j] += dv()*f[]*u.z[];
         if (f[] > 1e-6 && f[] < 1. - 1e-6) {
            coord n = mycs (point, f), p;
            double alpha = plane_alpha (f[], n);
            area[j] += pow(Delta, dimension - 1)*plane_area_center (n, alpha, &p);
         }
      }

#if _MPI
    MPI_Allreduce (MPI_IN_PLACE, v3, n3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce (MPI_IN_PLACE, b3, 3*n3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce (MPI_IN_PLACE, ux, n3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce (MPI_IN_PLACE, uy, n3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce (MPI_IN_PLACE, uz, n3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce (MPI_IN_PLACE, area, n3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif


  sprintf (name, "drop-stats-%04d.dat", i);
  FILE * fp = fopen (name, "w");
  fprintf (fp, "1:drop_id; 2:vol; 3:area; 4:xc; 5:yc; 6:zc; 7:uc; 8:vc; 9:wc\n"); 
  for (int j = 0; j < n3; j++)
    fprintf (fp, "%d %g %g %g %g %g %g %g %g\n", 
                    j, v3[j], area[j],
                         b3[j].x/v3[j], b3[j].y/v3[j], b3[j].z/v3[j],
               ux[j]/v3[j],  uy[j]/v3[j],  uz[j]/v3[j]);
  fclose (fp);

#endif   
  }
}

k=k+1;
i=0;
}//while loop ends here
