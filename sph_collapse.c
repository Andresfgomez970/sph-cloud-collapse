#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

//  Macros for positions and velocities
#define X 0
#define Y 1
#define Z 2

//  Canonical units of the system
#define UM 1.989e30  // Kg
#define UL 3.086e16   // meters
#define UT sqrt(UL*UL*UL / (UM*G) )  // Time unit
#define UTEMP 80.0  // Kelvins

//  Physical constants
#define K_B 1.38064852e-23*(UT*UT*UTEMP)/(UL*UL*UM)
#define K_BM_H 1.38064852e-23*(UT*UT*UTEMP)/(UL*UL)/1.6735575e-27
#define G 6.67430e-11

typedef struct
{
  int id;
  double pos[3];
  double vel[3];
  double accel[3];
  double mass;
  double rho;
  double h;
  double p;
  double c;
  double du;
  double u;
  int *nn;  // Pointer to savet the neighbors
  int nNeighbors;  // Number of neighbors
  double *r;
  double *W;
  double *dx;  // Pointer to save the derivate x
  double *dy;
  double *dz;
  double *dWx;  // Pointer to save the derivate of kernel respect x
  double *dWy;
  double *dWz;
  //  1 == fluid, 2 == sink particles
  int type;
}Particles;

// Pointer to structures of particles and auxiliar one in order to prevent
//   errors.
Particles *part, *auxPart;
int nFluid, nPart;  // SPH and gravity particles
double M_T = 1000.;
double R0;
int N_dead;
//  Mean distance between evenly distributed points inside a sphere of radius
//    R0
double h_0;


void ics(double R0, int N_tot);
double gaussian(double V,double sigma,double mu);
int conditions(double x, double y, double z,double vx, double vy, double vz,
  double R0, double sigma);


double W(double r, double h);
double dW(double r, double dx, double h);
void testKernel(void);
void NN(int i);
void test_NN(void);
void density(void);
void eos(void);
void navierStokes(void);
void viscosity(double dx);
void BorderInteraction(double dx);
void meanVelocity();
void acceleration(double dx);
void drift(double dt);
void kick(double dt);
void printState(char *outfile);
int gravitationalAcceleration(double dx);
void boundary();
void sinkParticles();

int main(int argc, char *argv[])
{
  /////////////////////////////////////////////////////////////////////////////
  //                    Initialization of variables                          //
  /////////////////////////////////////////////////////////////////////////////
  float N_iter;
  int i,counter;
  double dt;
  double t, tTotal;
  char outfiles[500];
  double dx;


  dt = 1e-2;
  N_iter = 1000;
  tTotal = N_iter*dt;
  t = 0;
  R0 = 10.0;
  h_0 = 148.0/(45.0*M_PI)*R0;

  nFluid = 1000;
  N_dead = 0;

  //  Allocating memory for global varibles
  part = (Particles *)malloc((size_t)nFluid*sizeof(Particles));

  if( part==NULL )
    {
      printf("Error alocando part\n");
      exit(0);
    }

  // Create the initial conditions
  ics(R0,nFluid);

  printf("------- Printing canonical units and valuable scalars ------ \n");

  printf("UM (Mass Unit)  = %e kg\n",UM);
  printf("UL (Lenght Unit)  = %e m\n",UL);
  printf("UTEMP  (Temperature Unit) = %e K\n",UTEMP);
  printf("UT (Time Unit) = %e s\n",UT/3.14e7);
  printf("M_T (Total Mass) = %e UM\n",M_T);
  printf("K_B  = %e\n",K_B);
  printf("K_BM_H  = %e\n",K_BM_H);
  //////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////
  //                          System evolution                               //
  /////////////////////////////////////////////////////////////////////////////

  counter = 0;
  //  printting system initial state
  sprintf(outfiles,"./output/state_%.4d",counter);
  printState(outfiles);


  //  main loop
  t = 0;
  while( t<=tTotal )
    {
      //  searching near neighbors for all fuid particles
      for( i=0; i<nFluid; i++ )
	      NN(i);

      //  Sink particles search and determination
      sinkParticles();
      //  computing density
      density();

      //  drift in leap-frog integration
      drift(dt);

      //  computing acceleration
      acceleration(dx);

      //  kick in leap-frog integration
      kick(dt);

      //  drift in leap-frog integration
      drift(dt);

      t = t + dt;
      counter++;
      printf("step = %d \n\n",counter);

      //  printting system state
      sprintf(outfiles,"./output/state_%.4d",counter);
      printState(outfiles);
    }

  // Print number of dead particles and free memory used
  printf("N_dead = %d\n",N_dead );
  free(part);

  return 0;
}


double gaussian(double v,double sigma,double mu)
{
  return exp(- (v-mu)*(v-mu)/(2*sigma*sigma) );
}

// Function that evaluates that randomly generated positions and
// velocities met the conditions of the Monte-Carlo sampling
int conditions(double x, double y, double z, double vx, double vy, double vz,  double R0, double sigma)
{
  double r;
  double dist_r,dist_vx,dist_vy,dist_vz;
  r = sqrt(x*x + y*y + z*z);
  dist_r = drand48();
  dist_vx = drand48();
  dist_vy = drand48();
  dist_vz = drand48();

  if( (r < R0) && (dist_vx < gaussian(vx, sigma/3.0, 0)) && (dist_vy < gaussian(vy, sigma/3.0, 0)) && (dist_vz < gaussian(vz, sigma/3.0, 0)) )
    return 1;
  else
    return 0;
}

void ics(double R0, int N_tot)
{
  int i,j,counter;
  double x,y,z,r,dist_N,sigma;
  double vx,vy,vz,vphi,L_specific;
  double ***Border;  // Border[theta][phi][x]  = x_i position
  FILE *fFluidIcs,*fBorder;

  fFluidIcs = fopen("fluid_ics.output","w");  //  fluid initital particles
  fBorder = fopen("border.output","w");

  //  Total specific angular momentum of a typical molecular cloud
  L_specific = (1e22*UT)/((UL*100)*(UL*100));

  // ics for fluid particles
  //srand48(time(NULL));
  counter = 0;
  //////////////////////////////////////////////////////////////////////////////
  //                       Generating fluid particles                         //
  //////////////////////////////////////////////////////////////////////////////
  while( counter < N_tot )
  {
    //  Generates random values between -R0 and R0
    x = (- 1 + 2.0*drand48())*R0;
    y = (- 1 + 2.0*drand48())*R0;
    z = (- 1 + 2.0*drand48())*R0;


    // assuming that the particles are independenly equally distributed the mean
    //  deviation decrases to dev/sqrt(N)

    sigma = 1e-3*sqrt((3.0/5.0) * (M_T*M_T/R0)); //Using virial theorem
    vx = (- 1 + 2.0*drand48())*sigma;
    vy = (- 1 + 2.0*drand48())*sigma;
    vz = (- 1 + 2.0*drand48())*sigma;

    vphi = L_specific/sqrt(x*x + y*y);  // Every particle has the same specific angular momentum


    if( conditions(x,y,z,vx,vy,vz,R0,sigma) )
    {
        r = sqrt(x*x + y*y + z*z);
        part[counter].id = counter;
        part[counter].pos[X] = x;
        part[counter].pos[Y] = y;
        part[counter].pos[Z] = z;
        //  Dot product of V_phi . X_unitary
        part[counter].vel[X] = vx - (y/r)*vphi;
        //  Dot product of V_phi . Y_unitary
        part[counter].vel[Y] = vy + (x/r)*vphi;
        part[counter].vel[Z] = vz;
        part[counter].accel[X] = 0.0;
        part[counter].accel[Y] = 0.0;
        part[counter].accel[Z] = 0.0;
        //  Estimated mean density
        part[counter].rho = 3*M_T/(4.0*M_PI*R0*R0*R0);
        //  Smoothing lenght
        part[counter].h = h_0;
        part[counter].mass = M_T/N_tot;
        part[counter].p = 0.0;
        part[counter].c = 0.0;
        part[counter].du = 0.0;
        part[counter].u = 3.0*K_B*UTEMP/2.0;
        part[counter].nn = NULL;
        part[counter].nNeighbors = 0;
        part[counter].r = NULL;
        part[counter].W = NULL;
        part[counter].dWx = NULL;
        part[counter].dWy = NULL;
        part[counter].type = 1;
        counter++;
      }
   }

   for( i=0 ; i < counter ; i++)
      fprintf(fFluidIcs, "%d %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf \n",
              i , part[i].pos[X], part[i].pos[Y], part[i].pos[Z],
              part[i].vel[X], part[i].vel[Y], part[i].vel[Z],
              part[i].accel[X], part[i].accel[Y], part[i].accel[Z],
              part[i].rho, part[i].mass, part[i].p, part[i].c, part[i].du, part[i].u);

  //////////////////////////////////////////////////////////////////////////////
  //                 Generando las particulas de las fronteras                //
  //////////////////////////////////////////////////////////////////////////////

  //printf("%d\n",counter);
  int kappa_f = 4;
  // Kappa_frontier is less than 2*kappa to avoid particles from scaping
  double theta_s = kappa_f*h_0/R0;
  int N_theta,N_phi;

  //  Evenly spaced particles in theta implies that in phi the spacing is less
  N_theta = (int) (M_PI/theta_s + 1); //  Theta range
  N_phi = (int) (2*M_PI/theta_s + 1); //  Phi range

  //  Setting auxPart to NULL for reallocation
  auxPart = NULL;

  //  Constructing border
  Border = (double ***)malloc( (size_t)(N_theta)*sizeof(double **) );
  auxPart = (Particles *)realloc(part,(size_t)(nFluid+N_theta*N_phi)*sizeof(Particles));
  if(auxPart==NULL)
    {
      printf("error en auxPart\n");
      exit(0);
    }
  else
    {
      part = auxPart;
      auxPart = NULL;
    }

  for (i=0; i<N_theta; i++){
    Border[i] = (double **)malloc( (size_t) (N_phi)*sizeof(double*) );
    for(j=0; j<N_phi; j++){
      Border[i][j] = (double *)malloc( (size_t) 3 * sizeof(double) );
      //  From spherical coordinates to cartesian
      Border[i][j][X] = R0*sin(i*theta_s)*cos(j*theta_s);
      Border[i][j][Y] = R0*sin(i*theta_s)*sin(j*theta_s);
      Border[i][j][Z] = R0*cos(i*theta_s);
      part[counter].id = counter;
      part[counter].pos[X] = Border[i][j][X];
      part[counter].pos[Y] = Border[i][j][Y];
      part[counter].pos[Z] = Border[i][j][Z];
      part[counter].vel[X] = 0;
      part[counter].vel[Y] = 0;
      part[counter].vel[Z] = 0;
      part[counter].accel[X] = 0.0;
      part[counter].accel[Y] = 0.0;
      part[counter].accel[Z] = 0.0;
      part[counter].rho = 0.0;
      part[counter].h = 0.0;
      part[counter].mass = M_T/N_tot;
      part[counter].p = 0.0;
      part[counter].c = 0.0;
      part[counter].du = 0.0;
      part[counter].u = 357.1;  ///Arreglar
      part[counter].nn = NULL;
      part[counter].nNeighbors = 0;
      part[counter].r = NULL;
      part[counter].W = NULL;
      part[counter].dWx = NULL;
      part[counter].dWy = NULL;
      part[counter].type = -1;
      counter++;
    }
  }


  for( i=nFluid ; i < counter; i++){
     fprintf(fFluidIcs, "%d %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf \n",
             i , part[i].pos[X], part[i].pos[Y], part[i].pos[Z],
             part[i].vel[X], part[i].vel[Y], part[i].vel[Z],
             part[i].accel[X], part[i].accel[Y], part[i].accel[Z],
             part[i].rho, part[i].mass, part[i].p, part[i].c, part[i].du, part[i].u);
    fprintf(fBorder, "%d %.10lf %.10lf %.10lf \n", i , part[i].pos[X], part[i].pos[Y], part[i].pos[Z]);
  }

  free(Border);
  fclose(fBorder);
  fclose(fFluidIcs);

}


double W(double r, double h)
{

  double R = r/h;

  double alpha = 15.0/(7.0*M_PI*h*h);

  if( (R >= 0.0) && (R < 1.0) )
    return alpha*((2.0/3.0) - R*R + 0.5*R*R*R);

  if( (R >= 1.0) && (R <= 2.0) )
    return alpha*((1.0/6.0)*(2.0-R)*(2.0-R)*(2.0-R));

  if( R>2.0)
    return 0.0;

  return 0.0;
}

double dW(double r, double dx, double h)
{
  double R = r/h;

  double alpha = 15.0/(7.0*M_PI*h*h);

  if( (R >= 0.0) && (R < 1.0) )
    return alpha*(-2.0 + 1.5*R)*dx/(h*h);

  if( (R >= 1.0) && (R <= 2.0) )
    return alpha*(-0.5*(2.0-R)*(2.0-R))*dx/(h*h*R);

  if( R>2.0)
    return 0.0;

  return 0.0;
}

void testKernel(void)
{
  double r, w, dw;

  FILE *fKernelTest;
  fKernelTest = fopen("kernel_test.output","w");

  for( r=-3.0; r<=3.0; r = r + 0.05)
  {
    w = W( fabs(r), 1.0);
    dw = dW( fabs(r), r/sqrt(3.0), 1.0);

    fprintf(fKernelTest,"%16.10lf %16.10lf %16.10lf\n",r,w,dw);

  }
  fclose(fKernelTest);
}

// Searching the near neighbors
void NN(int i)
{

  double kappa = R0/100.0;
  double xij, yij, zij, rij, hij;
  double *auxDouble;
  int j, *auxInt, nNeighbors;

  part[i].nn = NULL;
  part[i].dx = NULL;
  part[i].dy = NULL;
  part[i].dz = NULL;
  part[i].r = NULL;
  part[i].W = NULL;
  part[i].dWx = NULL;
  part[i].dWy = NULL;
  part[i].dWz = NULL;


  nNeighbors = 0;

  for( j=0; j<nFluid; j++ )
  {
    if(part[j].type == 1)
    {
      if( i!=j && part[i].type == 1)
	    {
    	  xij = part[i].pos[X] - part[j].pos[X];
    	  yij = part[i].pos[Y] - part[j].pos[Y];
        zij = part[i].pos[Z] - part[j].pos[Z];
        rij = sqrt( xij*xij + yij*yij + zij*zij );
    	  hij = 0.5*(part[i].h + part[j].h);


	      if( rij < kappa*hij  )
	      {
	      nNeighbors++;
        //printf("%d\n",nNeighbors);
	      // add neighbor id
	      auxInt = NULL;
	      auxInt = (int *)realloc(part[i].nn,(size_t)(nNeighbors)*sizeof(int));
	      part[i].nn = auxInt;

	      // add neighbor dx
	      auxDouble = NULL;
	      auxDouble = (double *)realloc(part[i].dx,(size_t)(nNeighbors)*sizeof(double));
	      part[i].dx = auxDouble;

	      // add neighbor dy
	      auxDouble = NULL;
	      auxDouble = (double *)realloc(part[i].dy,(size_t)(nNeighbors)*sizeof(double));
	      part[i].dy = auxDouble;

        // add neighbor dz
        auxDouble = NULL;
        auxDouble = (double *)realloc(part[i].dz,(size_t)(nNeighbors)*sizeof(double));
        part[i].dz = auxDouble;

	      // add neighbor r
	      auxDouble = NULL;
	      auxDouble = (double *)realloc(part[i].r,(size_t)(nNeighbors)*sizeof(double));
	      part[i].r = auxDouble;

	      // add neighbor W
	      auxDouble = NULL;
	      auxDouble = (double *)realloc(part[i].W,(size_t)(nNeighbors)*sizeof(double));
	      part[i].W = auxDouble;

	      // add neighbor dWx
	      auxDouble = NULL;
	      auxDouble = (double *)realloc(part[i].dWx,(size_t)(nNeighbors)*sizeof(double));
	      part[i].dWx = auxDouble;

	      // add neighbor dWy
	      auxDouble = NULL;
	      auxDouble = (double *)realloc(part[i].dWy,(size_t)(nNeighbors)*sizeof(double));
	      part[i].dWy = auxDouble;

        // add neighbor dWz
        auxDouble = NULL;
        auxDouble = (double *)realloc(part[i].dWz,(size_t)(nNeighbors)*sizeof(double));
        part[i].dWz = auxDouble;


	      part[i].nn[nNeighbors-1] = j;
	      part[i].dx[nNeighbors-1] = xij;
	      part[i].dy[nNeighbors-1] = yij;
        part[i].dz[nNeighbors-1] = zij;
	      part[i].r[nNeighbors-1] = rij;
	      part[i].W[nNeighbors-1] = W( rij, hij );
	      part[i].dWx[nNeighbors-1] = dW( rij, xij, hij);
	      part[i].dWy[nNeighbors-1] = dW( rij, yij, hij);
        part[i].dWz[nNeighbors-1] = dW( rij, zij, hij);

	     }
     }
   }
  }
  part[i].nNeighbors = nNeighbors;
}

void test_NN(void)
{
  int i,j,k;

  FILE *fTestNN, *fTestPart;
  fTestNN = fopen("NN_test_neighbors.output","w");
  fTestPart = fopen("NN_test_particle.output","w");

  //srand(time(NULL));
  srand(0);

  for(k=0 ; k<20 ; k++)
  {
    i = rand() % nFluid;

    //printf("Testing for particle %d \n",i);
    //printf("With %d neighbors \n",part[i].nNeighbors );

    fprintf(fTestPart, "%d %16.10lf %16.10lf %16.10lf\n",
            part[i].id,
            part[i].pos[X],
            part[i].pos[Y],
            part[i].pos[Z]);

    for(j = 0; j<part[i].nNeighbors; j++)
    {
      fprintf(fTestNN,"%d %16.10lf %16.10lf %16.10lf \n",
              part[i].nn[j],
              part[ part[i].nn[j] ].pos[X],
              part[ part[i].nn[j] ].pos[Y],
              part[ part[i].nn[j] ].pos[Z]);
    }
  }
  fclose(fTestNN);
  fclose(fTestPart);
}

void density(void)
{
  int i, j;
  double wii, norm;


  for( i=0; i<nFluid; i++ )
  {
    if(part[i].type != 0)
    {
      // self density
      wii = W( 0.0, part[i].h );
      //printf(" wii = %e\n",wii );
      part[i].rho = part[i].mass*wii;

      // computing density
      for( j=0; j<part[i].nNeighbors; j++ )
         part[i].rho = part[i].rho + part[part[i].nn[j]].mass*part[i].W[j];
      // normalizing the density

      norm = (part[i].mass/part[i].rho)*wii;
      for( j=0; j<part[i].nNeighbors; j++ )
         norm = norm + (part[part[i].nn[j]].mass/part[part[i].nn[j]].rho)*part[i].W[j];
      part[i].rho = part[i].rho/norm;
    }
  }
  printf("density computed\n");
}

void eos(void)
{
  int i;

  for( i=0; i<nFluid; i++ )
  {
    if(part[i].type != 0)
    {
      //  Equation of state for ideal gas
      part[i].p = (UTEMP*K_BM_H*part[i].rho)/(2.0) ;
    }
  }
}

void navierStokes(void)
{

  int i, j, k;
  double pij, vdw;
  // computing sound speed and pression
  eos();

  // computing acceleration
  for( i=0; i<nFluid; i++ )
  {
    if(part[i].type != 0)
    {
      part[i].accel[X] = part[i].accel[Y]  = part[i].accel[Z] = 0.0;
      part[i].du = 0.0;

      for( k=0; k<part[i].nNeighbors; k++ )
	    {
	     j = part[i].nn[k];
	     pij = ( part[i].p/(part[i].rho*part[i].rho) )
  	     + ( part[j].p/(part[j].rho*part[j].rho) );

	     part[i].accel[X] = part[i].accel[X] - part[j].mass*pij*part[i].dWx[k];
	     part[i].accel[Y] = part[i].accel[Y] - part[j].mass*pij*part[i].dWy[k];
       part[i].accel[Z] = part[i].accel[Z] - part[j].mass*pij*part[i].dWz[k];

	     vdw = (part[i].vel[X]-part[j].vel[X])*part[i].dWx[k]
	       + (part[i].vel[Y]-part[j].vel[Y])*part[i].dWy[k] + (part[i].vel[Z]-part[j].vel[Z])*part[i].dWz[k];
	     part[i].du = part[i].du + 0.5*part[j].mass*pij*vdw;
	    }
    }
  }
  printf("acceleration computed\n");
}

void viscosity(double dx)
{

  int i, j, k;

  double xij, yij, zij, vxij, vyij, vzij, vijrij, vdw;
  double hij, cij, phiij, rhoij, Piij;
  double alphapi = 1.0; // verify if this valus must be changed
  double betapi = 1.0;
  double eps = dx;
  double eps2 = 0.01*eps*eps;

  for( i=0; i<nFluid; i++ )
  {
    if(part[i].type != 0)
    {
    for( k=0; k<part[i].nNeighbors ; k++ )
	   {
	     j = part[i].nn[k];
	     xij = part[i].pos[X] - part[j].pos[X];
	     yij = part[i].pos[Y] - part[j].pos[Y];
       zij = part[i].pos[Z] - part[j].pos[Z];
  	   vxij = part[i].vel[X] - part[j].vel[X];
  	   vyij = part[i].vel[Y] - part[j].vel[Y];
       vzij = part[i].vel[Z] - part[j].vel[Z];

       vijrij = vxij*xij + vyij*yij + vzij*zij;

	     if( vijrij < 0.0 )
	     {
	      hij = 0.5*(part[i].h+part[j].h);
	      phiij = (hij*vijrij)/( xij*xij + yij*yij + zij*zij + eps2);
	      cij = 0.5*(part[i].c+part[j].c);
	      rhoij = 0.5*(part[i].rho+part[j].rho);

	      Piij = ( -alphapi*cij*phiij + betapi*phiij*phiij )/( rhoij );

	      part[i].accel[X] = part[i].accel[X] - part[j].mass*Piij*part[i].dWx[k];
	      part[i].accel[Y] = part[i].accel[Y] - part[j].mass*Piij*part[i].dWy[k];
        part[i].accel[Z] = part[i].accel[Z] - part[j].mass*Piij*part[i].dWz[k];

	      vdw = (part[i].vel[X]-part[j].vel[X])*part[i].dWx[k]
		      + (part[i].vel[Y]-part[j].vel[Y])*part[i].dWy[k] + (part[i].vel[Z]-part[j].vel[Z])*part[i].dWz[k];
	      part[i].du = part[i].du + 0.5*part[j].mass*Piij*vdw;
	     }
	   }
    }
  }
  printf("viscosity computed\n");
}

// this is by now commmented in the acceleration function
void BorderInteraction(double dx)
{

  int i, j;
  int n1 = 12, n2 = 4;
  double r0 = dx/2.0, D = 0.01;
  double xij, yij, zij, rij, PBxij, PByij, PBzij;

  for( i=0; i<nFluid; i++ )
  {
    for( j=0; j<part[i].nNeighbors; j++ )
  	{
	    if( part[part[i].nn[j]].type==-1 )
	    {
	      xij = part[i].pos[X] - part[part[i].nn[j]].pos[X];
	      yij = part[i].pos[Y] - part[part[i].nn[j]].pos[Y];
        zij = part[i].pos[Z] - part[part[i].nn[j]].pos[Z];
	      rij = sqrt( xij*xij + yij*yij +  zij*zij );

        if( rij<r0 )
		    {
		      PBxij = D*( pow((r0/rij),n1) - pow((r0/rij),n2) )*(xij/(rij*rij));
		      PByij = D*( pow((r0/rij),n1) - pow((r0/rij),n2) )*(yij/(rij*rij));
          PBzij = D*( pow((r0/rij),n1) - pow((r0/rij),n2) )*(zij/(rij*rij));

		      part[i].accel[X] = part[i].accel[X] + PBxij;
		      part[i].accel[Y] = part[i].accel[Y] + PByij;
          part[i].accel[Z] = part[i].accel[Z] + PBzij;
		    }
	    }
	  }
  }
  printf("interaction with Border computed\n");
}

void meanVelocity()
{

  int i , j;
  double epsilon = 0.3;
  double vxMean, vyMean , vzMean;
  double vxij, vyij, vzij, rhoij;

  for( i=0; i<nFluid; i++ )
  {
    if(part[i].type != 0)
    {
      vxMean = 0.0;
      vyMean = 0.0;

      for( j=0; j<part[i].nNeighbors; j++ )
    	{
  	    vxij = part[i].vel[X] - part[part[i].nn[j]].vel[X];
  	    vyij = part[i].vel[Y] - part[part[i].nn[j]].vel[Y];
        vzij = part[i].vel[Z] - part[part[i].nn[j]].vel[Z];
  	    rhoij = 0.5*(part[i].rho+part[part[i].nn[j]].rho);
  	    vxMean = vxMean + (part[part[i].nn[j]].mass/rhoij)*vxij*part[i].W[j];
  	    vyMean = vyMean + (part[part[i].nn[j]].mass/rhoij)*vyij*part[i].W[j];
        vzMean = vzMean + (part[part[i].nn[j]].mass/rhoij)*vzij*part[i].W[j];
  	  }

      part[i].vel[X] = part[i].vel[X] - epsilon*vxMean;
      part[i].vel[Y] = part[i].vel[Y] - epsilon*vyMean;
      part[i].vel[Z] = part[i].vel[Z] - epsilon*vzMean;
    }
  }
}

void acceleration(double dx)
{

  // computing acceleration and change of energy

  navierStokes();

  // computing viscosity contribution
  viscosity(dx);

  // computing interaction with Border -- This problem dos not require
  // this type of interaction
  //BorderInteraction(dx);

  // correction to mean velocity
  meanVelocity();

  gravitationalAcceleration(dx);

  boundary();

  printf("acceleration computed\n");

}

void drift(double dt)
{
  int i;
  for( i=0; i<nFluid; i++ )
  {
    part[i].pos[X] = part[i].pos[X] + 0.5*dt*part[i].vel[X];
    part[i].pos[Y] = part[i].pos[Y] + 0.5*dt*part[i].vel[Y];
    part[i].pos[Z] = part[i].pos[Z] + 0.5*dt*part[i].vel[Z];
    part[i].u = part[i].u + 0.5*dt*part[i].du;
  }
}

void kick(double dt)
{
  int i;
  for( i=0; i<nFluid; i++ )
  {
    part[i].vel[X] = part[i].vel[X] + dt*part[i].accel[X];
    part[i].vel[Y] = part[i].vel[Y] + dt*part[i].accel[Y];
    part[i].vel[Z] = part[i].vel[Z] + dt*part[i].accel[Z];
  }
}


int gravitationalAcceleration(double dx){
  int i, j;
  double dr[3];
  double r,r3;
  double eps2 = dx*dx;

  for(i = 0; i<nFluid; i++)
  {
    if(part[i].type != 0)
    {
      for (j = 0; j<nFluid; j++)
  	  {
	      if (j!=i)
       {
         dr[X] =  part[i].pos[X] - part[j].pos[X];
         dr[Y] =  part[i].pos[Y] - part[j].pos[Y];
         dr[Z] =  part[i].pos[Z] - part[j].pos[Z];
         r = sqrt( dr[X]*dr[X]+dr[Y]*dr[Y]+dr[Z]*dr[Z] + eps2);  // Distance( Particlesptr2[i].pos, Particlesptr2[j].pos);
         r3 = r*r*r;
         part[i].accel[X]  += -G*part[j].mass*dr[X]/r3;
         part[i].accel[Y]  += -G*part[j].mass*dr[Y]/r3;
         part[i].accel[Z]  += -G*part[j].mass*dr[Z]/r3;
       }
	    }
    }
  }
return 0;
}


void printState(char *outfile)
{

  int i;

  FILE *fState;
  fState = fopen(outfile,"w");

  for( i=0; i<nFluid; i++)
  {
    if( part[i].type != 0 )
    {
    fprintf(fState,"%d %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %d\n",
	      part[i].id,
	      part[i].pos[X],part[i].pos[Y], part[i].pos[Z],
	      part[i].vel[X],part[i].vel[Y], part[i].vel[Z],
	      part[i].accel[X],part[i].accel[Y], part[i].accel[Z],
	      part[i].rho,part[i].mass,
	      part[i].p,part[i].c,part[i].u,part[i].type);
    }
  }
  fclose(fState);
}


//Function to make the velocity of particles that try to scape equal to cero
void boundary()
{
  int i;
  double x,y,z;

  for (i = 0; i < nFluid; i++)
  {
    x = part[i].pos[X];
    y = part[i].pos[Y];
    z = part[i].pos[Z];
    if (x*x + y*y + z*z < R0*R0)
    {
        part[i].vel[X] = part[i].vel[Y] = part[i].vel[Z] = 0;
        part[i].pos[X] -= part[i].pos[X]*0.01;
        part[i].pos[Y] -= part[i].pos[Y]*0.01;
        part[i].pos[Z] -= part[i].pos[Z]*0.01;
    }
  }
}


void sinkParticles()
{
  int i,j,k;
  double auxMass,auxrho,m_jeans;

  for(i = 0; i<nFluid ; i++)
  {
    //  This is in case not particles are considered to become sink ones
    auxMass = 0;
    m_jeans = 1;

    //  Calculate local parameters to evaluate local collapse; only do this
    //    for gas-particles with neighbors.
    if( (part[i].type==1) && (part[i].nNeighbors>0))
    {
      auxMass = part[i].mass;
      auxrho  = part[i].rho;
      //  Loop over the neighborhood of i-th particle
      for(j = 0; j<part[i].nNeighbors; j++)
      {
        auxMass += part[ part[i].nn[j] ].mass;
        auxrho  += part[ part[i].nn[j] ].rho;
      }
      auxrho = auxrho/part[i].nNeighbors;
      m_jeans = 0.01*pow(5*UTEMP*K_BM_H/G,1.5)*pow(3.0/(4.0*M_PI*auxrho),0.5);
    }

    //  Evaluation of contition

    //  In this case the part[i] become sink particle, and its neighbors
    //    become dead.
    if(auxMass >= m_jeans)
    {
      part[i].mass = auxMass;
      part[i].type = 2;
      N_dead += part[i].nNeighbors;
      printf("N_dead = %d\n",N_dead );
      printf("sink particle id = %d\n",i);
      for(k = 0; k < part[i].nNeighbors; k++)
      {
        //  This is done to simplify the implementation given that
        //  the SFR efficiency is low and memory is not an issue for low
        //  resolution simulations
        part[ part[i].nn[k] ].mass = 0;
        part[ part[i].nn[j] ].type = 0; //  Mark down dead particles
      }
    }

  }
}
