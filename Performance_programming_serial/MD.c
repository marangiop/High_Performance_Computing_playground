/*
 *  Simple molecular dynamics code.
 *
 */
#include <stdio.h>
#include <math.h>
#include "coord.h"

void vis_force(int N,double *f, double *vis, double *vel);
void add_norm(int N,double *r, double *delta);
double force(double W, double delta, double r);
void wind_force(int N,double *f, double *vis, double vel);





void evolve(int count,double dt){
int  step;
int i,j,k,l;
int collided;
double Size;
double precalculated_force, precalculated_masses;
/*
 * Loop over timesteps.
 */
      for(step = 1;step<=count;step++){
        printf("timestep %d\n",step);
        printf("collisions %d\n",collisions);

/* set the viscosity term in the force calculation */
        for(j=0;j<Ndim;j++){   
          visc_force(Nbody,f[j],vis,velo[j]);
        }
/* add the wind term in the force calculation */
        for(j=0;j<Ndim;j++){
          wind_force(Nbody,f[j],vis,wind[j]);
        }
/* calculate distance from central mass */
       for(k=0;k<Nbody;k++){
          r[k] = 0.0;
        }
        for(i=0;i<Ndim;i++){
	  add_norm(Nbody,r,pos[i]);
        }
        for(k=0;k<Nbody;k++){
          r[k] = sqrt(r[k]);
        }
       /* calculate central force */
        for(i=0;i<Nbody;i++){
	  for(l=0;l<Ndim;l++){
                f[l][i] = f[l][i] - 
                   force(G*mass[i]*M_central,pos[l][i],r[i]);
	  }
	}

/* calculate pairwise separation of particles */
/* calculate norm of separation vector */

        for(i=0;i<Nbody;i++){
            for(j=i+1;j<Nbody;j++){
                delta_r[i][j] = 0.0;
                for(l=0;l<Ndim;l++){
                    delta_r[i][j] += ((pos[l][i]-pos[l][j]) * (pos[l][i]-pos[l][j]));
                    }
                    delta_r[i][j] = sqrt(delta_r[i][j]);
            }
        }

/*
 * add pairwise forces.
 */
        for(i=0;i<Nbody;i++){
          for(j=i+1;j<Nbody;j++){
            Size = radius[i] + radius[j];
            collided=0;
            precalculated_masses = G*mass[i]*mass[j];
/*  flip force if close in */
              if( delta_r[i][j] >= Size ){

                  for(l=0;l<Ndim;l++){

                      precalculated_force= force(precalculated_masses,(pos[l][i] -  pos[l][j]),delta_r[i][j]);

                      f[l][i] = f[l][i] - precalculated_force;
                      f[l][j] = f[l][j] + precalculated_force;
                  }
              }
              else{
                  for(l=0;l<Ndim;l++){
                      precalculated_force= force(precalculated_masses,(pos[l][i] -  pos[l][j]),delta_r[i][j]);

                      f[l][i] = f[l][i] +  precalculated_force;
                      f[l][j] = f[l][j] -  precalculated_force;
                      collided=1;
                  }
              }
                if( collided == 1 ){
                    collisions++;
                }
            }
        }

/* update positions */
        for(i=0;i<Nbody;i++){
          for(j=0;j<Ndim;j++){
            pos[j][i] = pos[j][i] + dt * velo[j][i];
          }
        }

/* update velocities */
        for(i=0;i<Nbody;i++){
          for(j=0;j<Ndim;j++){
            velo[j][i] = velo[j][i] + dt * (f[j][i]/mass[i]);
          }
        }

      }
}




