/*
 * This file defines static arrays that contains the primary coordinates
 * of the particles,
 *
 *  Nbody	Number of particles
 *  Npair	Number of particle pairs
 *  pos		Position of the particles
 *  r           distance of partice from central mass 
 * vel		velocity of the particles
 *  f		Forces acting on each particle
 *  visc		viscosity coefficient for each particle
 *  mass		mass of each particle
 *  delta_pos	seperation vector for each particle pair
 *  delta_r		seperation for each particle pair
 */

#ifdef DECL
#define DEF
#else
#define DEF extern
#endif
#define Nbody 4*1024
#define  Npair ((Nbody*(Nbody-1))/2)

enum{ Xcoord=0, Ycoord, Zcoord, Ndim };

DEF double pos[Ndim][Nbody] __attribute__((aligned(64)));
DEF double velo[Ndim][Nbody] __attribute__((aligned(64)));
DEF double f[Ndim][Nbody] __attribute__((aligned(64)));
DEF double vis[Nbody] __attribute__((aligned(64)));
DEF double mass[Nbody] __attribute__((aligned(64)));
DEF double radius[Nbody] __attribute__((aligned(64)));
DEF double r[Nbody] __attribute__((aligned(64)));
DEF double delta_r[Nbody][Nbody] __attribute__((aligned(64)));
DEF double wind[Ndim] __attribute__((aligned(64)));
DEF int collisions;


#define G 2.0
#define M_central 1000.0

void evolve(int Nstep, double dt);


