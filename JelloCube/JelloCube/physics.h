/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#ifndef _PHYSICS_H_
#define _PHYSICS_H_

double Vector_Dot_Product(struct point a, struct point b);
void computeAcceleration(struct world * jello, struct point a[8][8][8]);
struct point structForce(struct world * jello, int p_x, int p_y, int p_z);
struct point shearForce(struct world * jello, int p_x, int p_y, int p_z);
struct point bendForce(struct world * jello, int p_x, int p_y, int p_z);
struct point externalForce(struct world * jello, int p_x, int p_y, int p_z);
struct point collideForce(struct world * jello, int p_x, int p_y, int p_z);

// perform one step of Euler and Runge-Kutta-4th-order integrators
// updates the jello structure accordingly
void Euler(struct world * jello);
void RK4(struct world * jello);

#endif

