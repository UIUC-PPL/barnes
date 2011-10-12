/* adapted from barnes benchmark in SPLASH2 */

/*
 * TESTDATA: generate Plummer model initial conditions for test runs,
 * scaled to units such that M = -4E = G = 1 (Henon, Hegge, etc).
 * See Aarseth, SJ, Henon, M, & Wielen, R (1974) Astr & Ap, 37, 183.
 */

#include "Vector3D.h"
#include "Particle.h"

#define MFRAC  0.999                /* mass cut off at MFRAC of total */

Real xrand();
Real pow();

void testdata(int nbody)
{
   Real rsc, vsc, rsq, r, v, x, y;
   Vector3D<Real> cmr, cmv;
   Particle *p;
   int rejects = 0;
   int k;
   int halfnbody, i;
   Real offset;
   Particle *cp;
   Real tmp;

   Particle *bodytab = (bodyptr) G_MALLOC(nbody * sizeof(body));
   if (bodytab == NULL) {
      error("testdata: not enuf memory\n");
   }
   rsc = 9 * PI / 16;
   vsc = sqrt(1.0 / rsc);

   cmr = Vector3D<Real>(0.0);
   cmv = Vector3D<Real>(0.0);

   halfnbody = nbody / 2;
   if (nbody % 2 != 0) halfnbody++;
   for (p = bodytab; p < bodytab+halfnbody; p++) {
      p->mass = 1.0/nbody;
      r = 1 / sqrt(pow(xrand(0.0, MFRAC), -2.0/3.0) - 1);
      /*   reject radii greater than 10 */
      while (r > 9.0) {
	 rejects++;
	 r = 1 / sqrt(pow(xrand(0.0, MFRAC), -2.0/3.0) - 1);
      }        
      pickshell(p->position, rsc * r);
      cmr += p->position;
      do {
	 x = xrand(0.0, 1.0);
	 y = xrand(0.0, 0.1);

      } while (y > x*x * pow(1 - x*x, 3.5));

      v = sqrt(2.0) * x / pow(1 + r*r, 0.25);
      pickshell(p->velocity, vsc * v);
      cmv += p->velocity;
   }

   offset = 4.0;

   for (p = bodytab + halfnbody; p < bodytab+nbody; p++) {
      p->mass = 1.0 / nbody;
      cp = p - halfnbody;
      p->position = cp->position+offset;
      p->velocity = cp->velocity;
      cmr += p->position;
      cmv += p->position;
   }

   cmr /= nbody;
   cmv /= nbody;

   for (p = bodytab; p < bodytab+nbody; p++) {
     p->position -= cmr;
     p->velocity -= cmv;
   }
}

/*
 * PICKSHELL: pick a random point on a sphere of specified radius.
 */

void pickshell(Vector3D<Real> &vec, Real rad)
   //Real vec[];                     /* coordinate vector chosen */
   //Real rad;                       /* radius of chosen point */
{
   register int k;
   Real rsq, rsc;

   do {
     vec.x = xrand(-1.0,1.0);
     vec.y = xrand(-1.0,1.0);
     vec.z = xrand(-1.0,1.0);
     rsq = vec.lengthSquared();
   } while (rsq > 1.0);

   rsc = rad / sqrt(rsq);
   vec = rsc*vec;
}


int intpow(i,j)
  int i,j;
{   
    int k;
    int temp = 1;

    for (k = 0; k < j; k++)
        temp = temp*i;
    return temp;
}


