/******************* Constant Declarations *********************************/
#define L 40                    /* lattice linear size*/
#define N (L*L)                 /* total number of spins */

#define transient 10000         /* transient time */
#define samples 10000           /* number of samples */
#define tau 10                  /* decorrelation time */

#define initialState 0          /* Initial state: Choose '0' for an ordered state or '1' for a disordered state */      
#define Ti 3.00                 /* Initial temperature */
#define Tf 4.00                 /* Final temperature */
#define dT 0.05                 /* Temperature increment */

#define algorithm 'S'           /* Use 'S' for single-flip algorithm or 'W' for Wolff algorithm */

#ifdef DEBUG
#define SEED 123456
#else
#define SEED time(NULL)         /* Seed of the random number generator */
#endif

