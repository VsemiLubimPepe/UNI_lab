#ifndef THRESHOLD_H

#define THRESHOLD_H

/*
 * Threshold for double type numbers close to zero. Number which absolute 
 * value is lower than THRESHOLD is considered to be zero. 
 */
extern double threshold;
/*
 * Threshold for elements of MTRX_A subdiagonal.
 */
extern double epsilon;

#endif