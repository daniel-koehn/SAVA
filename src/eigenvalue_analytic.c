/*------------------------------------------------------------------------
 *   Determine eigenvalues of a real symmetric 3x3 matrix
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"


void eigenvalue(float * s, float * e){

	/* s ... elements of symmetric 3x3 matrix S as vector [s11,s22,s33,s23,s13,s12]	(in)
			S=[s11 s12 s13;s12 s22 s23;s13 s23 s33]
	   e ... eigenvalues of 3x3 matrix [l1,l2,l3]					(out) */


	/* local variables */
	int	i, j, l;
//	float	sphi, cphi, tphi;
	float	a, b, c, d, p, q, r, t, dis;
	float	syz2, sxz2, sxy2;


	/* compute eigenvalues of symmetric (3x3) matrix */
	syz2 = s[4]*s[4];
	sxz2 = s[5]*s[5];
	sxy2 = s[6]*s[6];
	a    = syz2 + sxz2 + sxy2;
	
	if (a == 0.0){
		/* matrix is diagonal */
		e[1] = s[1];
		e[2] = s[2];
		e[3] = s[3];
	}
	else{
		/* matrix is non-diagonal */
		b = (s[1] + s[2] + s[3])/3.0;
		c =  s[2]*s[3] + s[1]*s[3] + s[1]*s[2] - a;
		d = -(s[1]*(s[2]*s[3]-syz2) + 2.0*s[4]*s[5]*s[6] - s[2]*sxz2 - s[3]*sxy2);

		p = b*b - c/3.0;
		q = b*b*b - 0.5*(b*c + d);

		dis = q*q - p*p*p;

		if (dis < 0.0){
			r = q/sqrt(p*p*p);
			if (r >= 1.0){
				t = 0.0;
			}
			else{
				if (r <= -1.0){
					t = PI/3.0;
				}
				else{
					t = acos(r)/3.0;
				}
			}
			q = 2.0*sqrt(p);
			e[1] = b + q*cos(t);
			e[2] = b - q*cos(t + PI/3.0);
//			e[3] = b - q*cos(t - PI/3.0);
			e[3] = 3.0*b - e[1] - e[2];
		}
		else{
			if (p != 0.0){
				e[1] = b + 2.0*q/p;
				e[2] = b - q/p;
				e[3] = e[2];
			}
			else{
				e[1] = b;
				e[2] = b;
				e[3] = b;
			}
		}
	}

	/* sort eigenvalues */
	for (i=1;i<3;i++){
		l = i;
		for (j=i+1;j<=3;j++){
			if (e[j] > e[l]){
				l = j;
			}
		}
		if (i != l){
			a    = e[l];
			e[l] = e[i];
			e[i] = a;
		}
	}

}
