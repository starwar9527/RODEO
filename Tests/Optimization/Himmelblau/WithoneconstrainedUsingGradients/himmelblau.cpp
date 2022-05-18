#include<stdio.h>
#include<math.h>
#include <unistd.h>

double HimmelblauAdj(double *x, double *xb) {
	double tempb;
	double tempb0;
	tempb = 2.0*pow(x[0]*x[0]+x[1]-11.0, 2.0-1);
	tempb0 = 2.0*pow(x[0]+x[1]*x[1]-7.0, 2.0-1);
	xb[0] = tempb0 + 2*x[0]*tempb;
	xb[1] = 2*x[1]*tempb0 + tempb;

	return pow( (x[0]*x[0]+x[1]-11.0), 2.0 ) + pow( (x[0]+x[1]*x[1]-7.0), 2.0 );

}

double constraint(double *x, double *xb1){
		return x[0]+x[1];
		xb1[0] = 1;
		xb1[1] = 1;
}


int main(void){

double x[2];
double xb[2];
double xb1[2];

FILE *inp = fopen("dv.dat","r");
fscanf(inp,"%lf",&x[0]);
fscanf(inp,"%lf",&x[1]);
fclose(inp);

double result = HimmelblauAdj(x, xb);
double constraintValue = constraint(x,xb1);

FILE *outp = fopen("objFunVal.dat","w");
fprintf(outp,"himmelblau_function = %15.10f\n",result);
fprintf(outp,"himmelblau_gradient = %15.10f, %15.10f\n",xb[0],xb[1]);
fprintf(outp,"constraint_function = %15.10f\n",constraintValue);
fprintf(outp,"constraint_gradient = %15.10f\n",xb1[0],xb1[1]);
fclose(outp);

return 0;
}

