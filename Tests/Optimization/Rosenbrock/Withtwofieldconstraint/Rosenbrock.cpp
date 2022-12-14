#include<stdio.h>
#include<math.h>
#include <unistd.h>
#include <iostream>

#include <armadillo>
using namespace arma;
using namespace std;

double Rosenbrock10(double *x) {

	double temp = 0.0;

	for(int i=0; i<9; i++){

		temp+= 100.0*(x[i+1]-x[i]*x[i])*(x[i+1]-x[i]*x[i]) + (1-x[i])*(1-x[i]);
	}

	return temp;
}

vec constraint(double *x){          // g(x,t) = sum(x)*t; t =  [0,1]

	int m = 1000;
    double sum = 0;

	vec con = zeros(m);
	vec t = linspace(0,m-1,m);

	 for(int i=0; i<10; i++){
        sum = sum + x[i]*x[i];
	 }

	 sum = sqrt(sum);

	 vec output = zeros(m);

	 for(int i=0; i< m; i++){
		 output(i) = sum*t(i)/m + x[0]*x[1]*x[2];
     }

	 return output;
}


vec constraint1(double *x){          // g(x,t) = sum(x)*t; t =  [0,1]

	int m = 3000;
	double sum = 0;

	vec con = zeros(m);
    vec t = linspace(0,m-1,m);

	 for(int i=0; i<10; i++){
	   sum = sum + x[i]*x[i];
	 }

     sum = sqrt(sum);

	 vec output = zeros(m);

     for(int i=0; i< m; i++){
		 output(i) = sum*t(i)/m + x[0]*x[1]*x[2];
	 }

	return output;
}

int main(void){

  double x[10];

  FILE *inp = fopen("dv.dat","r");

  for (unsigned int k = 0; k < sizeof(x); k++ ){
		fscanf(inp,"%lf",&x[k]);
   }
  fclose(inp);


  double result = Rosenbrock10(x);
  vec con_result = constraint(x);
  vec con_result1 = constraint1(x);

  /*FILE *outp = fopen("objFunVal.dat","w");
  fprintf(outp,"Rosenbrock_function = %15.10f\n",result);
  fclose(outp);

  /* export constraint function value

   FILE *outp1 = fopen("conFunVal.dat","w");
   for (long unsigned int i=0; i< con_result.size(); i++){
        fprintf(outp1,"constraint_function = %15.10f\n",con_result[i]);
   }
   fclose(outp1);*/

   std::ofstream obj_value;
   obj_value.open("objFunVal.dat");
   obj_value << result << std::endl;
   obj_value.close();

   std::ofstream constraintValue;
   constraintValue.open("conFunVal.dat");
   for (long unsigned int i=0; i<con_result.size();i++){
	   constraintValue << con_result[i] << std::endl;
   }
   constraintValue.close();

   std::ofstream constraintValue1;
   constraintValue1.open("conFunVal1.dat");
   for (long unsigned int i=0; i<con_result1.size();i++){
   	   constraintValue1 << con_result1[i] << std::endl;
    }
    constraintValue1.close();

  return 0;
}

