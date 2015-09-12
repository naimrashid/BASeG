#include <stdio.h>
#include <stddef.h>
#include <time.h>
#include <string.h>
#include <R.h>
#include <Rmath.h>
#include <math.h>
#include "glm.h"
#include "utility.h"
#include "ase.h"
#include "asa121.h"
#include "asa103.h"
//#include "lbfgsb1.h"

//used for maximization of posterior mode of biv random effect via lbfgsb
double l1(int n, double* para, void* ex, SEXP x1){
			int i, k, N, P, h, h2, l;
		double *exPara, *b1, *b2, *mean0, *sigma0, *x_1, *x_2, *t,p, y1, y2, val;
		double l1, l2;
		
		exPara = (double *) ex;
		N      = exPara[0];
		P      = exPara[1];
		p    	 = exPara[2];
		y1   	 = exPara[3];
		y2   	 = exPara[4];
		h 	 	 = exPara[5];
		h2 	 	 = exPara[6];
		b1     = exPara + 7;
		b2     = b1 + P;
		mean0  = b2 + P; 
		sigma0 = mean0 + 2; 
		x_1    = sigma0 + 2;
		x_2    = x_1 +P;
		t      = para;

		//for(l=0;l<11+4*P;l++) Rprintf("h %f ", exPara[l]);


		if(h ==0 )  p = 0.0;
		if(h2 ==0)  b1[P-1] = b2[P-1] = 0.0;

		l1 = l2 = 0;
		for(l = 0; l<P; l++){
			l1 = l1 + x_1[l]*b1[l];
			l2 = l2 + x_2[l]*b2[l];
		}  
		l1 = exp(l1 + t[0]);
		l2 = exp(l2 + t[1]);

		val = dpois(y1, l1, 0)*dpois(y2, l2, 0)*dnorm(t[0], mean0[0] + p*sqrt(sigma0[0])*(t[1]- mean0[1])/sqrt(sigma0[1]), sqrt(sigma0[0]*(1-p*p)),0)*dnorm(t[1], mean0[1], sqrt(sigma0[1]),0);
		//Rprintf("val that %f %f\n", val, -log(val + pow(10,-100)));
 	//Rprintf("l1 val %e l1 %f l2 %f t1 %f t2 %f %f %f %f %f %f %e\n", -log(val + pow(10,-100)), l1 , l2, t[0], t[1], mean0[0], mean0[1],sigma0[0], sigma0[1],p, dnorm(t[0], mean0[0] + p*sqrt(sigma0[0])*(t[1]- mean0[1])/sqrt(sigma0[1]), sqrt(sigma0[0]*(1-p*p)),0)*dnorm(t[1], mean0[1], sqrt(sigma0[1]),0));
	if(val!=val){
		return(1000000);	
	}else{
		return(-log(val + pow(10,-100)));
	}
}

void gr(int n, double* para, double* gr0, void* ex, SEXP x1)
{
  	int i, k, N, P, h, h2, l;
		double *exPara, *b1, *b2, *mean0, *sigma0, *x_1, *x_2, *t,p, y1, y2;
		double l1, l2, dt1, dt2;
		
		exPara = (double *) ex;
		N      = exPara[0];
		P      = exPara[1];
		p    	 = exPara[2];
		y1   	 = exPara[3];
		y2   	 = exPara[4];
		h 	 	 = exPara[5];
		h2 	 	 = exPara[6];
		b1     = exPara + 7;
		b2     = b1 + P;
		mean0  = b2 + P; 
		sigma0 = mean0 + 2; 
		x_1    = sigma0 + 2;
		x_2    = x_1 +P;
		t      = para;
		
		if(h ==0 )  p = 0.0;
		if(h2 ==0)  b1[P-1] = b2[P-1] = 0.0;

		l1 = l2 = 0;
		for(l = 0; l<P; l++){
			l1 = l1 + x_1[l]*b1[l];
			l2 = l2 + x_2[l]*b2[l];
		}  
		l1 = exp(l1 + t[0]);
		l2 = exp(l2 + t[1]);

  	dt1 = -l1 + y1 - (1/(2*(1-p*p))) * (2*(t[0]-mean0[0])/sigma0[0] - 2*p*(t[1] - mean0[1])/sqrt(sigma0[0]*sigma0[1]));
		dt2 = -l2 + y2 - (1/(2*(1-p*p))) * (2*(t[1]-mean0[1])/sigma0[1] - 2*p*(t[0] - mean0[0])/sqrt(sigma0[0]*sigma0[1])); 
	
  	gr0[0] = -dt1;
		gr0[1] = -dt2;
		//Rprintf("gr l1 %f l2 %f t1 %f t2 %f dtq %f dt2 %f\n",l1 , l2, t[0], t[1], dt1, dt2);
		//error("");
}


void hess(double* t, double *xx, double y1, double y2, double *x_1, double *x_2, int h, int h2, int P, double *result){

	int l;
	double *b1, *b2, *sigma0, p, l1, l2, ddt1, ddt2, ddt12;
	if(h==0) xx[2*P+2] = 0;
	if(h2==0) xx[P-1] = xx[2*P-1] = 0; 

	b1 = xx;
	b2 = xx + P;
	sigma0 = b2 + P;  
	p  = xx[2*P+2];

	l1 = l2 = 0;
	for(l = 0; l<P; l++){
		l1 = l1 + x_1[l]*b1[l];
		l2 = l2 + x_2[l]*b2[l];
	}  
	l1 = exp(l1 + t[0]);
	l2 = exp(l2 + t[1]);
	
  ddt1 = l1 + 1/((1-p*p)*sigma0[0]);
	ddt2 = l2 + 1/((1-p*p)*sigma0[1]);
	ddt12 = -p/((1-p*p)*sqrt(sigma0[0]*sigma0[1]));
	//return vector representing 2x2 matrix with columns stacked
	result[0] = ddt1;
	result[1] = ddt12;
	result[2] = ddt12;
	result[3] = ddt2;	
} 

double adaQuad(double* xx, double* u, double* w, double* t,double* Qhat,double y1, double y2, double* x_1, double* x_2,int* dims){

	int l, j, k;
	double *b1, *b2, *sigma0, p, l1, l2, l01,l02,val, wj, wk, zstar1, zstar2;
	double Qhat2[4];
	double mean0[2];

	int N = dims[0];
	int h = dims[1];
	int h2 = dims[2];
	int P = dims[3];
	int nhrule = dims[4];

	if(h==0) xx[2*P+2] = 0;
	if(h2==0) xx[P-1] = xx[2*P-1] = 0; 

	b1 = xx;
	b2 = xx + P;
	sigma0 = b2 + P;
	p  = xx[2*P+2];

	mean0[0] = 0;///m -sigma0[0]/2;
	mean0[1] = 0;///m-sigma0[1]/2;

	l1 = l2 = 0;
	for(l = 0; l<P; l++){
		l1 = l1 + x_1[l]*b1[l];
		l2 = l2 + x_2[l]*b2[l];
	}  
	//Rprintf("adaquad l1 %e, l2 %e\n", l1, l2);

	//computes the expression zstar = that + sqrt(2)Qhat^1/2z for single observation
	chol(Qhat, Qhat2);
	
	//for(l = 0; l<4;l++) Rprintf("adaquad  qhat2 %e \n", Qhat2[l]);

	//Rprintf("%f %f %f %f %f %f\n", Qhat2[0],Qhat2[1],Qhat2[2],Qhat2[3], det(Qhat),mean0[1]);
	val=0.0;
	for(j=0; j<nhrule; j++){
		for(k=0; k<nhrule; k++){
			wj = w[j];
			wk = w[k];
			zstar1 = t[0] + sqrt(2)*(Qhat2[0]*u[j] + Qhat2[2]*u[k]);
			zstar2 = t[1] + sqrt(2)*(Qhat2[1]*u[j] + Qhat2[3]*u[k]);
			l01 = exp(l1 + zstar1);
			l02 = exp(l2 + zstar2);
			//Rprintf("%f %f %f %f %f %f %e\n", zstar1, zstar2, l1, l2, l01, l02,  2*sqrt(det(Qhat))*wj*wk*dpois(y1, l01, 0)*dpois(y2, l02, 0)*dnorm(zstar1, mean0[0] + p*sqrt(sigma0[0])*(zstar2- mean0[1])/sqrt(sigma0[1]), sqrt(sigma0[0]*(1-p*p)),0)*dnorm(zstar2, mean0[1], sqrt(sigma0[1]),0)*exp(u[j]*u[j])*exp(u[k]*u[k]));
			val = val + 2*sqrt(det(Qhat))*wj*wk*dpois(y1, l01, 0)*dpois(y2, l02, 0)*dnorm(zstar1, mean0[0] + p*sqrt(sigma0[0])*(zstar2- mean0[1])/sqrt(sigma0[1]), sqrt(sigma0[0]*(1-p*p)),0)*dnorm(zstar2, mean0[1], sqrt(sigma0[1]),0)*exp(u[j]*u[j])*exp(u[k]*u[k]);
	//Rprintf("2*det %e ww %e pois1 %e, pois2 %e, mvnorm %e, exps %e, val %e, cumsum val %e\n", 2*sqrt(det(Qhat)), wj*wk,dpois(y1, l01, 0), dpois(y2, l02, 0), dnorm(zstar1, mean0[0] + p*sqrt(sigma0[0])*(zstar2- mean0[1])/sqrt(sigma0[1]), sqrt(sigma0[0]*(1-p*p)),0)*dnorm(zstar2, mean0[1], sqrt(sigma0[1]),0),exp(u[j]*u[j])*exp(u[k]*u[k]),2*sqrt(det(Qhat))*wj*wk*dpois(y1, l01, 0)*dpois(y2, l02, 0)*dnorm(zstar1, mean0[0] + p*sqrt(sigma0[0])*(zstar2- mean0[1])/sqrt(sigma0[1]), sqrt(sigma0[0]*(1-p*p)),0)*dnorm(zstar2, mean0[1], sqrt(sigma0[1]),0)*exp(u[j]*u[j])*exp(u[k]*u[k]), val); 

		}		
	}
	//Rprintf("%e final val \n",val);
	return(-log(val+pow(10,-100)));
}

void adaQuadgr(double* xx, double* u, double* w, double* t,double* Qhat,double y1, double y2, double* x_1, double* x_2,int* dims, double* gr, int *trace){

	int l, j, k;
	double *b1, *b2, *sigma0, p, l1, l2, l01,l02,val, wj, wk, zstar1, zstar2,L;
	double Qhat2[4];
	double mean0[2];

	int N = dims[0];
	int h = dims[1];
	int h2 = dims[2];
	int P = dims[3];
	int nhrule = dims[4];

	if(h==0) xx[2*P+2] = 0;
	if(h2==0) xx[P-1] = xx[2*P-1] = 0; 

	b1 = xx;
	b2 = xx + P;
	sigma0 = b2 + P;
	p  = xx[2*P+2];

	mean0[0] = 0;///m-sigma0[0]/2;
	mean0[1] = 0;///m-sigma0[1]/2;

	l1 = l2 = 0;
	for(l = 0; l<P; l++){
		l1 = l1 + x_1[l]*b1[l];
		l2 = l2 + x_2[l]*b2[l];
	}  
	//if(*trace > 1) Rprintf("gr l1 %f l2 %f\n",l1,l2);   
	//computes the expression zstar = that + sqrt(2)Qhat^1/2z for single observation
	chol(Qhat, Qhat2);
	val=0.0;
	for(j=0; j<nhrule; j++){
		for(k=0; k<nhrule; k++){
			wj = w[j];
			wk = w[k];
			zstar1 = t[0] + sqrt(2)*(Qhat2[0]*u[j] + Qhat2[2]*u[k]);
			zstar2 = t[1] + sqrt(2)*(Qhat2[1]*u[j] + Qhat2[3]*u[k]);
			l01 = exp(l1 + zstar1);
			l02 = exp(l2 + zstar2);
			//if(*trace>2) Rprintf("%f %f %f %f %f %f\n", zstar1, zstar2, l1, l2, l01, l02);
			val = val + 2*sqrt(det(Qhat))*wj*wk*dpois(y1, l01, 0)*dpois(y2, l02, 0)*dnorm(zstar1, mean0[0] + p*sqrt(sigma0[0])*(zstar2- mean0[1])/sqrt(sigma0[1]), sqrt(sigma0[0]*(1-p*p)),0)*dnorm(zstar2, mean0[1], sqrt(sigma0[1]),0)*exp(u[j]*u[j])*exp(u[k]*u[k]);
		}		
	}

	val = val + pow(10,-100);
	for(l=0;l<2*P+3;l++){ 
		gr[l]		=	0.0;
	}

	for(j=0; j<nhrule; j++){
		for(k=0; k<nhrule; k++){
			wj = w[j];
			wk = w[k];
			zstar1 = t[0] + sqrt(2)*(Qhat2[0]*u[j] + Qhat2[2]*u[k]);
			zstar2 = t[1] + sqrt(2)*(Qhat2[1]*u[j] + Qhat2[3]*u[k]);
			l01 = exp(l1 + zstar1);
			l02 = exp(l2 + zstar2);
			L = 2*sqrt(det(Qhat))*wj*wk*dpois(y1, l01, 0)*dpois(y2, l02, 0)*dnorm(zstar1, mean0[0] + p*sqrt(sigma0[0])*(zstar2- mean0[1])/sqrt(sigma0[1]), sqrt(sigma0[0]*(1-p*p)),0)*dnorm(zstar2, mean0[1], sqrt(sigma0[1]),0)*exp(u[j]*u[j])*exp(u[k]*u[k]);
			//Rprintf("%f %f %f %f %f %f\n", zstar1, zstar2, l1, l2, l01, l02);
			for(l=0;l<P-1+h2;l++){ 
				gr[l]		=	gr[l]   + x_1[l]*(y1 - l01)*L/val;
				gr[l+P]	= gr[l+P] + x_2[l]*(y2 - l02)*L/val;
			}

			///m gr[2*P] = gr[2*P] + (-1/(2*sigma0[0]) -(1/(2*(1-p*p)))*( -pow(zstar1+sigma0[0]/2,2)/pow(sigma0[0],2) + (zstar1+sigma0[0]/2)/sigma0[0] + p*(zstar1+sigma0[0]/2)*(zstar2+sigma0[1]/2)*sigma0[1]/pow(sigma0[0]*sigma0[1],1.5) - p*( zstar2+sigma0[1]/2)/sqrt(sigma0[0]*sigma0[1])))*L/val;
			///m gr[2*P+1] = gr[2*P+1] + (-1/(2*sigma0[1]) -(1/(2*(1-p*p)))*( -pow(zstar2+sigma0[1]/2,2)/pow(sigma0[1],2) + (zstar2+sigma0[1]/2)/sigma0[1] + p*(zstar1+sigma0[0]/2)*(zstar2+sigma0[1]/2)*sigma0[0]/pow(sigma0[0]*sigma0[1],1.5) - p*( zstar1+sigma0[0]/2)/sqrt(sigma0[0]*sigma0[1])))*L/val;
			///m if(h==1) gr[2*P+2] = gr[2*P+2] + ((p/(1-p*p) +(zstar1+sigma0[0]/2)*(zstar2+sigma0[1]/2)/(sqrt(sigma0[0]*sigma0[1])*(1-p*p)) - p*(pow(zstar1+sigma0[0]/2,2)/sigma0[0] + pow(zstar2+sigma0[1]/2,2)/sigma0[1] - 2*p*(zstar1+sigma0[0]/2)*(zstar2+sigma0[1]/2)/sqrt(sigma0[0]*sigma0[1]))/pow(1-p*p,2))*L)/val;

			gr[2*P] = gr[2*P] + (-1/(2*sigma0[0]) -(1/(2*(1-p*p)))*( -pow(zstar1-mean0[0],2)/pow(sigma0[0],2) + p*(zstar1-mean0[0])*(zstar2-mean0[1])*sigma0[1]/pow(sigma0[0]*sigma0[1],1.5)))*L/val;
			gr[2*P+1] = gr[2*P+1] + (-1/(2*sigma0[1]) -(1/(2*(1-p*p)))*( -pow(zstar2-mean0[1],2)/pow(sigma0[1],2) + p*(zstar1-mean0[0])*(zstar2-mean0[1])*sigma0[0]/pow(sigma0[0]*sigma0[1],1.5)))*L/val;

			if(h==1) gr[2*P+2] = gr[2*P+2] + ((p/(1-p*p) +(zstar1-mean0[0])*(zstar2-mean0[1])/(sqrt(sigma0[0]*sigma0[1])*(1-p*p)) - p*(pow(zstar1-mean0[0],2)/sigma0[0] + pow(zstar2-mean0[1],2)/sigma0[1] - 2*p*(zstar1-mean0[0])*(zstar2-mean0[1])/sqrt(sigma0[0]*sigma0[1]))/pow(1-p*p,2))*L)/val;

//		if(*trace > 3){
//			for(l=0;l<2*P+3;l++){ 
//				Rprintf("%f\n", gr[l]);
//			}
//		}
		}		
	}
	
	for(l=0;l<2*P+3;l++){ 
		gr[l]		=	-gr[l];
	}

}


void adaQuad_vec(double *xx, double *y1, double *y2, double *x1, double *x2, double *u, double *w, int *dims, int *trace, double* sum){
	

	int i,l;
	int N = dims[0];
	int h = dims[1];
	int h2 = dims[2];
	int P = dims[3];
	int nhrule = dims[4];
	*sum=0.0;
	double that[2];
	//		Rprintf("ok3\n");
	double *exPara, *b1, *b2, *mean0, *sigma0, *x_1, *x_2, hessian[4], Qhat[4],p;
	exPara = (double *) Calloc(7+4*P+4, double); 


	//Rprintf("ada starting coef\n");
	//Rprint_v(xx, 0, 2*P+2);

	//be careful of conversion from ints to doubles
	exPara[0] = N;
	exPara[1] = P;
	exPara[2] = xx[2*P+2]; //p;
	exPara[3] = 0.0; //yi1 is updated in the loop below with xi1
	exPara[4] = 0.0; //same for second one
	exPara[5] = h;
	exPara[6] = h2;
	b1     = exPara + 7;
	b2     = b1 + P;
	mean0  = b2 + P; 
	sigma0 = mean0 + 2; 
	x_1    = sigma0 + 2;
	x_2    = x_1 +P;
	
  int npara   = 2; 
  int lmm     = 5; 
  int fail    = 0;
  int failA   = 0;
  int failB   = 0;
  int fncount = 0;//?
  int grcount = 0;//?
  int maxit   = 500;
  int nREPORT = 5;
			
	int nbd[npara];

	//technical parameters below:
		  double *wa, *g1;
		  int *iwa;
		  SEXP xx1;
		  PROTECT(xx1 = allocVector(REALSXP, npara));
		  //consider replacing with simple Calloc
		  //wa  = (double *) S_alloc(2*lmm*npara+4*npara+11*lmm*lmm+8*lmm,sizeof(double));
		  //iwa = (int*) R_alloc(3*npara,sizeof(int));
		  //g1 = (double *)R_alloc(npara, sizeof(double));
		  wa = (double *) Calloc(2*lmm*npara+4*npara+11*lmm*lmm+8*lmm, double);
                  iwa = (int*) Calloc(3*npara, int);
                  g1 = (double *) Calloc(npara, double);

		  double gr0[npara];
		  double initPara[npara];
		  double lower[npara];
		  double upper[npara];
		  double Fmin, factr, pgtol, th0, th1, th01;
		  
		  factr = 1e7;
			pgtol = 0.0;
		  
			char msg[1023];
		  
		  for(l=0; l < npara; l++){
					nbd[l] = 2;			  	
					initPara[l] = 0;
					lower[l] = -10;
					upper[l] = 10; 
			}
			//Rprintf("ok4\n");
	//set b1 b2 mean0 sigma to exPara	
	for(l=0;l<P;l++){
		b1[l] = xx[l];
		b2[l] = xx[l+P];
	}
	for(l=0; l<2;l++){
		mean0[l] = 0;///m -xx[2*P+l]/2;
		sigma0[l] = xx[2*P+l];
	}

	//loop over each observation
	*sum=0.0;
			//Rprintf("ok5\n");
	for(i=0;i<N;i++){
			// update xi and yi
			for(l=0;l<P;l++){
				x_1[l] = x1[l*N+i];
				x_2[l] = x2[l*N+i];
				//Rprintf("x1 %f x2 %f \n",x_1[l],x_2[l]); 
			}
			exPara[3] = y1[i];
			exPara[4] = y2[i];

		//for(l=0;l<11+4*P;l++) Rprintf("h %f", exPara[l]);
					//	Rprintf("ok6 %d\n", i);
			//l1(npara, initPara, (void*)exPara, xx1);
			//gr(npara, initPara, gr0, (void*)exPara, xx1);
      lbfgsb1(npara, lmm, initPara, lower, upper, nbd, &Fmin, 
           l1, gr, &failA, (void*)exPara, factr, pgtol,  
           &fncount, &grcount, maxit, msg, 0, nREPORT, wa, iwa, g1,xx1);
	    //if(*trace > 1) Rprintf("%d %s\n", failA, msg);
	    //if (failA) {
	      //if (*trace>0)
	      //  Rprintf("  i=%d, fail to fit t0 model in iteration\n", i);
	    
	      //continue;
	    //}else{
				//save final estimates of that
				that[0] = initPara[0];
				that[1] = initPara[1];
				//Rprint_v(that, 0, 1);
				//reset initPara to 0
				initPara[0]=initPara[1] = 0.0;
		//	}
			//calculate Qhat
						//Rprintf("ok6b %d\n", i);
			hess(that,xx,y1[i],y2[i],x_1,x_2, h, h2, P, hessian);
						//Rprintf("ok6c %d\n", i);
			inv(hessian,Qhat);
			//for(l=0;l<4;l++){
			//	Rprintf("hess %e Qhat %e\n", hessian[l], Qhat[l]);
			//}
			//calculate value and update sum
						//Rprintf("ok7 %d\n", i);
			*sum = *sum + adaQuad(xx, u, w, that,Qhat,y1[i], y2[i], x_1, x_2,dims);
			//nif(*trace > 0) Rprintf("%f %d\n",*sum,i);
	}
	//Rprintf("sum ada: %f\n", *sum);
  Free(exPara);
  Free(wa);
  Free(iwa);
  Free(g1);
  UNPROTECT(1);
}

void adaQuadgr_vec(double *xx, double *y1, double *y2, double *x1, double *x2, double *u, double *w, int *dims, int *trace, double* gr0){
	
	int i,l;
	int N = dims[0];
	int h = dims[1];
	int h2 = dims[2];
	int P = dims[3];
	int nhrule = dims[4];
	double that[2];

	double *exPara, *b1, *b2, *mean0, *sigma0, *x_1, *x_2, hessian[4], Qhat[4],p,grtemp[2*P+3];
	exPara = (double *) Calloc(7+4*P+4, double); 


	//printf("gr starting coef\n");
	//Rprint_v(xx, 0, 2*P+2);

	//be careful of conversion from ints to doubles
	exPara[0] = N;
	exPara[1] = P;
	exPara[2] = xx[2*P+2]; //p;
	exPara[3] = 0.0; //yi1 is updated in the loop below with xi1
	exPara[4] = 0.0; //same for second one
	exPara[5] = h;
	exPara[6] = h2;
	b1     = exPara + 7;
	b2     = b1 + P;
	mean0  = b2 + P; 
	sigma0 = mean0 + 2; 
	x_1    = sigma0 + 2;
	x_2    = x_1 +P;
	
	//Rprintf("gr starting expra\n");
	//Rprint_v(exPara, 0, 7);

  int npara   = 2; 
  int lmm     = 5; 
  int fail    = 0;
  int failA   = 0;
  int failB   = 0;
  int fncount = 0;//?
  int grcount = 0;//?
  int maxit   = 500;
  int nREPORT = 5;
			
	int nbd[npara];

	//technical parameters below:
		  double *wa, *g1;
		  int *iwa;
		  SEXP xx1;
		  PROTECT(xx1 = allocVector(REALSXP, npara));
		  //consider replacing with simple Calloc
		  //wa  = (double *) S_alloc(2*lmm*npara+4*npara+11*lmm*lmm+8*lmm,sizeof(double));
		  //iwa = (int*) R_alloc(3*npara,sizeof(int));
		  //g1 = (double *)R_alloc(npara, sizeof(double));
		
		  wa = (double *) Calloc(2*lmm*npara+4*npara+11*lmm*lmm+8*lmm, double);
		  iwa = (int*) Calloc(3*npara, int);
		  g1 = (double *) Calloc(npara, double);

		  //double gr0[npara];
		  double initPara[npara];
		  double lower[npara];
		  double upper[npara];
		  double Fmin, factr, pgtol, th0, th1, th01;
		  
		  factr = 1e7;
			pgtol = 0.0;
		  
			char msg[1023];
		  
		  for(l=0; l < npara; l++){
					nbd[l] = 2;			  	
					initPara[l] = 0;
					lower[l] = -10;
					upper[l] = 10; 
			}

	//set b1 b2 mean0 sigma to exPara	
	for(l=0;l<P;l++){
		b1[l] = xx[l];
		b2[l] = xx[l+P];
	}
	for(l=0; l<2;l++){
		mean0[l] = 0;///m -xx[2*P+l]/2;
		sigma0[l] = xx[2*P+l];
	}

	//loop over each observation
	for(i=0;i<N;i++){
			// update xi and yi
			for(l=0;l<P;l++){
				x_1[l] = x1[l*N+i];
				x_2[l] = x2[l*N+i];
			}
			exPara[3] = y1[i];
			exPara[4] = y2[i];

		//for(l=0;l<11+4*P;l++) Rprintf("h %f", exPara[l]);
			//get estimate of that
      lbfgsb1(npara, lmm, initPara, lower, upper, nbd, &Fmin, 
           l1, gr, &failA, (void*)exPara, factr, pgtol,  
           &fncount, &grcount, maxit, msg, 0, nREPORT, wa, iwa, g1,xx1);
	    //Rprintf("%d %s\n", failA, msg);
	    //if (failA) {
	      //if (*trace)
	      //  Rprintf("  i=%d, fail to fit t0 model in iteration\n", i);
	    
	      //continue;
	    //}else{
				//save final estimates of that
				that[0] = initPara[0];
				that[1] = initPara[1];
				//reset initPara to 0
				initPara[0]=initPara[1] = 0.0;
		//	}
			//calculate Qhat
			hess(that,xx,y1[i],y2[i],x_1,x_2, h, h2, P, hessian);
			inv(hessian,Qhat);
			//calculate value and update sum
			adaQuadgr(xx, u, w, that,Qhat,y1[i], y2[i], x_1, x_2,dims, grtemp, trace); 
			for(l=0;l<2*P+3;l++){
				gr0[l] = gr0[l] + grtemp[l];
			}
			//Rprintf("%f\n",*sum);
			//error("");
	}
	//Rprint_v(gr0, 0, 2*P+2);
  Free(exPara);
  Free(wa);
  Free(iwa);
  Free(g1);
  UNPROTECT(1);
}

// Wrapper for C-based bfgs
double adaQuad_vec_bfgs(int n, double* x, void* ex, SEXP x1){

	int i, j, k, l, N, nX, nX2,dims[5];	
	double quad, rho, sum=0.0;
	double mu1 = 0.0, mu2 = 0.0;
	double neglogquad = 0.0;
	double *exPara, *y1, *y2, *xx1, *xx2, *pXlast, *pXlast2, *ddimsnew, *u, *w, *trace, *ll;

	/*exPara[0] = N;
	exPara[1] = nX;
	exPara[2] = nX2;
  //YY     = exPara ; filled in later for each Y
  //YY2    = YY + N; filled in later for each Y2
  XX     = exPara + 2*N+3;
	pXXlast = XX + nX*N;
  XX2     = pXXlast + N; //add extra column for snp effect
	pXXlast2 = XX2 + nX2*N;
	ddimsNew = pXXlast2 + N; //add extra column for snp effect
	uu = ddimsNew + 10;
	ww = uu + pts;
	*/

  exPara = (double *) ex;
	N = (int) exPara[0]; //ceil(exPara[0]-0.5);
	nX = (int) exPara[1];//ceil(exPara[1]-0.5);
	nX2 = (int) exPara[2];//ceil(exPara[2]-0.5);
	y1  = exPara + 3;
	y2 = y1 + N; 
  xx1 = y2 + N;
	pXlast = xx1 + (nX-1)*N;
  xx2     = pXlast + N; //add extra column for snp effect
	pXlast2 = xx2 + (nX2-1)*N;
	ddimsnew = pXlast2 + N; 
	u = ddimsnew + 10;
	w = u + (int) ddimsnew[6];

	/*int i,l;
	int N = dims[0];
	int h = dims[1];
	int h2 = dims[2];
	int P = dims[3];
	int nhrule = dims[4];
	double that[2];
	*/
	dims[0] = (int) N;
	dims[1] = (int) ddimsnew[7];
	dims[2] = (int) ddimsnew[8];
	dims[3] = nX; //assume nX is the same as nX2, account for int and snp column
	dims[4] = (int) ddimsnew[6];

	//Rprint_v(x, 0, 5);
	//Rprint_v(y1, 0, 5);
	//Rprint_v(y2, 0, 5);
	//Rprint_v(xx1, 0, 5);
	//Rprint_v(xx2, 0, 5);
	//for(l=0;l<nX;l++){
	//	Rprintf("XX1 %f XX2 %f \n",xx1[N*l],xx2[N*l]); 
	//}

	//Rprint_v(u, 0, 4);
	//Rprint_v(w, 0, 4);
	//Rprint_vi(dims, 0, 4);


	sum = 0.0;
	adaQuad_vec(x, y1, y2, xx1, xx2, u, w, dims, 0, &sum);
	//Rprintf("sum %f\n", sum);
	if(sum!=sum){
		return(1000000);
	}else{
		return(sum);
	}
}

void adaQuadgr_vec_bfgs(int n, double* x, double* gr0, void* ex, SEXP x1){

		int i, j, k, l, N, nX, nX2,dims[5];	
	double quad, rho, sum=0.0;
	double mu1 = 0.0, mu2 = 0.0;
	double neglogquad = 0.0;
	double *exPara, *y1, *y2, *xx1, *xx2, *pXlast, *pXlast2, *ddimsnew, *u, *w, *trace, *ll;

	/*exPara[0] = N;
	exPara[1] = nX;
	exPara[2] = nX2;
  //YY     = exPara ; filled in later for each Y
  //YY2    = YY + N; filled in later for each Y2
  XX     = exPara + 2*N+3;
	pXXlast = XX + nX*N;
  XX2     = pXXlast + N; //add extra column for snp effect
	pXXlast2 = XX2 + nX2*N;
	ddimsNew = pXXlast2 + N; //add extra column for snp effect
	uu = ddimsNew + 10;
	ww = uu + pts;
	*/

  exPara = (double *) ex;
	N = (int) exPara[0]; //ceil(exPara[0]-0.5);
	nX = (int) exPara[1];//ceil(exPara[1]-0.5);
	nX2 = (int) exPara[2];//ceil(exPara[2]-0.5);
	y1  = exPara + 3;
	y2 = y1 + N; 
  xx1 = y2 + N;
	pXlast = xx1 + (nX-1)*N;
  xx2     = pXlast + N; //add extra column for snp effect
	pXlast2 = xx2 + (nX2-1)*N;
	ddimsnew = pXlast2 + N; 
	u = ddimsnew + 10;
	w = u + (int) ddimsnew[6];

	dims[0] = (int) N;
	dims[1] = (int) ddimsnew[7];
	dims[2] = (int) ddimsnew[8];
	dims[3] = nX; //assume nX is the same as nX2, account for int and snp column
	dims[4] = (int) ddimsnew[6];

	for(i=0;i<n;i++) gr0[i] = 0.0;

	/*
	Rprint_v(x, 0, 5);
	Rprint_v(y1, 0, 5);
	Rprint_v(y2, 0, 5);
	//Rprint_v(xx1, 0, 5);
	//Rprint_v(xx2, 0, 5);
	for(l=0;l<nX;l++){
		Rprintf("XX1 %f XX2 %f \n",xx1[N*l],xx2[N*l]); 
	}
	Rprint_v(u, 0, 4);
	Rprint_v(w, 0, 4);
	Rprint_vi(dims, 0, 4);
	*/

	for(i=0;i<n;i++) gr0[i] = 0.0;

	adaQuadgr_vec(x, y1, y2, xx1, xx2, u, w, dims, 0, gr0);

	for(i=0;i<n;i++){
		if(gr0[i] != gr0[i]) gr0[i]=0.0;
	}
}

void glmEQTL_joint (int* dims, double* Y, double* X, double* Z, double* z1, 
              int* link, double* offset, int* adjZ, char** output,  
              double* RP_cut, int* cis_only, int* cis_distance, 
              int* eChr, int* ePos, int* mChr, int* mPos, double* conv, 
              double* convGLM, int* yFailBaselineModel, double* scoreTestP, 
              int* trace, int* succeed,
							//new params
							double* Y2, double* X2, int* dChr, int* dPos, int* cis_distance2, 
							int* yFailBaselineModel2, double* z2, double* u, double* w /// add in quad points
)
{
  int i, j, k, c, m, nIter, df0, df1, convBase, convSNPj, p , pcount;
  double chisq[4], pval[4], phi0, phi1, scale, twologlik0, twologlik1, l0, l1, ll[4];
/////
	int l, c2, nIter2, df02, df12, convBase2, convSNPj2;
  double chisq2, pval2, phi02, phi12, scale2, twologlik02, twologlik12,  logLik0, logLik1;
///
	

  int family = NB;
  int linkR  = *link;
  int adjZR  = *adjZ;
  int npara, lmm, fail, failA, failB, fncount, grcount, nREPORT;
		
  /* pointers to Y, the last column of X, and Z */
  double *pY, *pXlast, *pZ;
/////
  double *pY2, *pXlast2,*pXXlast,*pXXlast2, *YY, *YY2, *XX, *XX2, *ddimsNew, *uu, *ww;
/////
  
  /* Work array */
  double *Xb, *fitted0, *fitted1, *fitted2, *resid, *weights, *offsetN;
/////
  double *Xb2, *fitted02, *fitted12, *fitted22, *resid2, *weights2, *offsetN2;
/////

  
  /* rank of model matrix */
  int rank0, rank1;
/////
  int rank02, rank12;
/////    

  /* position difference between gene and marker */
  int pos_diff;
/////
  int pos_diff2;
/////
  int g_start =0;
  /* grid used to output frequency */
  double grid;

  /* 
   * p-value frequencies
   * freqs[100] = #{ pval < P_cut }
   * i = 0:99
   * freqs[i]   = #{ pval < [i/100, (i+1)/100) }
   */
  
  unsigned long freqs[101];
  for(p=0; p<=100; p++){ freqs[p] = 0; }
  
  int nY = dims[0];
  int nX = dims[1];
  int nZ = dims[2];
  int N  = dims[3];
  int maxit = dims[4];
  int useOffset = dims[5];
///// adding to dimension of dims to accomodate new vars
	int nY2 = dims[6];
  int nX2 = dims[7];
	///need to add number of quad points to dims
	int pts		= dims[8];

  double P_cut = *RP_cut;

  /* dimsNew is passed to function glmNB, and then function glmFit */
  ///int dimsNew[7];
	int dimsNew[10]; //now has 10 parameters to hold l, h, h1
  dimsNew[0] = N;
  dimsNew[1] = nX;
  dimsNew[2] = maxit;
  dimsNew[3] = 0; /* whetehr to use initial values */
  dimsNew[4] = useOffset;
/////
  dimsNew[5] = nX2;
  dimsNew[6] = pts; //number of quadrature points
	

  /* allocate memory. The extra column in Xb is used to store 
   * one SNP, the same as for X
   */

  //Xb      = (double *) Calloc(N*(nX+1), double);
  Xb      = (double *) Calloc(N*(nX), double); //now we account for snp column in nX, also assume intercept it another cov
  fitted0 = (double *) Calloc(N, double);
  fitted1 = (double *) Calloc(N, double);
  fitted2 = (double *) Calloc(N, double);
  resid   = (double *) Calloc(N, double);
  weights = (double *) Calloc(N, double);
  offsetN = (double *) Calloc(N, double);
///// 
  //Xb2      = (double *) Calloc(N*(nX2+1), double);
	Xb2      = (double *) Calloc(N*(nX2), double);  
  fitted02 = (double *) Calloc(N, double);
  fitted12 = (double *) Calloc(N, double);
  fitted22 = (double *) Calloc(N, double);
  resid2   = (double *) Calloc(N, double);
  weights2 = (double *) Calloc(N, double);
  offsetN2 = (double *) Calloc(N, double);
/////

	/* define beta */ //00 and 01 have nX+NX2 length plus 2 for intercept.  Others have addition +2 for snp effect in each
	//double beta00[nX+nX2+2], beta01[nX+nX2+2], beta20[nX+nX2+2+2], beta21[nX+nX2+2+2];
	double beta00[nX+nX2-2], beta01[nX+nX2-2], beta20[nX+nX2], beta21[nX+nX2];
	

  /* point to the last column of X */
  pXlast  = X + N*(nX-1);
/////
  pXlast2  = X2 + N*(nX2-1);

	double *X_0;
	X_0 = X + N; //skip intercept

  double *X2_0;
	X2_0 = X2 + N; //skip intercept

	///temporory vector of length n
	double diff[N];

  if(*trace){
    Rprintf("\n--------------------------------------------------------------\n");
    Rprintf("(nY, nZ, nX, N, adjZ) = (%d, %d, %d, %d, %d)\t", nY, nZ, nX, N, adjZR);
    Rprintf("P_cut=%e\n", P_cut);
  }
  

	double *exPara;
	//X, X2, Snp col for X, Snp col for X2, Y, Y2, dimsnew (10) + N, nX, and nX2
	///now we also add 2*l for u (l quad points) and w (l weights)
  //exPara = (double *) Calloc(nX*N+nX2*N+ 2*N + 2*N + 10 + 3 , double); 
	
	/// new version	
	//exPara = (double *) Calloc(nX*N+nX2*N+ 2*N + 2*N + 10 + 3 + 2*pts, double); 
    exPara = (double *) Calloc(nX*N+nX2*N+ 2*N + 10 + 3 + 2*pts, double); 
    

	exPara[0] = N;
	exPara[1] = nX;
	exPara[2] = nX2;
  //YY     = exPara ; filled in later for each Y
  //YY2    = YY + N; filled in later for each Y2
  XX     = exPara + 2*N+3;
	pXXlast = XX + (nX-1)*N;
  XX2     = pXXlast + N; //add extra column for snp effect
	pXXlast2 = XX2 + (nX2-1)*N;
	ddimsNew = pXXlast2 + N; //add extra column for snp effect
	uu = ddimsNew + 10;
	ww = uu + pts;

	//copy X, X2 without snp effect now to locations	
	for(i=0; i<N; i++){
		 for(j=0; j<nX; j++){
			XX[N*j + i] = X[N*j + i];
			//if(i < 6) Rprintf("i %d, j %d, X %f, XX %f\n",i, j, X[N*j + i], XX[N*j + i]  );	
		}
	}

	for(i=0; i<N; i++){
  	for(j=0; j<nX2; j++){
			XX2[N*j + i] = X2[N*j + i];	
		}
	}

	//copy u and w to new location
	for(i=0; i<pts; i++){
  	uu[i] = u[i];
		ww[i] = w[i];
	}
	/*Rprint_v(u, 0, 4);
	Rprint_v(uu, 0, 4);
	Rprint_v(w, 0, 4);
	Rprint_v(ww, 0, 4);
	*/

j=0;

		  //npara   = nX+nX2+2+2+2+1; //ncovs for each, 2 intercepts, 2 snp covs, 2 phi, 1 lambda 
			npara   = nX+nX2+2+1; //ncovs for each, 2 phi, 1 lambda 
		  lmm     = 5; 
		  fail    = 0;
		  failA   = 0;
		  failB   = 0;
		  fncount = 0;//?
		  grcount = 0;//?
		  maxit   = 500;
		  nREPORT = 5;
			
			int nbd[npara];
			for(p=0; p < npara; p++){
				
					nbd[p] = 2;
				
			}

		  //technical parameters below:
		  double *wa, *g1;
		  int *iwa;
		  SEXP x1;
		  PROTECT(x1 = allocVector(REALSXP, npara));
		  //consider replacing with simple Calloc
		  wa  = (double *) S_alloc(2*lmm*npara+4*npara+11*lmm*lmm+8*lmm,sizeof(double));
		  iwa = (int*) R_alloc(3*npara,sizeof(int));
		  g1 = (double *)R_alloc(npara, sizeof(double));
		
		  double gr[npara];
		  double initPara[npara];
			double initPara2[npara];
		  double lower[npara];
		  double upper[npara];
		  double Fmin, factr, pgtol, th0, th1, th01;
		  
		  factr = 1e7;
			pgtol = 0.0;
		  
			char msg[1023];
		  
		  /**********************************************************/
		
	Rprintf("ok1 \n");
  /* output file handles *‫/
  

  /* time records */
  time_t sec_s;
  time_t sec_e;

   FILE *fo, *ff;

  /* starting time */
  sec_s = time(NULL);
  
  /* output file for the eQTL mapping results */
  fo = fopen (output[0], "w");
  
  /**
   * write out the header in the output file
   */
  fprintf(fo, "GeneRowID\tMarkerRowID\tChisq0\tPvalue0\tChisq1\tPvalue1\tChisq2\tPvalue2\tChisq3\tPvalue3\tlambda0\tlambda1\n");
  
		Rprint_v(Y, 0, 5);
		Rprint_v(Y2, 0, 5);
		Rprint_v(X, 0, 5);
		Rprint_v(X2, 0, 5); 

  /***
   * identifify eQTL gene by gene
   */
  
  /* pY is the pointer to gene expression data */
  pY = Y;
/////

	Rprintf("ok2 \n");
			
  for(i=0; i<nY; i++,pY+=N){

    if(*trace == 1){
      if(i%100 == 0){ Rprintf("\ni=%d\n", i); }
    }else if(*trace > 1){
      Rprintf("\ni=%d\n", i);
    }

	pY2 = Y2;
  g_start = 0;
  for(l=0; l<nY2; l++,pY2+=N){
	
	//if gene TSS - DNase distance too great, then skip    
		
				if(*cis_only){
      if(eChr[i] != dChr[l]) continue;
        
      pos_diff2 = fabs(ePos[i] - dPos[l]);
      if(pos_diff2 > *cis_distance2 ) continue; //this may need to be modified later
      Rprintf("starting, %d, %d %f\n", i ,l, fabs(ePos[i] - dPos[l]));  
    }
  	   
		//error("stop"); 
 
  
		/*Rprint_v(pY, 0, 5);
		Rprint_v(pY2, 0, 5);
		Rprint_v(X, 0, 5);
		Rprint_v(X2, 0, 5);
 */
    /* **********************************************************
     * fit a baseline model using only the confouding covariates 
     * family is assigned to a value of either Poisson or NB
     * **********************************************************/

		// in this version, just assume nb
    if(g_start==0){   
    dimsNew[1] = nX-2;
    /** 
     * no initial values, Note this value will be changed 
     * in glmNB after initail iteration 
     */
    dimsNew[3] = 0; 




		// *trace = 0;
		///here we are using the nb marginals to get initial estimates for beta, since only glmNB has procedure to get actually coefs
		///fit covariate matrix with intercept but without SNP effect. X has no intercept column.  nX also exlcudes snp effect
		///everything after this section must use BFGS since assuming lambda ~ MVN(0, Sigma) 

    convBase = glmNB(dimsNew, &nIter, pY, z1, &linkR, offset, X_0, convGLM,
                     &rank0, Xb, fitted0, resid, weights, &phi0, &scale, 
                     &df0, &family, &twologlik0, scoreTestP, trace, beta00);
   
    if (!convBase) {
      if(*trace>0){
        Rprintf("\n  glmEQTL: i=%d, initial model (H0) fail to converge\n", i);
      }

      yFailBaselineModel[i] = 1; 
      continue;
    }
    g_start=0;
    }
    dimsNew[3] = 0; 
		//setting offset to the same one used as for GE data
    convBase2 = glmNB(dimsNew, &nIter2, pY2, z2, &linkR, offset, X2_0, convGLM,
                     &rank02, Xb2, fitted02, resid2, weights2, &phi02, &scale2, 
                     &df02, &family, &twologlik02, scoreTestP, trace, beta01);
        
    if (!convBase2) {
      if(*trace>0){
        Rprintf("\n  glmEQTL2: i=%d, initial model (H0) fail to converge\n", i);
      }

      yFailBaselineModel2[i] = 1; 
      continue;
    }

		/*for(p=0;p<nX-1;p++){
			Rprintf("beta1 %d %f",p,beta00[p]);
			Rprintf("\n");
			Rprintf("beta2 %d %f",p,beta01[p]);
			Rprintf("\n");
		}*/
    
    if(*trace > 1){
      Rprintf("\n  glmEQTL: i=%d, finish fitting initial model: ");
      Rprintf("phi0 = %e, ", phi0);
      if (family==2) {
        Rprintf("family=Poisson\n");
      }else if (family==5) {
        Rprintf("family=NB\n");
      }else {
        Rprintf("family=%d\n", family);
      }
							//print_v(beta, 0, 6);
							Rprintf("%f\n",mean(pY, 0, 59)); 

    }
    
    /* **********************************************************
     * Now start to fit models with both confounding covariates
     * and each of the SNPs 
     * **********************************************************/
    
    pZ = Z;
    
    for(j=0; j<nZ; j++,pZ+=N){
			//Rprintf("i %d l %d j %d\n", i, l, j);
			//here we assume DNase is close to gene, so just check to see if snp is close to gene
      if(*cis_only){
        if(dChr[l] != mChr[j]) continue;
        
        pos_diff = fabs(dPos[l] - mPos[j]);
        
        if(pos_diff > *cis_distance) continue;
      }
      
      //if(*trace > 3){
        Rprintf("l=%d, i=%d, j=%d, pos_pair=%d pos_snp=%d %d %d %d\n", l+1, i+1, j+1, pos_diff2,pos_diff, ePos[i], dPos[l], mPos[j]);
      //}
            
      /* *
       * fill the last column of X by Z[j], genotype of one SNP
       * start with the fitted values of the confounder-only model
       */
      
      for (k=0; k<N; k++) {
				//only for now, update this later 
				if(pZ[k] == 3.0) pZ[k] = 1.0; 
				if(pZ[k] == 4.0) pZ[k] = 2.0;
        pXlast[k]  = pZ[k];
				pXXlast[k]  = pZ[k]; //for lbfgs
        fitted1[k] = fitted0[k];
/////
				pXlast2[k]  = pZ[k];
				pXXlast2[k]  = pZ[k]; //for lbfgs
        fitted12[k] = fitted02[k];
/////
      }
      
      //phi1  = phi0;
      scale = 1.0;
/////
      //phi12  = phi02;
      scale2 = 1.0;
/////
      
      //if(*trace > 2){
      //  Rprintf("copied covs\n");
      //}
			
			/// in this version we are skipping this whole section that follows

      /* family has been decided in the previous fitting of baseline model *
      //wrote code assuming adjZ = 0
      if (adjZR) {
        dimsNew[1] = nX;
        dimsNew[3] = 1; /* use initial values *
        
        convSNPj = glmNBlog(dimsNew, &nIter, pY, z1, &linkR, offset, X, conv, convGLM, 
                            &rank1, Xb, fitted1, resid, weights, &phi1, &scale, 
                            &df1, &family, &twologlik1, scoreTestP, trace, beta20,
                            fitted2, offsetN);        

        dimsNew[3] = 1; /* use initial values *       
				dimsNew[5] = nX2;
        convSNPj2 = glmNBlog(dimsNew, &nIter2, pY2, z2, &linkR, offset, X2, conv, convGLM, 
                            &rank12, Xb2, fitted12, resid2, weights2, &phi12, &scale2, 
                            &df12, &family, &twologlik12, scoreTestP, trace, beta21,
                            fitted22, offsetN2);        
      }else {
        dimsNew[1] = nX + 1;
        dimsNew[3] = 1; /* use initial values *
        
        convSNPj = glmNB(dimsNew, &nIter, pY, z1, &linkR, offset, X, convGLM, 
                         &rank1, Xb, fitted1, resid, weights, &phi1, &scale, 
                         &df1, &family, &twologlik1, scoreTestP, trace, beta20);        
       
				dimsNew[5] = nX2 + 1;
        convSNPj2 = glmNB(dimsNew, &nIter2, pY2, z2, &linkR, offset, X2, convGLM, 
                         &rank12, Xb2, fitted12, resid2, weights2, &phi12, &scale2, 
                         &df12, &family, &twologlik12, scoreTestP, trace, beta21);        

      }
			
      if(convSNPj == 0){
        if(*trace){
          Rprintf("\n  Fail GLM: i=%d, j=%d, adjZR=%d, ", i, j, adjZR);
          if (family==2) {
            Rprintf("family=Poisson\n");
          }else if (family==5) {
            Rprintf("family=NB\n");
          }else {
            Rprintf("family=%d\n", family);
          }          
        }
        continue;
      }

      if(convSNPj2 == 0){
        if(*trace){
          Rprintf("\n  Fail GLM: i=%d, j=%d, adjZR=%d, ", i, j, adjZR);
          if (family==2) {
            Rprintf("family=Poisson\n");
          }else if (family==5) {
            Rprintf("family=NB\n");
          }else {
            Rprintf("family=%d\n", family);
          }          
        }
        continue;
      }
      Rprintf("done marginals");
      /**
       * it is possible that df0 - df1 != 1
       * in glmFit, df = Nu - x_rank. It is possible that Nu is smaller than sample size
       * due to invalid fitted values for certain glm, e.g., negative values for Poisson
       *
      
      if(df0 - df1 != 1){
        Rprintf("i=%d, j=%d, df0=%d, df1=%d, rank0=%d, rank1=%d\n", i, j, df0, df1, rank0, rank1);
        
        if(df0 - df1 < 0.5) continue;
      }

      if(df0 - df1 != 1){
        Rprintf("2 i=%d, j=%d, df0=%d, df1=%d, rank0=%d, rank1=%d\n", l, j, df02, df12, rank02, rank12);
        
        if(df02 - df12 < 0.5) continue;
      }
      
			///// Now fit joint expression-dnase model:  with and without snp effect and compare to sum of marginals
			// load up C implementation of L-BFGS-B, negloglik, and neggradient functions for with and without snp effect
			// make negloglik and gradient functions flexible to handle with and without snp
/* **********************************************************
		   * parameters for function lbfgsb, which will be used to
		     obtain MLE for H1: with allelic imbalance
		   
		   void lbfgsb(int n, int lmm, double *x, double *lower,
		          double *upper, int *nbd, double *Fmin, optimfn fn,
		          optimgr gr, int *fail, void *ex, double factr,
		          double pgtol, int *fncount, int *grcount,
		          int maxit, char *msg, int trace, int nREPORT);
		   
		   n:       the number of parameters
		   
		   lmm:     is an integer giving the number of BFGS updates 
		            retained in the "L-BFGS-B" method, It defaults to 5.
		   
		   x:       starting parameters on entry and the final parameters on exit
		   
		   lower:   lower bounds
		   
		   upper:   upper bounds
		   
		   nbd:     specifies which bounds are to be used. 
		            nbd(i)=0 if x(i) is unbounded,
		            1 if x(i) has only a lower bound,
		            2 if x(i) has both lower and upper bounds, and
		            3 if x(i) has only an upper bound.
		            On exit nbd is unchanged.
		   
		   Fmin:    final value of the function
		   
		   fn:      the function to be minimized
		   
		   gr:      the gradient function
		   
		   fail:    integer code, 0 for success, 51 for warning and 52 for error
		   
		   ex:      extra parameters for the function to be minimized
		   
		   factr:   controls the convergence of the "L-BFGS-B" method. 
		            Convergence occurs when the reduction in the objective is 
		            within this factor of the machine tolerance. Default is 1e7, 
		            that is a tolerance of about 1e-8.
		   
		   pgtol:   helps control the convergence of the "L-BFGS-B" method. 
		            It is a tolerance on the projected gradient in the current 
		            search direction. This defaults to zero, when the check 
		            is suppressed.
		   
		   fncount: the number of calls to fn 
		   
		   grcount: the number of calls to gr
		   
		   maxit:   maximum of iterations
		   
		   msg:     A character string giving any additional information 
		            returned by the optimizer, or NULL
		   
		   trace:   Non-negative integer. If positive, tracing information 
		            on the progress of the optimization is produced. 
		            Higher values may produce more tracing information: 
		            for method "L-BFGS-B" there are six levels of tracing. 
		   
		   nREPORT: The frequency of reports for the "BFGS", "L-BFGS-B" 
		            and "SANN" methods if control$trace is positive. 
		            Defaults to every 10 iterations for "BFGS" and "L-BFGS-B"
		   
		   * **********************************************************/
		  
			
    	///update expara for the current values of Y and Y2 
	    YY = exPara+3;
	    YY2 = exPara + N+3;
			for(p=0; p < N; p++){
				YY[p] = pY[p];
				YY2[p] = pY2[p];
			}

			///update dimsNew
			dimsNew[1] = (double) nX; /* assuming adjZ = FALSE */
			dimsNew[5] = (double) nX2; /* assuming adjZ = FALSE */			
			dimsNew[7]   =  (double) 0; // we will update this later in ddimsnew
			dimsNew[8]   =  (double) 0; // we will update this later in ddimsnew

			for(p=0; p<10; p++) ddimsNew[p] = dimsNew[p];
 	
			//npara   = nX+nX2+2+2+2+1; //ncovs for each + 2 intercepts, 2 snp effects, 2 sigma, 1 rho 
		  npara   = nX+nX2+2+1; //ncovs for each + 2 sigma, 1 rho 
		  
   //if(*trace > 2){
   //     Rprintf("copied dims\n");
   //   }
		
        
			/* initPara = c(beta, beta2, phi, phi2, lambda) */
			// using beta and phi from null snp and lambda initialization

		  for(p=0; p < npara; p++){		
			  	if(p < nX-1){ 
					initPara[p] = beta00[p];
					lower[p] = -300; 
					nbd[p] = 0;
		  	}else if(p == nX-1){ 
					initPara[p] = 0; //snp effect is 0 under h0
					lower[p] = -300;  
					nbd[p] = 0; 
		  	}else if(p > nX-1 & p<nX+nX2-1){
					initPara[p] = beta01[p-nX];
					lower[p] = -300;
					nbd[p] = 0;
				}else if(p == nX+nX2-1){
					initPara[p] = 0; //snp effect is 0
					lower[p] = -300;  
					nbd[p] = 0;
				}else if(p == nX+nX2){
					///sigma1
					for(m=0;m<N;m++) diff[m] = log(YY[m]+1) - log(fitted0[m]+1);
					initPara[p] = var(diff, 0, N-1, mean(diff, 0, N-1)); ///this may be problematic
					lower[p] = 1e-3;
					nbd[p] = 1;
				}else if(p == nX+nX2+1){
					///sigma2
					for(m=0;m<N;m++) diff[m] = log(YY2[m]+1) - log(fitted02[m]+1);
					initPara[p] = var(diff, 0, N-1, mean(diff, 0, N-1)); ///this may be problematic
					lower[p] = 1e-3;
					nbd[p] = 1;
				}else{
					///rho
					initPara[p] = 0.0;  //moment estimates are insane due to0 small m1 and m2, need better way to initalize
					lower[p] = -1.0+1e-3; 
					nbd[p] = 2;
				}
				upper[p] = 300; //ignored because of nbd?
			   //if(*trace > 2){
    		    //Rprintf("p: %d %f\n", p, initPara[p]);
		     // }
			}	
			upper[npara-1] = 1-1e-3; //for rho

			//Rprintf("End init, start lbfgs, npara is %d\n", npara);

			/*Rprint_v(initPara, 0, 20);
			Rprint_v(exPara, 0, 5);
			Rprint_v(YY, 0, 5);       			
			Rprint_v(YY2, 0, 5);
	Rprint_v(u, 0, 4);
	Rprint_v(uu, 0, 4);
	Rprint_v(w, 0, 4);
	Rprint_v(ww, 0, 4);*/

			///This part is crucial. Need to pay attention to ddimsNew, where we update h and h1, nothing else changes
			//assumption is that P = nX + 1 = nX2 + 1 
			//for each round of estimation, everything stays the same except h and h1 which get updated 
			//also need to reset initPara after each round, as well as failA
			///this may cause issues here

			//Rprint_v(pXlast, 0, 20);
			//Rprint_v(pXlast2, 0, 20);
			double fail[4];
			for(m=0;m<4;m++){
				if(m==0){
					ddimsNew[7]   =  (double) 1;
					ddimsNew[8]   =  (double) 1;
				}else if(m==1){
					ddimsNew[7]   =  (double) 0;
					ddimsNew[8]   =  (double) 1;
				}else if(m==2){
					ddimsNew[7]   =  (double) 1;
					ddimsNew[8]   =  (double) 0;
				}else if(m==3){
					ddimsNew[7]   =  (double) 0;
					ddimsNew[8]   =  (double) 0;
				}
		
				for(p = 0; p<npara; p++) initPara2[p]=initPara[p];	


      	lbfgsb1(npara, lmm, initPara2, lower, upper, nbd, &Fmin, 
           adaQuad_vec_bfgs, adaQuadgr_vec_bfgs, &failA, (void*)exPara, factr, pgtol,  
           &fncount, &grcount, maxit, msg, 0, nREPORT, wa, iwa, g1,x1);
	      	ll[m] = Fmin;
		fail[m]=failA;
		if(fail[m] == 1 & strcmp(msg,"NEW_X")==0) fail[m] = -2; 
		for(p=0; p<npara;p++) if((initPara2[p] == lower[p] | initPara2[p] == upper[p]) & initPara2[p]!=0) fail[m] = -3;

					Rprintf("i %d ll %f coef ", m, Fmin);
					Rprint_v(initPara2, 0, npara-1);	
					Rprintf("%d %s\n", failA, msg);					
					if (failA) {
	      		if (*trace){
	        		Rprintf("  i=%d, fail to fit baseline ASE model @ situation A\n", i);	    
  						Rprintf("%d %s\n", failA, msg);
						}
	      			//break;
	    		}else{
	      		//th0 = initPara[0];
	      
	      		if(*trace > 1){
	        		Rprintf("\n  Obtained MLE for H0: twologlik0=%.4e", logLik0);
	        		Rprintf(" theta0A=%.4e", th0);
	      		}
	    		}
			}
			
			//if(failA == 1) continue;  //failA gets reset to zero by lbfgsb?  otherwise need to reset to 0 here

			/*
			//setup for h1, here we are including the snp effect            
			//need to move some of these initialionns up to the top
			// npara changes based on h0 and h1, may cause probs later
		  
			dimsNew[6]   = 1.0;  /* null hyp is true*
			for(p=0; p<7; p++) ddimsNew[p] = dimsNew[p];

	    //lbfgsb(npara, lmm, initPara, lower, upper, nbd, &Fmin, 
	    //       negLogH0, negGradLogH0, &failA, (void*)exPara, factr, pgtol,  
	    //       &fncount, &grcount, maxit, msg, 0, nREPORT);
      lbfgsb1(npara, lmm, initPara, lower, upper, nbd, &Fmin, 
             negLog, negGradLog, &failB, (void*)exPara, factr, pgtol,  
             &fncount, &grcount, maxit, msg, 0, nREPORT, wa, iwa, g1,x1);

	    if (failB) {
	      if (*trace)
	        Rprintf("  i=%d, fail to fit baseline ASE model @ situation A\n", i);
	    
	     // continue;
	    }else{
	      th01 = initPara[0];
	      logLik1 = -Fmin;
	      
	      if(*trace > 1){
	        Rprintf("\n  Obtained MLE for H0: twologlik0=%.4e", logLik0);
	        Rprintf(" theta0A=%.4e", th0);
	      }
	    }
			*/
			//Rprintf("%f %f \n%f %f \n%f %f\n", twologlik02 ,  twologlik0, twologlik12 ,  twologlik1 , logLik1, logLik0);
      //chisq = twologlik1 - twologlik0;
			//chisq2 = twologlik12 - twologlik02;

      /*if (chisq < -1e-5) {
				Rprintf("\n");
				//print_v(beta, 0, 7);
        error("wrong twologlik! i=%d, j=%d, twologlik=(%f, %f)\n", 
              i, j, twologlik0, twologlik1);
      }

      /*if (chisq < 1e-5) { 
        pval = 1.0; 
      }else{
        /* pchisq(double x, double df, int lower_tail, int give_log) *
        pval  = pchisq(chisq, (double)(df0 - df1), 0, 0);
      }
      
      k = (int) (pval / 0.01);
      freqs[k] += 1;
      */

      //chisq[0] = twologlik1 - twologlik0; //snp effect in GE only, 1df test
			chisq[0] = 2*(ll[1]-ll[0]); //cor in presence of snp
			chisq[1] = 2*(ll[2]-ll[0]); //snp effect in cor
			chisq[2] = 2*(ll[3]-ll[1]);//snp effect withour corr
			chisq[3] = 2*(ll[3]-ll[2]); //cor effect without snp
			//Rprintf(" 1 j %d\n",j);

			
      /*if (chisq < -1e-5) {
				Rprintf("\n");
				//print_v(beta, 0, 7);
        error("wrong twologlik! i=%d, j=%d, twologlik=(%f, %f)\n", 
              i, j, twologlik0, twologlik1);
      }*/

      /*if (chisq < 1e-5) { 
        pval = 1.0; 
      }else{
        /* pchisq(double x, double df, int lower_tail, int give_log) */
      //}

			//pval  = pchisq(chisq, (double)(df0 - df1), 0, 0);
      //pval[0]  = pchisq(chisq[0], 1.0, 0, 0);
			pval[0]  = 1-pchisq(chisq[0], 1.0, 1.0, 0);
			pval[1]  = 1-pchisq(chisq[1], 2.0, 1.0, 0);
			pval[2]  = 1-pchisq(chisq[2], 2.0, 1.0, 0);
			pval[3]  = 1-pchisq(chisq[3], 1.0, 1.0, 0);
			pcount = 0;  
			for(p=0; p<4; p++){  
				      
      	//if(pval[p] <= P_cut){
					//pcount = pcount + 1;
					//Rprintf("2 j %d\n",j);
       		if(p==0) fprintf(fo, "%d\t%d\t%d\t", i+1, l+1, j+1);
					fprintf(fo, "%.3f\t%.2e\t", fail[p], pval[p]);
				//}else{
				//	if(pcount>0) fprintf(fo, "%s\t%s\t", "NA", "NA");
				//}
				//Rprintf("chi %f p %f\n", chisq[p], pval[p]);
      }
			fprintf(fo, "\n");
			//error("done, i %d l %d j %d\n", i+1, l+1, j+1);
	    		
			//if(pcount>0) fprintf(fo, "%.3f\t%.3f\n", l0, l1);
		//Rprintf("3 j %d\n",j);
      /*if(pval < P_cut){
        freqs[100] += 1;
        
        /* gene ID and SNP ID *
        fprintf(fo, "%d\t%d\t", i+1, j+1);
        fprintf(fo, "%d\t%.3f\t%.2e\n", family, chisq, pval);
      }*/


      
    }
	}
  }
  
  fclose(fo);

  // print out the frequencies
  /*ff   = fopen(output[1], "w");
  grid = 0.0;

  fprintf(ff, "<%.2e", P_cut);
  for(i=0; i<100; i++){
    fprintf(ff, "\t%.2f-%.2f", grid, grid+0.01);
    grid += 0.01;
  }
  fprintf(ff, "\n");
  
  fprintf(ff, "%lu", freqs[100]);
  for(i=0; i<100; i++){
    fprintf(ff, "\t%lu", freqs[i]);
  }
  fprintf(ff, "\n");  
  fclose(ff);
  

  /* end time */
  sec_e  = time(NULL);
  if(*trace){
    Rprintf("total time spent is %ld secs\n", sec_e-sec_s);
    Rprintf("\n--------------------------------------------------------------\n");
  }
    
  Free(Xb);
  Free(fitted0);
  Free(fitted1);
  Free(fitted2);
  Free(resid);
  Free(weights);
  Free(offsetN);

  Free(Xb2);
  Free(fitted02);
  Free(fitted12);
  Free(fitted22);
  Free(resid2);
  Free(weights2);
  Free(offsetN2);
  Free(exPara);
  *succeed = 1;
 UNPROTECT(1);

//}
}

