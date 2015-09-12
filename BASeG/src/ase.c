#include "ase.h"
#include "asa103.h"
#include "utility.h"

//used for maximization of posterior mode of biv random effect via lbfgsb
double l1ASE(int n, double* para, void* ex, SEXP x1){
			int i, k, N, P, h, h2, l, nA1, nA2, n1, n2, geno;
		double *exPara, *mean0, *sigma0,*t,p,val, pi1, pi2, cov;
		double l1, l2;
		
		exPara = (double *) ex;
		N      = exPara[0];
		geno   = exPara[1];
		p      = exPara[2];
		nA1    = (int) exPara[3];
		nA2    = (int) exPara[4];
		h 	 	 = exPara[5];
		h2 	 	 = exPara[6];
		pi1    = exPara[7];
		pi2    = exPara[8];
		n1		 = (int) exPara[9];
		n2	 	 = (int) exPara[10];
		mean0  = exPara + 11; 
		sigma0 = mean0 + 2; 
		t      = para;

		//for(l=0;l<11+4*P;l++) Rprintf("h %f", exPara[l]);

		if(h ==0 )  p = 0.0;
		if(h2 ==0)  pi1 = pi2 = 0.0;

	        if(geno==1){
        	        cov = 1.0;
        	}else if(geno==3){
                	cov = -1.0;
        	}else{
                	cov = 0.0;
        	}

		l1 = exp(pi1*cov + t[0])/(1+exp(pi1*cov + t[0]));
		l2 = exp(pi2*cov + t[1])/(1+exp(pi2*cov + t[1]));

		if(l1!=l1){
			if(pi1*cov+t[0]>0.0){
				l1 = 1.0;
			}else{
				l1 = 0;
			}
		}

		if(l2!=l2){
			if(pi2*cov+t[1]>0.0){
				l2= 1.0;
			}else{
				l2 = 0;
			}
		}

		val = dbinom(nA1,n1,l1,0)*dbinom(nA2,n2,l2,0)*dnorm(t[0], mean0[0] + p*sqrt(sigma0[0])*(t[1]- mean0[1])/sqrt(sigma0[1]), sqrt(sigma0[0]*(1-p*p)),0)*dnorm(t[1], mean0[1], sqrt(sigma0[1]),0);
 	//Rprintf("l1 val %e l1 %f l2 %f t1 %f t2 %f %f %f %f %f %f %e\n", -log(val + pow(10,-100)), l1 , l2, t[0], t[1], mean0[0], mean0[1],sigma0[0], sigma0[1],p, dnorm(t[0], mean0[0] + p*sqrt(sigma0[0])*(t[1]- mean0[1])/sqrt(sigma0[1]), sqrt(sigma0[0]*(1-p*p)),0)*dnorm(t[1], mean0[1], sqrt(sigma0[1]),0));
  return(-log(val + pow(10,-100)));
}

void grASE(int n, double* para, double* gr0, void* ex, SEXP x1)
{
			int i, k, N, P, h, h2, l, nA1, nA2, n1, n2, geno;
		double *exPara, *mean0, *sigma0,*t,p,val, pi1, pi2, dt1, dt2, cov;
		double l1, l2;
		
		exPara = (double *) ex;
		N      = exPara[0];
		geno   = exPara[1];
		p    	 = exPara[2];
		nA1    = (int) exPara[3];
		nA2    = (int) exPara[4];
		h 	 	 = exPara[5];
		h2 	 	 = exPara[6];
		pi1    = exPara[7];
		pi2    = exPara[8];
		n1		 = (int) exPara[9];
		n2	 	 = (int) exPara[10];
		mean0  = exPara + 11; 
		sigma0 = mean0 + 2; 
		t      = para;

		//for(l=0;l<11+4*P;l++) Rprintf("h %f", exPara[l]);

		if(h ==0 )  p = 0.0;
		if(h2 ==0)  pi1 = pi2 = 0.0;

                if(geno==1){
                        cov = 1.0;
                }else if(geno==3){
                        cov = -1.0;
                }else{
                      	cov = 0.0;
                }

		l1 = exp(pi1*cov + t[0])/(1+exp(pi1*cov + t[0]));
		l2 = exp(pi2*cov + t[1])/(1+exp(pi2*cov + t[1]));

		if(l1!=l1){
			if(pi1*cov+t[0]>0.0){
				l1 = 1.0;
			}else{
				l1 = 0;
			}
		}

		if(l2!=l2){
			if(pi2*cov+t[1]>0.0){
				l2= 1.0;
			}else{
				l2 = 0;
			}
		}

  	dt1 = (nA1-n1*l1) - (1/(2*(1-p*p))) * (2*(t[0]-mean0[0])/sigma0[0] - 2*p*(t[1] - mean0[1])/sqrt(sigma0[0]*sigma0[1]));
		dt2 = (nA2-n2*l2) - (1/(2*(1-p*p))) * (2*(t[1]-mean0[1])/sigma0[1] - 2*p*(t[0] - mean0[0])/sqrt(sigma0[0]*sigma0[1])); 
	
  	gr0[0] = -dt1;
		gr0[1] = -dt2;
		//Rprintf("gr l1 %f l2 %f t1 %f t2 %f dtq %f dt2 %f\n",l1 , l2, t[0], t[1], dt1, dt2);
		//error("");
}

void hessASE(double* t, double *xx, double nA1, double nA2, double n1, double n2, int h, int h2, int P, double *result, int geno){

	//P = 1 here for ASE
	int l;
	double pi1, pi2, *sigma0, p, l1, l2, ddt1, ddt2, ddt12, cov;
	if(h==0) xx[2*P+2] = 0;
	if(h2==0) xx[P-1] = xx[2*P-1] = 0; 

        if(geno==1){
                cov = 1.0;
        }else if(geno==3){
                cov = -1.0;
        }else{
                cov = 0.0;
        }

	pi1 = xx[0];
	pi2 = xx[1];
	sigma0 = xx + 2;  
	p  = xx[2*P+2];

	l1 = exp(pi1*cov + t[0])/(1+exp(pi1*cov + t[0]));
	l2 = exp(pi2*cov + t[1])/(1+exp(pi2*cov + t[1]));
	
		if(l1!=l1){
			if(pi1*cov+t[0]>0.0){
				l1 = 1.0;
			}else{
				l1 = 0;
			}
		}

		if(l2!=l2){
			if(pi2*cov+t[1]>0.0){
				l2= 1.0;
			}else{
				l2 = 0;
			}
		}

  ddt1 = n1*l1*(1-l1) + 1/((1-p*p)*sigma0[0]);
	ddt2 = n2*l2*(1-l2) + 1/((1-p*p)*sigma0[1]);
	ddt12 = -p/((1-p*p)*sqrt(sigma0[0]*sigma0[1]));
	//return vector representing 2x2 matrix with columns stacked
	result[0] = ddt1;
	result[1] = ddt12;
	result[2] = ddt12;
	result[3] = ddt2;	
} 

double adaQuadASE(double* xx, double* u, double* w, double* t,double* Qhat,double nA1, double nA2, double n1, double n2,int* dims, int geno){

	int l, j, k;
	double pi1, pi2, *sigma0, p, l1, l2, l01,l02,val, wj, wk, zstar1, zstar2;
	double Qhat2[4];
	double mean0[2];
	double cov;

	int N = dims[0];
	int h = dims[1];
	int h2 = dims[2];
	int P = dims[3];
	int nhrule = dims[4];

	if(h==0) xx[2*P+2] = 0;
	if(h2==0) xx[P-1] = xx[2*P-1] = 0; 

        if(geno==1){
                cov = 1.0;
        }else if(geno==3){
                cov = -1.0;
        }else{
              	cov = 0.0;
        }
	

	pi1 = xx[0];
	pi2 = xx[1];
	sigma0 = xx + 2;
	p  = xx[2*P+2];

	mean0[0] = 0;///m -sigma0[0]/2;
	mean0[1] = 0;///m-sigma0[1]/2;

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
			l01 = exp(pi1*cov + zstar1)/(1+exp(pi1*cov + zstar1));
			l02 = exp(pi2*cov + zstar2)/(1+exp(pi2*cov + zstar2));

			if(l01!=l01){
				if(pi1*cov+zstar1>0.0){
					l01 = 1.0;
				}else{
					l01 = 0;
				}
			}

			if(l02!=l02){
				if(pi2*cov+zstar2>0.0){
					l02 = 1.0;
				}else{
					l02 = 0;
				}
			}

			//Rprintf("%f %f %f %f %f %f %e\n", zstar1, zstar2, l1, l2, l01, l02,  2*sqrt(det(Qhat))*wj*wk*dpois(y1, l01, 0)*dpois(y2, l02, 0)*dnorm(zstar1, mean0[0] + p*sqrt(sigma0[0])*(zstar2- mean0[1])/sqrt(sigma0[1]), sqrt(sigma0[0]*(1-p*p)),0)*dnorm(zstar2, mean0[1], sqrt(sigma0[1]),0)*exp(u[j]*u[j])*exp(u[k]*u[k]));
			val = val + 2*sqrt(det(Qhat))*wj*wk*dbinom(nA1,n1,l01,0)*dbinom(nA2,n2,l02,0)*dnorm(zstar1, mean0[0] + p*sqrt(sigma0[0])*(zstar2- mean0[1])/sqrt(sigma0[1]), sqrt(sigma0[0]*(1-p*p)),0)*dnorm(zstar2, mean0[1], sqrt(sigma0[1]),0)*exp(u[j]*u[j])*exp(u[k]*u[k]);
	//Rprintf("2*det %e ww %e pois1 %e, pois2 %e, mvnorm %e, exps %e, val %e, cumsum val %e\n", 2*sqrt(det(Qhat)), wj*wk,dpois(y1, l01, 0), dpois(y2, l02, 0), dnorm(zstar1, mean0[0] + p*sqrt(sigma0[0])*(zstar2- mean0[1])/sqrt(sigma0[1]), sqrt(sigma0[0]*(1-p*p)),0)*dnorm(zstar2, mean0[1], sqrt(sigma0[1]),0),exp(u[j]*u[j])*exp(u[k]*u[k]),2*sqrt(det(Qhat))*wj*wk*dpois(y1, l01, 0)*dpois(y2, l02, 0)*dnorm(zstar1, mean0[0] + p*sqrt(sigma0[0])*(zstar2- mean0[1])/sqrt(sigma0[1]), sqrt(sigma0[0]*(1-p*p)),0)*dnorm(zstar2, mean0[1], sqrt(sigma0[1]),0)*exp(u[j]*u[j])*exp(u[k]*u[k]), val); 
		}		
	}
	//Rprintf("%e final val \n",val);
	return(-log(val+pow(10,-100)));
}

void adaQuadgrASE(double* xx, double* u, double* w, double* t,double* Qhat,double nA1, double nA2, double n1, double n2,int* dims, double* gr, int *trace, int geno){

	int l, j, k;
	double pi1,pi2,*sigma0, p, l1, l2, l01,l02,val, wj, wk, zstar1, zstar2,L;
	double Qhat2[4];
	double mean0[2];

	int N = dims[0];
	int h = dims[1];
	int h2 = dims[2];
	int P = dims[3];
	int nhrule = dims[4];
	double cov;

	if(h==0) xx[2*P+2] = 0;
	if(h2==0) xx[P-1] = xx[2*P-1] = 0; 

	if(geno==1){
		cov = 1.0;
	}else if(geno==3){
		cov = -1.0;
	}else{
		cov = 0.0;
	}
	
	pi1 = xx[0];
	pi2 = xx[1];
	sigma0 = xx + 2;
	p  = xx[2*P+2];

	mean0[0] = 0;///m -sigma0[0]/2;
	mean0[1] = 0;///m-sigma0[1]/2;

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
			l01 = exp(pi1*cov + zstar1)/(1+exp(pi1*cov + zstar1));
			l02 = exp(pi2*cov + zstar2)/(1+exp(pi2*cov + zstar2));

			if(l01!=l01){
				if(pi1*cov+zstar1>0.0){
					l01 = 1.0;
				}else{
					l01 = 0;
				}
			}

			if(l02!=l02){
				if(pi2*cov+zstar2>0.0){
					l02 = 1.0;
				}else{
					l02 = 0;
				}
			}

			//Rprintf("%f %f %f %f %f %f %e\n", zstar1, zstar2, l1, l2, l01, l02,  2*sqrt(det(Qhat))*wj*wk*dpois(y1, l01, 0)*dpois(y2, l02, 0)*dnorm(zstar1, mean0[0] + p*sqrt(sigma0[0])*(zstar2- mean0[1])/sqrt(sigma0[1]), sqrt(sigma0[0]*(1-p*p)),0)*dnorm(zstar2, mean0[1], sqrt(sigma0[1]),0)*exp(u[j]*u[j])*exp(u[k]*u[k]));
			val = val + 2*sqrt(det(Qhat))*wj*wk*dbinom(nA1,n1,l01,0)*dbinom(nA2,n2,l02,0)*dnorm(zstar1, mean0[0] + p*sqrt(sigma0[0])*(zstar2- mean0[1])/sqrt(sigma0[1]), sqrt(sigma0[0]*(1-p*p)),0)*dnorm(zstar2, mean0[1], sqrt(sigma0[1]),0)*exp(u[j]*u[j])*exp(u[k]*u[k]);
	//Rprintf("2*det %e ww %e pois1 %e, pois2 %e, mvnorm %e, exps %e, val %e, cumsum val %e\n", 2*sqrt(det(Qhat)), wj*wk,dpois(y1, l01, 0), dpois(y2, l02, 0), dnorm(zstar1, mean0[0] + p*sqrt(sigma0[0])*(zstar2- mean0[1])/sqrt(sigma0[1]), sqrt(sigma0[0]*(1-p*p)),0)*dnorm(zstar2, mean0[1], sqrt(sigma0[1]),0),exp(u[j]*u[j])*exp(u[k]*u[k]),2*sqrt(det(Qhat))*wj*wk*dpois(y1, l01, 0)*dpois(y2, l02, 0)*dnorm(zstar1, mean0[0] + p*sqrt(sigma0[0])*(zstar2- mean0[1])/sqrt(sigma0[1]), sqrt(sigma0[0]*(1-p*p)),0)*dnorm(zstar2, mean0[1], sqrt(sigma0[1]),0)*exp(u[j]*u[j])*exp(u[k]*u[k]), val); 
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
			l01 = exp(pi1*cov + zstar1)/(1+exp(pi1*cov + zstar1));
			l02 = exp(pi2*cov + zstar2)/(1+exp(pi2*cov + zstar2));

			if(l01!=l01){
				if(pi1*cov+zstar1>0.0){
					l01 = 1.0;
				}else{
					l01 = 0;
				}
			}

			if(l02!=l02){
				if(pi2*cov+zstar2>0.0){
					l02 = 1.0;
				}else{
					l02 = 0;
				}
			}
			L = 2*sqrt(det(Qhat))*wj*wk*dbinom(nA1,n1,l01,0)*dbinom(nA2,n2,l02,0)*dnorm(zstar1, mean0[0] + p*sqrt(sigma0[0])*(zstar2- mean0[1])/sqrt(sigma0[1]), sqrt(sigma0[0]*(1-p*p)),0)*dnorm(zstar2, mean0[1], sqrt(sigma0[1]),0)*exp(u[j]*u[j])*exp(u[k]*u[k]);
			//Rprintf("%f %f %f %f %f %f\n", zstar1, zstar2, l1, l2, l01, l02);

			if(h2 ==1 & (geno == 1 | geno == 3)){
				gr[0] = gr[0] + cov*(nA1 - n1*l01)*L/val;
				gr[1] = gr[1] + cov*(nA2 - n2*l02)*L/val;
			}

			gr[2*P] = gr[2*P] + (-1/(2*sigma0[0]) -(1/(2*(1-p*p)))*( -pow(zstar1-mean0[0],2)/pow(sigma0[0],2) + p*(zstar1-mean0[0])*(zstar2-mean0[1])*sigma0[1]/pow(sigma0[0]*sigma0[1],1.5)))*L/val;
			gr[2*P+1] = gr[2*P+1] + (-1/(2*sigma0[1]) -(1/(2*(1-p*p)))*( -pow(zstar2-mean0[1],2)/pow(sigma0[1],2) + p*(zstar1-mean0[0])*(zstar2-mean0[1])*sigma0[0]/pow(sigma0[0]*sigma0[1],1.5)))*L/val;

			if(h==1) gr[2*P+2] = gr[2*P+2] + ((p/(1-p*p) +(zstar1-mean0[0])*(zstar2-mean0[1])/(sqrt(sigma0[0]*sigma0[1])*(1-p*p)) - p*(pow(zstar1-mean0[0],2)/sigma0[0] + pow(zstar2-mean0[1],2)/sigma0[1] - 2*p*(zstar1-mean0[0])*(zstar2-mean0[1])/sqrt(sigma0[0]*sigma0[1]))/pow(1-p*p,2))*L)/val;

		/*if(*trace > 3){
			for(l=0;l<2*P+3;l++){ 
				Rprintf("%f\n", gr[l]);
			}
		}*/
		}		
	}
	
	for(l=0;l<2*P+3;l++){ 
		if(gr[l] != gr[l]){
                        Rprintf("%f %f %f %f %f %f %f %f\n", gr[l], L, val, xx[0], xx[1], xx[2], xx[3], xx[4]);
                }
		gr[l]		=	-gr[l];
	}

}

void adaQuad_vecASE(double *xx, double *nA1, double *nA2, double *n1, double *n2, double *u, double *w, int *dims, int *trace, double* sum, double *geno){
	

	int i,l;
	int N = dims[0];
	int h = dims[1];
	int h2 = dims[2];
	int P = dims[3];
	int nhrule = dims[4];
	*sum=0.0;
	double that[2], xx0[5];


	//Rprintf("vecASE: N %d h %d h2 %d P %d nhrule\n",N, h, h2, P, nhrule);
	//Rprint_v(xx,0,4);
	

	double *exPara, *b1, *b2, *mean0, *sigma0, *x_1, *x_2, hessian[4], Qhat[4],p,grtemp[2*P+3], pinit[2];
	exPara = (double *) Calloc(11+4, double); 

	//be careful of conversion from ints to doubles
	exPara[0] = N;
	exPara[1] = 0;//placeholder, updated in loop below  //P;
	exPara[2] = xx[2*P+2]; //p;
	exPara[3] = 0.0; //nA1 is updated in the loop below with n1
	exPara[4] = 0.0; //same for second one
	exPara[5] = h;
	exPara[6] = h2;
	exPara[7] = 0.0; //pi1 and pi2 are set in loop below
	exPara[8] = 0.0; 
	exPara[9] = 0.0; //n1 is updated in the loop below with nA1
	exPara[10] = 0.0; //same for second one
	mean0  = exPara + 11; 
	sigma0 = mean0 + 2; 

  int npara   = 2; 
  int lmm     = 5; 
  int fail    = 0;
  int failA   = 0;
  int failB   = 0;
  int fncount = 0;//?
  int grcount = 0;//?
  int maxit   = 1000;
  int nREPORT = 5;
			
	int nbd[npara];

	//technical parameters below:
		  double *wa, *g1;
		  int *iwa;
		  SEXP xx1;
		  PROTECT(xx1 = allocVector(REALSXP, npara));
		  //consider replacing with simple Calloc
		  wa  = (double *) Calloc(2*lmm*npara+4*npara+11*lmm*lmm+8*lmm, double);
		  iwa = (int*) Calloc(3*npara, int);
		  g1 = (double *) Calloc(npara, double);
		
		  double gr0[npara];
		  double initPara[npara];
		  double lower[npara];
		  double upper[npara];
		  double Fmin, factr, pgtol, th0, th1, th01;
		  double cov;	  
		  factr = 1e7;
			pgtol = 0.0;
		  
			char msg[1023];
		  
		  for(l=0; l < npara; l++){
					nbd[l] = 2;			  	
					initPara[l] = 0;
					lower[l] = -10;
					upper[l] = 10; 
			}

	for(l=0; l<2;l++){
		mean0[l] = 0;///m -xx[2*P+l]/2;
		sigma0[l] = xx[2*P+l];
	}

	//loop over each observation, but first save values of pi1 and pi2
	pinit[0] = xx[0];
	pinit[1] = xx[1];
	

	for(i=0; i<5;i++) xx0[i] = xx[i];

	for(i=0;i<N;i++){
			//if no AS reads at this site, then skip
			if(n1[i] <1 | n2[i] < 1) continue;
		
			//h2 may need to be two dimensional so can test the hyp separately for each genotype
			//pi1 and pi2 are also set to 0.5 if genotype is not heterozygous, use nA1, nA2 computed from other haplotypes
			//in this case we do not estimate pi1 and pi2, only sigma1, sigma2, and rho
			if(h2==0){
				exPara[7] =  0;
				exPara[8] =  0;
				xx0[0] =  0;
				xx0[1] =  0;
			}else{
				exPara[7] =  pinit[0];
				exPara[8] =  pinit[1];
				xx0[0] =  pinit[0];
				xx0[1] =  pinit[1];
			}

			// update nA, n
        		exPara[1] = geno[i];
			exPara[3] = nA1[i];
			exPara[4] = nA2[i];
			exPara[9] = n1[i];
			exPara[10] = n2[i];

		//for(l=0;l<11+4*P;l++) Rprintf("h %f", exPara[l]);
			
			//l1(npara, initPara, (void*)exPara, xx1);
			//gr(npara, initPara, gr0, (void*)exPara, xx1);
	        //Rprint_v(exPara, 0, 14);
      lbfgsb1(npara, lmm, initPara, lower, upper, nbd, &Fmin, 
           l1ASE, grASE, &failA, (void*)exPara, factr, pgtol,  
           &fncount, &grcount, maxit, msg, 0, nREPORT, wa, iwa, g1,xx1);
	    //if(*trace > 1) Rprintf("%d %s\n", failA, msg);
	    //if (failA) {
	    //  if (*trace)
	        //Rprintf("  i=%d, fail to fit t0 model in iteration\n", i);
	    	//                Rprintf("  i=%d, n1=%f, n2=%f, nA1=%f, nA2=%f fail to fit t0 model in iteration\n", i, n1[i], n2[i], nA1[i], nA2[i]);
	      //continue;
	    //}else{
				//save final estimates of that
				that[0] = initPara[0];
				that[1] = initPara[1];
				//reset initPara to 0
				initPara[0]=initPara[1] = 0.0;
		//	}
			//calculate Qhat
			hessASE(that,xx0,nA1[i],nA2[i],n1[i],n2[i], h, h2, P, hessian, (int) geno[i]);
			inv(hessian,Qhat);
			//for(l=0;l<4;l++){
			//	Rprintf("hess %e Qhat %e\n", hessian[l], Qhat[l]);
			//}
			//calculate value and update sum
			*sum = *sum + adaQuadASE(xx0, u, w, that,Qhat,nA1[i],nA2[i],n1[i],n2[i],dims, (int) geno[i]);
			//if(*trace > 0) Rprintf("%f %d\n",*sum,i);

	}

  Free(exPara);
  Free(wa);
  Free(iwa);
  Free(g1);
  UNPROTECT(1);
}


void adaQuadgr_vecASE(double *xx, double *nA1, double *nA2, double *n1, double *n2, double *u, double *w, int *dims, int *trace, double* gr0, double *geno){
	
	int i,l;
	int N = dims[0];
	int h = dims[1];
	int h2 = dims[2];
	int P = dims[3];
	int nhrule = dims[4];
	double that[2], xx0[5];
	

	//Rprintf("vecgrASE: N %d h %d h2 %d P %d nhrule\n",N, h, h2, P, nhrule);
	//Rprint_v(xx,0,4);
	double *exPara, *b1, *b2, *mean0, *sigma0, *x_1, *x_2, hessian[4], Qhat[4],p,grtemp[2*P+3], pinit[2];
	exPara = (double *) Calloc(11+4, double); 

	//be careful of conversion from ints to doubles
	exPara[0] = N;
	exPara[1] = 0;
	exPara[2] = xx[2*P+2]; //p;
	exPara[3] = 0.0; //nA1 is updated in the loop below with n1
	exPara[4] = 0.0; //same for second one
	exPara[5] = h;
	exPara[6] = h2;
	exPara[7] = xx[0]; //pi1
	exPara[8] = xx[1]; //pi2
	exPara[9] = 0.0; //n1 is updated in the loop below with nA1
	exPara[10] = 0.0; //same for second one
	mean0  = exPara + 11; 
	sigma0 = mean0 + 2; 

  int npara   = 2; 
  int lmm     = 5; 
  int fail    = 0;
  int failA   = 0;
  int failB   = 0;
  int fncount = 0;//?
  int grcount = 0;//?
  int maxit   = 1000;
  int nREPORT = 5;
			
	int nbd[npara];

	//technical parameters below:
		  double *wa, *g1;
		  int *iwa;
		  SEXP xx1;
		  PROTECT(xx1 = allocVector(REALSXP, npara));
		  //consider replacing with simple Calloc
		  wa  = (double *) Calloc(2*lmm*npara+4*npara+11*lmm*lmm+8*lmm, double);
		  iwa = (int*) Calloc(3*npara, int);
		  g1 = (double *) Calloc(npara, double);
		
		  //double gr0[npara];
		  double initPara[npara];
		  double lower[npara];
		  double upper[npara];
		  double Fmin, factr, pgtol, th0, th1, th01;
		  double cov;		  
		  factr = 1e7;
			pgtol = 0.0;
		  
			char msg[1023];
		  
		  for(l=0; l < npara; l++){
					nbd[l] = 0;			  	
					initPara[l] = 2;
					lower[l] = -10;
					upper[l] = 10; 
			}

	for(l=0; l<2;l++){
		mean0[l] = 0;///m -xx[2*P+l]/2;
		sigma0[l] = xx[2*P+l];
	}

	//loop over each observation
	pinit[0] = xx[0];
	pinit[1] = xx[1];

	for(i=0; i<5;i++) xx0[i] = xx[i];

	for(i=0;i<N;i++){
			//if no AS reads at this site, then skip
			if(n1[i] < 1 | n2[i] < 1) continue;
		
                        //h2 may need to be two dimensional so can test the hyp separately for each genotype
                        //pi1 and pi2 are also set to 0.5 if genotype is not heterozygous, use nA1, nA2 computed from other haplotypes
                        //in this case we do not estimate pi1 and pi2, only sigma1, sigma2, and rho
                        if(h2==0 | (geno[i] != 1.0 & geno[i] != 3.0)){
                                exPara[7] =  0;
                                exPara[8] =  0;
                                xx0[0] =  0;
                                xx0[1] =  0;
                        }else{
                              	exPara[7] =  pinit[0];
                                exPara[8] =  pinit[1];
                                xx0[0] =  pinit[0];
                                xx0[1] =  pinit[1];
                        }	

			// update nA, n
                        exPara[1] = geno[i];
			exPara[3] = nA1[i];
			exPara[4] = nA2[i];
			exPara[9] = n1[i];
			exPara[10] = n2[i];

		//for(l=0;l<11+4*P;l++) Rprintf("h %f", exPara[l]);
			//get estimate of that
      lbfgsb1(npara, lmm, initPara, lower, upper, nbd, &Fmin, 
           l1ASE, grASE, &failA, (void*)exPara, factr, pgtol,  
           &fncount, &grcount, maxit, msg, 0, nREPORT, wa, iwa, g1,xx1);
	    //Rprintf("%d %s\n", failA, msg);
	    //if (failA) {
	      //if (*trace)
	        //Rprintf("  i=%d, fail to fit t0 model in iteration\n", i);
               // Rprintf("  i=%d, n1=%d, n2=%d, nA1=%d, nA2=%d fail to fit t0 model in iteration\n", i, n1[i], n2[i], nA1[i], nA2[i]);	    
	      //continue;
	    //}else{
				//save final estimates of that
				that[0] = initPara[0];
				that[1] = initPara[1];
				//reset initPara to 0
				initPara[0]=initPara[1] = 0.0;
	//		}
			//calculate Qhat
			hessASE(that,xx0,nA1[i],nA2[i],n1[i],n2[i], h, h2, P, hessian, (int) geno[i]);
			inv(hessian,Qhat);
			//calculate value and update sum
			adaQuadgrASE(xx0, u, w, that,Qhat,nA1[i],nA2[i],n1[i],n2[i],dims, grtemp, trace, (int) geno[i]); 
			for(l=0;l<2*P+3;l++){
				gr0[l] = gr0[l] + grtemp[l];
			}
			//Rprintf("%f\n",*sum);
			//error("");

		}
	//Rprintf("end adaQuadgr\n");
	//Rprint_v(xx, 0, 4);
	//Rprintf("\n");
  Free(exPara);
  Free(wa);
  Free(iwa);
  Free(g1);

  UNPROTECT(1);
}

//double adaQuad_vec_bfgs_ase(int n, double* x2, void* ex, SEXP x1){
double adaQuad_vec_bfgs_ase(int n, double* x, void* ex, SEXP x1){
	/* double x[5];
	x[0] = x2[0];
	x[1] = x2[1];
	x[2] = exp(x2[2]);
	x[3] = exp(x2[3]);
	x[4] = (exp(2*x2[4])-1)/(exp(2*x2[4])+1);
	*/
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

	/*Rprintf("start adaQuad_vecASE\n");
	Rprint_v(x, 0, 4);
	Rprintf("\n");*/

  exPara = (double *) ex;
	N = (int) exPara[0]; //ceil(exPara[0]-0.5);
	nX = (int) exPara[1];//ceil(exPara[1]-0.5);
	nX2 = (int) exPara[2];//ceil(exPara[2]-0.5);
	y1  = exPara + 3;
	y2 = y1 + N; 
  xx1 = y2 + N;
	pXlast = xx1 + N;
  xx2     = pXlast + N; //add extra column for snp effect
	pXlast2 = xx2 + N;
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
	dims[3] = 1; //assume nX is the same as nX2, account for int and snp column
	dims[4] = (int) ddimsnew[6];

	/*Rprintf("bfgs_ase: \n");
	Rprint_v(x, 0, 4);
	Rprint_v(y1, 0, 5);
	Rprint_v(y2, 0, 5);
	Rprint_v(xx1, 0, 5);
	Rprint_v(xx2, 0, 5);*/
	//for(l=0;l<nX;l++){
	//	Rprintf("XX1 %f XX2 %f \n",xx1[N*l],xx2[N*l]); 
	//}

	//Rprint_v(u, 0, 4);
	//Rprint_v(w, 0, 4);
	//Rprint_vi(dims, 0, 4);


	sum = 0.0;
	//adaQuad_vecASE(double *xx, double *nA1, double *nA2, double *n1, double *n2, double *u, double *w, int *dims, int *trace, double* sum, int *geno)

	adaQuad_vecASE(x, y2, xx2, y1, xx1, u, w, dims, 0, &sum, pXlast);
	//Rprintf("end adaQuad_vecASE\n");
	//Rprint_v(x, 0, 4);
	//Rprintf("\n");
	//Rprintf("sum %f\n", sum);
	//Rprintf("sum %f\n", sum);
	if(sum!=sum){
		return(1000000);
	}else{
		return(sum);
	}
}


//void adaQuadgr_vec_bfgs_ase(int n, double* x2, double* gr0, void* ex, SEXP x1){
void adaQuadgr_vec_bfgs_ase(int n, double* x, double* gr0, void* ex, SEXP x1){
 	/* double x[5];
        x[0] = x2[0];
        x[1] = x2[1];
        x[2] = exp(x2[2]);
        x[3] = exp(x2[3]);
        x[4] = (exp(2*x2[4])-1)/(exp(2*x2[4])+1);
	*/
 
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
	//Rprintf("start adaQuadgr_vecASE\n");
	//Rprint_v(x, 0, 4);
	//Rprintf("\n");
  exPara = (double *) ex;
	N = (int) exPara[0]; //ceil(exPara[0]-0.5);
	nX = (int) exPara[1];//ceil(exPara[1]-0.5);
	nX2 = (int) exPara[2];//ceil(exPara[2]-0.5);
	y1  = exPara + 3;
	y2 = y1 + N; 
  xx1 = y2 + N;
	pXlast = xx1 + N;
  xx2     = pXlast + N; //add extra column for snp effect
	pXlast2 = xx2 + N;
	ddimsnew = pXlast2 + N; 
	u = ddimsnew + 10;
	w = u + (int) ddimsnew[6];

	dims[0] = (int) N;
	dims[1] = (int) ddimsnew[7];
	dims[2] = (int) ddimsnew[8];
	dims[3] = 1; //assume nX is the same as nX2, account for int and snp column
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
//adaQuadgr_vecASE(double *xx, double *nA1, double *nA2, double *n1, double *n2, double *u, double *w, int *dims, int *trace, double* gr0, int *geno)
	adaQuadgr_vecASE(x, y2, xx2, y1, xx1, u, w, dims, 0, gr0, pXlast);

	//Rprintf("end adaQuadgr_vecASE\n");
	//Rprint_v(gr0, 0, 4);
	//Rprint_v(x,0,4);
	//Rprintf("\n");
	for(i=0;i<n;i++){
		if(gr0[i] != gr0[i]){
			Rprintf("na\n");
			gr0[i]=0.0;
		}
	}
}


void glmEQTL_joint_ase (int* dims, double* Y, double* X, double* Z, double* z1, 
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
  int i, j, k, c, m, m2,nIter, df0, df1, convBase, convSNPj, p , pcount, m2FLAG=0;
  double chisq[4], pval[4], phi0, phi1, scale, twologlik0, twologlik1, l0, l1, ll[4];
/////
	int l, c2, nIter2, df02, df12, convBase2, convSNPj2;
  double chisq2, pval2, phi02, phi12, scale2, twologlik02, twologlik12,  logLik0, logLik1;
///
	

  int family = 0;
  int linkR  = *link;
  int adjZR  = *adjZ;
  int npara, lmm, fail, failA, failB, fncount, grcount, nREPORT;
		
  /* pointers to Y, the last column of X, and Z */
  double *pY, *pXlast, *pZ;
/////
  double *pY2, *pXlast2,*pXXlast,*pXXlast2, *YY, *YY2, *XX, *XX2, *ddimsNew, *uu, *ww, *pX, *pX2;
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
  
  int nY = dims[0]; //ntot ASE GE
  int nX = dims[1]; //ntot ASE dnase
  int nZ = dims[2]; //n genetic effects
  int N  = dims[3]; //N samples
  int maxit = dims[4];
  int useOffset = dims[5];
///// adding to dimension of dims to accomodate new vars
	int nY2 = dims[6]; //hap1 GE
  int nX2 = dims[7]; //hap2 dnase
	///need to add number of quad points to dims
	int pts		= dims[8];

  double P_cut = *RP_cut;
  /* dimsNew is passed to function glmNB, and then function glmFit */
  ///int dimsNew[7]; we are not using glmNB here anymore
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

	
	//since glmnb is not used anymore, this can be ignored
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

	int failv[4];

	/* define beta */ //00 and 01 have nX+NX2 length plus 2 for intercept.  Others have addition +2 for snp effect in each
	//double beta00[nX+nX2+2], beta01[nX+nX2+2], beta20[nX+nX2+2+2], beta21[nX+nX2+2+2];

	// not used here
	///double beta00[nX+nX2-2], beta01[nX+nX2-2], beta20[nX+nX2], beta21[nX+nX2];
	

	// here we are simply tacking on the snp column at the end of the X and X2 dnase columns in expara

  /* point to the last column of X */  
	/// these pointers are still used here, but only point to snp effect
	//  X does not correspond to covs anymore
  //pXlast  = X + N*(nX-1);  
	pXlast  = X + N;

/////
  //pXlast2  = X2 + N*(nX2-1);
	pXlast2  = X2 + N;

	/*  These are also not used
	double *X_0;
	X_0 = X + N; //skip intercept

  double *X2_0;
	X2_0 = X2 + N; //skip intercept

	*/

	///temporory vector of length n
	double diff[N];
	int count=0;
  if(*trace){
    Rprintf("\n--------------------------------------------------------------\n");
    Rprintf("(nY, nZ, nX, N, adjZ) = (%d, %d, %d, %d, %d)\t", nY, nZ, nX, N, adjZR);
    Rprintf("P_cut=%e\n", P_cut);
  }
  

	double *exPara;
	//X, X2, 2 Snp cols, Y, Y2, dimsnew (10) + N, nX, and nX2
	///now we also add 2*l for u (l quad points) and w (l weights)
  //exPara = (double *) Calloc(nX*N+nX2*N+ 2*N + 2*N + 10 + 3 , double); 
	
	/// new version	
	//exPara = (double *) Calloc(nX*N+nX2*N+ 2*N + 2*N + 10 + 3 + 2*pts, double); 

	// ase version:  Column for Y, Y1, X, X2, two columns for snp
    exPara = (double *) Calloc(2*N + 2*N + 2*N + 10 + 3 + 2*pts, double); 
    

	exPara[0] = N; //this is going to be primarily used since in innermost loop we only have single columned data for each type
	exPara[1] = nX; //nX is irrelevant down here, need to make use it is used properly as number of DNase sites
	exPara[2] = nX2; //nX2 is the number of Dnase sites
  //YY     = exPara ; filled in later for each Y
  //YY2    = YY + N; filled in later for each Y2
  XX     = exPara + 2*N+3;
	pXXlast = XX + N;  
  XX2     = pXXlast + N; //add extra column for snp effect
	pXXlast2 = XX2 + N;
	ddimsNew = pXXlast2 + N; //add extra column for snp effect
	uu = ddimsNew + 10;
	ww = uu + pts;

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
			npara   = 5; //pi1, pi2, 2 phi, 1 lambda 
		  lmm     = 5; 
		  fail    = 0;
		  failA   = 0;
		  failB   = 0;
		  fncount = 0;//?
		  grcount = 0;//?
		  maxit   = 1000;
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
  fprintf(fo, "GeneRowID\tMarkerRowID\tChisq0\tPvalue0\tChisq1\tPvalue1\tChisq2\tPvalue2\tChisq3\tPvalue3\n");
  
		/*Rprint_v(Y, 0, 5);
		Rprint_v(Y2, 0, 5);
		Rprint_v(X, 0, 5);
		Rprint_v(X2, 0, 5); */

  /***
   * identifify eQTL gene by gene
   */
  
  /* pY is the pointer to total ase gene expression data */
  pY = Y;
	pY2 = Y2;
/////

	Rprintf("ok2 \n");
			
  for(i=0; i<nY; i++,pY+=N,pY2+=N){

    /* Rprintf("just y, %d\n", i);  

		Rprint_v(pY, 0, 5);
		Rprint_v(pY2, 0, 5);*/

    if(*trace == 1){
      if(i%100 == 0){ Rprintf("\ni=%d\n", i); }
    }else if(*trace > 1){
      Rprintf("\ni=%d\n", i);
    }

	pX = X;    
	pX2 = X2;
  for(l=0; l<nX; l++,pX+=N,pX2+=N){


   Rprintf("running, %d, %d %f\n", i ,l, fabs(ePos[i] - dPos[l]));  

		Rprint_v(pY, 0, 5);
		Rprint_v(pY2, 0, 5);
		Rprint_v(pX, 0, 5);
		Rprint_v(pX2, 0, 5); 
 
	//if(l==75 & i==65) error("stopped");
	
	//if gene TSS - DNase distance too great, then skip    
		
				if(*cis_only){
      if(eChr[i] != dChr[l]) continue;
        
      pos_diff2 = fabs(ePos[i] - dPos[l]);
      if(pos_diff2 > *cis_distance2 ) continue; //this may need to be modified later
      Rprintf("starting, %d, %d %f\n", i ,l, fabs(ePos[i] - dPos[l]));  
    }
  	   
      	count = 0;
        for(p=0; p<N;p++){
                if(pY[p] >= 10 & pX[p] >= 10) count++;
        }
        Rprintf("count: %d: %d %d\n", count, i, l); //j should be 0 here
        if(count < 10){
                for(p=0; p<4; p++){
                        if(p==0) fprintf(fo, "%d\t%d\t%d\t", i+1, l+1, 0);
                        fprintf(fo, "%d\t%d\t", -1, -1);
                }
                fprintf(fo, "\n");
                continue;
        }
        Rprintf("Proceeding\n");

		//error("stop"); 
 
  /*
		Rprint_v(pY, 0, 5);
		Rprint_v(pY2, 0, 5);
		Rprint_v(pX, 0, 5);
		Rprint_v(pX2, 0, 5);
 */
    /* **********************************************************
     * fit a baseline model using only the confouding covariates 
     * family is assigned to a value of either Poisson or NB
     * **********************************************************/

		// in this version, just assume nb
    
    //dimsNew[1] = nX-2;
    /** 
     * no initial values, Note this value will be changed 
     * in glmNB after initail iteration 
     */
    //dimsNew[3] = 0; 




		//*trace = 0;
		///here we are using the nb marginals to get initial estimates for beta, since only glmNB has procedure to get actually coefs
		///fit covariate matrix with intercept but without SNP effect. X has no intercept column.  nX also exlcudes snp effect
		///everything after this section must use BFGS since assuming lambda ~ MVN(0, Sigma) 

		/* not initialization is used for ASE
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
		}*
    
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
    */
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
   
   Rprintf("running, %d, %d, %d, ePos %d, dPos %d, mPos %d posdiff2 %d posdiff %d \n", i ,l, j, ePos[i], dPos[l], mPos[j],pos_diff2, pos_diff );

                Rprint_v(pY, 0, 5);
                Rprint_v(pY2, 0, 5);
                Rprint_v(pX, 0, 5);
                Rprint_v(pX2, 0, 5);
		Rprint_v(pZ, 0, 5);
   
      if(*trace > 3){
        Rprintf("\nl=%d i=%d, j=%d\n",l, i, j);
      }
            
      /* *
       * fill the last column of X by Z[j], genotype of one SNP
       * start with the fitted values of the confounder-only model
       */
      
      for (k=0; k<N; k++) {
				//only for now, update this later 
				pXlast[k] = pZ[k];
				//if(pZ[k] == 3.0) pXlast[k] = 1.0; 
				//if(pZ[k] == 4.0) pXlast[k] = 2.0;
				pXXlast[k]  = pXlast[k]; //for lbfgs
        //fitted1[k] = fitted0[k];
/////
				pXlast2[k]  = pXlast[k];
				pXXlast2[k]  = pXlast[k]; //for lbfgs
        //fitted12[k] = fitted02[k];
/////
      }
      
      //phi1  = phi0;
      scale = 1.0;
/////
      //phi12  = phi02;
      scale2 = 1.0;
/////
   	count = 0;
	for (k=0; k<N; k++) {
		if(pXlast[k]==1 | pXlast[k]==3) count++;
	}
	Rprintf("hetcount: %d ", count);
   	count =	0;
        for (k=0; k<N; k++) {
                if((pXlast[k]==1 | pXlast[k]==3) & pY[k]>=10 & pX[k]>=10) count++;
        }
	if(count<0){
		fprintf(fo, "%d\t%d\t%d\t", i+1, l+1, j+1);
                fprintf(fo, "-2\t-2\t-2\t-2\t-2\t-2\t-2\t-2\n");
		continue;
	}
	Rprintf("hetand10count: %d\n", count);
   
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
                        	if(pZ[p] == 3.0){
                                	//YY2[p] = YY[p] - YY2[p];
                        	}
			}

		//copy X, X2 without snp effect now to locations	
		for(p=0; p<N; p++){
			XX[p] = pX[p];
			XX2[p] = pX2[p];
			if(pZ[p] == 3.0){
                                //XX2[p] = XX[p] - XX2[p];
                        }	
		}


			///update dimsNew
			dimsNew[1] = (double) nX; /* assuming adjZ = FALSE */
			dimsNew[5] = (double) nX2; /* assuming adjZ = FALSE */			
			dimsNew[7]   =  (double) 0; // we will update this later in ddimsnew
			dimsNew[8]   =  (double) 0; // we will update this later in ddimsnew

			for(p=0; p<10; p++) ddimsNew[p] = dimsNew[p];
 	
			//npara   = nX+nX2+2+2+2+1; //ncovs for each + 2 intercepts, 2 snp effects, 2 sigma, 1 rho 
		  npara   = 5; //ncovs for each + 2 sigma, 1 rho 
		  
   //if(*trace > 2){
   //     Rprintf("copied dims\n");
   //   }
		
        
			/* initPara = c(beta, beta2, phi, phi2, lambda) */
			// using beta and phi from null snp and lambda initialization
		  for(m=0;m<N;m++) diff[m] = 0.0;
		  for(p=0; p < npara; p++){		
			//reset diff
   		if(p == 0){
					initPara[p] = 0;
					lower[p] = -10;
					upper[p] = 10;
					nbd[p] = 0;
				}else if(p == 1){
					initPara[p] = 0; //snp effect is 0
					lower[p] = -10;
					upper[p] = 10;
					nbd[p] = 0;
				}else if(p == 2){
					///sigma1 log((nA1[geno==1]/n1[geno==1])/(1-nA1[geno==1]/n1[geno==1]))
					m2=0;
					for(m=0;m<N;m++){
						if(YY[m] > 0.0 & YY2[m] > 0.0 & YY2[m] !=YY[m] ){
							//diff[m2] = log( (YY2[m]/YY[m])/(1-(YY2[m]/YY[m])) );
							diff[m2] = log((YY2[m]/YY[m])/(1-(YY2[m]/YY[m]))) ;
							//Rprintf("YY2 %f YY %f diff %f m %d m2 %d\n", YY2[m], YY[m], diff[m2], m, m2); 
							m2++;
						}
					}				
					if(m2>1){
						initPara[p] = var(diff, 0, m2-1, mean(diff, 0, m2-1)); ///this may be problematic
					}else{
						m2FLAG = 1;					
					}					
					nbd[p] = 2;
					lower[p] = .0001;
					upper[p] = 5;
					for(m=0;m<N;m++) diff[m] = 0.0;
				}else if(p == 3){
					///sigma2
					m2=0;
					for(m=0;m<N;m++){
						if( XX[m] > 0.0 & XX2[m] > 0.0 & XX[m] != XX2[m]){
							//diff[m2] = log( (XX2[m]/XX[m])/(1-(XX2[m]/XX[m])) );
							diff[m2] = log((XX2[m]/XX[m])/(1-(XX2[m]/XX[m]))) ;
							//Rprintf("XX2 %f XX %f diff %f m %d m %d\n", XX2[m], XX[m], diff[m2], m, m2);
							m2++;
						}
					}				
					if(m2>1){
						initPara[p] = var(diff, 0, m2-1, mean(diff, 0, m2-1)); ///this may be problematic
					}else{
						m2FLAG = 1;					
					}					
					nbd[p] = 2;
					lower[p] = .0001;
					upper[p] = 5;
					for(m=0;m<N;m++) diff[m] = 0.0;
				}else{
					///rho
					initPara[p] = 0.0;  //moment estimates are insane due to0 small m1 and m2, need better way to initalize
					nbd[p] = 2;
					lower[p] = -1+.0001; 
					upper[p] = 1-.0001; //for rho		
				}
			  //if(*trace > 2){
    		Rprintf("p: %d %f\n", p, initPara[p]);
		    //}
			}	
			
			//Rprintf("End init, start lbfgs, npara is %d\n", npara);

			if(m2FLAG == 1){
				m2FLAG = 0;
				continue;
			}

			Rprint_v(initPara, 0, 4);
			Rprint_v(exPara, 0, 5);
			Rprint_v(YY, 0, 5);       			
			Rprint_v(YY2, 0, 5);
			Rprint_v(pXlast, 0, 5);
			Rprint_v(XX, 0, 5);       			
			Rprint_v(XX2, 0, 5);
			Rprint_v(pXlast2, 0, 5);
			Rprint_v(u, 0, 4);
			Rprint_v(uu, 0, 4);
			Rprint_v(w, 0, 4);
			Rprint_v(ww, 0, 4);

			///This part is crucial. Need to pay attention to ddimsNew, where we update h and h1, nothing else changes
			//assumption is that P = nX + 1 = nX2 + 1 
			//for each round of estimation, everything stays the same except h and h1 which get updated 
			//also need to reset initPara after each round, as well as failA
			///this may cause issues here

			//Rprint_v(pXlast2, 0, 20);
		
			for(m=3;m>-1;m--){
				if(m==0){
					ddimsNew[7]   =  (double) 1;
					ddimsNew[8]   =  (double) 1;
					//use initpara2 form previous step
				}else if(m==1){
					ddimsNew[7]   =  (double) 0;
					ddimsNew[8]   =  (double) 1;
					//for(p = 0; p<npara; p++) initPara2[p]=initPara[p]; 
				}else if(m==2){
					ddimsNew[7]   =  (double) 1;
					ddimsNew[8]   =  (double) 0;
					//use initpara2 from m=3
				}else if(m==3){
					ddimsNew[7]   =  (double) 0;
					ddimsNew[8]   =  (double) 0;
                                	//for(p = 0; p<npara; p++) initPara2[p]=initPara[p];
				}
		
for(p = 0; p<npara; p++) initPara2[p]=initPara[p];
      	lbfgsb1(npara, lmm, initPara2, lower, upper, nbd, &Fmin, 
           adaQuad_vec_bfgs_ase, adaQuadgr_vec_bfgs_ase, &failA, (void*)exPara, factr, pgtol,  
           &fncount, &grcount, maxit, msg, 0, nREPORT, wa, iwa, g1,x1);
	      	ll[m] = Fmin;
        	failv[m] = failA;
		for(p=0; p<npara;p++){
			if((initPara2[p] == lower[p] | initPara2[p]== upper[p]) & initPara2[p] !=0.0) failv[m] = -3;
		}
					Rprintf("m %d ll %f coef ", m, Fmin);
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
	failv[i] = failB;
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
					fprintf(fo, "%d\t%.2e\t", failv[p], pval[p]);
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
  //Free(exPara);
  *succeed = 1;
//}
}


///// Retest


void glmEQTL_joint_ase_test (int* dims, double* Y, double* X, double* Z,  
              int* trace, int* succeed,double* Y2, double* X2, 
							 double* z2, double* u, double* w , double* ll/// add in quad points
)
{
  int i, j, k, c, m, m2,nIter, df0, df1, convBase, convSNPj, p , pcount, m2FLAG=0;
  double chisq[4], pval[4], phi0, phi1, scale, twologlik0, twologlik1, l0, l1;
/////
	int l, c2, nIter2, df02, df12, convBase2, convSNPj2;
  double chisq2, pval2, phi02, phi12, scale2, twologlik02, twologlik12,  logLik0, logLik1;
///
	

  int family = 0;
  int npara, lmm, fail, failA, failB, fncount, grcount, nREPORT;
		
  /* pointers to Y, the last column of X, and Z */
  double *pY, *pXlast, *pZ;
/////
  double *pY2, *pXlast2,*pXXlast,*pXXlast2, *YY, *YY2, *XX, *XX2, *ddimsNew, *uu, *ww, *pX, *pX2;
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
  
  int nY = dims[0]; //ntot ASE GE
  int nX = dims[1]; //ntot ASE dnase
  int nZ = dims[2]; //n genetic effects
  int N  = dims[3]; //N samples
  int maxit = dims[4];
  int useOffset = dims[5];
///// adding to dimension of dims to accomodate new vars
	int nY2 = dims[6]; //hap1 GE
  int nX2 = dims[7]; //hap2 dnase
	///need to add number of quad points to dims
	int pts		= dims[8];

  /* dimsNew is passed to function glmNB, and then function glmFit */
  ///int dimsNew[7]; we are not using glmNB here anymore
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

	
	int failv[4];
	pXlast  = X + N;
 	pXlast2  = X2 + N;

		///temporory vector of length n
	double diff[N];
	int count=0;
  if(*trace){
    Rprintf("\n--------------------------------------------------------------\n");
    Rprintf("(nY, nZ, nX, N, adjZ) = (%d, %d, %d, %d, %d)\t", nY, nZ, nX, N, 0);
    Rprintf("P_cut=%e\n", 0);
  }
  

	double *exPara;
  exPara = (double *) Calloc(2*N + 2*N + 2*N + 10 + 3 + 2*pts, double); 
    

	exPara[0] = N; //this is going to be primarily used since in innermost loop we only have single columned data for each type
	exPara[1] = nX; //nX is irrelevant down here, need to make use it is used properly as number of DNase sites
	exPara[2] = nX2; //nX2 is the number of Dnase sites
  //YY     = exPara ; filled in later for each Y
  //YY2    = YY + N; filled in later for each Y2
  XX     = exPara + 2*N+3;
	pXXlast = XX + N;  
  XX2     = pXXlast + N; //add extra column for snp effect
	pXXlast2 = XX2 + N;
	ddimsNew = pXXlast2 + N; //add extra column for snp effect
	uu = ddimsNew + 10;
	ww = uu + pts;

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
			npara   = 5; //pi1, pi2, 2 phi, 1 lambda 
		  lmm     = 5; 
		  fail    = 0;
		  failA   = 0;
		  failB   = 0;
		  fncount = 0;//?
		  grcount = 0;//?
		  maxit   = 1000;
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
  pY = Y;
	pY2 = Y2;
	i=0;
 	pX = X;    
	pX2 = X2;
	l=0;  
  pZ = Z;
  j=0;
                Rprint_v(pY, 0, 5);
                Rprint_v(pY2, 0, 5);
                Rprint_v(pX, 0, 5);
                Rprint_v(pX2, 0, 5);
		Rprint_v(pZ, 0, 5);
   
      
      for (k=0; k<N; k++) {
				pXlast[k] = pXlast2[k] = pZ[k];
				//only for now, update this later 
				if(pZ[k] == 3.0) pXlast[k] = pXlast2[k]  = 1.0; 
				if(pZ[k] == 4.0) pXlast[k] = pXlast2[k]  = 2.0;
				pXXlast[k]  = pXlast[k]; //for lbfgs
        //fitted1[k] = fitted0[k];
/////
				pXXlast[k]  = pXlast[k];
				pXXlast2[k]  = pXlast2[k]; //for lbfgs
        //fitted12[k] = fitted02[k];
/////
      }
      		Rprint_v(pXlast, 0, 5);
				Rprint_v(pXXlast, 0, 5);
      //phi1  = phi0;
      scale = 1.0;
/////
      //phi12  = phi02;
      scale2 = 1.0;
/////
 		
    	///update expara for the current values of Y and Y2 
	    YY = exPara+3;
	    YY2 = exPara + N+3;
			for(p=0; p < N; p++){
				YY[p] = pY[p];
				YY2[p] = pY2[p];
                        	if(pZ[p] == 3.0){
                                	//YY2[p] = YY[p] - YY2[p];
                        	}
			}

		//copy X, X2 without snp effect now to locations	
		for(p=0; p<N; p++){
			XX[p] = pX[p];
			XX2[p] = pX2[p];
			if(pZ[p] == 3.0){
                                //XX2[p] = XX[p] - XX2[p];
                        }	
		}


			///update dimsNew
			dimsNew[1] = (double) nX; /* assuming adjZ = FALSE */
			dimsNew[5] = (double) nX2; /* assuming adjZ = FALSE */			
			dimsNew[7]   =  (double) 0; // we will update this later in ddimsnew
			dimsNew[8]   =  (double) 0; // we will update this later in ddimsnew

			for(p=0; p<10; p++) ddimsNew[p] = dimsNew[p];
 	
		  npara   = 5; //ncovs for each + 2 sigma, 1 rho 
		  
      for(m=0;m<N;m++) diff[m] = 0.0;
		  for(p=0; p < npara; p++){		
			//reset diff
   		if(p == 0){
					initPara[p] = 0;
					lower[p] = -300;
					nbd[p] = 0;
				}else if(p == 1){
					initPara[p] = 0; //snp effect is 0
					lower[p] = -300;
					nbd[p] = 0;
				}else if(p == 2){
					///sigma1 log((nA1[geno==1]/n1[geno==1])/(1-nA1[geno==1]/n1[geno==1]))
					m2=0;
					for(m=0;m<N;m++){
						if(YY[m] > 0.0 & YY2[m] > 0.0 & YY2[m] !=YY[m] ){
							diff[m2] = log( (YY2[m]/YY[m])/(1-(YY2[m]/YY[m])) );
							//Rprintf("YY2 %f YY %f diff %f m %d m2 %d\n", YY2[m], YY[m], diff[m2], m, m2); 
							m2++;
						}
					}				
					if(m2>1){
						initPara[p] = var(diff, 0, m2-1, mean(diff, 0, m2-1)); ///this may be problematic
					}else{
						m2FLAG = 1;					
					}					
					nbd[p] = 2;
					lower[p] = 1e-10;
					upper[p] = 3;
					for(m=0;m<N;m++) diff[m] = 0.0;
				}else if(p == 3){
					///sigma2
					m2=0;
					for(m=0;m<N;m++){
						if( XX[m] > 0.0 & XX2[m] > 0.0 & XX[m] != XX2[m]){
							diff[m2] = log( (XX2[m]/XX[m])/(1-(XX2[m]/XX[m])) );
							//Rprintf("XX2 %f XX %f diff %f m %d m %d\n", XX2[m], XX[m], diff[m2], m, m2);
							m2++;
						}
					}				
					if(m2>1){
						initPara[p] = var(diff, 0, m2-1, mean(diff, 0, m2-1)); ///this may be problematic
					}else{
						m2FLAG = 1;					
					}					
					nbd[p] = 1;
					lower[p] = 1e-10;
					upper[p] = 3;
					for(m=0;m<N;m++) diff[m] = 0.0;
				}else{
					///rho
					initPara[p] = 0.0;  //moment estimates are insane due to0 small m1 and m2, need better way to initalize
					nbd[p] = 2;
					lower[p] = -1.0+1e-10; 
					upper[npara-1] = 1-1e-10; //for rho		
				}
			  //if(*trace > 2){
    		Rprintf("p: %d %f\n", p, initPara[p]);
		    //}
			}	
			
			//Rprintf("End init, start lbfgs, npara is %d\n", npara);
			Rprint_v(initPara, 0, 4);
			Rprint_v(exPara, 0, 5);
			Rprint_v(YY, 0, 5);       			
			Rprint_v(YY2, 0, 5);
			Rprint_v(pXlast, 0, 5);
			Rprint_v(XX, 0, 5);       			
			Rprint_v(XX2, 0, 5);
			Rprint_v(pXlast2, 0, 5);
			Rprint_v(u, 0, 4);
			Rprint_v(uu, 0, 4);
			Rprint_v(w, 0, 4);
			Rprint_v(ww, 0, 4);

			///This part is crucial. Need to pay attention to ddimsNew, where we update h and h1, nothing else changes
			//assumption is that P = nX + 1 = nX2 + 1 
			//for each round of estimation, everything stays the same except h and h1 which get updated 
			//also need to reset initPara after each round, as well as failA
			///this may cause issues here

			//Rprint_v(pXlast2, 0, 20);
		
			for(m=3;m>-1;m--){
				if(m==0){
					ddimsNew[7]   =  (double) 1;
					ddimsNew[8]   =  (double) 1;
					//use initpara2 form previous step
				}else if(m==1){
					ddimsNew[7]   =  (double) 0;
					ddimsNew[8]   =  (double) 1;
					for(p = 0; p<npara; p++) initPara2[p]=initPara[p]; 
				}else if(m==2){
					ddimsNew[7]   =  (double) 1;
					ddimsNew[8]   =  (double) 0;
					//use initpara2 from m=3
				}else if(m==3){
					ddimsNew[7]   =  (double) 0;
					ddimsNew[8]   =  (double) 0;
                                	for(p = 0; p<npara; p++) initPara2[p]=initPara[p];
				}
		
for(p = 0; p<npara; p++) initPara2[p]=initPara[p];
      	lbfgsb1(npara, lmm, initPara2, lower, upper, nbd, &Fmin, 
           adaQuad_vec_bfgs_ase, adaQuadgr_vec_bfgs_ase, &failA, (void*)exPara, factr, pgtol,  
           &fncount, &grcount, maxit, msg, 0, nREPORT, wa, iwa, g1,x1);
	      	ll[m] = Fmin;
        	failv[m] = failA;
		for(p=0; p<npara;p++){
			if((initPara2[p] == lower[p] | initPara2[p]== upper[p]) & initPara2[p] !=0.0) failv[m] = -1;
		}
					Rprintf("m %d ll %f coef ", m, Fmin);
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
			
			chisq[0] = 2*(ll[1]-ll[0]); //cor in presence of snp
			chisq[1] = 2*(ll[2]-ll[0]); //snp effect in cor
			chisq[2] = 2*(ll[3]-ll[1]);//snp effect withour corr
			chisq[3] = 2*(ll[3]-ll[2]); //cor effect without snp
			pval[0]  = 1-pchisq(chisq[0], 1.0, 1.0, 0);
			pval[1]  = 1-pchisq(chisq[1], 2.0, 1.0, 0);
			pval[2]  = 1-pchisq(chisq[2], 2.0, 1.0, 0);
			pval[3]  = 1-pchisq(chisq[3], 1.0, 1.0, 0);
		
  Free(exPara);
  *succeed = 1;
//}
}



void adaQuad_vecASE_noquad(double *xx, double *nA1, double *nA2, double *n1, double *n2, double *u, double *w, int *dims, int *trace, double* sum, double *geno, double *uhat){
	

	int i,l;
	int N = dims[0];
	int h = dims[1];
	int h2 = dims[2];
	int P = dims[3];
	int nhrule = dims[4];
	*sum=0.0;
	double that[2], xx0[5];


	//Rprintf("vecASE: N %d h %d h2 %d P %d nhrule\n",N, h, h2, P, nhrule);
	//Rprint_v(xx,0,4);
	

	double *exPara, *b1, *b2, *mean0, *sigma0, *x_1, *x_2, hessian[4], Qhat[4],p,grtemp[2*P+3], pinit[2];
	exPara = (double *) Calloc(11+4, double); 

	//be careful of conversion from ints to doubles
	exPara[0] = N;
	exPara[1] = P;
	exPara[2] = xx[2*P+2]; //p;
	exPara[3] = 0.0; //nA1 is updated in the loop below with n1
	exPara[4] = 0.0; //same for second one
	exPara[5] = h;
	exPara[6] = h2;
	exPara[7] = 0.0; //pi1 and pi2 are set in loop below
	exPara[8] = 0.0; 
	exPara[9] = 0.0; //n1 is updated in the loop below with nA1
	exPara[10] = 0.0; //same for second one
	mean0  = exPara + 11; 
	sigma0 = mean0 + 2; 

  int npara   = 2; 
  int lmm     = 5; 
  int fail    = 0;
  int failA   = 0;
  int failB   = 0;
  int fncount = 0;//?
  int grcount = 0;//?
  int maxit   = 1000;
  int nREPORT = 5;
			
	int nbd[npara];

	//technical parameters below:
		  double *wa, *g1;
		  int *iwa;
		  SEXP xx1;
		  PROTECT(xx1 = allocVector(REALSXP, npara));
		  //consider replacing with simple Calloc
		  wa  = (double *) S_alloc(2*lmm*npara+4*npara+11*lmm*lmm+8*lmm,sizeof(double));
		  iwa = (int*) R_alloc(3*npara,sizeof(int));
		  g1 = (double *)R_alloc(npara, sizeof(double));
		
		  double gr0[npara];
		  double initPara[npara];
		  double lower[npara];
		  double upper[npara];
		  double Fmin, factr, pgtol, th0, th1, th01;
		  
		  factr = 1e7;
			pgtol = 0.0;
		  
			char msg[1023];
		  
		  for(l=0; l < npara; l++){
					nbd[l] = 0;			  	
					initPara[l] = 0;
					lower[l] = -1000;
					upper[l] = 1000; 
			}

	for(l=0; l<2;l++){
		mean0[l] = 0;///m -xx[2*P+l]/2;
		sigma0[l] = xx[2*P+l];
	}

	//loop over each observation, but first save values of pi1 and pi2
	pinit[0] = xx[0];
	pinit[1] = xx[1];


	for(i=0; i<5;i++) xx0[i] = xx[i];

	for(i=0;i<N;i++){
			//if no AS reads at this site, then skip
			if(n1[i] <1 | n2[i] < 1) continue;
		
			//h2 may need to be two dimensional so can test the hyp separately for each genotype
			//pi1 and pi2 are also set to 0.5 if genotype is not heterozygous, use nA1, nA2 computed from other haplotypes
			//in this case we do not estimate pi1 and pi2, only sigma1, sigma2, and rho
			if(h2==0 | geno[i] != 1.0){
				exPara[7] =  0;
				exPara[8] =  0;
				xx0[0] =  0;
				xx0[1] =  0;
			}else{
				exPara[7] =  pinit[0];
				exPara[8] =  pinit[1];
				xx0[0] =  pinit[0];
				xx0[1] =  pinit[1];
			}

			// update nA, n
			exPara[3] = nA1[i];
			exPara[4] = nA2[i];
			exPara[9] = n1[i];
			exPara[10] = n2[i];

		//for(l=0;l<11+4*P;l++) Rprintf("h %f", exPara[l]);
			
			//l1(npara, initPara, (void*)exPara, xx1);
			//gr(npara, initPara, gr0, (void*)exPara, xx1);
	        //Rprint_v(exPara, 0, 14);
      //lbfgsb1(npara, lmm, initPara, lower, upper, nbd, &Fmin, 
      //     l1ASE, grASE, &failA, (void*)exPara, factr, pgtol,  
      //     &fncount, &grcount, maxit, msg, 0, nREPORT, wa, iwa, g1,xx1);
	    //if(*trace > 1) Rprintf("%d %s\n", failA, msg);
	    //if (failA) {
	    //  if (*trace)
	        //Rprintf("  i=%d, fail to fit t0 model in iteration\n", i);
	    	//                Rprintf("  i=%d, n1=%f, n2=%f, nA1=%f, nA2=%f fail to fit t0 model in iteration\n", i, n1[i], n2[i], nA1[i], nA2[i]);
	      //continue;
	    //}else{
				//save final estimates of that
				that[0] = uhat[i];//initPara[0];
				that[1] = uhat[i+N];//initPara[1];
				//reset initPara to 0
				//initPara[0]=initPara[1] = 0.0;
			//}
			//calculate Qhat
			hessASE(that,xx0,nA1[i],nA2[i],n1[i],n2[i], h, h2, P, hessian, (int) geno[i]);
			inv(hessian,Qhat);
			//for(l=0;l<4;l++){
			//	Rprintf("hess %e Qhat %e\n", hessian[l], Qhat[l]);
			//}
			//calculate value and update sum
			*sum = *sum + adaQuadASE(xx0, u, w, that,Qhat,nA1[i],nA2[i],n1[i],n2[i],dims, (int) geno[i]);
			//if(*trace > 0) Rprintf("%f %d\n",*sum,i);

	}

  Free(exPara);
  UNPROTECT(1);
}


void adaQuad_max(double *xx, double *nA1, double *nA2, double *n1, double *n2, double *u, double *w, int *dims, int *trace, double* sum, double *geno, double *uhat){
	

	int i,l;
	int N = dims[0];
	int h = dims[1];
	int h2 = dims[2];
	int P = dims[3];
	int nhrule = dims[4];
	*sum=0.0;
	double that[2], xx0[5];


	//Rprintf("vecASE: N %d h %d h2 %d P %d nhrule\n",N, h, h2, P, nhrule);
	//Rprint_v(xx,0,4);
	

	double *exPara, *b1, *b2, *mean0, *sigma0, *x_1, *x_2, hessian[4], Qhat[4],p,grtemp[2*P+3], pinit[2];
	exPara = (double *) Calloc(11+4, double); 

	//be careful of conversion from ints to doubles
	exPara[0] = N;
	exPara[1] = P;
	exPara[2] = xx[2*P+2]; //p;
	exPara[3] = 0.0; //nA1 is updated in the loop below with n1
	exPara[4] = 0.0; //same for second one
	exPara[5] = h;
	exPara[6] = h2;
	exPara[7] = 0.0; //pi1 and pi2 are set in loop below
	exPara[8] = 0.0; 
	exPara[9] = 0.0; //n1 is updated in the loop below with nA1
	exPara[10] = 0.0; //same for second one
	mean0  = exPara + 11; 
	sigma0 = mean0 + 2; 

  int npara   = 2; 
  int lmm     = 5; 
  int fail    = 0;
  int failA   = 0;
  int failB   = 0;
  int fncount = 0;//?
  int grcount = 0;//?
  int maxit   = 1000;
  int nREPORT = 5;
			
	int nbd[npara];

	//technical parameters below:
		  double *wa, *g1;
		  int *iwa;
		  SEXP xx1;
		  PROTECT(xx1 = allocVector(REALSXP, npara));
		  //consider replacing with simple Calloc
		  wa  = (double *) S_alloc(2*lmm*npara+4*npara+11*lmm*lmm+8*lmm,sizeof(double));
		  iwa = (int*) R_alloc(3*npara,sizeof(int));
		  g1 = (double *)R_alloc(npara, sizeof(double));
		
		  double gr0[npara];
		  double initPara[npara];
		  double lower[npara];
		  double upper[npara];
		  double Fmin, factr, pgtol, th0, th1, th01;
		  
		  factr = 1e7;
			pgtol = 0.0;
		  
			char msg[1023];
		  
		  for(l=0; l < npara; l++){
					nbd[l] = 0;			  	
					initPara[l] = 0;
					lower[l] = -1000;
					upper[l] = 1000; 
			}

	for(l=0; l<2;l++){
		mean0[l] = 0;///m -xx[2*P+l]/2;
		sigma0[l] = xx[2*P+l];
	}

	//loop over each observation, but first save values of pi1 and pi2
	pinit[0] = xx[0];
	pinit[1] = xx[1];


	for(i=0; i<5;i++) xx0[i] = xx[i];

	for(i=0;i<N;i++){
			//if no AS reads at this site, then skip
			if(n1[i] <1 | n2[i] < 1) continue;
		
			//h2 may need to be two dimensional so can test the hyp separately for each genotype
			//pi1 and pi2 are also set to 0.5 if genotype is not heterozygous, use nA1, nA2 computed from other haplotypes
			//in this case we do not estimate pi1 and pi2, only sigma1, sigma2, and rho
			if(h2==0 | geno[i] != 1.0){
				exPara[7] =  0;
				exPara[8] =  0;
				xx0[0] =  0;
				xx0[1] =  0;
			}else{
				exPara[7] =  pinit[0];
				exPara[8] =  pinit[1];
				xx0[0] =  pinit[0];
				xx0[1] =  pinit[1];
			}

			// update nA, n
			exPara[3] = nA1[i];
			exPara[4] = nA2[i];
			exPara[9] = n1[i];
			exPara[10] = n2[i];

		//for(l=0;l<11+4*P;l++) Rprintf("h %f", exPara[l]);
			
			//l1(npara, initPara, (void*)exPara, xx1);
			//gr(npara, initPara, gr0, (void*)exPara, xx1);
	        //Rprint_v(exPara, 0, 14);
      lbfgsb1(npara, lmm, initPara, lower, upper, nbd, &Fmin, 
           l1ASE, grASE, &failA, (void*)exPara, factr, pgtol,  
           &fncount, &grcount, maxit, msg, 0, nREPORT, wa, iwa, g1,xx1);
	    //if(*trace > 1) Rprintf("%d %s\n", failA, msg);
	    if (failA) {
	    //  if (*trace)
	        Rprintf("  i=%d, fail to fit t0 model in iteration\n", i);
	    	                Rprintf("  i=%d, n1=%f, n2=%f, nA1=%f, nA2=%f fail to fit t0 model in iteration\n", i, n1[i], n2[i], nA1[i], nA2[i]);
	      //continue;
	    }else{
				//save final estimates of that
				uhat[i] = initPara[0];
				uhat[i+N] = initPara[1];
				//reset initPara to 0
				initPara[0]=initPara[1] = 0.0;
			}
			//calculate Qhat
			//hessASE(that,xx0,nA1[i],nA2[i],n1[i],n2[i], h, h2, P, hessian);
			//inv(hessian,Qhat);
			//for(l=0;l<4;l++){
			//	Rprintf("hess %e Qhat %e\n", hessian[l], Qhat[l]);
			//}
			//calculate value and update sum
			//*sum = *sum + adaQuadASE(xx0, u, w, that,Qhat,nA1[i],nA2[i],n1[i],n2[i],dims);
			//if(*trace > 0) Rprintf("%f %d\n",*sum,i);

	}

  Free(exPara);
  UNPROTECT(1);
}



double l1ASE2(int n, double* para, void* ex, SEXP x1){
			int i, k, N, P, h, h2, l, nA1, nA2, n1, n2, geno;
		double *exPara, *mean0, sigma0[2], *sigma00, *t,p,val, pi1, pi2, cov;
		double l1, l2;
		
    //sigma is on log scale here and is not variance, big difference
    
		exPara = (double *) ex;
		N      = exPara[0];
		geno   = exPara[1];
		p      = exPara[2];
		nA1    = (int) exPara[3];
		nA2    = (int) exPara[4];
		h 	 	 = exPara[5];
		h2 	 	 = exPara[6];
		pi1    = exPara[7];
		pi2    = exPara[8];
		n1		 = (int) exPara[9];
		n2	 	 = (int) exPara[10];
		mean0  = exPara + 11; 
		sigma00 = mean0 + 2; 
		t      = para;

		//for(l=0;l<11+4*P;l++) Rprintf("h %f", exPara[l]);

    sigma0[0] = exp(sigma00[0]);
    sigma0[1] = exp(sigma00[1]);    
    p = tanh(p);

		if(h ==0 )  p = 0.0;
		if(h2 ==0)  pi1 = pi2 = 0.0;

	  if(geno==1){
      cov = 1.0;
    }else if(geno==3){
      cov = -1.0;
    }else{
      cov = 0.0;
    }

		l1 = exp(pi1*cov + t[0]*sigma0[0])/(1+exp(pi1*cov + t[0]*sigma0[0]));
		l2 = exp(pi2*cov + sigma0[1]*(t[0]*p + t[1]*sqrt(1-p*p)))/(1+exp(pi2*cov + sigma0[1]*(t[0]*p + t[1]*sqrt(1-p*p))));

    //check for non-finite values in l1
		if(l1!=l1){
			if(pi1*cov+t[0]*sigma0[0] > 0.0){
				l1 = 1.0;
			}else{
				l1 = 0;
			}
		}

  //check for non-finite values in l2
		if(l2!=l2){
			if(pi2*cov + sigma0[1]*(t[0]*p + t[1]*sqrt(1-p*p)) > 0.0 ){
				l2= 1.0;
			}else{
				l2 = 0;
			}
		}

		val = dbinom(nA1,n1,l1,0)*dbinom(nA2,n2,l2,0)*dnorm(t[0], 0.0, 1.0, 0)*dnorm(t[1], 0, 1, 0);
 	//Rprintf("l1 val %e l1 %f l2 %f t1 %f t2 %f %f %f %f %f %f %e\n", -log(val + pow(10,-100)), l1 , l2, t[0], t[1], mean0[0], mean0[1],sigma0[0], sigma0[1],p, dnorm(t[0], mean0[0] + p*sqrt(sigma0[0])*(t[1]- mean0[1])/sqrt(sigma0[1]), sqrt(sigma0[0]*(1-p*p)),0)*dnorm(t[1], mean0[1], sqrt(sigma0[1]),0));
  return(-log(val + pow(10,-100)));
}

void grASE2(int n, double* para, double* gr0, void* ex, SEXP x1)
{
		int i, k, N, P, h, h2, l, nA1, nA2, n1, n2, geno;
		double *exPara, *mean0, *sigma00,sigma0[2], *t,p,val, pi1, pi2, dt1, dt2, cov;
		double l1, l2;
		
		exPara = (double *) ex;
		N      = exPara[0];
		geno   = exPara[1];
		p    	 = exPara[2];
		nA1    = (int) exPara[3];
		nA2    = (int) exPara[4];
		h 	 	 = exPara[5];
		h2 	 	 = exPara[6];
		pi1    = exPara[7];
		pi2    = exPara[8];
		n1		 = (int) exPara[9];
		n2	 	 = (int) exPara[10];
		mean0  = exPara + 11; 
		sigma00 = mean0 + 2; 
		t      = para;

		//for(l=0;l<11+4*P;l++) Rprintf("h %f", exPara[l]);
    sigma0[0] = exp(sigma00[0]);
    sigma0[1] = exp(sigma00[1]);    
    p = tanh(p);
    
		if(h ==0 )  p = 0.0;
  	if(h2 ==0)  pi1 = pi2 = 0.0;

	  if(geno==1){
      cov = 1.0;
    }else if(geno==3){
      cov = -1.0;
    }else{
      cov = 0.0;
    }

		l1 = exp(pi1*cov + t[0]*sigma0[0])/(1+exp(pi1*cov + t[0]*sigma0[0]));
		l2 = exp(pi2*cov + sigma0[1]*(t[0]*p + t[1]*sqrt(1-p*p)))/(1+exp(pi2*cov + sigma0[1]*(t[0]*p + t[1]*sqrt(1-p*p))));

    //check for non-finite values in l1
		if(l1!=l1){
			if(pi1*cov+t[0]*sigma0[0] > 0.0){
				l1 = 1.0;
			}else{
				l1 = 0.0;
			}
		}

  //check for non-finite values in l2
		if(l2!=l2){
			if(pi2*cov + sigma0[1]*(t[0]*p + t[1]*sqrt(1-p*p)) > 0.0 ){
				l2= 1.0;
			}else{
				l2 = 0.0;
			}
		}
         
  	dt1 = sigma0[0]*(nA1-n1*l1) + sigma0[1]*p*(nA2 - n2*l2) - t[0];
		dt2 = sigma0[1]*sqrt(1-p*p)*(nA2 - n2*l2) - t[1];
  	gr0[0] = -dt1;
		gr0[1] = -dt2;
		//Rprintf("gr l1 %f l2 %f t1 %f t2 %f dtq %f dt2 %f\n",l1 , l2, t[0], t[1], dt1, dt2);
		//error("");
}

void hessASE2(double* t, double *xx, double nA1, double nA2, double n1, double n2, int h, int h2, int P, double *result, int geno){

	//P = 1 here for ASE
	int l;
	double pi1, pi2, sigma0[2], p, l1, l2, ddt1, ddt2, ddt12, cov;

  
  pi1 = xx[P-1];
  pi2 = xx[2*P-1];
	sigma0[0] = exp(xx[2*P]);
  sigma0[1] = exp(xx[2*P+1]);
	p  = tanh(xx[2*P+2]);

  if(h ==0 )  p = 0.0;
  if(h2 ==0)  pi1 = pi2 = 0.0;

	  if(geno==1){
      cov = 1.0;
    }else if(geno==3){
      cov = -1.0;
    }else{
      cov = 0.0;
    }

		l1 = exp(pi1*cov + t[0]*sigma0[0])/(1+exp(pi1*cov + t[0]*sigma0[0]));
		l2 = exp(pi2*cov + sigma0[1]*(t[0]*p + t[1]*sqrt(1-p*p)))/(1+exp(pi2*cov + sigma0[1]*(t[0]*p + t[1]*sqrt(1-p*p))));

    //check for non-finite values in l1
		if(l1!=l1){
			if(pi1*cov+t[0]*sigma0[0] > 0.0){
				l1 = 1.0;
			}else{
				l1 = 0;
			}
		}

  //check for non-finite values in l2
		if(l2!=l2){
			if(pi2*cov + sigma0[1]*(t[0]*p + t[1]*sqrt(1-p*p)) > 0.0 ){
				l2= 1.0;
			}else{
				l2 = 0;
			}
		}
    
  ddt1 = n1*sigma0[0]*sigma0[0]*l1*(1-l1) + n2*sigma0[1]*sigma0[1]*p*p*l2*(1-l2) + 1;
	ddt2 = n2*sigma0[1]*sigma0[1]*(1-p*p)*l2*(1-l2)  + 1;
	ddt12 = n2*sigma0[1]*sigma0[1]*p*sqrt(1-p*p)*l2*(1-l2);
	//return vector representing 2x2 matrix with columns stacked
	result[0] = ddt1;
	result[1] = ddt12;
	result[2] = ddt12;
	result[3] = ddt2;	
} 

double adaQuadASE2(double* xx, double* u, double* w, double* t,double* Qhat,double nA1, double nA2, double n1, double n2,int* dims, int geno){

	int l, j, k;
	double pi1, pi2, *sigma00, sigma0[2], p, l1, l2, l01,l02,val, wj, wk, zstar1, zstar2;
	double Qhat2[4];
	double mean0[2];
	double cov;

	int N = dims[0];
	int h = dims[1];
	int h2 = dims[2];
	int P = dims[3];
	int nhrule = dims[4];

  pi1 = xx[P-1];
	pi2 = xx[2*P-1];
	sigma00 = xx + 2;
	p  = xx[2*P+2];

  	sigma0[0] = exp(sigma00[0])+pow(10,-7);
  	sigma0[1] = exp(sigma00[1])+pow(10,-7);
  	p = tanh(p)-copysign(1,p)*pow(10,-6);
  
	if(h ==0 )  p = 0.0;
  if(h2 ==0)  pi1 = pi2 = 0.0;

	  if(geno==1){
      cov = 1.0;
    }else if(geno==3){
      cov = -1.0;
    }else{
      cov = 0.0;
    }

		
    //check for non-finite values in l1
		if(l1!=l1){
			if(pi1*cov+t[0]*sigma0[0] > 0.0){
				l1 = 1.0;
			}else{
				l1 = 0;
			}
		}

  //check for non-finite values in l2
		if(l2!=l2){
			if(pi2*cov + sigma0[1]*(t[0]*p + t[1]*sqrt(1-p*p)) > 0.0 ){
				l2= 1.0;
			}else{
				l2 = 0;
			}
		}
    
	//computes the expression zstar = that + sqrt(2)Qhat^1/2z for single observation
	chol(Qhat, Qhat2);
	
	//for(l = 0; l<4;l++) Rprintf("adaquad  qhat2 %e \n", Qhat2[l]);

	//Rprintf("%f %f %f %f %f %f\n", Qhat2[0],Qhat2[1],Qhat2[2],Qhat2[3], det(Qhat),mean0[1]);
	val=0.0;
	for(j=0; j<nhrule; j++){
		for(k=0; k<nhrule; k++){
			wj = w[j];
			wk = w[k];
			zstar1 = t[0] + sqrt(2)*(Qhat2[0]*u[j] + Qhat2[1]*u[k]); //this may be another issue in old code
			zstar2 = t[1] + sqrt(2)*(Qhat2[2]*u[j] + Qhat2[3]*u[k]); //should be chol(qhat)%*% z, this doesnt match
			//l01 = exp(pi1*cov + zstar1)/(1+exp(pi1*cov + zstar1));
			//l02 = exp(pi2*cov + zstar2)/(1+exp(pi2*cov + zstar2));
      l01 = exp(pi1*cov + zstar1*sigma0[0])/(1+exp(pi1*cov + zstar1*sigma0[0]));
      l02 = exp(pi2*cov + sigma0[1]*(zstar1*p + zstar2*sqrt(1-p*p)))/(1+exp(pi2*cov + sigma0[1]*(zstar1*p + zstar2*sqrt(1-p*p))));

			if(l01!=l01){
        if(pi1*cov+zstar1*sigma0[0] > 0.0){
          l01 = 1.0;
        }else{
          l01 = 0.0;
        }
      }

      //check for non-finite values in l02
      if(l02!=l02){
        if(pi2*cov + sigma0[1]*(zstar1*p + zstar2*sqrt(1-p*p)) > 0.0 ){
          l02= 1.0;
        }else{
          l02 = 0.0;
        }
      }

	val = val + 2*sqrt(det(Qhat))*wj*wk*dbinom(nA1,n1,l01,0)*dbinom(nA2,n2,l02,0)*dnorm(zstar1, 0.0, 1.0, 0)*dnorm(zstar2, 0.0, 1.0, 0.0)*exp(u[j]*u[j])*exp(u[k]*u[k]);
	}		
}
	//Rprintf("%e final val \n",val);
	return(-log(val+pow(10,-100)));
}


void adaQuad_vecASE2(double *xx, double *nA1, double *nA2, double *n1, double *n2, double *u, double *w, int *dims, int *trace, double* sum, double 
*geno){
	

	int i,l;
	int N = dims[0];
	int h = dims[1];
	int h2 = dims[2];
	int P = dims[3];
	int nhrule = dims[4];
	*sum=0.0;
	double that[2], xx0[5];


	//Rprintf("vecASE: N %d h %d h2 %d P %d nhrule\n",N, h, h2, P, nhrule);
	//Rprint_v(xx,0,4);
	

	double *exPara, *b1, *b2, *mean0, *sigma0, *x_1, *x_2, hessian[4], Qhat[4],p,grtemp[2*P+3], pinit[2];
	exPara = (double *) Calloc(11+4, double); 

	//be careful of conversion from ints to doubles
	exPara[0] = N;
	exPara[1] = 0;//placeholder, updated in loop below  //P;
	exPara[2] = xx[2*P+2]; //p;
	exPara[3] = 0.0; //nA1 is updated in the loop below with n1
	exPara[4] = 0.0; //same for second one
	exPara[5] = h;
	exPara[6] = h2;
	exPara[7] = 0.0; //pi1 and pi2 are set in loop below
	exPara[8] = 0.0; 
	exPara[9] = 0.0; //n1 is updated in the loop below with nA1
	exPara[10] = 0.0; //same for second one
	mean0  = exPara + 11; 
	sigma0 = mean0 + 2; 

  int npara   = 2; 
  int lmm     = 5; 
  int fail    = 0;
  int failA   = 0;
  int failB   = 0;
  int fncount = 0;//?
  int grcount = 0;//?
  int maxit   = 1000;
  int nREPORT = 5;
			
	int nbd[npara];

	//technical parameters below:
		  double *wa, *g1;
		  int *iwa;
		  SEXP xx1;
		  PROTECT(xx1 = allocVector(REALSXP, npara));
		  //consider replacing with simple Calloc
		  wa  = (double *) Calloc(2*lmm*npara+4*npara+11*lmm*lmm+8*lmm, double);
		  iwa = (int*) Calloc(3*npara, int);
		  g1 = (double *) Calloc(npara, double);
		
		  double gr0[npara];
		  double initPara[npara];
		  double lower[npara];
		  double upper[npara];
		  double Fmin, factr, pgtol, th0, th1, th01;
		  double cov;	  
		  factr = 1e7;
			pgtol = 0.0;
		  
			char msg[1023];
		  
		  for(l=0; l < npara; l++){
					nbd[l] = 2;			  	
					initPara[l] = 0;
					lower[l] = -100;
					upper[l] = 100; 
			}

	for(l=0; l<2;l++){
		mean0[l] = 0;///m -xx[2*P+l]/2;
		sigma0[l] = xx[2*P+l];
	}

	//loop over each observation, but first save values of pi1 and pi2
	pinit[0] = xx[0];
	pinit[1] = xx[1];
	

	for(i=0; i<5;i++) xx0[i] = xx[i];

	for(i=0;i<N;i++){
			//if no AS reads at this site, then skip
			if(n1[i] <1 | n2[i] < 1) continue; // may also skip if not het
		
			//h2 may need to be two dimensional so can test the hyp separately for each genotype
			//pi1 and pi2 are also set to 0.5 if genotype is not heterozygous, use nA1, nA2 computed from other haplotypes
			//in this case we do not estimate pi1 and pi2, only sigma1, sigma2, and rho
			if(h2==0){
				exPara[7] =  0;
				exPara[8] =  0;
				xx0[0] =  0;
				xx0[1] =  0;
			}else{
				exPara[7] =  pinit[0];
				exPara[8] =  pinit[1];
				xx0[0] =  pinit[0];
				xx0[1] =  pinit[1];
			}

			// update nA, n
      exPara[1] = geno[i];
			exPara[3] = nA1[i];
			exPara[4] = nA2[i];
			exPara[9] = n1[i];
			exPara[10] = n2[i];

		//for(l=0;l<11+4*P;l++) Rprintf("h %f", exPara[l]);
			
			//l1(npara, initPara, (void*)exPara, xx1);
			//gr(npara, initPara, gr0, (void*)exPara, xx1);
	        //Rprint_v(exPara, 0, 14);
      lbfgsb1(npara, lmm, initPara, lower, upper, nbd, &Fmin, 
           l1ASE2, grASE2, &failA, (void*)exPara, factr, pgtol,  
           &fncount, &grcount, maxit, msg, 0, nREPORT, wa, iwa, g1,xx1);
	    //if(*trace > 1) Rprintf("%d %s\n", failA, msg);
	    //if (failA) {
	    //  if (*trace)
	        //Rprintf("  i=%d, fail to fit t0 model in iteration\n", i);
	    	//                Rprintf("  i=%d, n1=%f, n2=%f, nA1=%f, nA2=%f fail to fit t0 model in iteration\n", i, n1[i], n2[i], nA1[i], nA2[i]);
	      //continue;
	    //}else{
				//save final estimates of that
				that[0] = initPara[0];
				that[1] = initPara[1];
				//reset initPara to 0
				initPara[0]=initPara[1] = 0.0;
		//	}
			//calculate Qhat
			hessASE2(that,xx0,nA1[i],nA2[i],n1[i],n2[i], h, h2, P, hessian, (int) geno[i]);
			inv(hessian,Qhat);
			//for(l=0;l<4;l++){
			//	Rprintf("hess %e Qhat %e\n", hessian[l], Qhat[l]);
			//}
			//calculate value and update sum
			*sum = *sum + adaQuadASE2(xx0, u, w, that,Qhat,nA1[i],nA2[i],n1[i],n2[i],dims, (int) geno[i]);
			//if(*sum != *sum){
			//	Rprintf("%f %d\n",*sum,i);
			//	Rprint_v(xx0, 0, 4);
			//}
			//if(*trace > 0) Rprintf("%f %d\n",*sum,i);

	}

  Free(exPara);
  Free(wa);
  Free(iwa);
  Free(g1);
  UNPROTECT(1);
}

void adaQuadgrASE2(double* xx, double* u, double* w, double* t,double* Qhat,double nA1, double nA2, double n1, double n2,int* dims, double* gr, int 
*trace, int geno){

	
  int l, j, k;
	double pi1, pi2, *sigma00, sigma0[2], p, l1, l2, l01,l02,val, wj, wk, zstar1, zstar2, L;
	double Qhat2[4];
	double mean0[2];
	double cov;

	int N = dims[0];
	int h = dims[1];
	int h2 = dims[2];
	int P = dims[3];
	int nhrule = dims[4];

  	pi1 = xx[P-1];
	pi2 = xx[2*P-1];
	sigma00 = xx + 2;
	p  = xx[2*P+2];

  sigma0[0] = exp(sigma00[0])+pow(10, -7);
  sigma0[1] = exp(sigma00[1])+pow(10, -7);    
  p = tanh(p)-copysign(1, p)*pow(10, -6);
  
	if(h ==0 )  p = 0.0;
  if(h2 ==0)  pi1 = pi2 = 0.0;

	  if(geno==1){
      cov = 1.0;
    }else if(geno==3){
      cov = -1.0;
    }else{
      cov = 0.0;
    }

		
    //check for non-finite values in l1
		if(l1!=l1){
			if(pi1*cov+t[0]*sigma0[0] > 0.0){
				l1 = 1.0;
			}else{
				l1 = 0;
			}
		}

  //check for non-finite values in l2
		if(l2!=l2){
			if(pi2*cov + sigma0[1]*(t[0]*p + t[1]*sqrt(1-p*p)) > 0.0 ){
				l2= 1.0;
			}else{
				l2 = 0;
			}
		}
    
	//computes the expression zstar = that + sqrt(2)Qhat^1/2z for single observation
	chol(Qhat, Qhat2);
	
	//for(l = 0; l<4;l++) Rprintf("adaquad  qhat2 %e \n", Qhat2[l]);
        //Rprintf("pre0: %f %f %f %f %f %f %f\n", t[0], t[1], Qhat[0],Qhat[1],Qhat[2],Qhat[3], det(Qhat));
	//Rprintf("pre: %f %f %f %f %f %f %f\n", t[0], t[1], Qhat2[0],Qhat2[1],Qhat2[2],Qhat2[3], det(Qhat2));
	val=0.0;
	for(j=0; j<nhrule; j++){
		for(k=0; k<nhrule; k++){
			wj = w[j];
			wk = w[k];
			zstar1 = t[0] + sqrt(2)*(Qhat2[0]*u[j] + Qhat2[1]*u[k]); //this may be another issue in old code
			zstar2 = t[1] + sqrt(2)*(Qhat2[2]*u[j] + Qhat2[3]*u[k]); //should be chol(qhat)%*% z, this doesnt match
			//l01 = exp(pi1*cov + zstar1)/(1+exp(pi1*cov + zstar1));
			//l02 = exp(pi2*cov + zstar2)/(1+exp(pi2*cov + zstar2));
      l01 = exp(pi1*cov + zstar1*sigma0[0])/(1+exp(pi1*cov + zstar1*sigma0[0]));
      l02 = exp(pi2*cov + sigma0[1]*(zstar1*p + zstar2*sqrt(1-p*p)))/(1+exp(pi2*cov + sigma0[1]*(zstar1*p + zstar2*sqrt(1-p*p))));

			if(l01!=l01){
        if(pi1*cov+zstar1*sigma0[0] > 0.0){
          l01 = 1.0;
        }else{
          l01 = 0.0;
        }
      }

      //check for non-finite values in l02
      if(l02!=l02){
        if(pi2*cov + sigma0[1]*(zstar1*p + zstar2*sqrt(1-p*p)) > 0.0 ){
          l02= 1.0;
        }else{
          l02 = 0.0;
        }
      }

			val = val + 2*sqrt(det(Qhat))*wj*wk*dbinom(nA1,n1,l01,0)*dbinom(nA2,n2,l02,0)*dnorm(zstar1, 0.0, 1.0, 0)*dnorm(zstar2, 0.0, 
1.0, 0.0)*exp(u[j]*u[j])*exp(u[k]*u[k]);
		}
		if(val!=val) Rprintf("naval: %f %f %f %f %f %f %f %f\n", nA1, n1, l01,nA2,n2, l02,zstar1, zstar2);		
	}
  
	for(l=0;l<2*P+3;l++){ 
		gr[l]		=	0.0;
	}

	for(j=0; j<nhrule; j++){
		for(k=0; k<nhrule; k++){
			wj = w[j];
  		wk = w[k];
			zstar1 = t[0] + sqrt(2)*(Qhat2[0]*u[j] + Qhat2[1]*u[k]); //this may be another issue in old code
			zstar2 = t[1] + sqrt(2)*(Qhat2[2]*u[j] + Qhat2[3]*u[k]); //should be chol(qhat)%*% z, this doesnt match
			//l01 = exp(pi1*cov + zstar1)/(1+exp(pi1*cov + zstar1));
			//l02 = exp(pi2*cov + zstar2)/(1+exp(pi2*cov + zstar2));
      l01 = exp(pi1*cov + zstar1*sigma0[0])/(1+exp(pi1*cov + zstar1*sigma0[0]));
      l02 = exp(pi2*cov + sigma0[1]*(zstar1*p + zstar2*sqrt(1-p*p)))/(1+exp(pi2*cov + sigma0[1]*(zstar1*p + zstar2*sqrt(1-p*p))));

			if(l01!=l01){
        if(pi1*cov+zstar1*sigma0[0] > 0.0){
          l01 = 1.0;
        }else{
          l01 = 0.0;
        }
      }

      //check for non-finite values in l02
      if(l02!=l02){
        if(pi2*cov + sigma0[1]*(zstar1*p + zstar2*sqrt(1-p*p)) > 0.0 ){
          l02= 1.0;
        }else{
          l02 = 0.0;
        }
      }
      
			L = 2*sqrt(det(Qhat))*wj*wk*dbinom(nA1,n1,l01,0)*dbinom(nA2,n2,l02,0)*dnorm(zstar1, 0.0, 1.0, 0)*dnorm(zstar2, 0.0, 1.0, 0.0)*exp(u[j]*u[j])*exp(u[k]*u[k]);
			if(L !=L ) Rprintf("naL: %f %f %f %f %f %f %f %f\n", nA1, n1, l01,nA2,n2, l02,zstar1, zstar2);
			if(h2 ==1 & (geno == 1 | geno == 3)){
				gr[0] = gr[0] + cov*(nA1 - n1*l01)*L/val;
				gr[1] = gr[1] + cov*(nA2 - n2*l02)*L/val;
			}


			gr[2*P] = gr[2*P] + zstar1*(nA1 - n1*l01)*L/val;
      gr[2*P+1] = gr[2*P+1] + (p*zstar1 + sqrt(1-p*p)*zstar2)*(nA2 - n2*l02)*L/val;
      
			if(h==1) gr[2*P+2] = gr[2*P+2] + (sigma0[1]*(zstar1 - zstar2*p/sqrt(1-p*p)))*(nA2 - n2*l02)*L/val;
		//if(*trace > 3){
		
			for(l=0;l<2*P+3;l++){ 
				if(gr[l] != gr[l]){
					Rprintf("nagr: %d %f %f %f %f %f %f %f %f %f \n",l, sigma0[0], sigma0[1], zstar1, zstar2, p, l01, l02, L, val);
				}
			}
			
		//}
		}		
	}
	
	for(l=0;l<2*P+3;l++){
                if(gr[l] != gr[l]){
                        Rprintf("badgr: %f %f %f %f %f %f %f %f\n", gr[l], L, val, xx[0], xx[1], xx[2], xx[3], xx[4]);
                } 
		gr[l]		=	-gr[l];
	}

}

void adaQuadgr_vecASE2(double *xx, double *nA1, double *nA2, double *n1, double *n2, double *u, double *w, int *dims, int *trace, double* gr0, double *geno){
	
	int i,l;
	int N = dims[0];
	int h = dims[1];
	int h2 = dims[2];
	int P = dims[3];
	int nhrule = dims[4];
	double that[2], xx0[5];
	

	//Rprintf("vecgrASE: N %d h %d h2 %d P %d nhrule\n",N, h, h2, P, nhrule);
	//Rprint_v(xx,0,4);
	double *exPara, *b1, *b2, *mean0, *sigma0, *x_1, *x_2, hessian[4], Qhat[4],p,grtemp[2*P+3], pinit[2];
	exPara = (double *) Calloc(11+4, double); 

	//be careful of conversion from ints to doubles
	exPara[0] = N;
	exPara[1] = 0;
	exPara[2] = xx[2*P+2]; //p;
	exPara[3] = 0.0; //nA1 is updated in the loop below with n1
	exPara[4] = 0.0; //same for second one
	exPara[5] = h;
	exPara[6] = h2;
	exPara[7] = xx[0]; //pi1
	exPara[8] = xx[1]; //pi2
	exPara[9] = 0.0; //n1 is updated in the loop below with nA1
	exPara[10] = 0.0; //same for second one
	mean0  = exPara + 11; 
	sigma0 = mean0 + 2; 

  int npara   = 2; 
  int lmm     = 5; 
  int fail    = 0;
  int failA   = 0;
  int failB   = 0;
  int fncount = 0;//?
  int grcount = 0;//?
  int maxit   = 1000;
  int nREPORT = 5;
			
	int nbd[npara];

	//technical parameters below:
		  double *wa, *g1;
		  int *iwa;
		  SEXP xx1;
		  PROTECT(xx1 = allocVector(REALSXP, npara));
		  //consider replacing with simple Calloc
		  wa  = (double *) Calloc(2*lmm*npara+4*npara+11*lmm*lmm+8*lmm, double);
		  iwa = (int*) Calloc(3*npara, int);
		  g1 = (double *) Calloc(npara, double);
		
		  //double gr0[npara];
		  double initPara[npara];
		  double lower[npara];
		  double upper[npara];
		  double Fmin, factr, pgtol, th0, th1, th01;
		  double cov;		  
		  factr = 1e7;
			pgtol = 0.0;
		  
			char msg[1023];
		  
		  for(l=0; l < npara; l++){
					nbd[l] = 0;			  	
					initPara[l] = 2;
					lower[l] = -100;
					upper[l] = 100; 
			}

	for(l=0; l<2;l++){
		mean0[l] = 0;///m -xx[2*P+l]/2;
		sigma0[l] = xx[2*P+l];
	}

	//loop over each observation
	pinit[0] = xx[0];
	pinit[1] = xx[1];

	for(i=0; i<5;i++) xx0[i] = xx[i];

	for(i=0;i<N;i++){
			//if no AS reads at this site, then skip
			if(n1[i] < 1 | n2[i] < 1) continue;
		
                        //h2 may need to be two dimensional so can test the hyp separately for each genotype
                        //pi1 and pi2 are also set to 0.5 if genotype is not heterozygous, use nA1, nA2 computed from other haplotypes
                        //in this case we do not estimate pi1 and pi2, only sigma1, sigma2, and rho
                        if(h2==0){
                                exPara[7] =  0;
                                exPara[8] =  0;
                                xx0[0] =  0;
                                xx0[1] =  0;
                        }else{
                              	exPara[7] =  pinit[0];
                                exPara[8] =  pinit[1];
                                xx0[0] =  pinit[0];
                                xx0[1] =  pinit[1];
                        }	

			// update nA, n
                        exPara[1] = geno[i];
			exPara[3] = nA1[i];
			exPara[4] = nA2[i];
			exPara[9] = n1[i];
			exPara[10] = n2[i];

		//for(l=0;l<11+4*P;l++) Rprintf("h %f", exPara[l]);
			//get estimate of that
      lbfgsb1(npara, lmm, initPara, lower, upper, nbd, &Fmin, 
           l1ASE2, grASE2, &failA, (void*)exPara, factr, pgtol,  
           &fncount, &grcount, maxit, msg, 0, nREPORT, wa, iwa, g1,xx1);
	    //Rprintf("%d %s\n", failA, msg);
	    //if (failA) {
	      //if (*trace)
	        //Rprintf("  i=%d, fail to fit t0 model in iteration\n", i);
               // Rprintf("  i=%d, n1=%d, n2=%d, nA1=%d, nA2=%d fail to fit t0 model in iteration\n", i, n1[i], n2[i], nA1[i], nA2[i]);	    
	      //continue;
	    //}else{
				//save final estimates of that
				that[0] = initPara[0];
				that[1] = initPara[1];
				//reset initPara to 0
				initPara[0]=initPara[1] = 0.0;
	//		}
			//calculate Qhat
			hessASE2(that,xx0,nA1[i],nA2[i],n1[i],n2[i], h, h2, P, hessian, (int) geno[i]);
			inv(hessian,Qhat);
                        //for(l=0;l<4;l++){
                        //	Rprintf("hess %e Qhat %e\n", hessian[l], Qhat[l]);
                        //}
			//calculate value and update sum
			adaQuadgrASE2(xx0, u, w, that,Qhat,nA1[i],nA2[i],n1[i],n2[i],dims, grtemp, trace, (int) geno[i]); 
			for(l=0;l<2*P+3;l++){
				gr0[l] = gr0[l] + grtemp[l];
			}
			//Rprintf("%f\n",*sum);
			//error("");

		}
	//Rprintf("end adaQuadgr\n");
	//Rprint_v(xx, 0, 4);
	//Rprintf("\n");
  Free(exPara);
  Free(wa);
  Free(iwa);
  Free(g1);

  UNPROTECT(1);
}


//double adaQuad_vec_bfgs_ase(int n, double* x2, void* ex, SEXP x1){
double adaQuad_vec_bfgs_ase2(int n, double* x, void* ex, SEXP x1){
	/* double x[5];
	x[0] = x2[0];
	x[1] = x2[1];
	x[2] = exp(x2[2]);
	x[3] = exp(x2[3]);
	x[4] = (exp(2*x2[4])-1)/(exp(2*x2[4])+1);
	*/
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

	/*Rprintf("start adaQuad_vecASE\n");
	Rprint_v(x, 0, 4);
	Rprintf("\n");*/

  exPara = (double *) ex;
	N = (int) exPara[0]; //ceil(exPara[0]-0.5);
	nX = (int) exPara[1];//ceil(exPara[1]-0.5);
	nX2 = (int) exPara[2];//ceil(exPara[2]-0.5);
	y1  = exPara + 3;
	y2 = y1 + N; 
  xx1 = y2 + N;
	pXlast = xx1 + N;
  xx2     = pXlast + N; //add extra column for snp effect
	pXlast2 = xx2 + N;
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
	dims[3] = 1; //assume nX is the same as nX2, account for int and snp column
	dims[4] = (int) ddimsnew[6];

	/*Rprintf("bfgs_ase: \n");
	Rprint_v(x, 0, 4);
	Rprint_v(y1, 0, 5);
	Rprint_v(y2, 0, 5);
	Rprint_v(xx1, 0, 5);
	Rprint_v(xx2, 0, 5);*/
	//for(l=0;l<nX;l++){
	//	Rprintf("XX1 %f XX2 %f \n",xx1[N*l],xx2[N*l]); 
	//}

	//Rprint_v(u, 0, 4);
	//Rprint_v(w, 0, 4);
	//Rprint_vi(dims, 0, 4);


	sum = 0.0;
	//adaQuad_vecASE(double *xx, double *nA1, double *nA2, double *n1, double *n2, double *u, double *w, int *dims, int *trace, double* sum, int *geno)

	adaQuad_vecASE2(x, y2, xx2, y1, xx1, u, w, dims, 0, &sum, pXlast);
	//Rprintf("end adaQuad_vecASE\n");
	//Rprint_v(x, 0, 4);
	//Rprintf("\n");
	//Rprintf("sum %f\n", sum);
	//Rprintf("sum %f\n", sum);
	if(sum!=sum){
		return(1000000);
	}else{
		return(sum);
	}
}


//void adaQuadgr_vec_bfgs_ase(int n, double* x2, double* gr0, void* ex, SEXP x1){
void adaQuadgr_vec_bfgs_ase2(int n, double* x, double* gr0, void* ex, SEXP x1){
 	/* double x[5];
        x[0] = x2[0];
        x[1] = x2[1];
        x[2] = exp(x2[2]);
        x[3] = exp(x2[3]);
        x[4] = (exp(2*x2[4])-1)/(exp(2*x2[4])+1);
	*/
 
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
	//Rprintf("start adaQuadgr_vecASE\n");
	//Rprint_v(x, 0, 4);
	//Rprintf("\n");
  exPara = (double *) ex;
	N = (int) exPara[0]; //ceil(exPara[0]-0.5);
	nX = (int) exPara[1];//ceil(exPara[1]-0.5);
	nX2 = (int) exPara[2];//ceil(exPara[2]-0.5);
	y1  = exPara + 3;
	y2 = y1 + N; 
  xx1 = y2 + N;
	pXlast = xx1 + N;
  xx2     = pXlast + N; //add extra column for snp effect
	pXlast2 = xx2 + N;
	ddimsnew = pXlast2 + N; 
	u = ddimsnew + 10;
	w = u + (int) ddimsnew[6];

	dims[0] = (int) N;
	dims[1] = (int) ddimsnew[7];
	dims[2] = (int) ddimsnew[8];
	dims[3] = 1; //assume nX is the same as nX2, account for int and snp column
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
//adaQuadgr_vecASE(double *xx, double *nA1, double *nA2, double *n1, double *n2, double *u, double *w, int *dims, int *trace, double* gr0, int *geno)
	adaQuadgr_vecASE2(x, y2, xx2, y1, xx1, u, w, dims, 0, gr0, pXlast);

	//Rprintf("end adaQuadgr_vecASE\n");
	//Rprint_v(gr0, 0, 4);
	//Rprint_v(x,0,4);
	//Rprintf("\n");
	for(i=0;i<n;i++){
		if(gr0[i] != gr0[i]){
			Rprintf("na\n");
			gr0[i]=0.0;
		}
	}
}


void glmEQTL_joint_ase2 (int* dims, double* Y, double* X, double* Z, double* z1, 
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
  int i, j, k, c, m, m2,nIter, df0, df1, convBase, convSNPj, p , pcount, m2FLAG=0;
  double chisq[4], pval[4], phi0, phi1, scale, twologlik0, twologlik1, l0, l1, ll[4];
/////
	int l, c2, nIter2, df02, df12, convBase2, convSNPj2;
  double chisq2, pval2, phi02, phi12, scale2, twologlik02, twologlik12,  logLik0, logLik1;
///
	

  int family = 0;
  int linkR  = *link;
  int adjZR  = *adjZ;
  int npara, lmm, fail, failA, failB, fncount, grcount, nREPORT;
		
  /* pointers to Y, the last column of X, and Z */
  double *pY, *pXlast, *pZ;
/////
  double *pY2, *pXlast2,*pXXlast,*pXXlast2, *YY, *YY2, *XX, *XX2, *ddimsNew, *uu, *ww, *pX, *pX2;
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
  
  int nY = dims[0]; //ntot ASE GE
  int nX = dims[1]; //ntot ASE dnase
  int nZ = dims[2]; //n genetic effects
  int N  = dims[3]; //N samples
  int maxit = dims[4];
  int useOffset = dims[5];
///// adding to dimension of dims to accomodate new vars
	int nY2 = dims[6]; //hap1 GE
  int nX2 = dims[7]; //hap2 dnase
	///need to add number of quad points to dims
	int pts		= dims[8];

  double P_cut = *RP_cut;
  /* dimsNew is passed to function glmNB, and then function glmFit */
  ///int dimsNew[7]; we are not using glmNB here anymore
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

	
	//since glmnb is not used anymore, this can be ignored
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

	int failv[4];

	/* define beta */ //00 and 01 have nX+NX2 length plus 2 for intercept.  Others have addition +2 for snp effect in each
	//double beta00[nX+nX2+2], beta01[nX+nX2+2], beta20[nX+nX2+2+2], beta21[nX+nX2+2+2];

	// not used here
	///double beta00[nX+nX2-2], beta01[nX+nX2-2], beta20[nX+nX2], beta21[nX+nX2];
	

	// here we are simply tacking on the snp column at the end of the X and X2 dnase columns in expara

  /* point to the last column of X */  
	/// these pointers are still used here, but only point to snp effect
	//  X does not correspond to covs anymore
  //pXlast  = X + N*(nX-1);  
	pXlast  = X + N;

/////
  //pXlast2  = X2 + N*(nX2-1);
	pXlast2  = X2 + N;

	/*  These are also not used
	double *X_0;
	X_0 = X + N; //skip intercept

  double *X2_0;
	X2_0 = X2 + N; //skip intercept

	*/

	///temporory vector of length n
	double diff[N];
	int count=0;
  if(*trace){
    Rprintf("\n--------------------------------------------------------------\n");
    Rprintf("(nY, nZ, nX, N, adjZ) = (%d, %d, %d, %d, %d)\t", nY, nZ, nX, N, adjZR);
    Rprintf("P_cut=%e\n", P_cut);
  }
  

	double *exPara;
	//X, X2, 2 Snp cols, Y, Y2, dimsnew (10) + N, nX, and nX2
	///now we also add 2*l for u (l quad points) and w (l weights)
  //exPara = (double *) Calloc(nX*N+nX2*N+ 2*N + 2*N + 10 + 3 , double); 
	
	/// new version	
	//exPara = (double *) Calloc(nX*N+nX2*N+ 2*N + 2*N + 10 + 3 + 2*pts, double); 

	// ase version:  Column for Y, Y1, X, X2, two columns for snp
    exPara = (double *) Calloc(2*N + 2*N + 2*N + 10 + 3 + 2*pts, double); 
    

	exPara[0] = N; //this is going to be primarily used since in innermost loop we only have single columned data for each type
	exPara[1] = nX; //nX is irrelevant down here, need to make use it is used properly as number of DNase sites
	exPara[2] = nX2; //nX2 is the number of Dnase sites
  //YY     = exPara ; filled in later for each Y
  //YY2    = YY + N; filled in later for each Y2
  XX     = exPara + 2*N+3;
	pXXlast = XX + N;  
  XX2     = pXXlast + N; //add extra column for snp effect
	pXXlast2 = XX2 + N;
	ddimsNew = pXXlast2 + N; //add extra column for snp effect
	uu = ddimsNew + 10;
	ww = uu + pts;

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
			npara   = 5; //pi1, pi2, 2 phi, 1 lambda 
		  lmm     = 5; 
		  fail    = 0;
		  failA   = 0;
		  failB   = 0;
		  fncount = 0;//?
		  grcount = 0;//?
		  maxit   = 1000;
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
		
	//Rprintf("ok1 \n");
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
  fprintf(fo, "GeneRowID\tMarkerRowID\tChisq0\tPvalue0\tChisq1\tPvalue1\tChisq2\tPvalue2\tChisq3\tPvalue3\n");
  
		/*Rprint_v(Y, 0, 5);
		Rprint_v(Y2, 0, 5);
		Rprint_v(X, 0, 5);
		Rprint_v(X2, 0, 5); */

  /***
   * identifify eQTL gene by gene
   */
  
  /* pY is the pointer to total ase gene expression data */
  pY = Y;
	pY2 = Y2;
/////

	Rprintf("ok2 \n");
			
  for(i=0; i<nY; i++,pY+=N,pY2+=N){

    /* Rprintf("just y, %d\n", i);  

		Rprint_v(pY, 0, 5);
		Rprint_v(pY2, 0, 5);*/

    if(*trace == 1){
      if(i%100 == 0){ Rprintf("\ni=%d\n", i); }
    }else if(*trace > 1){
      Rprintf("\ni=%d\n", i);
    }

	pX = X;    
	pX2 = X2;
  for(l=0; l<nX; l++,pX+=N,pX2+=N){


   Rprintf("running, %d, %d %f\n", i ,l, fabs(ePos[i] - dPos[l]));  

		//Rprint_v(pY, 0, 5);
		//Rprint_v(pY2, 0, 5);
		//Rprint_v(pX, 0, 5);
		//Rprint_v(pX2, 0, 5); 
 
	//if(l==75 & i==65) error("stopped");
	
	//if gene TSS - DNase distance too great, then skip    
		
				if(*cis_only){
      if(eChr[i] != dChr[l]) continue;
        
      pos_diff2 = fabs(ePos[i] - dPos[l]);
      if(pos_diff2 > *cis_distance2 ) continue; //this may need to be modified later
      Rprintf("starting, %d, %d %f\n", i ,l, fabs(ePos[i] - dPos[l]));  
    }
  	   
      	count = 0;
        for(p=0; p<N;p++){
                if(pY[p] >= 10 & pX[p] >= 10) count++;
        }
        Rprintf("count: %d: %d %d\n", count, i, l); //j should be 0 here
        if(count < 10){
                for(p=0; p<4; p++){
                        if(p==0) fprintf(fo, "%d\t%d\t%d\t", i+1, l+1, 0);
                        fprintf(fo, "%d\t%d\t", -1, -1);
                }
                fprintf(fo, "\n");
                continue;
        }
        Rprintf("Proceeding\n");

		//error("stop"); 
 
  /*
		Rprint_v(pY, 0, 5);
		Rprint_v(pY2, 0, 5);
		Rprint_v(pX, 0, 5);
		Rprint_v(pX2, 0, 5);
 */
    /* **********************************************************
     * fit a baseline model using only the confouding covariates 
     * family is assigned to a value of either Poisson or NB
     * **********************************************************/

		// in this version, just assume nb
    
    //dimsNew[1] = nX-2;
    /** 
     * no initial values, Note this value will be changed 
     * in glmNB after initail iteration 
     */
    //dimsNew[3] = 0; 




		//*trace = 0;
		///here we are using the nb marginals to get initial estimates for beta, since only glmNB has procedure to get actually coefs
		///fit covariate matrix with intercept but without SNP effect. X has no intercept column.  nX also exlcudes snp effect
		///everything after this section must use BFGS since assuming lambda ~ MVN(0, Sigma) 

		/* not initialization is used for ASE
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
		}*
    
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
							//Rprintf("%f\n",mean(pY, 0, 59)); 

    }
    */
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
   
   Rprintf("running, %d, %d, %d, ePos %d, dPos %d, mPos %d posdiff2 %d posdiff %d \n", i ,l, j, ePos[i], dPos[l], mPos[j],pos_diff2, pos_diff );

                /*Rprint_v(pY, 0, 5);
                Rprint_v(pY2, 0, 5);
                Rprint_v(pX, 0, 5);
                Rprint_v(pX2, 0, 5);
		Rprint_v(pZ, 0, 5);
   */
      if(*trace > 3){
        Rprintf("\nl=%d i=%d, j=%d\n",l, i, j);
      }
            
      /* *
       * fill the last column of X by Z[j], genotype of one SNP
       * start with the fitted values of the confounder-only model
       */
      
      for (k=0; k<N; k++) {
				//only for now, update this later 
				pXlast[k] = pZ[k];
				//if(pZ[k] == 3.0) pXlast[k] = 1.0; 
				//if(pZ[k] == 4.0) pXlast[k] = 2.0;
				pXXlast[k]  = pXlast[k]; //for lbfgs
        //fitted1[k] = fitted0[k];
/////
				pXlast2[k]  = pXlast[k];
				pXXlast2[k]  = pXlast[k]; //for lbfgs
        //fitted12[k] = fitted02[k];
/////
      }
      
      //phi1  = phi0;
      scale = 1.0;
/////
      //phi12  = phi02;
      scale2 = 1.0;
/////
   	count = 0;
	for (k=0; k<N; k++) {
		if(pXlast[k]==1 | pXlast[k]==3) count++;
	}
	Rprintf("hetcount: %d ", count);
   	count =	0;
        for (k=0; k<N; k++) {
                if((pXlast[k]==1 | pXlast[k]==3) & pY[k]>=10 & pX[k]>=10) count++;
        }
	if(count<0){
		fprintf(fo, "%d\t%d\t%d\t", i+1, l+1, j+1);
                fprintf(fo, "-2\t-2\t-2\t-2\t-2\t-2\t-2\t-2\n");
		continue;
	}
	Rprintf("hetand10count: %d\n", count);
   
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
                        	if(pZ[p] == 3.0){
                                	//YY2[p] = YY[p] - YY2[p];
                        	}
			}

		//copy X, X2 without snp effect now to locations	
		for(p=0; p<N; p++){
			XX[p] = pX[p];
			XX2[p] = pX2[p];
			if(pZ[p] == 3.0){
                                //XX2[p] = XX[p] - XX2[p];
                        }	
		}


			///update dimsNew
			dimsNew[1] = (double) nX; /* assuming adjZ = FALSE */
			dimsNew[5] = (double) nX2; /* assuming adjZ = FALSE */			
			dimsNew[7]   =  (double) 0; // we will update this later in ddimsnew
			dimsNew[8]   =  (double) 0; // we will update this later in ddimsnew

			for(p=0; p<10; p++) ddimsNew[p] = dimsNew[p];
 	
			//npara   = nX+nX2+2+2+2+1; //ncovs for each + 2 intercepts, 2 snp effects, 2 sigma, 1 rho 
		  npara   = 5; //ncovs for each + 2 sigma, 1 rho 
		  
   //if(*trace > 2){
   //     Rprintf("copied dims\n");
   //   }
		
        
			/* initPara = c(beta, beta2, phi, phi2, lambda) */
			// using beta and phi from null snp and lambda initialization
		  for(m=0;m<N;m++) diff[m] = 0.0;
		  for(p=0; p < npara; p++){		
			//reset diff
   		if(p == 0){
					initPara[p] = 0;
					lower[p] = -100;
					upper[p] = 100;
					nbd[p] = 0;
				}else if(p == 1){
					initPara[p] = 0; //snp effect is 0
					lower[p] = -100;
					upper[p] = 100;
					nbd[p] = 0;
				}else if(p == 2){
					///sigma1 log((nA1[geno==1]/n1[geno==1])/(1-nA1[geno==1]/n1[geno==1]))
					m2=0;
					for(m=0;m<N;m++){
						if(YY[m] > 0.0 & YY2[m] > 0.0 & YY2[m] !=YY[m] ){
							//diff[m2] = log( (YY2[m]/YY[m])/(1-(YY2[m]/YY[m])) );
							diff[m2] = log((YY2[m]/YY[m])/(1-(YY2[m]/YY[m]))) ;
							//Rprintf("YY2 %f YY %f diff %f m %d m2 %d\n", YY2[m], YY[m], diff[m2], m, m2); 
							m2++;
						}
					}				
					if(m2>1){
						initPara[p] = log(sqrt(var(diff, 0, m2-1, mean(diff, 0, m2-1)))); ///this may be problematic
					}else{
						m2FLAG = 1;					
					}					
					nbd[p] = 2;
					lower[p] = -100;
					upper[p] = 5;
					for(m=0;m<N;m++) diff[m] = 0.0;
				}else if(p == 3){
					///sigma2
					m2=0;
					for(m=0;m<N;m++){
						if( XX[m] > 0.0 & XX2[m] > 0.0 & XX[m] != XX2[m]){
							//diff[m2] = log( (XX2[m]/XX[m])/(1-(XX2[m]/XX[m])) );
							diff[m2] = log((XX2[m]/XX[m])/(1-(XX2[m]/XX[m]))) ;
							//Rprintf("XX2 %f XX %f diff %f m %d m %d\n", XX2[m], XX[m], diff[m2], m, m2);
							m2++;
						}
					}				
					if(m2>1){
						initPara[p] = log(sqrt(var(diff, 0, m2-1, mean(diff, 0, m2-1)))); ///this may be problematic
					}else{
						m2FLAG = 1;					
					}					
					nbd[p] = 2;
					lower[p] = -100;
					upper[p] = 5;
					for(m=0;m<N;m++) diff[m] = 0.0;
				}else{
					///rho
					initPara[p] = 0.0;  //moment estimates are insane due to0 small m1 and m2, need better way to initalize
					nbd[p] = 2;
					lower[p] = -6; 
					upper[p] = 6; //for rho		
				}
			  //if(*trace > 2){
    		Rprintf("p: %d %f\n", p, initPara[p]);
		    //}
			}	
			
			//Rprintf("End init, start lbfgs, npara is %d\n", npara);

			if(m2FLAG == 1){
				m2FLAG = 0;
				continue;
			}

			/*Rprint_v(initPara, 0, 4);
			Rprint_v(exPara, 0, 5);
			Rprint_v(YY, 0, 5);       			
			Rprint_v(YY2, 0, 5);
			Rprint_v(pXlast, 0, 5);
			Rprint_v(XX, 0, 5);       			
			Rprint_v(XX2, 0, 5);
			Rprint_v(pXlast2, 0, 5);
			Rprint_v(u, 0, 4);
			Rprint_v(uu, 0, 4);
			Rprint_v(w, 0, 4);
			Rprint_v(ww, 0, 4);
*/
			///This part is crucial. Need to pay attention to ddimsNew, where we update h and h1, nothing else changes
			//assumption is that P = nX + 1 = nX2 + 1 
			//for each round of estimation, everything stays the same except h and h1 which get updated 
			//also need to reset initPara after each round, as well as failA
			///this may cause issues here

			//Rprint_v(pXlast2, 0, 20);
		
			for(m=3;m>-1;m--){
				if(m==0){
					ddimsNew[7]   =  (double) 1;
					ddimsNew[8]   =  (double) 1;
					//use initpara2 form previous step
				}else if(m==1){
					ddimsNew[7]   =  (double) 0;
					ddimsNew[8]   =  (double) 1;
					//for(p = 0; p<npara; p++) initPara2[p]=initPara[p]; 
				}else if(m==2){
					ddimsNew[7]   =  (double) 1;
					ddimsNew[8]   =  (double) 0;
					//use initpara2 from m=3
				}else if(m==3){
					ddimsNew[7]   =  (double) 0;
					ddimsNew[8]   =  (double) 0;
                                	//for(p = 0; p<npara; p++) initPara2[p]=initPara[p];
				}
		
for(p = 0; p<npara; p++) initPara2[p]=initPara[p];
      	lbfgsb1(npara, lmm, initPara2, lower, upper, nbd, &Fmin, 
           adaQuad_vec_bfgs_ase2, adaQuadgr_vec_bfgs_ase2, &failA, (void*)exPara, factr, pgtol,  
           &fncount, &grcount, maxit, msg, 0, nREPORT, wa, iwa, g1,x1);
	      	ll[m] = Fmin;
        	failv[m] = failA;
		for(p=0; p<npara;p++){
			if((initPara2[p] == lower[p] | initPara2[p]== upper[p]) & initPara2[p] !=0.0) failv[m] = -3;
		}
					Rprintf("m %d ll %f coef ", m, Fmin);
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
	failv[i] = failB;
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
					fprintf(fo, "%d\t%.2e\t", failv[p], pval[p]);
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
  //Free(exPara);
  *succeed = 1;
//}
}

