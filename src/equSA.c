#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <R.h>
#include <Rmath.h>


/* over definition */

/*** parameter settings for corsel.c and pcorsel.c ***/
static int  MAX_NUM=2000000, GRID=2, mingrid=2;
static int  WARM=1, stepscale=10000, samsize=10000;
static double rho=0.01;


/* pre_function*/


int nrerror(error_text)
char error_text[];
{
    //void exit();
    
    Rprintf("Numerical Recipes run-time error...\n");
    Rprintf("%s\n",error_text);
    Rprintf("...now exiting to system...\n");
    return 0;
}



float *vector(nl,nh)
int nl,nh;
{
    float *v;
    
    v=(float *)malloc((unsigned) (nh-nl+1)*sizeof(float));
    if (!v) nrerror("allocation failure in vector()");
        return v-nl;
}

int *ivector(nl,nh)
int nl,nh;
{
    int *v;
    
    v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
    if (!v) nrerror("allocation failure in ivector()");
        return v-nl;
}

double *dvector(nl,nh)
int nl,nh;
{
    double *v;
    
    v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
    if (!v) nrerror("allocation failure in dvector()");
        return v-nl;
}



float **matrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
    int i;
    float **m;
    
    m=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*));
    if (!m) nrerror("allocation failure 1 in matrix()");
        m -= nrl;
        
        for(i=nrl;i<=nrh;i++) {
            m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float));
            if (!m[i]) nrerror("allocation failure 2 in matrix()");
            m[i] -= ncl;
        }
    return m;
}

double **dmatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
    int i;
    double **m;
    
    m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
    if (!m) nrerror("allocation failure 1 in dmatrix()");
        m -= nrl;
        
        for(i=nrl;i<=nrh;i++) {
            m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
            if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
            m[i] -= ncl;
        }
    return m;
}

int **imatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
    int i,**m;
    
    m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
    if (!m) nrerror("allocation failure 1 in imatrix()");
        m -= nrl;
        
        for(i=nrl;i<=nrh;i++) {
            m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
            if (!m[i]) nrerror("allocation failure 2 in imatrix()");
            m[i] -= ncl;
        }
    return m;
}



float **submatrix(a,oldrl,oldrh,oldcl,oldch,newrl,newcl)
float **a;
int oldrl,oldrh,oldcl,oldch,newrl,newcl;
{
    int i,j;
    float **m;
    
    m=(float **) malloc((unsigned) (oldrh-oldrl+1)*sizeof(float*));
    if (!m) nrerror("allocation failure in submatrix()");
        m -= newrl;
        
        for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+oldcl-newcl;
            
            return m;
}



void free_vector(v,nl,nh)
float *v;
int nl,nh;
{
    free((char*) (v+nl));
}

void free_ivector(v,nl,nh)
int *v,nl,nh;
{
    free((char*) (v+nl));
}

void free_dvector(v,nl,nh)
double *v;
int nl,nh;
{
    free((char*) (v+nl));
}



void free_matrix(m,nrl,nrh,ncl,nch)
float **m;
int nrl,nrh,ncl,nch;
{
    int i;
    
    for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
        free((char*) (m+nrl));
        }

void free_dmatrix(m,nrl,nrh,ncl,nch)
double **m;
int nrl,nrh,ncl,nch;
{
    int i;
    
    for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
        free((char*) (m+nrl));
        }

void free_imatrix(m,nrl,nrh,ncl,nch)
int **m;
int nrl,nrh,ncl,nch;
{
    int i;
    
    for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
        free((char*) (m+nrl));
        }



void free_submatrix(b,nrl,nrh,ncl,nch)
float **b;
int nrl,nrh,ncl,nch;
{
    free((char*) (b+nrl));
}



float **convert_matrix(a,nrl,nrh,ncl,nch)
float *a;
int nrl,nrh,ncl,nch;
{
    int i,j,nrow,ncol;
    float **m;
    
    nrow=nrh-nrl+1;
    ncol=nch-ncl+1;
    m = (float **) malloc((unsigned) (nrow)*sizeof(float*));
    if (!m) nrerror("allocation failure in convert_matrix()");
        m -= nrl;
        for(i=0,j=nrl;i<=nrow-1;i++,j++) m[j]=a+ncol*i-ncl;
            return m;
}



void free_convert_matrix(b,nrl,nrh,ncl,nch)
float **b;
int nrl,nrh,ncl,nch;
{
    free((char*) (b+nrl));
}


#define TINY 1.0e-20;
#define pi 3.14159265

void ludcmp(a,n,indx,d)
int n,*indx;
double **a,*d;
{
    int i,imax,j,k;
    double big,dum,sum,temp;
    double *vv;
    
    vv=dvector(1,n);
    *d=1.0;
    for (i=1;i<=n;i++) {
        big=0.0;
        for (j=1;j<=n;j++)
            if ((temp=fabs(a[i][j])) > big) big=temp;
        if (big <= 0.0) nrerror("Singular matrix in routine LUDCMP");
        vv[i]=1.0/big;
    }
    for (j=1;j<=n;j++) {
        for (i=1;i<j;i++) {
            sum=a[i][j];
            for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
            a[i][j]=sum;
        }
        big=0.0;
        for (i=j;i<=n;i++) {
            sum=a[i][j];
            for (k=1;k<j;k++)
                sum -= a[i][k]*a[k][j];
            a[i][j]=sum;
            if ( (dum=vv[i]*fabs(sum)) >= big) {
                big=dum;
                imax=i;
            }
        }
        if (j != imax) {
            for (k=1;k<=n;k++) {
                dum=a[imax][k];
                a[imax][k]=a[j][k];
                a[j][k]=dum;
            }
            *d = -(*d);
            vv[imax]=vv[j];
        }
        indx[j]=imax;
        if (a[j][j] == 0.0) a[j][j]=TINY;
        if (j != n) {
            dum=1.0/(a[j][j]);
            for (i=j+1;i<=n;i++) a[i][j] *= dum;
        }
    }
    free_dvector(vv,1,n);
}

#undef TINY


void lubksb(a,n,indx,b)
double **a,*b;
int n,*indx;
{
    int i,ii=0,ip,j;
    double sum;
    
    for (i=1;i<=n;i++) {
        ip=indx[i];
        sum=b[ip];
        b[ip]=b[i];
        if (ii)
            for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
        else if (sum) ii=i;
        b[i]=sum;
    }
    for (i=n;i>=1;i--) {
        sum=b[i];
        for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
        b[i]=sum/a[i][i];
    }
}



double matrix_logdet(X, n)
double **X;
int n;
{
    int j, *indx;
    double d, logdet;
    
    indx=ivector(1,n);
    ludcmp(X, n, indx, &d);
    for(logdet=0.0,j=1; j<=n; j++) logdet+=log(fabs(X[j][j]));
        
        free_ivector(indx,1,n);
        return logdet;
}


/* Y=inv(X), return d=log(det(X)) */
double matrix_inverse(X, Y, n)
double **X, **Y;
int n;
{
    double d, *col;
    int i, j, *indx;
    double logdet;
    
    col=dvector(1,n);
    indx=ivector(1,n);
    
    ludcmp(X, n, indx, &d);
    
    for(logdet=0.0,j=1; j<=n; j++) logdet+=log(fabs(X[j][j]));
    /*
     if(d==0.0){ Rprintf("Singular matrix\n");  return; }
     */
        
        for(j=1; j<=n; j++){
            for(i=1; i<=n; i++) col[i]=0.0;
            col[j]=1.0;
            lubksb(X, n, indx, col);
            for(i=1; i<=n; i++) Y[i][j]=col[i];
        }
    
    for(i=1; i<=n; i++)
        for(j=1; j<=n; j++) { Y[i][j]=(Y[i][j]+Y[j][i])*0.5; Y[j][i]=Y[i][j]; }
    
    free_dvector(col,1,n); free_ivector(indx,1,n);
    
    return logdet;
}


/* Y=inv(X), return d=log(det(X)) */
int matrix_inverse_diag(X, Y, diag, n)
double **X, **Y, *diag;
int n;
{
    double d, *col;
    int i, j, *indx;
    
    col=dvector(1,n);
    indx=ivector(1,n);
    
    ludcmp(X, n, indx, &d);
    for(j=1; j<=n; j++) diag[j]=X[j][j];
        
        for(j=1; j<=n; j++){
            for(i=1; i<=n; i++) col[i]=0.0;
            col[j]=1.0;
            lubksb(X, n, indx, col);
            for(i=1; i<=n; i++) Y[i][j]=col[i];
        }
    
    for(i=1; i<=n; i++)
        for(j=1; j<=n; j++) { Y[i][j]=(Y[i][j]+Y[j][i])*0.5; Y[j][i]=Y[i][j]; }
    
    free_dvector(col,1,n); free_ivector(indx,1,n);
    
    return 0;
}



double matrix_trace(A,p)
double **A;
int p;
{
    int i;
    double sum;
    
    for(sum=0.0, i=1; i<=p; i++) sum+=A[i][i];
        
        return sum;
}


int matrix_sum(A,B,C,n,p)
double **A, **B, **C;
int n, p;
{
    int i, j;
    for(i=1; i<=n; i++)
        for(j=1; j<=p; j++) C[i][j]=A[i][j]+B[i][j];
            return 0;
}


/* Matrix: A: n by p; B: p by m;  C: n by m */
int matrix_multiply(A,B,C,n,p,m)
double **A, **B, **C;
int n, p, m;
{
    int i, j, k;
    for(i=1; i<=n; i++)
        for(j=1; j<=m; j++){
            C[i][j]=.0;
            for(k=1; k<=p; k++) C[i][j]+=A[i][k]*B[k][j];
        }
    return 0;
}


int matrix_vector_prod(A,b,d,n,p)
double **A, *b, *d;
int n,p;
{
    int i,j;
    for(i=1; i<=n; i++){
        d[i]=0.0;
        for(j=1; j<=p; j++) d[i]+=A[i][j]*b[j];
    }
    return 0;
}

double vector_matrix_vector(a,X,b,m,n)
double *a,**X,*b;
int m, n;
{
    double sum;
    int i, j;
    
    for(sum=0.0, i=1; i<=m; i++)
        for(j=1; j<=n; j++) sum+=a[i]*X[i][j]*b[j];
            
            return sum;
}


void copy_vector(a,b,p)
double *a, *b;
int p;
{
    int i;
    for(i=1; i<=p; i++) b[i]=a[i];
        }

void copy_matrix(a,b,n,p)
double **a,**b;
int n, p;
{
    int i, j;
    for(i=1; i<=n; i++)
        for(j=1; j<=p; j++) b[i][j]=a[i][j];
            }


int choldc(a, n, D)
double **a;
int n;
double **D;
{
    int i, j, k;
    double sum, *p;
    
    p=dvector(1,n);
    for (i=1; i<=n; i++) {
        for (j=i; j<=n; j++) {
            for (sum=a[i][j], k=i-1; k>=1; k--) sum -= a[i][k]*a[j][k];
            if (i==j) {
                if (sum <=0.0){ Rprintf("choldc failed"); return 1; }
                p[i]=sqrt(sum);
            } else a[j][i]=sum/p[i];
        }
    }
    /* Transfer the lower triangular part of A to the lower triangular matrix D */
    for(i=1; i<=n; i++){
        D[i][i]=p[i];
        for(j=1; j<i; j++){ D[i][j]=a[i][j]; D[j][i]=0.0; }
    }
    free_dvector(p,1,n);
    return 0;
}


/* calculate log(Gamma(s))  */
double loggamma(xx)
double xx;
{
    double x,tmp,ser;
    static double cof[6]={76.18009173,-86.50532033,24.01409822,
        -1.231739516,0.120858003e-2,-0.536382e-5};
    int j;
    
    x=xx-1.0;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.0;
    for (j=0;j<=5;j++) {
        x += 1.0;
        ser += cof[j]/x;
    }
    return -tmp+log(2.50662827465*ser);
}


/* calculate log(k!) */
double logpum(k)
int k;
{
    double value;
    int i;
    
    for(value=0.0, i=1; i<=k; i++) value+=log(1.0*i);
    
    return value;
}


/* generate the random variable form Gamma(a,b) */
double Rgamma(a,b)
double a, b;
{
    int ok;
    double d, q, un, u1, y, z;
    
    if(a<=0.0 || b<=0.0) { Rprintf("Gamma parameter error (<0.0)\n"); return 0; }
    
    if(a<1.0){  /* Ahrens, P.213 */
        ok=0;
        while(ok==0){
            un=0.0;
            GetRNGstate();
            while(un<=0.0 || un>=1.0) un=unif_rand();PutRNGstate();
            d=(2.718282+a)/2.718282;
            q=d*un;
            
            if(q<=1.0){
                z=exp(1.0/a*log(q));
                GetRNGstate();
                u1=unif_rand();PutRNGstate();
                if(u1<exp(-z)) ok=1;
            }
            else{
                z=-1.0*log((d-q)/a);
                GetRNGstate();
                u1=unif_rand();PutRNGstate();
                if(u1<exp((a-1)*log(z))) ok=1;
            }
        } /* end ok */
    }
    else {  /* a>=1.0 Fishman, P.214 */
        ok=0;
        while(ok==0){
            un=0.0;
            GetRNGstate();
            while(un<=0.0 || un>=1.0) un=unif_rand();PutRNGstate();
            y=-1.0*log(un);
            GetRNGstate();
            u1=unif_rand();PutRNGstate();
            if(u1<exp((a-1)*(log(y)-(y-1)))) { z=a*y; ok=1; }
        }
    }
    
    return z/b;
}


/* Generate a random variable from Beta(1,k), where
 the first parameter is 1, the second parameter is b */
double Rbeta(b)
double b;
{
    double un;
    un=0.0;
    GetRNGstate();
    
    while(un<=0.0 || un>=1.0) un=unif_rand();PutRNGstate();
    return 1.0-exp(1.0/b*log(un));
}


/* Generate deviates from Dirichlet(a1,a2,\ldots,a_k) */
int RDirichlet(w,a,k)
double *w,*a;
int k;
{
    double sum;
    int i;
    for(sum=0.0,i=1; i<=k; i++){
        w[i]=Rgamma(a[i],1.0);
        sum+=w[i];
    }
    for(i=1; i<=k; i++) w[i]/=sum;
        return 0;
}


double gasdev()
{
    static int iset=0;
    static double gset;
    double fac,r,v1,v2;
    
    if  (iset == 0) {
        do {    GetRNGstate();
            v1=unif_rand()*2.0-1.0;
            PutRNGstate();
            GetRNGstate();
            v2=unif_rand()*2.0-1.0;
            PutRNGstate();
            r=v1*v1+v2*v2;
        } while (r >= 1.0);
        fac=sqrt(-2.0*log(r)/r);
        gset=v1*fac;
        iset=1;
        return v2*fac;
    } else {
        iset=0;
        return gset;
    }
}


double Rgasdev(mean,variance)
double mean,variance;
{
    static int iset=0;
    static double gset;
    double fac,r,v1,v2;
    
    if  (iset == 0) {
        do {    GetRNGstate();
            v1=unif_rand()*2.0-1.0;
            PutRNGstate();
            GetRNGstate();
            v2=unif_rand()*2.0-1.0;
            PutRNGstate();
            r=v1*v1+v2*v2;
        } while (r >= 1.0);
        fac=sqrt(-2.0*log(r)/r);
        gset=v1*fac;
        iset=1;
        return v2*fac*sqrt(variance)+mean;
    } else {
        iset=0;
        return gset*sqrt(variance)+mean;
    }
}


int RNORM(x,mu,Sigma,p)
double *x,*mu,**Sigma;
int p;
{
    int i, j;
    double **D, *z;
    
    D=dmatrix(1,p,1,p);
    z=dvector(1,p);
    
    choldc(Sigma,p,D);
    
    for(i=1; i<=p; i++) z[i]=gasdev();
        for(i=1; i<=p; i++){
            x[i]=mu[i];
            for(j=1; j<=i; j++) x[i]+=D[i][j]*z[j];
        }
    
    free_dmatrix(D,1,p,1,p);
    free_dvector(z,1,p);
    
    return 0;
}


int Rwishart(B,df,Sigma,p)
double **B,df,**Sigma;
int p;
{
    double **Z, **A, *Y;
    int i, j, k;
    
    Z=dmatrix(1,p,1,p);
    A=dmatrix(1,p,1,p);
    Y=dvector(1,p);
    
    for(i=1; i<=p; i++)
        for(j=1; j<=p; j++) Z[i][j]=A[i][j]=B[i][j]=0.0;
            
            choldc(Sigma, p, A);
            
            for(j=1; j<=p; j++)
                for(i=1; i<=p; i++) Z[i][j]=gasdev();
                    for(i=1; i<=p; i++) Y[i]=Rgamma(0.5*df,0.5);
                        
                        B[1][1]=Y[1];
                        for(j=2; j<=p; j++){
                            B[j][j]=Y[j];
                            for(i=1; i<j; i++) B[j][j]+=Z[i][j]*Z[i][j];
                            B[1][j]=Z[1][j]*sqrt(Y[1]);
                        }
    for(j=2; j<=p; j++)
        for(i=2; i<j; i++){
            B[i][j]=Z[i][j]*sqrt(Y[i]);
            for(k=1; k<=i-1; k++) B[i][j]+=Z[k][i]*Z[k][j];
        }
    for(i=1; i<=p; i++)
        for(j=1; j<i; j++) B[i][j]=B[j][i];
            
            matrix_multiply(A,B,Z,p,p,p);
            
            for(i=1; i<=p; i++)
                for(j=1; j<i; j++){ A[j][i]=A[i][j]; A[i][j]=0.0; }
    
    matrix_multiply(Z,A,B,p,p,p);
    
    free_dmatrix(Z,1,p,1,p);
    free_dmatrix(A,1,p,1,p);
    free_dvector(Y,1,p);
    
    return 0;
}


/* calculated the log-density of  z~gamma(a,b) */
double dloggamma(x,a,b)
double x,a,b;
{
    double logcon, den;
    logcon=loggamma(a);
    den=log(b)-b*x+(a-1)*log(b*x)-logcon;
    return den;
}


double dloggauss(z,mean,variance)
double z,mean,variance;
{
    double sum;
    sum=-0.5*log(2.0*pi*variance);
    sum+=-0.5*(z-mean)*(z-mean)/variance;
    return sum;
}

double dlogstudent(z,k)
double z;
int k;  /* the degree of freedom */
{
    double logprob;
    logprob=-0.5*(k+1)*log(1.0+z*z/k);
    return logprob;
}



double DLOGGAUSS(z,mean,variance,p)
double *z, *mean, **variance;
int p;
{
    int i;
    double logdet, sum, **mat,*vect, *mu;
    
    mat=dmatrix(1,p,1,p);
    vect=dvector(1,p);
    mu=dvector(1,p);
    
    for(i=1; i<=p; i++) mu[i]=z[i]-mean[i];
        logdet=matrix_inverse(variance,mat,p);
        matrix_vector_prod(mat,mu,vect,p,p);
        
        for(sum=0.0,i=1; i<=p; i++) sum+=mu[i]*vect[i];
            sum*=-0.5;
            sum+=-0.5*logdet-0.5*p*log(2.0*pi);
            
            free_dmatrix(mat,1,p,1,p);
            free_dvector(vect,1,p);
            free_dvector(mu,1,p);
            
            return sum;
}


double Dlogwishart(D,df,Sigma,p)
double **D,df,**Sigma;
int p;
{
    int i, j;
    double a, sum, logdet1, logdet2, **mt1, **mt2;
    
    mt1=dmatrix(1,p,1,p);
    mt2=dmatrix(1,p,1,p);
    
    for(i=1; i<=p; i++)
        for(j=1; j<=p; j++) mt1[i][j]=D[i][j];
            logdet1=matrix_inverse(mt1,mt2,p);
            
            for(i=1; i<=p; i++)
                for(j=1; j<=p; j++) mt1[i][j]=Sigma[i][j];
                    logdet2=matrix_inverse(mt1,mt2,p);
                    
                    matrix_multiply(mt2,D,mt1,p,p,p);
                    
                    for(sum=0.0,i=1; i<=p; i++) sum+=mt1[i][i];
                        sum*=-0.5;
                        sum+=0.5*(df-p-1)*logdet1;
                        sum+=-0.5*df*logdet2;
                        sum+=-0.5*df*p*log(2.0);
                        sum+=-0.25*p*(p-1)*log(pi);
                        for(i=1; i<=p; i++){
                            a=df-i+1;
                            if(a<0.0001) a=0.0001;
                            sum+=-loggamma(0.5*a);
                        }
    
    free_dmatrix(mt1,1,p,1,p);
    free_dmatrix(mt2,1,p,1,p);
    
    return sum;
}


int uniform_direction(d, n)
double *d;
int n;
{
    double sum;
    int k;
    
    for(sum=0, k=1; k<=n; k++){
        d[k]=gasdev();
        sum+=d[k]*d[k];
    }
    for(k=1; k<=n; k++)
        d[k]=d[k]/sqrt(sum);
        
        return 0;
}


int dmaxclass(z,n)
double *z;
int n;
{
    int i, maxi;
    double maxd;
    
    maxd=z[1]; maxi=1;
    for(i=2; i<=n; i++)
        if(z[i]>maxd){ maxd=z[i]; maxi=i; }
    
    return maxi;
}

double dmaximum(z,n)
double *z;
int n;
{
    int i;
    double maxd;
    
    maxd=z[1];
    for(i=2; i<=n; i++)
        if(z[i]>maxd){ maxd=z[i];}
    
    return maxd;
}


int imaxclass(z,n)
int *z, n;
{
    int i, maxi;
    int maxd;
    
    maxd=z[1]; maxi=1;
    for(i=2; i<=n; i++)
        if(z[i]>maxd){ maxd=z[i]; maxi=i; }
    
    return maxi;
}


int binary_trans(k,l,d)
int k, l,*d;
{
    int i, j;
    for(i=1; i<=l; i++) d[i]=0;
        j=l;
        while(k>0){
            d[j]=k%2;
            k=k/2;
            j--;
        }
    return 0;
}


double logsum(a,b)
double a, b;
{
    double sum;
    if(a>b) sum=a+log(1.0+exp(b-a));
        else sum=b+log(1.0+exp(a-b));
            return sum;
}


double maxvector(x,n)
double *x;
int n;
{
    double max;
    int i;
    max=x[1];
    for(i=2; i<=n; i++)
        if(x[i]>max) max=x[i];
            return max;
}


double minvector(x,n)
double *x;
int n;
{
    double min;
    int i;
    min=x[1];
    for(i=2; i<=n; i++)
        if(x[i]<min) min=x[i];
            return min;
}



double sample_variance(x,n)
double *x;
int n;
{
    int i;
    double sum1, sum2, mean, var;
    
    sum1=sum2=0.0;
    for(i=1; i<=n; i++){
        sum1+=x[i];
        sum2+=x[i]*x[i];
    }
    mean=sum1/n;
    var=(sum2-n*mean*mean)/(n-1);
    
    return var;
}



/* Return the value ln[Gamma(xx)] for xx>0 */
double gammln(xx)
double xx;
{
    double x,tmp,ser;
    static double cof[6]={76.18009173,-86.50532033,24.01409822,
        -1.231739516,0.120858003e-2,-0.536382e-5};
    int j;
    
    x=xx-1.0;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.0;
    for (j=0;j<=5;j++) {
        x += 1.0;
        ser += cof[j]/x;
    }
    return -tmp+log(2.50662827465*ser);
}


#define ITMAX 100
#define EPS 3.0e-7
void gser(gamser,a,x,gln)
double a,x,*gamser,*gln;
{
    int n;
    double sum,del,ap;
    
    *gln=gammln(a);
    if (x <= 0.0) {
        if (x < 0.0) nrerror("x less than 0 in routine GSER");
        *gamser=0.0;
        return;
    } else {
        ap=a;
        del=sum=1.0/a;
        for (n=1;n<=ITMAX;n++) {
            ap += 1.0;
            del *= x/ap;
            sum += del;
            if (fabs(del) < fabs(sum)*EPS) {
                *gamser=sum*exp(-x+a*log(x)-(*gln));
                return;
            }
        }
        nrerror("a too large, ITMAX too small in routine GSER");
        return;
    }
}

#undef ITMAX
#undef EPS

#define ITMAX 500
#define EPS 1.0e-9
void gserln(gamser,a,x,gln)
double a,x,*gamser,*gln;
{
    int n;
    double sum,del,ap;
    
    *gln=gammln(a);
    if (x <= 0.0) {
        if (x < 0.0) nrerror("x less than 0 in routine GSERLN");
        *gamser=0.0;
        return;
    } else {
        ap=a;
        del=sum=1.0/a;
        for (n=1;n<=ITMAX;n++) {
            ap += 1.0;
            del *= x/ap;
            sum += del;
            if (fabs(del) < fabs(sum)*EPS) {
                *gamser=log(sum)-x+a*log(x)-(*gln);
                return;
            }
        }
        nrerror("a too large, ITMAX too small in routine GSERLN");
        return;
    }
}
#undef ITMAX
#undef EPS



#define ITMAX 100
#define EPS 3.0e-7
void gcf(gammcf,a,x,gln)
double a,x,*gammcf,*gln;
{
    int n;
    double gold=0.0,g,fac=1.0,b1=1.0;
    double b0=0.0,anf,ana,an,a1,a0=1.0;
    /* float gammln();
     void nrerror(); */
    
    *gln=gammln(a);
    a1=x;
    for (n=1;n<=ITMAX;n++) {
        an=(double) n;
        ana=an-a;
        a0=(a1+a0*ana)*fac;
        b0=(b1+b0*ana)*fac;
        anf=an*fac;
        a1=x*a0+anf*a1;
        b1=x*b0+anf*b1;
        if (a1) {
            fac=1.0/a1;
            g=b1*fac;
            if (fabs((g-gold)/g) < EPS) {
                *gammcf=exp(-x+a*log(x)-(*gln))*g;
                return;
            }
            gold=g;
        }
    }
    nrerror("a too large, ITMAX too small in routine GCF");
}
#undef ITMAX
#undef EPS

#define ITMAX 500
#define EPS 1.0e-9
void gcfln(gammcf,a,x,gln)
double a,x,*gammcf,*gln;
{
    int n;
    double gold=0.0,g,fac=1.0,b1=1.0;
    double b0=0.0,anf,ana,an,a1,a0=1.0;
    /* float gammln();
     void nrerror(); */
    
    *gln=gammln(a);
    a1=x;
    for (n=1;n<=ITMAX;n++) {
        an=(double) n;
        ana=an-a;
        a0=(a1+a0*ana)*fac;
        b0=(b1+b0*ana)*fac;
        anf=an*fac;
        a1=x*a0+anf*a1;
        b1=x*b0+anf*b1;
        if (a1) {
            fac=1.0/a1;
            g=b1*fac;
            if (fabs((g-gold)/g) < EPS) {
                *gammcf=-x+a*log(x)-(*gln)+log(g);
                return;
            }
            gold=g;
        }
    }
    nrerror("a too large, ITMAX too small in routine GCFLN");
}
#undef ITMAX
#undef EPS



double gammp(a,x)
double a,x;
{
    double gamser,gammcf,gln;
    /* void gser(),gcf(),nrerror(); */
    
    if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine GAMMP");
        if (x < (a+1.0)) {
            gser(&gamser,a,x,&gln);
            return gamser;
        } else {
            gcf(&gammcf,a,x,&gln);
            return 1.0-gammcf;
        }
}


double gammq(a,x)
double a,x;
{
    double gamser,gammcf,gln;
    // void gcf(),gser(),nrerror();
    
    if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine GAMMQ");
        if (x < (a+1.0)) {
            gser(&gamser,a,x,&gln);
            return 1.0-gamser;
        } else {
            gcf(&gammcf,a,x,&gln);
            return gammcf;
        }
}



double gammpln(a,x,tail)
double a,x;
int *tail;
{
    double gamserln,gammcfln,gln;
    /* void gser(),gcf(),nrerror(); */
    
    
    if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine GAMMP");
        if (x < (a+1.0)) {
            *tail=1;
            gserln(&gamserln,a,x,&gln);
            return gamserln;
        } else {
            gcfln(&gammcfln,a,x,&gln);
            *tail=-1;
            return gammcfln;
        }
}


/* Return the CDF of the standard normal distribution */
double gauss_cdf(x)
double x;
{
    double s, prob;
    
    s=0.5*x*x;
    if(x>0) prob=0.5+0.5*gammp(0.5,s);
    else prob=0.5-0.5*gammp(0.5,s);
    
    return prob;
}


/* Return the log(CDF) of the standard normal distribution */
double gauss_cdf_ln(x)
double x;
{
    double s, a, b, logprob;
    int tail;
    
    s=0.5*x*x;
    a=log(0.5);
    b=gammpln(0.5,s,&tail)+a;
    
    if(tail==1){
        if(x>0.0) logprob=a+log(1+exp(b-a));
        else logprob=a+log(1-exp(b-a));
    }
    else{
        if(x<0.0) logprob=b;
        else logprob=b+log(exp(-b)-1.0);
    }
    
    return logprob;
}



/* return Gamma'(z)/Gamma(z)   */
/* Refer to "mathematics handbook pp.287" */
double diGamma(z)
double z;
{
    int i;
    double sum, delta, epsilon=3.0e-7;
    
    sum=-1.0/z-0.5772156649;
    
    delta=1.0-1.0/(1+z);
    sum+=delta;
    i=1;
    while(delta>epsilon){
        i++;
        delta=1.0/i-1.0/(i+z);
        sum+=delta;
    }
    
    return sum;
}


/* return Gamma'(z)/Gamma(z)   */
/* Refer to "mathematics handbook pp.287" */
double derivative_gamma(z)
double z;
{
    int i;
    double sum, delta, epsilon=3.0e-7;
    
    sum=-1.0/z-0.5772156649;
    
    delta=1.0-1.0/(1+z);
    sum+=delta;
    i=1;
    while(delta>epsilon){
        i++;
        delta=1.0/i-1.0/(i+z);
        sum+=delta;
    }
    
    return sum;
}

double dlogGnormal(z,mu,alpha,beta)
double z, mu, alpha,beta;
{
    double sum;
    
    if(fabs(z-mu)<1.0e-10) sum=log(beta)-log(2.0*alpha)-gammln(1.0/beta);
        else sum=log(beta)-log(2.0*alpha)-gammln(1.0/beta)-exp(beta*log(fabs(z-mu)/alpha));
            
            return sum;
}

double correlation(z1,z2,p)
double *z1, *z2;
int p;
{
    double ave1, ave2, sq1, sq2, sum;
    int i;
    
    ave1=ave2=0.0; sq1=sq2=0.0;
    for(i=1; i<=p; i++){
        ave1+=z1[i]; ave2+=z2[i];
        sq1+=z1[i]*z1[i]; sq2+=z2[i]*z2[i];
    }
    ave1/=p; ave2/=p;
    sq1=(sq1-p*ave1*ave1)/(p-1);
    sq2=(sq2-p*ave2*ave2)/(p-1);
    
    if(sq1<=0.0 || sq2<=0.0) return 0.0;
    else{
        for(sum=0.0,i=1; i<=p; i++) sum+=(z1[i]-ave1)*(z2[i]-ave2);
        sum/=p;
        sum/=sqrt(sq1*sq2);
        return sum;
    }
}


int permut_sample(sam,n)
int *sam, n;
{
    int j,k,u,v,*b;
    
    b=ivector(1,n);
    
    for(j=1; j<=n; j++) b[j]=j;
        k=0;
        while(k<n){
            u=0;
            GetRNGstate();
            while(u<=0 || u>n-k) u=floor(unif_rand()*(n-k))+1;PutRNGstate();
            sam[k+1]=b[u];
            for(v=u; v<n-k; v++) b[v]=b[v+1];
            k++;
        }
    
    return 0;
}



int random_order(x,n)
int *x, n;
{
    int i, j, k, m, *y;
    
    y=ivector(1,n);
    
    m=n;
    for(i=1; i<=m; i++) y[i]=i;
        for(i=1; i<=n; i++){
            j=0;
            GetRNGstate();
            while(j<1 || j>m) j=(int)(unif_rand()*m)+1;PutRNGstate();
            x[i]=y[j];
            for(k=j+1; k<=m; k++) y[k-1]=y[k];
            m--;
        }
    
    free_ivector(y,1,n);
    
    return 0;
}


/* Generate a subset sample of size M from the set 1:N */
int subset_sample(x,z,M,N)
int *x, *z,M, N;
{
    int i, j, k, m, *y;
    
    y=ivector(1,N);
    
    m=N;
    for(i=1; i<=N; i++) y[i]=i;
        for(i=1; i<=M; i++){
            j=0;
            GetRNGstate();
            while(j<1 || j>m) j=floor(unif_rand()*m)+1;PutRNGstate();
            x[i]=y[j];
            for(k=j+1; k<=m; k++) y[k-1]=y[k];
            m--;
        }
    
    for(i=1; i<=N-M; i++) z[i]=y[i];
        
        free_ivector(y,1,N);
        
        return 0;
}



void indexx(n,arrin,indx)
int n,*indx;
double *arrin;
{
    int l,j,ir,indxt,i;
    double q;
    
    for (j=1;j<=n;j++) indx[j]=j;
        l=(n >> 1) + 1;
        ir=n;
        for (;;) {
            if (l > 1)
                q=arrin[(indxt=indx[--l])];
            else {
                q=arrin[(indxt=indx[ir])];
                indx[ir]=indx[1];
                if (--ir == 1) {
                    indx[1]=indxt;
                    return;
                }
            }
            i=l;
            j=l << 1;
            while (j <= ir) {
                if (j < ir && arrin[indx[j]] < arrin[indx[j+1]]) j++;
                if (q < arrin[indx[j]]) {
                    indx[i]=indx[j];
                    j += (i=j);
                }
                else j=ir+1;
            }
            indx[i]=indxt;
        }
}


void indexx_integer(n,arrin,indx)
int n,*indx;
int *arrin;
{
    int l,j,ir,indxt,i;
    double q;
    
    for (j=1;j<=n;j++) indx[j]=j;
        l=(n >> 1) + 1;
        ir=n;
        for (;;) {
            if (l > 1)
                q=arrin[(indxt=indx[--l])];
            else {
                q=arrin[(indxt=indx[ir])];
                indx[ir]=indx[1];
                if (--ir == 1) {
                    indx[1]=indxt;
                    return;
                }
            }
            i=l;
            j=l << 1;
            while (j <= ir) {
                if (j < ir && arrin[indx[j]] < arrin[indx[j+1]]) j++;
                if (q < arrin[indx[j]]) {
                    indx[i]=indx[j];
                    j += (i=j);
                }
                else j=ir+1;
            }
            indx[i]=indxt;
        }
}


void indexx_convert_double(n,indx,x,y)
int n, *indx;
double *x, *y;
{
    int i;
    for(i=1; i<=n; i++) y[indx[i]]=x[i];
        }


void indexx_convert_integer(n,indx,x,y)
int n, *indx;
int *x, *y;
{
    int i;
    for(i=1; i<=n; i++) y[indx[i]]=x[i];
        }

/* coefficients for the rational approximants for the normal probit: */
#define a1	(-3.969683028665376e+01)
#define a2	( 2.209460984245205e+02)
#define a3	(-2.759285104469687e+02)
#define a4	( 1.383577518672690e+02)
#define a5	(-3.066479806614716e+01)
#define a6	( 2.506628277459239e+00)
#define b1	(-5.447609879822406e+01)
#define b2	( 1.615858368580409e+02)
#define b3	(-1.556989798598866e+02)
#define b4	( 6.680131188771972e+01)
#define b5	(-1.328068155288572e+01)
#define c1	(-7.784894002430293e-03)
#define c2	(-3.223964580411365e-01)
#define c3	(-2.400758277161838e+00)
#define c4	(-2.549732539343734e+00)
#define c5	( 4.374664141464968e+00)
#define c6	( 2.938163982698783e+00)
#define d1	( 7.784695709041462e-03)
#define d2	( 3.224671290700398e-01)
#define d3	( 2.445134137142996e+00)
#define d4	( 3.754408661907416e+00)
#define p_low	0.02425
#define logp_low -3.719339
#define p_high	(1.0 - p_low)
#define logp_high -0.02454887

/**
 * Returns the probit value of the normal distribution CDF.  This is
 * an implementation of the algorithm published at
 * http://home.online.no/~pjacklam/notes/invnorm/
 */
double inverse_normal_cdf(double p) {
    double q, x;
    
    if(0.0 < p && p < p_low) {
        /* rational approximation for the lower region */
        q = sqrt(-2.0*log(p));
        x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1);
    } else if(p_low <= p && p <= p_high) {
        double r;
        /* rational approximation for the central region */
        q = p - 0.5;
        r = q*q;
        x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q / (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0);
    } else /* if(p_high < p && p < 1.0) */ {
        /* rational approximation for the upper region */
        q = sqrt(-2.0*log(1.0-p));
        x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1.0);
    }
    
    if(0.0 < p && p < 1.0) {
        double u, e;
        e = 0.5 * erfc(-x/sqrt(2.0)) - p;
        u = e * sqrt(2.0*M_PI) * exp(x*x/2.0);
        x = x - u/(1.0 + x*u/2.0);
    }
    
    return x;
}

/* the input p is given in logairithm */
double inverse_normal_cdf_log(double logp) {
    double q, x;
    
    if(logp < logp_low) {
        /* rational approximation for the lower region */
        q = sqrt(-2.0*logp);
        x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1);
    } else if(logp_low <= logp && logp <= logp_high) {
        double r;
        /* rational approximation for the central region */
        q = exp(logp) - 0.5;
        r = q*q;
        x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q / (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0);
    } else /* if(p_high < p && p < 1.0) */ {
        /* rational approximation for the upper region */
        q = sqrt(-2.0*(logp+log(exp(-logp)-1)));
        x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1.0);
    }
    
    if(logp >= logp_low && logp < 0.0) {
        double u, e;
        e = 0.5 * erfc(-x/sqrt(2.0)) - exp(logp);
        u = e * sqrt(2.0*M_PI) * exp(x*x/2.0);
        x = x - u/(1.0 + x*u/2.0);
    }
    
    return x;
}

#undef a1
#undef a2
#undef a3
#undef a4
#undef a5
#undef a6
#undef b1
#undef b2
#undef b3
#undef b4
#undef b5
#undef c1
#undef c2
#undef c3
#undef c4
#undef c5
#undef c6
#undef d1
#undef d2
#undef d3
#undef d4
#undef p_low
#undef p_high
#undef logp_low
#undef logp_high

/**
 * Returns the quantile for the two-tailed student's t distribution.
 * This is an implementation of the algorithm in
 * G. W. Hill. "Algorithm 396: Student's t-Quantiles." Communications
 * of the ACM 13(10):619--620.  ACM Press, October, 1970.
 */
double inverse_t_cdf(double p, long n) {
    double a, b, c, d, x, y;
    
    if(n < 1) {
        /* you can put your own error handling here */
        Rprintf("tquantile(%f, %d): error: second argument must be >= 1 !", p, n);
        return 0.0;
    } else if(p > 1.0 || p <= 0.0) {
        /* you can put your own error handling here */
        Rprintf("tquantile(%f, %d): error: first argument must be in (0.0, 1.0] !", p, n);
        return 0.0;
    }
    
    if(n == 1) {
        /* special case */
        p *= M_PI_2;
        return cos(p) / sin(p);
    }
    
    a = 1.0 / (n-0.5);
    b = pow(48.0 / a, 2.0);
    c = ((20700.0 * a / b - 98.0) * a - 16.0) * a + 96.36;
    d = ((94.5 / (b + c) - 3.0) / b + 1.0) * sqrt(a * M_PI_2) * (double)n;
    x = d * p;
    y = pow(x, 2.0/(double)n);
    if(y > 0.05 + a) {
        /* asymptotic inverse expansion about the normal */
        x = inverse_normal_cdf(p * 0.5);
        y = x * x;
        if(n < 5) {
            c += 0.3 * ((double)n - 4.5) * (x + 0.6);
            c = (((0.5 * d * x - 0.5) * x - 7.0) * x - 2.0) * x + b + c;
            y = (((((0.4 * y + 6.3) * y + 36.0) * y + 94.5) / c - y - 3.0) / b + 1.0) * x;
            y *= a * y;
            if(y > 0.002)
                y = exp(y) - 1.0;
            else
                y += 0.5 * y * y;
        }
    } else
        y = ((1.0/((((double)n + 6.0)/((double)n * y) - 0.089 * d - 0.822) * ((double)n+2.0) * 3.0) + 0.5 / ((double)n+4.0))*y - 1.0) * ((double)n + 1.0) / ((double)n + 2.0) + 1.0 / y;
    
    return sqrt((double)n * y);
}

int iminimum(d,n)
int *d, n;
{
    int min, i;
    
    min=d[1];
    for(i=2; i<=n; i++){
        if(d[i]<min) min=d[i];
    }
    
    return min;
}

int imaximum(d,n)
int *d, n;
{
    int max, i;
    
    max=d[1];
    for(i=2; i<=n; i++){
        if(d[i]>max) max=d[i];
    }
    
    return max;
}

int iminmax(d,n,min,max)
int *d, n,*min,*max;
{
    int i;
    
    *min=d[1]; *max=d[1];
    for(i=2; i<=n; i++){
        if(d[i]<*min) *min=d[i];
        else if(d[i]>*max) *max=d[i];
    }
    
    return 0;
}

/* Generate n MCMC samples from generalized normal(mu,alpha,beta) and
 return the samples z  */
double MCMCGnormal(mu,alpha,beta,n,z)
double mu, alpha, beta, *z;
int n;
{
    double x,y, fx, fy, un, r;
    int i, accept, warm=50;
    
    x=mu+0.5*alpha;
    fx=dlogGnormal(x,mu,alpha,beta);
    
    for(i=1; i<=n+warm; i++){
        
        y=x+gasdev()*alpha;
        fy=dlogGnormal(y,mu,alpha,beta);
        
        r=fy-fx;
        if(r>0.0) accept=1;
        else{
            un=0.0;
            GetRNGstate();
            while(un<=0.0) un=unif_rand();PutRNGstate();
            if(un<exp(r)) accept=1;
            else accept=0;
        }
        
        if(accept==1){ x=y; fx=fy; }
        if(i>warm) z[i-warm]=x;
    }
    
    return 0;
}

// x is a vector of z-values, gamma is a parameter with a default value of 0.1,
// the output is the estimated mean mu and standard deviation sigma
//
int EstNull_fdr(x,n,gamma,mu,sigma)
double *x,gamma,*mu,*sigma;
int n;
{
    double *t, *phiplus, *phiminus, *dphiplus,*dphiminus, *phi, *dphi;
    double gan, shat, uhat;
    double s, tt, a, b, c, da, db;
    int i, j;
    
    t=dvector(1,1000);
    phiplus=dvector(1,1000);
    phiminus=dvector(1,1000);
    dphiplus=dvector(1,1000);
    dphiminus=dvector(1,1000);
    phi=dvector(1,1000);
    dphi=dvector(1,1000);
    
    for(i=1; i<=1000; i++){
        t[i]=1.0*i/200;
        phiplus[i]=1;
        phiminus[i]=1;
        dphiplus[i]=1;
        phi[i]=1;
        dphi[i]=1;
    }
    gan=exp(-gamma*log(1.0*n));
    
    for(i=1; i<=1000; i++){
        s=t[i];
        for(phiplus[i]=0, j=1; j<=n; j++) phiplus[i]+=cos(s*x[j]);
        phiplus[i]/=n;
        
        for(phiminus[i]=0, j=1; j<=n; j++) phiminus[i]+=sin(s*x[j]);
        phiminus[i]/=n;
        
        for(dphiplus[i]=0, j=1; j<=n; j++) dphiplus[i]+=-x[j]*sin(s*x[j]);
        dphiplus[i]/=n;
        
        for(dphiminus[i]=0, j=1; j<=n; j++) dphiminus[i]+=x[j]*cos(s*x[j]);
        dphiminus[i]/=n;
        
        phi[i]=sqrt(phiplus[i]*phiplus[i]+phiminus[i]*phiminus[i]);
    }
    
    
    i=1;
    while(phi[i]>gan && i<1000) i++;
    
    // Rprintf("i=%d gan=%g phi=%g\n",i, gan, phi[i]);
    
    tt=t[i];
    a=phiplus[i];
    b=phiminus[i];
    da=dphiplus[i];
    db=dphiminus[i];
    c=phi[i];
    
    shat=sqrt(-(a*da+b*db)/(tt*c*c));
    uhat=-(da*b-db*a)/(c*c);
   // epshat=1-c*exp((tt*shat)*(tt*shat)/2);
    
    *mu=uhat;
    *sigma=shat;
    
    
    free_dvector(t,1,1000);
    free_dvector(phiplus,1,1000);
    free_dvector(phiminus,1,1000);
    free_dvector(dphiplus,1,1000);
    free_dvector(dphiminus,1,1000);
    free_dvector(phi,1,1000);
    free_dvector(dphi,1,1000);
    
    return 0;
}


double EpsEst_fdr(x,n,mu,sigma)
double *x, mu, sigma;
int n;
{
    double *z, tmax, epshat, epsest, *xi,*f,*w,*co,*tt;
    double t, sum1, sum2;
    int i, j, KK, l; 
    
    z=dvector(1,n);
    xi=dvector(0,100);
    f=dvector(0,100);
    w=dvector(0,100);
    co=dvector(0,100);
    
    
    for(i=1; i<=n; i++) z[i]=(x[i]-mu)/sigma;
        for(i=0; i<=100; i++) xi[i]=1.0*i/100;
            
            tmax=sqrt(log(1.0*n));
            
            KK=floor(tmax/0.1); 
            tt=dvector(0,KK);
            for(i=0; i<=KK; i++) tt[i]=0.1*i;
                
                epsest=0.0;
                for(j=0; j<=KK; j++){
                    
                    t=tt[j];
                    for(i=0; i<=100; i++){
                        f[i]=exp((t*xi[i])*(t*xi[i])/2);
                        w[i]=1-fabs(xi[i]);
                    }
                    
                    for(i=0; i<=100; i++){
                        co[i]=0.0;
                        for(l=1; l<=n; l++) co[i]+=cos(t*xi[i]*z[l]);
                        co[i]/=n;
                    }
                    
                    for(sum1=sum2=0.0, i=0; i<=100; i++){ 
                        sum1+=w[i]*f[i]*co[i];
                        sum2+=w[i];
                    }
                    
                    epshat=1.0-sum1/sum2;
                    if(epshat>epsest) epsest=epshat;
                }
    
    free_dvector(z,1,n);
    free_dvector(xi,0,100);
    free_dvector(f,0,100);
    free_dvector(w,0,100);
    free_dvector(co,0,100);
    free_dvector(tt,0,KK);
    
    return epsest;
}


int corsel(datanum,row,col,datax, Qmat)
int datanum,*row, *col;
double *datax, **Qmat;
{
    int  grid, bestgrid, rep;
    int *indx, *indx2,  *sam, i, j, k, k1, k2, m, com, iter, ok;
    double **resam,*rez,*dis,*dif,*mu,*alpha,*beta,*prob,*logprob,*wei,*FDR,*Q, entropy;
    double delta,z,un,max,min,sum,tep,tep1,tep2,eta,dmu,dalpha,dbeta;
    double bestentropy, *bestprob, *bestmu, *bestalpha, *bestbeta;
    //double **DataX;
    
    
    /* initialize the random number generator */
    /*
     ltime=time(NULL);
     stime=(unsigned int)ltime/2;
     srand(stime);
     
     ins=fopen("ee.log","a");
     fprintf(ins, "random seed=%d\n", stime);
     fclose(ins);
     */
    
    mu=dvector(1,GRID);
    alpha=dvector(1,GRID);
    beta=dvector(1,GRID);
    prob=dvector(1,GRID);
    wei=dvector(1,GRID);
    sam=ivector(1,MAX_NUM);
    indx=ivector(1,MAX_NUM);
    indx2=ivector(1,GRID);
    FDR=dvector(1,MAX_NUM);
    Q=dvector(1,MAX_NUM);
    bestmu=dvector(1,GRID);
    bestalpha=dvector(1,GRID);
    bestbeta=dvector(1,GRID);
    bestprob=dvector(1,GRID);
    logprob=dvector(1,GRID);
    resam=dmatrix(1,GRID,1,samsize);
    rez=dvector(1,samsize);
    dis=dvector(1,GRID);
    dif=dvector(1,GRID);
    
    if(datanum>MAX_NUM) Rprintf("Data size is TOO LARGE in corsel.c\n");
        else Rprintf("working subdata size=%d\n", datanum);
            
            indexx(datanum,datax,indx);
            
            grid=GRID+1;  bestentropy=1.0e+100;
            while(grid>mingrid){
                
                grid--;
                
                /* Initialization */
                mu[1]=0.0; alpha[1]=0.0; k=0;
                for(i=1; i<=(int)(datanum*0.99); i++){
                    mu[1]+=datax[indx[i]];
                    alpha[1]+=datax[indx[i]]*datax[indx[i]];
                    k++;
                }
                mu[1]=mu[1]/k;
                alpha[1]=sqrt(alpha[1]/k-mu[1]*mu[1])*sqrt(2.0);
                beta[1]=2.0; prob[1]=0.99;
                
                for(j=2; j<=grid; j++){
                    tep=(1-prob[1])/(grid-1);
                    mu[j]=0.0; alpha[j]=0.0; k=0;
                    for(i=(int)(datanum*(prob[1]+(j-2)*tep)); i<=(int)(datanum*(prob[1]+(j-1)*tep)); i++){
                        mu[j]+=datax[indx[i]];
                        alpha[j]+=datax[indx[i]]*datax[indx[i]];
                        k++;
                    }
                    mu[j]=mu[j]/k;
                    alpha[j]=sqrt(alpha[j]/k-mu[j]*mu[j])*sqrt(2.0);
                    beta[j]=2.0; prob[j]=tep;
                }
                for(j=1; j<=grid; j++){
                    logprob[j]=log(prob[j])-log(prob[grid]);
                }
                
                
                for(iter=1; iter<=100; iter++){
                    
                    random_order(sam,datanum);
                    for(i=1; i<=datanum; i++){
                        
                        rep=datanum*(iter-1)+i;
                        if(rep<=WARM*stepscale) delta=rho;
                        else delta=rho*exp(-log(1.0*(rep-(WARM-1)*stepscale)/stepscale));
                        
                        for(sum=0.0,j=1; j<=grid; j++) sum+=exp(logprob[j]);
                        for(j=1; j<=grid; j++) prob[j]=exp(logprob[j])/sum;
                        
                        z=datax[sam[i]];
                        max=-1.0e+100;
                        for(j=1; j<=grid; j++){
                            wei[j]=dlogGnormal(z,mu[j],alpha[j],beta[j])+log(prob[j]);
                            if(wei[j]>max){ max=wei[j]; }
                        }
                        for(sum=0.0, j=1; j<=grid; j++) sum+=exp(wei[j]-max);
                        
                        indexx(grid,mu,indx2);
                        for(k=1; k<=grid; k++){
                            j=indx2[k];
                            eta=exp(wei[j]-max)/sum;
                            
                            if(z>mu[j]) dmu=beta[j]/alpha[j]*exp((beta[j]-1.0)*log(fabs(z-mu[j])/alpha[j]));
                            else dmu=-beta[j]/alpha[j]*exp((beta[j]-1.0)*log(fabs(z-mu[j])/alpha[j]));
                            
                            dalpha=-1.0/alpha[j]+beta[j]/alpha[j]*exp(beta[j]*log(fabs(z-mu[j])/alpha[j]));
                            dbeta=1.0/beta[j]+1.0/beta[j]/beta[j]*derivative_gamma(1.0/beta[j])
                            -exp(beta[j]*log(fabs(z-mu[j])/alpha[j]))*log(fabs(z-mu[j])/alpha[j]);
                            
                            mu[j]+=delta*eta*dmu;
                            alpha[j]+=delta*eta*dalpha;
                            beta[j]+=delta*eta*dbeta;
                            logprob[j]+=delta*(eta-prob[j]);
                            
                            if(mu[j]<-20.0) mu[j]=-20.0;
                            if(mu[j]>40.0) mu[j]=40.0;
                            if(alpha[j]<0.1) alpha[j]=0.1;
                            if(alpha[j]>10.0) alpha[j]=10.0;
                            if(beta[j]<0.1) beta[j]=0.1;
                            if(beta[j]>10.0) beta[j]=10.0;
                        }
                    } /* end one epoach */
                    
                    if(iter%1==0) Rprintf("iter=%d delta=%g %g %g %g %g: %g %g %g %g\n",
                                          iter,delta,prob[1],mu[1],alpha[1],beta[1],prob[2],mu[2],alpha[2],beta[2]);
                    
                }      /* end for iterations */
                
                
                for(sum=0.0,j=1; j<=grid; j++) sum+=exp(logprob[j]);
                for(j=1; j<=grid; j++) prob[j]=exp(logprob[j])/sum;
                
                for(entropy=0.0,i=1; i<=datanum; i++){
                    max=-1.0e+100;
                    for(j=1; j<=grid; j++){
                        wei[j]=dlogGnormal(z,mu[j],alpha[j],beta[j])+log(prob[j]);
                        if(wei[j]>max){ max=wei[j]; }
                    }
                    for(sum=0.0, j=1; j<=grid; j++) sum+=exp(wei[j]-max);
                    entropy+=log(sum)+max;
                }
                entropy*=-1.0/datanum;
                entropy+=0.5*(grid*3+grid-1)*log(1.0*datanum)/datanum;
                
                
                
                if(entropy<=bestentropy){
                    bestentropy=entropy;
                    bestgrid=grid;
                    for(j=1; j<=grid; j++){ bestprob[j]=prob[j]; bestmu[j]=mu[j];
                        bestalpha[j]=alpha[j]; bestbeta[j]=beta[j]; }
                }
            }  /* end for different grid size */
    
    /* calculating the distance between the mixture components */
    for(k=1; k<=bestgrid; k++){
        MCMCGnormal(bestmu[k],bestalpha[k],bestbeta[k],samsize,rez);
        for(j=1; j<=samsize; j++) resam[k][j]=rez[j];
    }
    indexx(bestgrid, bestmu, indx2);
    for(j=1; j<bestgrid; j++){
        k1=indx2[j]; k2=indx2[j+1];
        tep1=tep2=0.0;
        for(i=1; i<=samsize; i++){
            tep1+=dlogGnormal(resam[k1][i],bestmu[k1],bestalpha[k1],bestbeta[k1])-
            dlogGnormal(resam[k1][i],bestmu[k2],bestalpha[k2],bestbeta[k2]);
            tep2+=dlogGnormal(resam[k2][i],bestmu[k2],bestalpha[k2],bestbeta[k2])-
            dlogGnormal(resam[k2][i],bestmu[k1],bestalpha[k1],bestbeta[k1]);
        }
        if(tep1<tep2) dis[j]=tep1/samsize;
        else dis[j]=tep2/samsize;
    }
    
    
    /* determine the number of components for f_0 */
    dif[1]=1.0;
    for(j=2; j<bestgrid; j++) dif[j]=(dis[j]-dis[j-1])/dis[j-1];
        dif[bestgrid]=-1.0;
        com=1; ok=0;
        while(ok==0 && com<bestgrid-1){
            if(dif[com]*dif[com+1]<0.0) ok=1;
            else com++;
        }
    /*
     max=0.0; k=1;
     for(j=2; j<bestgrid; j++)
     if(dif[j]>max){ max=dif[j]; k=j; }
     if(k<com) com=k;
     */
    
    
    
    // com=bestgrid-1;
    com=1;
    
    /* calculating FDR */
    for(i=1; i<=datanum; i++){
        z=datax[indx[i]];
        for(sum=0.0,m=1; m<=com; m++){
            k=indx2[m];
            if(z>bestmu[k]){
                un=exp(bestbeta[k]*log((z-bestmu[k])/bestalpha[k]));
                tep=(0.5-0.5*gammp(1.0/bestbeta[k],un))*bestprob[k];
            } else if(z<bestmu[k]){
                un=exp(bestbeta[k]*log((bestmu[k]-z)/bestalpha[k]));
                tep=(0.5+0.5*gammp(1.0/bestbeta[k],un))*bestprob[k];
            } else tep=0.5*bestprob[k];
            sum+=tep;
        }
        FDR[i]=sum;
        
        for(sum=0.0, k=1; k<=bestgrid; k++){
            if(z>bestmu[k]){
                un=exp(bestbeta[k]*log((z-bestmu[k])/bestalpha[k]));
                tep=(0.5-0.5*gammp(1.0/bestbeta[k],un))*bestprob[k];
            } else if(z<bestmu[k]){
                un=exp(bestbeta[k]*log((bestmu[k]-z)/bestalpha[k]));
                tep=(0.5+0.5*gammp(1.0/bestbeta[k],un))*bestprob[k];
            } else tep=0.5*bestprob[k];
            sum+=tep;
        }
        FDR[i]/=sum;
    }
    
    min=FDR[1]; Q[1]=min;
    for(i=2; i<=datanum; i++){
        if(FDR[i]<min) min=FDR[i];
        Q[i]=min;
    }
    
    
    
    
    
    /* alternative method of FDR calculation */
    for(i=1; i<=datanum; i++){
        z=datax[indx[i]];
        for(sum=0.0,m=1; m<=com; m++){
            k=indx2[m];
            if(z>bestmu[k]){
                un=exp(bestbeta[k]*log((z-bestmu[k])/bestalpha[k]));
                tep=(0.5-0.5*gammp(1.0/bestbeta[k],un))*bestprob[k];
            } else if(z<bestmu[k]){
                un=exp(bestbeta[k]*log((bestmu[k]-z)/bestalpha[k]));
                tep=(0.5+0.5*gammp(1.0/bestbeta[k],un))*bestprob[k];
            } else tep=0.5*bestprob[k];
            sum+=tep;
        }
        FDR[i]=sum*datanum/(datanum-i+1);
        if(FDR[i]>1.0) FDR[i]=1.0;
    }
    
    min=FDR[1]; Q[1]=min;
    for(i=2; i<=datanum; i++){
        if(FDR[i]<min) min=FDR[i];
        Q[i]=min;
    }
    
    for(i=1; i<=datanum; i++){
        Qmat[i][1]=row[indx[i]]; Qmat[i][2]=col[indx[i]]; Qmat[i][3]=datax[indx[i]]; Qmat[i][4]=Q[i];
    }
    
    
    free_dvector(mu,1,GRID);
    free_dvector(alpha,1,GRID);
    free_dvector(beta,1,GRID);
    free_dvector(prob,1,GRID);
    free_dvector(wei,1,GRID);
    free_ivector(sam,1,MAX_NUM);
    free_ivector(indx,1,MAX_NUM);
    free_ivector(indx2,1,GRID);
    free_dvector(FDR,1,MAX_NUM);
    free_dvector(Q,1,MAX_NUM);
    free_dvector(bestmu,1,GRID);
    free_dvector(bestalpha,1,GRID);
    free_dvector(bestbeta,1,GRID);
    free_dvector(bestprob,1,GRID);
    free_dvector(logprob,1,GRID);
    free_dmatrix(resam,1,GRID,1,samsize);
    free_dvector(rez,1,samsize);
    free_dvector(dis,1,GRID);
    free_dvector(dif,1,GRID);
    
    return 0;
}

int psi_calculation(MaxNei, DataNum, DataP, M, Qmat, COV, qthreshold, row, col, score, psi)
int M, *row, *col, MaxNei, DataNum, DataP;
double **Qmat, **COV, qthreshold, *score, *psi;
{

    int *x,*y,*indx, *indx0, **sep, *S;
    int i,j,k,k1,k2,l,m,v, dim;
    double **A,**IA,*z,eta;
    double un, qvalue, qscore, zscore, **SS;
    FILE *ins;
    
    x=ivector(1,DataP);
    y=ivector(1,DataP);
    A=dmatrix(1,DataP,1,DataP);
    IA=dmatrix(1,DataP,1,DataP);
    z=dvector(1,DataP);
    indx=ivector(1,DataP);
    indx0=ivector(1,DataP);
    sep=imatrix(1,DataP,0,DataP);
    SS=dmatrix(1,DataP,1,DataP);
    S=ivector(1,DataP);
    
    k=M; qvalue=Qmat[k][4];
    while(qvalue<qthreshold && k>1){
        k--;
        qvalue=Qmat[k][4];
    }
    
    if(k==M) { Rprintf("The given q-threshold value is too small\n");  goto loop; }
    else if(k==1){ Rprintf("The given q-threshold value is too big\n");  goto loop; }
    else {
        qscore=0.5*(fabs(Qmat[k][3])+fabs(Qmat[k+1][3]));
        // qscore=fabs(Qmat[k+1][3]);
        Rprintf("q-score=%g\n", qscore);
    }
    
    Rprintf("Thresholding is done\n");
    
    /*** determine significant neighbors  ****/
    for(i=1; i<=DataP; i++)
        for(j=0; j<=DataP; j++) sep[i][j]=0;
            
            
            for(i=1; i<=DataP*(DataP-1)/2; i++){
                
                k1=row[i]; k2=col[i]; un=score[i];
                
                if(un>qscore){
                    sep[k1][0]+=1;
                    sep[k1][sep[k1][0]]=k2;
                    sep[k2][0]+=1;
                    sep[k2][sep[k2][0]]=k1;
                    
                    SS[k1][sep[k1][0]]=un;
                    SS[k2][sep[k2][0]]=un;
                }
            }
    
    
    for(i=1; i<=DataP; i++){
        
        if(sep[i][0]>MaxNei){
            
            m=sep[i][0];
            indexx(m,SS[i],indx);
            for(j=1; j<=m; j++) S[j]=sep[i][j];
            
            k=m-MaxNei;
            for(j=1; j<=MaxNei; j++) sep[i][j]=S[indx[k+j]];
            sep[i][0]=MaxNei;
        }
    }
    
    
    v=0;
    ins=fopen("sim0.pcor.est.ini","w");
    for(i=1; i<=DataP; i++){
        for(j=i+1; j<=DataP; j++){
            
            if(sep[i][0]<sep[j][0]){
                m=sep[i][0];
                for(k=1; k<=m; k++) x[k]=sep[i][k];
            }
            else{
                m=sep[j][0];
                for(k=1; k<=m; k++) x[k]=sep[j][k];
            }
            
            if(m>1) indexx_integer(m,x,indx);
            else{
                if(m==1) indx[1]=1;
            }
            
            y[1]=i; y[2]=j; dim=2;
            
            if(m>=1){
                k=1;
                if(x[indx[k]]!=i && x[indx[k]]!=j){ dim++; y[dim]=x[indx[k]];}
                while(k<m){
                    k++;
                    if(x[indx[k]]!=x[indx[k-1]] && x[indx[k]]!=i && x[indx[k]]!=j){ dim++; y[dim]=x[indx[k]]; }
                }
            }
            
            for(k=1; k<=dim; k++)
                for(l=1; l<=dim; l++) A[k][l]=COV[y[k]][y[l]];
            
            for(k=1; k<=dim; k++) A[k][k]+=1.0e-3;
            matrix_inverse(A,IA,dim);
            
            eta=-IA[1][2]/sqrt(IA[1][1]*IA[2][2]);
            zscore=0.5*log((1+eta)/(1-eta))*sqrt(DataNum-dim-1);
            
            /*
             un=-fabs(zscore);
             q=gauss_cdf_ln(un)+log(2.0);
             pscore=-1.0*inverse_normal_cdf_log(q);
             */
            
            fprintf(ins, " %d %d %g\n", i,j,zscore);
            
            v++;
            psi[v]=zscore;
        }
        if(i%100==0) Rprintf("psi-score i=%d\n", i);
    }
    
    fclose(ins);
    
    
loop:free_ivector(x,1,DataP);
    free_ivector(y,1,DataP);
    free_dmatrix(A,1,DataP,1,DataP);
    free_dmatrix(IA,1,DataP,1,DataP);
    free_dvector(z,1,DataP);
    free_ivector(indx,1,DataP);
    free_ivector(indx0,1,DataP);
    free_imatrix(sep,1,DataP,0,DataP);
    free_dmatrix(SS,1,DataP,1,DataP);
    free_ivector(S,1,DataP);
    
    return 0;
}

int pcorsel(datanum,row,col,datax, Qmat)
int datanum,*row, *col;
double *datax, **Qmat;
{
    int  grid, bestgrid, rep;
    int *indx, *indx2,  *sam, i, j, k, k1, k2, m, com, iter, ok;
    double **resam,*rez,*dis,*dif,*mu,*alpha,*beta,*prob,*logprob,*wei,*FDR,*Q, entropy;
    double delta,z,un,max,min,sum,tep,tep1,tep2,eta,dmu,dalpha,dbeta;
    double bestentropy, *bestprob, *bestmu, *bestalpha, *bestbeta;

    
    
    
    
    
    /* initialize the random number generator */
    /*
     ltime=time(NULL);
     stime=(unsigned int)ltime/2;
     srand(stime);
     
     ins=fopen("aa.log","a");
     fprintf(ins, "random seed=%d\n", stime);
     fclose(ins);
     */
    
    mu=dvector(1,GRID);
    alpha=dvector(1,GRID);
    beta=dvector(1,GRID);
    prob=dvector(1,GRID);
    wei=dvector(1,GRID);
    sam=ivector(1,MAX_NUM);
    indx=ivector(1,MAX_NUM);
    indx2=ivector(1,GRID);
    FDR=dvector(1,MAX_NUM);
    Q=dvector(1,MAX_NUM);
    bestmu=dvector(1,GRID);
    bestalpha=dvector(1,GRID);
    bestbeta=dvector(1,GRID);
    bestprob=dvector(1,GRID);
    logprob=dvector(1,GRID);
    resam=dmatrix(1,GRID,1,samsize);
    rez=dvector(1,samsize);
    dis=dvector(1,GRID);
    dif=dvector(1,GRID);
    
    
    if(datanum>MAX_NUM) Rprintf("Data size is TOO LARGE in pcorsel.c\n");
        else Rprintf("working subdata size=%d\n", datanum);
            
            indexx(datanum,datax,indx);
            
            grid=GRID+1;  bestentropy=1.0e+100;
            while(grid>mingrid){
                
                grid--;
                
                /* Initialization */
                mu[1]=0.0; alpha[1]=0.0; k=0;
                for(i=1; i<=(int)(datanum*0.99); i++){
                    mu[1]+=datax[indx[i]];
                    alpha[1]+=datax[indx[i]]*datax[indx[i]];
                    k++;
                }
                mu[1]=mu[1]/k;
                alpha[1]=sqrt(alpha[1]/k-mu[1]*mu[1])*sqrt(2.0);
                beta[1]=2.0; prob[1]=0.99;
                
                for(j=2; j<=grid; j++){
                    tep=(1-prob[1])/(grid-1);
                    mu[j]=0.0; alpha[j]=0.0; k=0;
                    for(i=(int)(datanum*(prob[1]+(j-2)*tep)); i<=(int)(datanum*(prob[1]+(j-1)*tep)); i++){
                        mu[j]+=datax[indx[i]];
                        alpha[j]+=datax[indx[i]]*datax[indx[i]];
                        k++;
                    }
                    mu[j]=mu[j]/k;
                    alpha[j]=sqrt(alpha[j]/k-mu[j]*mu[j])*sqrt(2.0);
                    beta[j]=2.0; prob[j]=tep;
                }
                for(j=1; j<=grid; j++){
                    logprob[j]=log(prob[j])-log(prob[grid]);
                }
                
                
                for(iter=1; iter<=100; iter++){
                    
                    random_order(sam,datanum);
                    for(i=1; i<=datanum; i++){
                        
                        rep=datanum*(iter-1)+i;
                        if(rep<=WARM*stepscale) delta=rho;
                        else delta=rho*exp(-log(1.0*(rep-(WARM-1)*stepscale)/stepscale));
                        
                        for(sum=0.0,j=1; j<=grid; j++) sum+=exp(logprob[j]);
                        for(j=1; j<=grid; j++) prob[j]=exp(logprob[j])/sum;
                        
                        z=datax[sam[i]];
                        max=-1.0e+100;
                        for(j=1; j<=grid; j++){
                            wei[j]=dlogGnormal(z,mu[j],alpha[j],beta[j])+log(prob[j]);
                            if(wei[j]>max){ max=wei[j];}
                        }
                        for(sum=0.0, j=1; j<=grid; j++) sum+=exp(wei[j]-max);
                        
                        indexx(grid,mu,indx2);
                        for(k=1; k<=grid; k++){
                            j=indx2[k];
                            eta=exp(wei[j]-max)/sum;
                            
                            if(z>mu[j]) dmu=beta[j]/alpha[j]*exp((beta[j]-1.0)*log(fabs(z-mu[j])/alpha[j]));
                            else dmu=-beta[j]/alpha[j]*exp((beta[j]-1.0)*log(fabs(z-mu[j])/alpha[j]));
                            
                            dalpha=-1.0/alpha[j]+beta[j]/alpha[j]*exp(beta[j]*log(fabs(z-mu[j])/alpha[j]));
                            dbeta=1.0/beta[j]+1.0/beta[j]/beta[j]*derivative_gamma(1.0/beta[j])
                            -exp(beta[j]*log(fabs(z-mu[j])/alpha[j]))*log(fabs(z-mu[j])/alpha[j]);
                            
                            mu[j]+=delta*eta*dmu;
                            alpha[j]+=delta*eta*dalpha;
                            beta[j]+=delta*eta*dbeta;
                            logprob[j]+=delta*(eta-prob[j]);
                            
                            if(mu[j]<-20.0) mu[j]=-20.0;
                            if(mu[j]>40.0) mu[j]=40.0;
                            if(alpha[j]<0.1) alpha[j]=0.1;
                            if(alpha[j]>10.0) alpha[j]=10.0;
                            if(beta[j]<0.1) beta[j]=0.1;
                            if(beta[j]>10.0) beta[j]=10.0;
                        }
                    } /* end one epoach */
                    
                    if(iter%1==0) Rprintf("iter=%d delta=%g %g %g %g %g: %g %g %g %g\n",
                                          iter,delta,prob[1],mu[1],alpha[1],beta[1],prob[2],mu[2],alpha[2],beta[2]);
                    
                }      /* end for iterations */
                
                
                for(sum=0.0,j=1; j<=grid; j++) sum+=exp(logprob[j]);
                for(j=1; j<=grid; j++) prob[j]=exp(logprob[j])/sum;
                
                for(entropy=0.0,i=1; i<=datanum; i++){
                    max=-1.0e+100;
                    for(j=1; j<=grid; j++){
                        wei[j]=dlogGnormal(z,mu[j],alpha[j],beta[j])+log(prob[j]);
                        if(wei[j]>max){ max=wei[j];  }
                    }
                    for(sum=0.0, j=1; j<=grid; j++) sum+=exp(wei[j]-max);
                    entropy+=log(sum)+max;
                }
                entropy*=-1.0/datanum;
                entropy+=0.5*(grid*3+grid-1)*log(1.0*datanum)/datanum;
                
                
                
                if(entropy<=bestentropy){
                    bestentropy=entropy;
                    bestgrid=grid;
                    for(j=1; j<=grid; j++){ bestprob[j]=prob[j]; bestmu[j]=mu[j];
                        bestalpha[j]=alpha[j]; bestbeta[j]=beta[j]; }
                }
            }  /* end for different grid size */
    
    /* calculating the distance between the mixture components */
    for(k=1; k<=bestgrid; k++){
        MCMCGnormal(bestmu[k],bestalpha[k],bestbeta[k],samsize,rez);
        for(j=1; j<=samsize; j++) resam[k][j]=rez[j];
    }
    indexx(bestgrid, bestmu, indx2);
    for(j=1; j<bestgrid; j++){
        k1=indx2[j]; k2=indx2[j+1];
        tep1=tep2=0.0;
        for(i=1; i<=samsize; i++){
            tep1+=dlogGnormal(resam[k1][i],bestmu[k1],bestalpha[k1],bestbeta[k1])-
            dlogGnormal(resam[k1][i],bestmu[k2],bestalpha[k2],bestbeta[k2]);
            tep2+=dlogGnormal(resam[k2][i],bestmu[k2],bestalpha[k2],bestbeta[k2])-
            dlogGnormal(resam[k2][i],bestmu[k1],bestalpha[k1],bestbeta[k1]);
        }
        if(tep1<tep2) dis[j]=tep1/samsize;
        else dis[j]=tep2/samsize;
    }
    
    
    /* determine the number of components for f_0 */
    dif[1]=1.0;
    for(j=2; j<bestgrid; j++) dif[j]=(dis[j]-dis[j-1])/dis[j-1];
        dif[bestgrid]=-1.0;
        com=1; ok=0;
        while(ok==0 && com<bestgrid-1){
            if(dif[com]*dif[com+1]<0.0) ok=1;
            else com++;
        }
    /*
     max=0.0; k=1;
     for(j=2; j<bestgrid; j++)
     if(dif[j]>max){ max=dif[j]; k=j; }
     if(k<com) com=k;
     */
    
    
    
    // com=bestgrid-1;
    com=1;
    
    /* calculating FDR */
    for(i=1; i<=datanum; i++){
        z=datax[indx[i]];
        for(sum=0.0,m=1; m<=com; m++){
            k=indx2[m];
            if(z>bestmu[k]){
                un=exp(bestbeta[k]*log((z-bestmu[k])/bestalpha[k]));
                tep=(0.5-0.5*gammp(1.0/bestbeta[k],un))*bestprob[k];
            } else if(z<bestmu[k]){
                un=exp(bestbeta[k]*log((bestmu[k]-z)/bestalpha[k]));
                tep=(0.5+0.5*gammp(1.0/bestbeta[k],un))*bestprob[k];
            } else tep=0.5*bestprob[k];
            sum+=tep;
        }
        FDR[i]=sum;
        
        for(sum=0.0, k=1; k<=bestgrid; k++){
            if(z>bestmu[k]){
                un=exp(bestbeta[k]*log((z-bestmu[k])/bestalpha[k]));
                tep=(0.5-0.5*gammp(1.0/bestbeta[k],un))*bestprob[k];
            } else if(z<bestmu[k]){
                un=exp(bestbeta[k]*log((bestmu[k]-z)/bestalpha[k]));
                tep=(0.5+0.5*gammp(1.0/bestbeta[k],un))*bestprob[k];
            } else tep=0.5*bestprob[k];
            sum+=tep;
        }
        FDR[i]/=sum;
    }
    
    min=FDR[1]; Q[1]=min;
    for(i=2; i<=datanum; i++){
        if(FDR[i]<min) min=FDR[i];
        Q[i]=min;
    }
    
    
    
    
    /* alternative method of FDR calculation */
    for(i=1; i<=datanum; i++){
        z=datax[indx[i]];
        for(sum=0.0,m=1; m<=com; m++){
            k=indx2[m];
            if(z>bestmu[k]){
                un=exp(bestbeta[k]*log((z-bestmu[k])/bestalpha[k]));
                tep=(0.5-0.5*gammp(1.0/bestbeta[k],un))*bestprob[k];
            } else if(z<bestmu[k]){
                un=exp(bestbeta[k]*log((bestmu[k]-z)/bestalpha[k]));
                tep=(0.5+0.5*gammp(1.0/bestbeta[k],un))*bestprob[k];
            } else tep=0.5*bestprob[k];
            sum+=tep;
        }
        FDR[i]=sum*datanum/(datanum-i+1);
        if(FDR[i]>1.0) FDR[i]=1.0;
    }
    
    min=FDR[1]; Q[1]=min;
    for(i=2; i<=datanum; i++){
        if(FDR[i]<min) min=FDR[i];
        Q[i]=min;
    }
    
    
    for(i=1; i<=datanum; i++){
        Qmat[i][1]=row[indx[i]]; Qmat[i][2]=col[indx[i]]; Qmat[i][3]=datax[indx[i]]; Qmat[i][4]=Q[i];
    }
    
    
    
    free_dvector(mu,1,GRID);
    free_dvector(alpha,1,GRID);
    free_dvector(beta,1,GRID);
    free_dvector(prob,1,GRID);
    free_dvector(wei,1,GRID);
    free_ivector(sam,1,MAX_NUM);
    free_ivector(indx,1,MAX_NUM);
    free_ivector(indx2,1,GRID);
    free_dvector(FDR,1,MAX_NUM);
    free_dvector(Q,1,MAX_NUM);
    free_dvector(bestmu,1,GRID);
    free_dvector(bestalpha,1,GRID);
    free_dvector(bestbeta,1,GRID);
    free_dvector(bestprob,1,GRID);
    free_dvector(logprob,1,GRID);
    free_dmatrix(resam,1,GRID,1,samsize);
    free_dvector(rez,1,samsize);
    free_dvector(dis,1,GRID);
    free_dvector(dif,1,GRID);
    
    return 0;
}



void equSA_sub(double DataX[],
            int *iMaxNei,
            int *iDataNum,
            int *iDataP,
            double *ALPHA1)
{
#define pi 3.14159265
    int  DATA_NUM, DATA_NUM_PLUS, *no, *row, *col, *subrow, *subcol,*indx, **E;
    int i, j, k, m, M, m0;
    double *datax, **iDataX, **COV, **CORR, *score, *subscore, *psi, **Qmat;
    double q, sum, un;
    static int ratio;
    // static double ALPHA1=0.05, ALPHA2=0.05;
    

    
    

    
    
    
    

    
    DATA_NUM=*iDataP*(*iDataP-1)/2;
    DATA_NUM_PLUS=*iDataP*(*iDataP+1)/2;
    ratio=ceil(DATA_NUM/100000.0);
    
    
    datax=dvector(1,DATA_NUM_PLUS);
    score=dvector(1,DATA_NUM_PLUS);
    row=ivector(1,DATA_NUM_PLUS);
    col=ivector(1,DATA_NUM_PLUS);
    no=ivector(1,DATA_NUM);
    indx=ivector(1,DATA_NUM);
    iDataX=dmatrix(1,*iDataNum,1,*iDataP);
    COV=dmatrix(1,*iDataP,1,*iDataP);
    CORR=dmatrix(1,*iDataP,1,*iDataP);
    subrow=ivector(1,DATA_NUM);
    subcol=ivector(1,DATA_NUM);
    subscore=dvector(1,DATA_NUM);
    Qmat=dmatrix(1,DATA_NUM,1,4);
    psi=dvector(1,DATA_NUM);
    E=imatrix(1,*iDataP,1,*iDataP);
    
    
    
    for(i=1; i<=*iDataNum; i++){
        for(j=1; j<=*iDataP; j++)
            iDataX[i][j]=DataX[(i-1)*(*iDataP)+(j-1)];
    }
    
    
    
    
    
    
    
    /*** centralization ****/
    for(j=1; j<=*iDataP; j++){
        for(sum=0.0,i=1; i<=*iDataNum; i++) sum+=iDataX[i][j];
        un=sum/(*iDataNum);
        for(i=1; i<=*iDataNum; i++) iDataX[i][j]-=un;
    }
    
    /*** calculate covariance matrix ***/
    m=0;
    for(i=1; i<=*iDataP; i++){
        for(j=i; j<=*iDataP; j++){
            COV[i][j]=0.0;
            for(k=1; k<=*iDataNum; k++) COV[i][j]+=iDataX[k][i]*iDataX[k][j];
            COV[i][j]/=*iDataNum;
            if(i!=j) COV[j][i]=COV[i][j];
            
            m++;
            row[m]=i; col[m]=j; datax[m]=COV[i][j];
        }
        if(i%100==0) Rprintf("i=%d\n", i);
    }
    if(m!=DATA_NUM_PLUS){ Rprintf("Matrix error I in corz.c\n"); }
    

    
    
    
    for(i=1; i<=*iDataP; i++)
        for(j=i; j<=*iDataP; j++){
            CORR[i][j]=COV[i][j]/sqrt(COV[i][i]*COV[j][j]);
            if(i!=j) CORR[j][i]=CORR[i][j];
        }
    
    sum=0.0; k=0;
    for(i=1; i<=*iDataP; i++)
        for(j=i+1; j<=*iDataP; j++){
            k++;
            un=0.5*log((1+CORR[i][j])/(1-CORR[i][j]))*sqrt(*iDataNum-3.0);
            
            row[k]=i; col[k]=j; datax[k]=un;
            no[k]=k;
        }
    if(k!=DATA_NUM) Rprintf("Matrix size error\n");
    

    for(k=1; k<=DATA_NUM; k++){
        
        un=-fabs(datax[k]);
        q=gauss_cdf_ln(un);
        sum=q+log(2.0);
        score[k]=-1.0*inverse_normal_cdf_log(sum);
        
    }

    
    
    /*** subsampling  **************/
    if(ratio>1){
        indexx(DATA_NUM, score, indx);
        M=DATA_NUM/ratio;
        m0=DATA_NUM-M*ratio;
        
        for(i=1; i<=M; i++){
            
            k=0;
            GetRNGstate();
            while(k<=0 || k>ratio) k=floor(unif_rand()*ratio)+1;
            PutRNGstate();
            j=(i-1)*ratio+k;
            subrow[i]=row[indx[j]]; subcol[i]=col[indx[j]]; subscore[i]=score[indx[j]];
        }
        if(m0>0){
            k=0;
            GetRNGstate();
            while(k<=0 || k>m0) k=floor(unif_rand()*m0)+1;
            PutRNGstate();
            j=M*ratio+k;
            M++;
            subrow[M]=row[indx[j]]; subcol[M]=col[indx[j]]; subscore[M]=score[indx[j]];
        }
    }
    else{
        M=DATA_NUM;
        for(i=1; i<=M; i++){ subrow[i]=row[i]; subcol[i]=col[i]; subscore[i]=score[i]; }
    }
    

    
    corsel(M,subrow,subcol,subscore,Qmat);
    
    /*** psi-score calculation *****/
    
   psi_calculation(*iMaxNei, *iDataNum, *iDataP,M,Qmat,COV,*ALPHA1,row,col,score,psi);
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    free_dvector(datax,1,DATA_NUM_PLUS);
    free_dvector(score,1,DATA_NUM_PLUS);
    free_ivector(row,1,DATA_NUM_PLUS);
    free_ivector(col,1,DATA_NUM_PLUS);
    free_ivector(no,1,DATA_NUM);
    free_dmatrix(iDataX,1,*iDataNum,1,*iDataP);
    free_dmatrix(COV,1,*iDataP,1,*iDataP);
    free_dmatrix(CORR,1,*iDataP,1,*iDataP);
    free_ivector(indx,1,DATA_NUM);
    free_ivector(subrow,1,DATA_NUM);
    free_ivector(subcol,1,DATA_NUM);
    free_dvector(subscore,1,DATA_NUM);
    free_dmatrix(Qmat,1,DATA_NUM,1,4);
    free_dvector(psi,1,DATA_NUM);
    free_imatrix(E,1,*iDataP,1,*iDataP);
    
    //return 0;
}













void equSA1(double DataX[],
           int *iMaxNei,
           int *iDataNum,
           int *iDataP,
           double *ALPHA1,
           double *ALPHA2)
{
   #define pi 3.14159265
    FILE *ins;
    int  DATA_NUM, DATA_NUM_PLUS, *no, *row, *col, *subrow, *subcol,*indx, **E;
    int i, j, k, m, M, m0;
    double *datax, **iDataX, **COV, **CORR, *score, *subscore, *psi, **Qmat;
    double q, sum, un, qscore;
    static int ratio;
   // static double ALPHA1=0.05, ALPHA2=0.05;
    


    

    

    
    DATA_NUM=*iDataP*(*iDataP-1)/2;
    DATA_NUM_PLUS=*iDataP*(*iDataP+1)/2;
    ratio=ceil(DATA_NUM/100000.0);
    
    
    datax=dvector(1,DATA_NUM_PLUS);
    score=dvector(1,DATA_NUM_PLUS);
    row=ivector(1,DATA_NUM_PLUS);
    col=ivector(1,DATA_NUM_PLUS);
    no=ivector(1,DATA_NUM);
    indx=ivector(1,DATA_NUM);
    iDataX=dmatrix(1,*iDataNum,1,*iDataP);
    COV=dmatrix(1,*iDataP,1,*iDataP);
    CORR=dmatrix(1,*iDataP,1,*iDataP);
    subrow=ivector(1,DATA_NUM);
    subcol=ivector(1,DATA_NUM);
    subscore=dvector(1,DATA_NUM);
    Qmat=dmatrix(1,DATA_NUM,1,4);
    psi=dvector(1,DATA_NUM);
    E=imatrix(1,*iDataP,1,*iDataP);
    
    
    
    for(i=1; i<=*iDataNum; i++){
        for(j=1; j<=*iDataP; j++)
        iDataX[i][j]=DataX[(i-1)*(*iDataP)+(j-1)];
    }
    
    
    
    
    
    
    
    /*** centralization ****/
    for(j=1; j<=*iDataP; j++){
        for(sum=0.0,i=1; i<=*iDataNum; i++) sum+=iDataX[i][j];
        un=sum/(*iDataNum);
        for(i=1; i<=*iDataNum; i++) iDataX[i][j]-=un;
    }
    
    /*** calculate covariance matrix ***/
    m=0;
    for(i=1; i<=*iDataP; i++){
        for(j=i; j<=*iDataP; j++){
            COV[i][j]=0.0;
            for(k=1; k<=*iDataNum; k++) COV[i][j]+=iDataX[k][i]*iDataX[k][j];
            COV[i][j]/=*iDataNum;
            if(i!=j) COV[j][i]=COV[i][j];
            
            m++;
            row[m]=i; col[m]=j; datax[m]=COV[i][j];
        }
        if(i%100==0) Rprintf("i=%d\n", i);
    }
    if(m!=DATA_NUM_PLUS){ Rprintf("Matrix error I in corz.c\n"); }
    
    
    
    
    for(i=1; i<=*iDataP; i++)
        for(j=i; j<=*iDataP; j++){
            CORR[i][j]=COV[i][j]/sqrt(COV[i][i]*COV[j][j]);
            if(i!=j) CORR[j][i]=CORR[i][j];
        }
    
    sum=0.0; k=0;
    for(i=1; i<=*iDataP; i++)
        for(j=i+1; j<=*iDataP; j++){
            k++;
            un=0.5*log((1+CORR[i][j])/(1-CORR[i][j]))*sqrt(*iDataNum-3.0);
            
            row[k]=i; col[k]=j; datax[k]=un;
            no[k]=k;
        }
    if(k!=DATA_NUM) Rprintf("Matrix size error\n");
    
    for(k=1; k<=DATA_NUM; k++){
        
        un=-fabs(datax[k]);
        q=gauss_cdf_ln(un);
        sum=q+log(2.0);
        score[k]=-1.0*inverse_normal_cdf_log(sum);
        
    }
 
    
    /*** subsampling  **************/
    if(ratio>1){
        indexx(DATA_NUM, score, indx);
        M=DATA_NUM/ratio;
        m0=DATA_NUM-M*ratio;
        
        for(i=1; i<=M; i++){
            
            k=0;
            GetRNGstate();
            while(k<=0 || k>ratio) k=floor(unif_rand()*ratio)+1;
            PutRNGstate();
            j=(i-1)*ratio+k;
            subrow[i]=row[indx[j]]; subcol[i]=col[indx[j]]; subscore[i]=score[indx[j]];
        }
        if(m0>0){
            k=0;
            GetRNGstate();
            while(k<=0 || k>m0) k=floor(unif_rand()*m0)+1;
            PutRNGstate();
            j=M*ratio+k;
            M++;
            subrow[M]=row[indx[j]]; subcol[M]=col[indx[j]]; subscore[M]=score[indx[j]];
        }
    }
    else{
        M=DATA_NUM;
        for(i=1; i<=M; i++){ subrow[i]=row[i]; subcol[i]=col[i]; subscore[i]=score[i]; }
    }
    
    corsel(M,subrow,subcol,subscore,Qmat);
    
    /*** psi-score calculation *****/
    
 psi_calculation(*iMaxNei, *iDataNum, *iDataP,M,Qmat,COV,*ALPHA1,row,col,score,psi);
    
    /*** psi-score transformation  ***/
    for(i=1; i<=DATA_NUM; i++){
        un=-fabs(psi[i]);
        q=gauss_cdf_ln(un)+log(2.0);
        psi[i]=-1.0*inverse_normal_cdf_log(q);
    }
    
    
    
    /*** subsampling  for psi-score **************/
    if(ratio>1){
        indexx(DATA_NUM, psi, indx);
        M=DATA_NUM/ratio;
        m0=DATA_NUM-M*ratio;
        
        for(i=1; i<=M; i++){
            
            k=0;
            GetRNGstate();
            while(k<=0 || k>ratio) k=floor(unif_rand()*ratio)+1;
            PutRNGstate();
            j=(i-1)*ratio+k;
            subrow[i]=row[indx[j]]; subcol[i]=col[indx[j]]; subscore[i]=psi[indx[j]];
        }
        if(m0>0){
            k=0;
            GetRNGstate();
            while(k<=0 || k>m0) k=floor(unif_rand()*m0)+1;
            PutRNGstate();
            j=M*ratio+k;
            M++;
            subrow[M]=row[indx[j]]; subcol[M]=col[indx[j]]; subscore[M]=psi[indx[j]];
        }
    }
    else{
        M=DATA_NUM;
        for(i=1; i<=M; i++){ subrow[i]=row[i]; subcol[i]=col[i]; subscore[i]=psi[i]; }
    }
    
    
    pcorsel(M,subrow,subcol,subscore,Qmat);
    
    k=M; q=Qmat[k][4];
    while(q<*ALPHA2 && k>1){
        k--;
        q=Qmat[k][4];
    }
    
    if(k==M) { Rprintf("The given q-threshold value is too small\n");  goto loop;}
    else if(k==1){ Rprintf("The given q-threshold value is too big\n"); goto loop; }
    else {
        qscore=0.5*(fabs(Qmat[k][3])+fabs(Qmat[k+1][3]));
        // qscore=fabs(Qmat[k+1][3]);
        Rprintf("threshold psi-score=%g\n", qscore);
    }
    
    for(i=1; i<=*iDataP; i++)
        for(j=1; j<=*iDataP; j++) E[i][j]=0;
    
    for(i=1; i<=DATA_NUM; i++)
        if(psi[i]>qscore){
            E[row[i]][col[i]]=1;
            E[col[i]][row[i]]=1;
        }
    
    ins=fopen("sim0.mat","w");
    for(i=1; i<=*iDataP; i++){
        for(j=1; j<=*iDataP; j++) fprintf(ins, " %d", E[i][j]);
        fprintf(ins, "\n");
    }
    fclose(ins);
    
    
loop: free_dvector(datax,1,DATA_NUM_PLUS);
    free_dvector(score,1,DATA_NUM_PLUS);
    free_ivector(row,1,DATA_NUM_PLUS);
    free_ivector(col,1,DATA_NUM_PLUS);
    free_ivector(no,1,DATA_NUM);
    free_dmatrix(iDataX,1,*iDataNum,1,*iDataP);
    free_dmatrix(COV,1,*iDataP,1,*iDataP);
    free_dmatrix(CORR,1,*iDataP,1,*iDataP);
    free_ivector(indx,1,DATA_NUM);
    free_ivector(subrow,1,DATA_NUM);
    free_ivector(subcol,1,DATA_NUM);
    free_dvector(subscore,1,DATA_NUM);
    free_dmatrix(Qmat,1,DATA_NUM,1,4);
    free_dvector(psi,1,DATA_NUM);
    free_imatrix(E,1,*iDataP,1,*iDataP);

    //return 0;
}


void allR(int DataX[],
          int *DataNum,
          int *GeneNum,
          int *total_iteration,
          double *stepsize)

{
    
#define pi 3.14159265
    FILE *ins;
    int i, j, **x, accept, iter, *count;
    double newalpha, *alpha, *logalpha, newlogalpha, *beta, **theta, **ave,  *accept_loc, *total_loc, un, r, sum, newsum;
    
    
    
    /* parameters  */
    
    static int warm=500;
    static double a1=0.001, b1=1E8, a2=0.001, b2=1E8;
    
    

    
    /*
     ins=fopen("aa.log","a");
     fprintf(ins, "seed=%d\n", stime);
     fclose(ins);
     */
    x=imatrix(1, *GeneNum, 1, *DataNum);
    count = ivector(1, *GeneNum);
    theta=dmatrix(1, *GeneNum, 1, *DataNum);
    ave=dmatrix(1, *GeneNum, 1, *DataNum);
    alpha=dvector(1, *GeneNum);
    logalpha=dvector(1, *GeneNum);
    beta=dvector(1, *GeneNum);
    accept_loc=dvector(1, *GeneNum);
    total_loc=dvector(1, *GeneNum);
    
    
    
    for (i = 1; i <= *GeneNum; i++) {
        for (j = 1; j <= *DataNum; j++) {
            x[i][j]=DataX[(i-1)*(*DataNum)+(j-1)];
        }
    }
    
    
    
    for (i = 1; i <= *GeneNum; i++) {
        logalpha[i] = 0.0,
        alpha[i] = exp(logalpha[i]);
        beta[i] = 1.0;
        accept_loc[i] = 0.1;
        total_loc[i] = 0.1;
        count[i] = 0;
    }
    
    for (i = 1; i <= *GeneNum; i++) {
        for (j=1; j <= *DataNum; j++) {
            theta[i][j] = (x[i][j] + alpha[i]) / (1 + beta[i]);
            ave[i][j] = 0.0;
        }
    }
    
    
    
    
    
    
    
    for (iter=1; iter<=warm+*total_iteration; iter++) {
        for (i = 1; i <= *GeneNum; i++) {
            // update alpha
            newlogalpha=logalpha[i]+gasdev()*(*stepsize);
            alpha[i]=exp(logalpha[i]);
            //printf("alpha=%f\n",alpha[i]);
            newalpha=exp(newlogalpha);
            
            for(newsum=0.0, j=1; j<=*DataNum; j++) newsum+=(x[i][j]+newalpha-1)*log(theta[i][j]);
            newsum+=*DataNum*(newalpha*log(beta[i])+gammln(newalpha));
            newsum+=a1*newlogalpha-b1*newalpha;
            
            for(sum=0.0, j=1; j<=*DataNum; j++) sum+=(x[i][j]+alpha[i]-1)*log(theta[i][j]);
            sum+=*DataNum*(alpha[i]*log(beta[i])+gammln(alpha[i]));
            sum+=a1*logalpha[i]-b1*alpha[i];
            
            r=newsum-sum;
            
            if(r>0.0) accept=1;
            else{
                un=unif_rand();
                if(un<exp(r)) accept=1;
                else accept=0;
            }
            
            if(accept==1){ logalpha[i]=newlogalpha; alpha[i]=newalpha; accept_loc[i]+=1.0; }
            total_loc[i]+=1.0;
            
            
            // update beta
            for(sum=0.0, j=1; j<=*DataNum; j++) sum+=theta[i][j];
            beta[i]=Rgamma(alpha[i]**DataNum+a2, sum+b2);
            // printf("beta=%f\n",beta[i]);
            // update theta
            for(j=1; j<=*DataNum; j++) theta[i][j]=Rgamma(x[i][j]+alpha[i],1.0+beta[i]);
            
            
            if(iter>warm && iter%10==0){
                count[i] ++;
                for(j=1; j<=*DataNum; j++) ave[i][j]+=theta[i][j];
            }
            
        }
        b1=b1+1/iter;
        b2=b1;
        
        
    }
    
    
    ins=fopen("continuzed.dat","w");
    
    for(i = 1; i <= *DataNum; i++) {
        for(j = 1; j <= *GeneNum; j++) {
            fprintf(ins, " %g", ave[j][i]/count[j]);
        }
        fprintf(ins, "\n");
    }
    fclose(ins);
    
    
    
 /*   for (i = 1; i <= *GeneNum; i++) {
        Rprintf("acceptance rate=%g\n", 1.0*accept_loc[i]/total_loc[i]);
    }
*/
    
    free_imatrix(x,1,*GeneNum, 1, *DataNum);
    free_dmatrix(theta,1,*GeneNum, 1, *DataNum);
    free_dmatrix(ave,1,*GeneNum, 1, *DataNum);
    free_dvector(alpha, 1, *GeneNum);
    free_dvector(logalpha, 1, *GeneNum);
    free_dvector(beta, 1, *GeneNum);
    free_dvector(accept_loc, 1, *GeneNum);
    free_dvector(total_loc, 1, *GeneNum);
    free_ivector(count, 1, *GeneNum);
    
}



void pcorsel1(int irow[], int icol[], double idatax[], int *DATA_NUM )
{
#define pi 3.14159265
    static int  MAX_NUM=200000, GRID=2, mingrid=2;
    static int  WARM=1, stepscale=10000, samsize=10000;
    static double rho=0.01;
    
    //#define len(array,len){len=(sizeof(array)/sizeof(array[0]));}
    
    
    
    FILE *ins;
    int  grid, bestgrid, rep;
    int *indx, *indx2,  *sam, *row, *col, i, j, k, k1, k2, m, com, iter, ok;
    double **resam,*rez,*dis,*dif,*mu,*alpha,*beta,*prob,*logprob,*wei,*datax,*FDR,*Q;
    double entropy;
    double delta,z,un,max,min,sum,tep,tep1,tep2,eta,dmu,dalpha,dbeta;
    double bestentropy, *bestprob, *bestmu, *bestalpha, *bestbeta;
    

    

    
    
    datax=dvector(1,MAX_NUM);
    row=ivector(1,MAX_NUM);
    col=ivector(1,MAX_NUM);
    mu=dvector(1,GRID);
    alpha=dvector(1,GRID);
    beta=dvector(1,GRID);
    prob=dvector(1,GRID);
    wei=dvector(1,GRID);
    sam=ivector(1,MAX_NUM);
    indx=ivector(1,MAX_NUM);
    indx2=ivector(1,GRID);
    FDR=dvector(1,MAX_NUM);
    Q=dvector(1,MAX_NUM);
    bestmu=dvector(1,GRID);
    bestalpha=dvector(1,GRID);
    bestbeta=dvector(1,GRID);
    bestprob=dvector(1,GRID);
    logprob=dvector(1,GRID);
    resam=dmatrix(1,GRID,1,samsize);
    rez=dvector(1,samsize);
    dis=dvector(1,GRID);
    dif=dvector(1,GRID);
    
    
    
    if(*DATA_NUM > MAX_NUM){ Rprintf("The dataset is too big\n"); }
    
    
    Rprintf("working subdata size=%d\n", *DATA_NUM);
    
    // col is col2, row is col1, datax is col3.
    
    for( i=0; i<*DATA_NUM; i++)
    {
        row[i+1]=irow[i];
        col[i+1]=icol[i];
        datax[i+1]=idatax[i];
        
    }
    
    indexx(*DATA_NUM,datax,indx);
    
    grid=GRID+1;  bestentropy=1.0e+100;
    while(grid>mingrid){
        
        grid--;
        
        /* Initialization */
        mu[1]=0.0; alpha[1]=0.0; k=0;
        for(i=1; i<=(int)(*DATA_NUM*0.99); i++){
            mu[1]+=datax[indx[i]];
            alpha[1]+=datax[indx[i]]*datax[indx[i]];
            k++;
        }
        mu[1]=mu[1]/k;
        alpha[1]=sqrt(alpha[1]/k-mu[1]*mu[1])*sqrt(2.0);
        beta[1]=2.0; prob[1]=0.99;
        
        for(j=2; j<=grid; j++){
            tep=(1-prob[1])/(grid-1);
            mu[j]=0.0; alpha[j]=0.0; k=0;
            for(i=(int)(*DATA_NUM*(prob[1]+(j-2)*tep)); i<=(int)(*DATA_NUM*(prob[1]+(j-1)*tep)); i++){
                mu[j]+=datax[indx[i]];
                alpha[j]+=datax[indx[i]]*datax[indx[i]];
                k++;
            }
            mu[j]=mu[j]/k;
            alpha[j]=sqrt(alpha[j]/k-mu[j]*mu[j])*sqrt(2.0);
            beta[j]=2.0; prob[j]=tep;
        }
        for(j=1; j<=grid; j++){
            logprob[j]=log(prob[j])-log(prob[grid]);
        }
        
        
        for(iter=1; iter<=100; iter++){
            
            random_order(sam,*DATA_NUM);
            for(i=1; i<=*DATA_NUM; i++){
                
                rep=*DATA_NUM*(iter-1)+i;
                if(rep<=WARM*stepscale) delta=rho;
                else delta=rho*exp(-log(1.0*(rep-(WARM-1)*stepscale)/stepscale));
                
                for(sum=0.0,j=1; j<=grid; j++) sum+=exp(logprob[j]);
                for(j=1; j<=grid; j++) prob[j]=exp(logprob[j])/sum;
                
                z=datax[sam[i]];
                max=-1.0e+100;
                for(j=1; j<=grid; j++){
                    wei[j]=dlogGnormal(z,mu[j],alpha[j],beta[j])+log(prob[j]);
                    if(wei[j]>max){ max=wei[j];}
                }
                for(sum=0.0, j=1; j<=grid; j++) sum+=exp(wei[j]-max);
                
                indexx(grid,mu,indx2);
                for(k=1; k<=grid; k++){
                    j=indx2[k];
                    eta=exp(wei[j]-max)/sum;
                    
                    if(z>mu[j]) dmu=beta[j]/alpha[j]*exp((beta[j]-1.0)*log(fabs(z-mu[j])/alpha[j]));
                    else dmu=-beta[j]/alpha[j]*exp((beta[j]-1.0)*log(fabs(z-mu[j])/alpha[j]));
                    
                    dalpha=-1.0/alpha[j]+beta[j]/alpha[j]*exp(beta[j]*log(fabs(z-mu[j])/alpha[j]));
                    dbeta=1.0/beta[j]+1.0/beta[j]/beta[j]*derivative_gamma(1.0/beta[j])
                    -exp(beta[j]*log(fabs(z-mu[j])/alpha[j]))*log(fabs(z-mu[j])/alpha[j]);
                    
                    mu[j]+=delta*eta*dmu;
                    alpha[j]+=delta*eta*dalpha;
                    beta[j]+=delta*eta*dbeta;
                    logprob[j]+=delta*(eta-prob[j]);
                    
                    if(mu[j]<-20.0) mu[j]=-20.0;
                    if(mu[j]>40.0) mu[j]=40.0;
                    if(alpha[j]<0.1) alpha[j]=0.1;
                    if(alpha[j]>10.0) alpha[j]=10.0;
                    if(beta[j]<0.1) beta[j]=0.1;
                    if(beta[j]>10.0) beta[j]=10.0;
                }
            } /* end one epoach */
            
            if(iter%1==0) Rprintf("iter=%d delta=%g %g %g %g %g: %g %g %g %g\n",
                                 iter,delta,prob[1],mu[1],alpha[1],beta[1],prob[2],mu[2],alpha[2],beta[2]);
            
        }      /* end for iterations */
        
        
        for(sum=0.0,j=1; j<=grid; j++) sum+=exp(logprob[j]);
        for(j=1; j<=grid; j++) prob[j]=exp(logprob[j])/sum;
        
        for(entropy=0.0,i=1; i<=*DATA_NUM; i++){
            max=-1.0e+100;
            for(j=1; j<=grid; j++){
                wei[j]=dlogGnormal(z,mu[j],alpha[j],beta[j])+log(prob[j]);
                if(wei[j]>max){ max=wei[j]; }
            }
            for(sum=0.0, j=1; j<=grid; j++) sum+=exp(wei[j]-max);
            entropy+=log(sum)+max;
        }
        entropy*=-1.0/(*DATA_NUM);
        entropy+=0.5*(grid*3+grid-1)*log(1.0**DATA_NUM)/(*DATA_NUM);
        
        
        
        if(entropy<=bestentropy){
            bestentropy=entropy;
            bestgrid=grid;
            for(j=1; j<=grid; j++){ bestprob[j]=prob[j]; bestmu[j]=mu[j];
                bestalpha[j]=alpha[j]; bestbeta[j]=beta[j]; }
        }
    }  /* end for different grid size */
    
    /* calculating the distance between the mixture components */
    for(k=1; k<=bestgrid; k++){
        MCMCGnormal(bestmu[k],bestalpha[k],bestbeta[k],samsize,rez);
        for(j=1; j<=samsize; j++) resam[k][j]=rez[j];
    }
    indexx(bestgrid, bestmu, indx2);
    for(j=1; j<bestgrid; j++){
        k1=indx2[j]; k2=indx2[j+1];
        tep1=tep2=0.0;
        for(i=1; i<=samsize; i++){
            tep1+=dlogGnormal(resam[k1][i],bestmu[k1],bestalpha[k1],bestbeta[k1])-
            dlogGnormal(resam[k1][i],bestmu[k2],bestalpha[k2],bestbeta[k2]);
            tep2+=dlogGnormal(resam[k2][i],bestmu[k2],bestalpha[k2],bestbeta[k2])-
            dlogGnormal(resam[k2][i],bestmu[k1],bestalpha[k1],bestbeta[k1]);
        }
        if(tep1<tep2) dis[j]=tep1/samsize;
        else dis[j]=tep2/samsize;
    }
    
    
    /* determine the number of components for f_0 */
    dif[1]=1.0;
    for(j=2; j<bestgrid; j++) dif[j]=(dis[j]-dis[j-1])/dis[j-1];
    dif[bestgrid]=-1.0;
    com=1; ok=0;
    while(ok==0 && com<bestgrid-1){
        if(dif[com]*dif[com+1]<0.0) ok=1;
        else com++;
    }
    /*
     max=0.0; k=1;
     for(j=2; j<bestgrid; j++)
     if(dif[j]>max){ max=dif[j]; k=j; }
     if(k<com) com=k;
     */
    
   
    
    // com=bestgrid-1;
    com=1;
    
    
    /* calculating FDR */
    
    for(i=1; i<=*DATA_NUM; i++){
        z=datax[indx[i]];
        for(sum=0.0,m=1; m<=com; m++){
            k=indx2[m];
            if(z>bestmu[k]){
                un=exp(bestbeta[k]*log((z-bestmu[k])/bestalpha[k]));
                tep=(0.5-0.5*gammp(1.0/bestbeta[k],un))*bestprob[k];
            } else if(z<bestmu[k]){
                un=exp(bestbeta[k]*log((bestmu[k]-z)/bestalpha[k]));
                tep=(0.5+0.5*gammp(1.0/bestbeta[k],un))*bestprob[k];
            } else tep=0.5*bestprob[k];
            sum+=tep;
        }
        FDR[i]=sum;
        
        for(sum=0.0, k=1; k<=bestgrid; k++){
            if(z>bestmu[k]){
                un=exp(bestbeta[k]*log((z-bestmu[k])/bestalpha[k]));
                tep=(0.5-0.5*gammp(1.0/bestbeta[k],un))*bestprob[k];
            } else if(z<bestmu[k]){
                un=exp(bestbeta[k]*log((bestmu[k]-z)/bestalpha[k]));
                tep=(0.5+0.5*gammp(1.0/bestbeta[k],un))*bestprob[k];
            } else tep=0.5*bestprob[k];
            sum+=tep;
        }
        FDR[i]/=sum;
    }
    
    min=FDR[1]; Q[1]=min;
    for(i=2; i<=*DATA_NUM; i++){
        if(FDR[i]<min) min=FDR[i];
        Q[i]=min;
    }
    

    
    
    /* alternative method of FDR calculation */
    for(i=1; i<=*DATA_NUM; i++){
        z=datax[indx[i]];
        for(sum=0.0,m=1; m<=com; m++){
            k=indx2[m];
            if(z>bestmu[k]){
                un=exp(bestbeta[k]*log((z-bestmu[k])/bestalpha[k]));
                tep=(0.5-0.5*gammp(1.0/bestbeta[k],un))*bestprob[k];
            } else if(z<bestmu[k]){
                un=exp(bestbeta[k]*log((bestmu[k]-z)/bestalpha[k]));
                tep=(0.5+0.5*gammp(1.0/bestbeta[k],un))*bestprob[k];
            } else tep=0.5*bestprob[k];
            sum+=tep;
        }
        FDR[i]=sum**DATA_NUM/(*DATA_NUM-i+1);
        if(FDR[i]>1.0) FDR[i]=1.0;
    }
    
    min=FDR[1]; Q[1]=min;
    for(i=2; i<=*DATA_NUM; i++){
        if(FDR[i]<min) min=FDR[i];
        Q[i]=min;
    }
    
    ins=fopen("aaa.fdr", "w");
    for(i=1; i<=*DATA_NUM; i++)
        fprintf(ins, " %5d %5d %g %g %g\n", row[indx[i]],col[indx[i]],datax[indx[i]],FDR[i],Q[i]);
    fprintf(ins, "\n");
    fclose(ins);


    
    
    free_dvector(datax,1,MAX_NUM);
    free_ivector(row,1,MAX_NUM);
    free_ivector(col,1,MAX_NUM);
    free_dvector(mu,1,GRID);
    free_dvector(alpha,1,GRID);
    free_dvector(beta,1,GRID);
    free_dvector(prob,1,GRID);
    free_dvector(wei,1,GRID);
    free_ivector(sam,1,MAX_NUM);
    free_ivector(indx,1,MAX_NUM);
    free_ivector(indx2,1,GRID);
    free_dvector(FDR,1,MAX_NUM);
    free_dvector(Q,1,MAX_NUM);
    free_dvector(bestmu,1,GRID);
    free_dvector(bestalpha,1,GRID);
    free_dvector(bestbeta,1,GRID);
    free_dvector(bestprob,1,GRID);
    free_dvector(logprob,1,GRID);
    free_dmatrix(resam,1,GRID,1,samsize);
    free_dvector(rez,1,samsize);
    free_dvector(dis,1,GRID);
    free_dvector(dif,1,GRID);
}







