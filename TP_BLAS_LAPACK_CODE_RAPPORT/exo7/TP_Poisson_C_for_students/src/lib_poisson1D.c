/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_rowMajor_poisson1D(double* AB, int *lab,int *la){

  //TODO
  int ii_1,ii_2,ii_3,ii_4,ii_5;
  int jj;
  
  for(jj=0;jj<(*la);jj++){
  ii_1=jj;
  ii_2=(*la)+jj;
  ii_3=2*(*la)+jj;
  ii_4=3*(*la)+jj;
  AB[ii_1]=0.0;
  AB[ii_2]=-1.0;
  AB[ii_3]=2.0;
  AB[ii_4]=-1.0;
  
  }
  AB[*la]=0.0;
  AB[4*(*la)-1]=0.0;
}
void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    if (*kv>=0){
      for (ii=0;ii< *kv;ii++){
	AB[kk+ii]=0.0;
      }
    }
    AB[kk+ *kv]=-1.0;
    AB[kk+ *kv+1]=2.0;
    AB[kk+ *kv+2]=-1.0;
  }
  AB[0]=0.0;
  if (*kv == 1) {AB[1]=0;}
  
  AB[(*lab)*(*la)-1]=0.0;
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    if (*kv>=0){
      for (ii=0;ii< *kv;ii++){
	AB[kk+ii]=0.0;
      }
    }
    AB[kk+ *kv]=0.0;
    AB[kk+ *kv+1]=1.0;
    AB[kk+ *kv+2]=0.0;
  }
  AB[1]=0.0;
  AB[(*lab)*(*la)-1]=0.0;
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  int jj;
  RHS[0]= *BC0;
  RHS[(*la)-1]= *BC1;
  for (jj=1;jj<(*la)-1;jj++){
    RHS[jj]=0.0;
  }
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  int jj;
  double h, DELTA_T;
  DELTA_T=(*BC1)-(*BC0);
  for (jj=0;jj<(*la);jj++){
    EX_SOL[jj] = (*BC0) + X[jj]*DELTA_T;
  }
}  

void set_grid_points_1D(double* x, int* la){
  int jj;
  double h;
  h=1.0/(1.0*((*la)+1));
  for (jj=0;jj<(*la);jj++){
    x[jj]=(jj+1)*h;
  }
}

void write_GB_operator_rowMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*lab);ii++){
      for (jj=0;jj<(*la);jj++){
	fprintf(file,"%lf\t",AB[ii*(*la)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  //TODO
  int ii;
  int jj;
  FILE *file;
  // On parocours notre 
  file=fopen(filename,"w");
  if(file!=NULL){
    // On parocours jsqu'a notre dimension n dans notre cas c'est la 
    for(ii=0;ii<(*la);ii++){
    // on parcours jusqu'a la leading dimension
     for(jj=0;jj<(*lab);jj++){
    	 fprintf(file,"%lf\t",AB[ii*(*lab)+jj]);
     }
     fprintf(file,"\n");
   }
   fclose(file);
  }
  else{
  
  perror(filename);
  
  }
 
}

void write_vec(double* vec, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\n",vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  } 
}  

void write_xy(double* vec, double* x, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\t%lf\n",x[jj],vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  } 
}  

void eig_poisson1D(double* eigval, int *la){
  int ii;
  double scal;
  for (ii=0; ii< *la; ii++){
    scal=(1.0*ii+1.0)*M_PI_2*(1.0/(*la+1));
    eigval[ii]=sin(scal);
    eigval[ii]=4*eigval[ii]*eigval[ii];
  } 
}

double eigmax_poisson1D(int *la){
  double eigmax;
  eigmax=sin(*la *M_PI_2*(1.0/(*la+1)));
  eigmax=4*eigmax*eigmax;
  return eigmax;
}

double eigmin_poisson1D(int *la){
  double eigmin;
  eigmin=sin(M_PI_2*(1.0/(*la+1)));
  eigmin=4*eigmin*eigmin;
  return eigmin;
}

double richardson_alpha_opt(int *la){
  //TODO
  double opt;
  opt=eigmin_poisson1D(la);
  opt+=eigmax_poisson1D(la);
  opt=2.0/alpha_opt;
  
  return opt;
  
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit){
  //TODO
  double alpha,beta,res,normb,normes;
  double *Y;
  int incx,incy;
  int it;
  
  it=0;
  /*set incx and incy for bls orpeations*/
  
  incx=1;
  incy=1;
  
  Y=(double*)calloc(*la,sizeof(double));
  /*compute first residuel b-AX */
  cblas_dcopy(*la,RHS,incx,Y,incy);
  normb=cblas_ddot(*la,Y,1,Y,1);
  normb=sqrt(normb);
  
  /*set alpha and beta for dgmv*/
  alpha=-1.0;
  beta=1.0;
  
  //dgmv('N',la,la,kl,ku,alpha,AB,lab,x,incx,beta,RHS,incy);
  cblas_dgbmv(CblasColMajor,CblasNoTrans,*la,*la,*kl,*ku,alpha,
  AB,*lab,X,incx,beta,Y,incy);
  normes=cblas_ddot(*la,Y,1,Y,1);
  normes=sqrt(normes);
  res=normes/normb;
  printf("\n res0=%lf,normabs=%lf\n",res,normes);
  
  /*Richardson iterative process*/
  
  while ((res>(*tol))&(it<(*maxit-1)))
  {
    /*compute iterate*/
    cblas_daxpy(*la,*alpha_rich,Y,incx,X,incy);
    /*new residual*/
    cblas_dcopy(*la,RHS,incx,Y,incy);
    cblas_dgbmv(CblasColMajor,CblasNoTrans,*la,*la,*kl
  ,*ku,alpha,AB,*lab,X,incx,beta,Y,incy);
  
    /*compute relative residual*/
    normes=cblas_ddot(*la,Y,1,Y,1);
    normes=sqrt(normes);
    res=normes/normb;
    printf("\n resvec[%d]=%e \n",it,res);
    it++;
    
   }
   
   free(Y);
}

