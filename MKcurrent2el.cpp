//---------------------------------------------------------------------------

#include "fileutils.h"
#include "stringutils.h"
#include "CGanalysis.h"
#include "systemclass.h"
#include "treeutils.h"
#include "MKcurrent2el.h"

MKcurrent2el::MKcurrent2el(int steps)
  {
  from = new int[steps];
  to = new int[steps];
  energy = new double[steps];
  MCsteps = new int[steps];
  dxs = new int[steps];
  dys = new int[steps];
  ts = new double[steps];
  des = new double[steps];
  }

MKcurrent2el::~MKcurrent2el()
{

  if (from!=NULL) delete[] from;
  if (to!=NULL) delete[] to;
  if (jl!=NULL) delete[] jl;
  if (GammaT!=NULL) delete[] GammaT;
  if (energy!=NULL) delete[] energy;
  if (MCsteps!=NULL) delete[] MCsteps;
  if (dxs!=NULL) delete[] dxs;
  if (dys!=NULL) delete[] dys;
  if (ts!=NULL) delete[] ts;
  if (des!=NULL) delete[] des;
}


void MKcurrent2el::setE(double Ex1, double Ey1, double Ez1)
{
  Ex=Ex1; Ey=Ey1; Ez=Ez1;
}

void MKcurrent2el::setES(ESystem *e)
{
  es=e;
}

void MKcurrent2el::init(int D1, int L1, int N1, int maxJL1, double a1, double beta1,double maxProbability1, int maxJL2el1,int d2el1)
{
  int dx, dy, dx1,dy1,dx2,dy2,dxt,dyt, n,Nmem2el;
  double gamma,r12,r13,r14,r23,r24,r34,e1,e2,v1,v2;

  D=D1; L=L1; N=N1; maxJL=maxJL1; JL2=2*maxJL+1; A=a1; beta = beta1;
  maxProbability=maxProbability1/beta,maxJL2el=maxJL2el1,d2el=d2el1;

  printf("L: %d N: %d maxJL: %d  maxJL2el:  %d  d2el:  %d\n", L, N, maxJL,maxJL2el,d2el);
  if (D==2)
    {
      Nmem=JL2*JL2;
      jl = new double[Nmem];
      GammaT = new double[Nmem];
      for (dx=-maxJL;dx<maxJL+1;dx++)
	for (dy=-maxJL;dy<maxJL+1;dy++)
	  {
	    jl[(dx+maxJL)+JL2*(dy+maxJL)]=sqrt(dx*dx+dy*dy);
	  }
      jl[maxJL+JL2*maxJL]=0.0;
      gamma=0.0;
      for (n=0;n<Nmem/2;n++)
	{
	  gamma+=exp(-A*jl[n]); GammaT[n]=gamma;
	  //printf("%le ",gamma); if ((n+1)%JL2==0) printf("\n");
	}
      GammaT[Nmem/2]=gamma; //printf("0 ");
      for (n=Nmem/2+1;n<Nmem;n++)
	{
	  gamma+=exp(-A*jl[n]); GammaT[n]=gamma;
	  //printf("%le ",gamma); if ((n+1)%JL2==0) printf("\n");
	}

      GammaTtot=gamma;
      printf("\n GammaTtot =  %le \n",GammaTtot);

      JL22el = 2*maxJL2el+1;
      Nmem2el = JL22el*JL22el*JL22el*JL22el*(2*d2el+1)*(2*d2el+1);
      printf("Nmem2el= %d\n",Nmem2el);
      GammaT2el =  new double[Nmem2el];
      dx2el =  new int[Nmem2el];
      dy2el =  new int[Nmem2el];
      dx12el =  new int[Nmem2el];
      dy12el =  new int[Nmem2el];
      dx22el =  new int[Nmem2el];
      dy22el =  new int[Nmem2el];

      // NB: antar her at Coulomb ikke maa regnes rundt systemet (dvs at lengden paa hopp + avstanden mellom par < halve L)
      gamma = 0.0;
      n = 0;
      for (dx=-maxJL2el;dx<maxJL2el+1;dx++)
	for (dy=-maxJL2el;dy<maxJL2el+1;dy++)
	  if( dx != 0 || dy !=0 )
	    for (dx1=-d2el;dx1<d2el+1;dx1++)
	      for (dy1=-d2el;dy1<d2el+1;dy1++)
		if( (dx1 != 0 || dy1 !=0) && ( dx1 != dx || dy1 != dy))
		  for (dx2=-maxJL2el;dx2<maxJL2el+1;dx2++)
		    for (dy2=-maxJL2el;dy2<maxJL2el+1;dy2++)
		      if( (dx2 != 0 || dy2 !=0) && ((dx1+dx2) != 0 || (dy1+dy2) !=0) && ((dx1+dx2) != dx || (dy1+dy2) !=dy)) 
			{
			  r13 = jl[(dx+maxJL)+JL2*(dy+maxJL)];
			  r24 = jl[(dx2+maxJL)+JL2*(dy2+maxJL)];
			  r12 = sqrt(dx1*dx1+dy1*dy1);
			  dxt = dx1+dx2; dyt = dy1+dy2;
			  r14 = sqrt(dxt*dxt+dyt*dyt);
			  dxt = dx1-dx; dyt = dy1-dy;
			  r23 = sqrt(dxt*dxt+dyt*dyt);
			  dxt = dx1+dx2-dx; dyt = dy1+dy2-dy;
			  r34 = sqrt(dxt*dxt+dyt*dyt);
			  
			  e1 = exp(-A*(r13+r24));
			  e2 = exp(-A*(r23+r14));
			  v1 = 1/r12-1/r14-1/r23+1/r34;
			  v2 = 1/r12-1/r13-1/r24+1/r34;
		      
			  gamma += e1*e1*v1*v1 + e2*e2*v2*v2 + e1*e2*v1*v2;
			  GammaT2el[n] = gamma;
			  dx2el[n] = dx;
			  dy2el[n] = dy;
			  dx12el[n] = dx1;
			  dy12el[n] = dy1;
			  dx22el[n] = dx2;
			  dy22el[n] = dy2;
			  n++;
			  // printf("dx=%d dy=%d dx1=%d dy1=%d dx2=%d dy2=%d gamma= %le \n",dx,dy,dx1,dy1,dx2,dy2,gamma); 
			  //printf("%d ",n);
			  //printf("%le \n",e1*e1*v1*v1 + e2*e2*v2*v2 + e1*e2*v1*v2);
			}
      nGamma2el = n;
      GammaTtot2el = gamma;
      totGammaTtot = GammaTtot + GammaTtot2el;
      printf("Nmem= %d nGamma2el= %d GammaTtot= %le GammaTtot2el= %le\n",Nmem,nGamma2el,GammaTtot,GammaTtot2el);
    }
}

/*
double inline MKcurrent2el::calcRate2d(int n1, int n2, int dx, int dy)
{
  double dE, dU, r;

  r=jl[(dx+maxJL)+JL2*(dy+maxJL)];
  dU=dx*Ex+dy*Ey;
  dE=es->hoppEdiff(n1,-1,0,n2,0,0)+dU;
  if (dE>0)
    return exp(-A*r-beta*dE);
  else
    return exp(-A*r);
}
*/

void inline MKcurrent2el::getjump(int &i, int &j, int &dx, int &dy)
{
  double r2;
  int h,l,step, x, y, x2, y2;

  i=(int)(N*(double)es->ran2(0));

  while (es->getocci(i) == 0){
    i=(int)(N*(double)es->ran2(0));}

  r2=es->ran2(0)*GammaTtot;
  l=-1;h=Nmem-1; step=Nmem/2;

  while (step>0)
    {
      if (GammaT[l+step]>=r2)
	{ h=l+step;}
      else
	{ l=l+step;}
      step=(h-l)/2;
    }
  //at end h=correct jump;

  x=i%L;
  y=i/L;

  dx=h%JL2-maxJL;
  dy=h/JL2-maxJL;

  x2=x+dx; if (x2>=L) x2-=L; else if (x2<0) x2+=L;
  y2=y+dy; if (y2>=L) y2-=L; else if (y2<0) y2+=L;

  j=x2+y2*L;
  //printf("gj: %d %d %d %d h: %d r2: %le\n",i,j,dx,dy,h,r2);
}

void inline MKcurrent2el::getjump2el(int &i, int &j,int &k, int &l, int &dx, int &dy,int &dx2, int &dy2  )
{
  double r2;
  int h,ll,step, x, y, xj, yj,xk,yk,xl,yl;

  i=(int)(N*(double)es->ran2(0));

  while (es->getocci(i) == 0){
    i=(int)(N*(double)es->ran2(0));}

  r2=es->ran2(0)*GammaTtot2el;
  ll=-1;h=nGamma2el-1; step=nGamma2el/2;

  while (step>0)
    {
      if (GammaT2el[ll+step]>=r2)
	{ h=ll+step;}
      else
	{ ll=ll+step;}
      step=(h-ll)/2;
    }
  //at end h=correct jump;

  x=i%L;
  y=i/L;

  dx=dx2el[h];
  dy=dy2el[h];

  xj=x+dx; if (xj>=L) xj-=L; else if (xj<0) xj+=L;
  yj=y+dy; if (yj>=L) yj-=L; else if (yj<0) yj+=L;

  j=xj+yj*L;

  xk=x+dx12el[h]; if (xk>=L) xk-=L; else if (xk<0) xk+=L;
  yk=y+dy12el[h]; if (yk>=L) yk-=L; else if (yk<0) yk+=L;

  k=xk+yk*L;

  dx2 = dx22el[h];
  dy2 = dy22el[h];

  xl=xk+dx2; if (xl>=L) xl-=L; else if (xl<0) xl+=L;
  yl=yk+dy2; if (yl>=L) yl-=L; else if (yl<0) yl+=L;

  l=xl+yl*L;

 
  //printf("gj: %d %d %d %d h: %d r2: %le\n",i,j,dx,dy,h,r2);
}

bool inline MKcurrent2el::testjump(int i, int j, int dx, int dy)
{
  double dE1, dE2, rate, r2, exponent;

   if (es->getocci(j)!=0)
    {
      //printf("ij: %d %d %d %d  n[j]!=0\n",i,j,currnptr[i],currnptr[j]);
      return false;
    }

  dE1=es->hoppEdiffij(i,j);
  dE2=dx*Ex+dy*Ey+dE1;
  exponent=beta*dE2;
 
  if (dE2<0) rate = 1; 
  else rate= exp(-exponent);

  // if (dE2<0) rate = -dE2*(1+1/(exp(-exponent)-1)); 
  //else if (dE2==0) rate = 1/beta;
  //else rate= dE2/(exp(exponent)-1);
  //  if (rate > 3/beta) printf("%f\n",rate); 


  r2=es->ran2(0);
  //printf("r2: %le rate: %le \n",r2,rate);
  //if (r2<rate) printf("\nij: %d %d r2: %le rate: %le \n",i,j,r2,rate);
  // if (r2<rate) printf("dE1 = %le ",dE1);
  return (r2<rate);
}

bool inline MKcurrent2el::testjump2el(int i, int j, int k, int l,  int dx, int dy, int dx2, int dy2)
{
  double dE1, dE2, rate, r2, exponent;

  if (es->getocci(j)!=0 || es->getocci(k)==0 || es->getocci(l)!=0)
    {
      //printf("ij: %d %d %d %d  n[j]!=0\n",i,j,currnptr[i],currnptr[j]);
      return false;
    }
  
  dE1=es->hoppEdiffijkl(i,j,k,l,dx,dy,dx2,dy2);
  dE2=(dx+dx2)*Ex+(dy+dy2)*Ey+dE1;
  exponent=beta*dE2;
  if (dE2<0) rate = 1; 
  else rate= exp(-exponent);
  
  r2=es->ran2(0);
  if (r2<rate) printf("dE1 = %le ",dE1);
  return (r2<rate);
}


void MKcurrent2el::calcAllRates(string filename)
{
  int i,h,j,k,l,c1,c2,x,y,dx,dy,x2,y2,xj,yj,xk,yk,xl,yl,dx2,dy2,n,n1,limit;
  double dE1,dE2,rate,exponent,*rates,*rates2el,tot,*en,*en2el;
  int *indexes;
  string st;
  FILE* f;
  printf("Nmem= %i N = %i nGamma2el = %i\n",Nmem,N,nGamma2el);
  rates = new double[Nmem*N];
  en = new double[Nmem*N];

  printf("Har allokert 1el\n");
  rates2el = new double[nGamma2el*N/2];
  en2el = new double[nGamma2el*N/2];
  printf("Har allokert 2el\n");
  
  c1 = 0; c2 = 0;
  for(i=0;i<N;i++)
    if (es->getocci(i)==1){
      for(h=0;h<Nmem;h++){
	x=i%L;
	y=i/L;
	
	dx=h%JL2-maxJL;
	dy=h/JL2-maxJL;
	
	x2=x+dx; if (x2>=L) x2-=L; else if (x2<0) x2+=L;
	y2=y+dy; if (y2>=L) y2-=L; else if (y2<0) y2+=L;
	
	j=x2+y2*L;
	
	if (es->getocci(j)==0){
	  dE1=es->hoppEdiffij(i,j);
	  dE2=dx*Ex+dy*Ey+dE1;
	  exponent=beta*dE2;
	  
	  if (dE2<0) rate = 1; 
	  else rate= exp(-exponent);
	  if(h==0)
	    rate *= GammaT[0];
	  else
	    rate *= GammaT[h]-GammaT[h-1];
	  if(rate>0){
	    rates[c1] = rate;
	    en[c1] = dE2;
	    c1++;
	  }
	  // printf("1el: i = %i h = %i  rate = %e \n",i,h,rate);
	}	
      }
      for(h=0;h<nGamma2el;h++){
	x=i%L;
	y=i/L;
	
	dx=dx2el[h];
	dy=dy2el[h];
	
	xj=x+dx; if (xj>=L) xj-=L; else if (xj<0) xj+=L;
	yj=y+dy; if (yj>=L) yj-=L; else if (yj<0) yj+=L;
	
	j=xj+yj*L;
	
	xk=x+dx12el[h]; if (xk>=L) xk-=L; else if (xk<0) xk+=L;
	yk=y+dy12el[h]; if (yk>=L) yk-=L; else if (yk<0) yk+=L;
	
	k=xk+yk*L;
	
	dx2 = dx22el[h];
	dy2 = dy22el[h];
	
	xl=xk+dx2; if (xl>=L) xl-=L; else if (xl<0) xl+=L;
	yl=yk+dy2; if (yl>=L) yl-=L; else if (yl<0) yl+=L;
	
	l=xl+yl*L;
	if (es->getocci(j)==0 && es->getocci(k)==1 && es->getocci(l)==0){
	  //  printf("%i %i %i %i \n",i,j,k,l);
	  dE1=es->hoppEdiffijkl(i,j,k,l,dx,dy,dx2,dy2);
	  dE2=(dx+dx2)*Ex+(dy+dy2)*Ey+dE1;
	  exponent=beta*dE2;
	  if (dE2<0) rate = 1; 
	  else rate= exp(-exponent);
	  if(h==0)
	    rate *= GammaT2el[0];
	  else
	    rate *= GammaT2el[h]-GammaT2el[h-1];
	  if(rate>0){
	    rates2el[c2] = rate;
	    en2el[c2] = dE2;
	    c2++;
	  }
	  //  printf("2el: i = %i h = %i  rate = %e \n",i,h,rate);
	}	
      }
    }
  printf("c1 = %i  c2 = %i \n",c1,c2);

  n = c1;
  indexes = new int[n];

  for (n1=0;n1<n;n1++)
    indexes[n1]=n1;

  qsort(rates,indexes,0,n);
  
  

  f=FileCreate(filename+".1el");

  limit=20000;

  if (n<limit)
    {
      tot = 0;
      for(n1=0;n1<n;n1++)
	{
	  rate = rates[indexes[n1]];
	  tot+=rate;
	  st=IntToStr(n-n1) + "\t" + FloatToStr(rate)+ "\t" + FloatToStr(tot)+ "\t" + FloatToStr(en[indexes[n1]])+"\n";
	  FileWrite(f,st.c_str(),st.length());
     
	}
      FileClose(f);
    }
  else
    {
      tot = 0;
      for(n1=0;n1<n-limit;n1++)
	{
	  rate = rates[indexes[n1]];
	  tot+=rate;
	  //	  st=IntToStr(n-n1) + "\t" + FloatToStr(rate)+ "\t" + FloatToStr(tot)+"\n";
	  //      FileWrite(f,st.c_str(),st.length());
     
	}
       for(n1=n-limit;n1<n;n1++)
	{
	  rate = rates[indexes[n1]];
	  tot+=rate;
	  st=IntToStr(n-n1) + "\t" + FloatToStr(rate)+ "\t" + FloatToStr(tot)+ "\t" + FloatToStr(en[indexes[n1]])+"\n";
	  FileWrite(f,st.c_str(),st.length());
     
	}
     FileClose(f);
    }

 delete[] indexes;


  n = c2;
  indexes = new int[n];

  for (n1=0;n1<n;n1++)
    indexes[n1]=n1;
  printf("starter sortering\n");
  qsort(rates2el,indexes,0,n);
 printf("slutter sortering\n");
   
  

  f=FileCreate(filename+".2el");

  limit=20000;

  if (n<limit)
    {
      tot = 0;
      for(n1=0;n1<n;n1++)
	{
	  rate = rates2el[indexes[n1]];
	  tot+=rate;
	  st=IntToStr(n-n1) + "\t" + FloatToStr(rate)+ "\t" + FloatToStr(tot)+ "\t" + FloatToStr(en2el[indexes[n1]])+"\n";
	  FileWrite(f,st.c_str(),st.length());
     
	}
      FileClose(f);
    }
  else
    {
      tot = 0;
      for(n1=0;n1<n-limit;n1++)
	{
	  rate = rates2el[indexes[n1]];
	  tot+=rate;
	  //	  st=IntToStr(n-n1) + "\t" + FloatToStr(rate)+ "\t" + FloatToStr(tot)+"\n";
	  //      FileWrite(f,st.c_str(),st.length());
     
	}
       for(n1=n-limit;n1<n;n1++)
	{
	  rate = rates2el[indexes[n1]];
	  tot+=rate;
	  st=IntToStr(n-n1) + "\t" + FloatToStr(rate)+ "\t" + FloatToStr(tot)+ "\t" + FloatToStr(en2el[indexes[n1]])+"\n";
	  FileWrite(f,st.c_str(),st.length());
     
	}
     FileClose(f);
    }


  delete[] rates;
  delete[] rates2el;
  delete[] en;
  delete[] en2el;
  delete[] indexes;
	
}

void MKcurrent2el::qsort(double *val,int *pa, int beg, int end)
{
	if (end > beg + 1) {
		double piv = val[pa[beg]];
		int l = beg + 1, r = end;
		while (l < r) {
			if (val[pa[l]] <= piv)
				l++;
			else
				iswap(pa[l], pa[--r]);
		}
		iswap(pa[--l], pa[beg]);
		qsort(val,pa, beg, l);
		qsort(val,pa, r, end);
	}
};


void MKcurrent2el::runCurrent(int steps,double &E, double &t)
{
  int s, MCs, i, j, k,l,dx, dy,dx2,dy2;
  double de,dt,prob1e,de1,e1;
  bool jump1el,jump2el;

  
  tMC = 1/(N*es->nu()*totGammaTtot);
  
  prob1e = GammaTtot/totGammaTtot;

  for (s=0;s<steps;s++)
    {
      //printf(":");
      // currnptr=es->getnptr();
      MCs=0;
      jump1el=false; jump2el=false;
      while (!jump1el && !jump2el)
	{
	  //printf(".");
	  MCs++;
	  if ((es->ran2(0))<prob1e){
	    getjump(i,j,dx,dy);
	    jump1el=testjump(i,j,dx,dy);
	  }
	  else{
	    getjump2el(i,j, k,l,dx,dy,dx2,dy2);
	    jump2el=testjump2el(i,j, k,l,dx, dy,dx2,dy2);
	  }
	  //  if (MCs%10000==0) printf(".");
	}
      //de=es->hoppEdiffij(i,j);

      de = es->hopp(i,-1,0,j,0,0);
      from[s]=i;
      to[s]=j;
      dt = MCs*tMC; 
      E += de;
      t += dt;
      ts[s] = t;
      energy[s]=E;
      des[s]=de;
      dxs[s]=dx;
      dys[s]=dy;
      MCsteps[s]=MCs;
      //	  printf("de = %le \n",de);
      
      if (jump2el)
	{
	  s++;
	  de1 = de;
	  de = es->hopp(k,-1,0,l,0,0);
	  printf("de1 = %le de = %le de1+de = %le \n",de1, de,de1+de);

	  from[s]=k;
	  to[s]=l;
	  E += de;
	  ts[s] = t;
	  energy[s]=E;
	  des[s]=de;
	  dxs[s]=dx2;
	  dys[s]=dy2;
	  MCsteps[s]=-1;
	}
      
      //      e1 = es->calcenergy();
      //      printf("E = %le  e1 = %le \n",E,e1);
      if (s%1000==0) printf("%d E=%le MCs=%d\n",s,E,MCs);
      // printf("Step: %d: %d %d dE: %le I:%le MCs: %d\n", s, i, j,energy[s],dx,MCsteps[s]);
    }
}

/*void MKcurrent2el::toFile(string filename, int steps)
{
  int s;
  string st;
  FILE *f;

  f=FileCreate(filename);
  for(s=0;s<steps;s++)
    {
      st=IntToStr(from[s])+"\t"+IntToStr(to[s])+"\t"+FloatToStr(energy[s])+"\t"+FloatToStr(current[s])+"\t"+IntToStr(MCsteps[s])+"\n";
      FileWrite(f,st.c_str(),st.length());
    }
  FileClose(f);

}*/


void MKcurrent2el::jumpsToFileSmall(string filename, int steps)
{
  int s, x1, x2, y1, y2, dx, dy, cumdx;
  double cumc, r, dt, t, e, de;
  string st;
  FILE *f,*f2;

  cumdx=0;
  f=FileCreate(filename);
  f2 = FileCreate(filename+".2el");
  
  st="t  \t E \t from \t to \t dx \t MCs \t dE  tMC="+DoubleToStr(tMC)+"\n";
  FileWrite(f,st.c_str(),st.length());
  FileWrite(f2,st.c_str(),st.length());
  
  for(s=0;s<steps;s++)
    {
       cumdx+=dxs[s];
       if(MCsteps[s] == -1)
	 {
           st=FloatToStr(ts[s-1])+"\t"+DoubleToStr(energy[s-1])+"\t"+IntToStr(from[s-1])+"\t"+IntToStr(to[s-1])+"\t"+IntToStr(cumdx)+"\t"+IntToStr(MCsteps[s-1])+"\t"+FloatToStr(des[s-1])+"\n";
           FileWrite(f2,st.c_str(),st.length());
           st=FloatToStr(ts[s])+"\t"+DoubleToStr(energy[s])+"\t"+IntToStr(from[s])+"\t"+IntToStr(to[s])+"\t"+IntToStr(cumdx)+"\t"+IntToStr(MCsteps[s])+"\t"+FloatToStr(des[s])+"\n";
           FileWrite(f2,st.c_str(),st.length());
	 }
      if (s%10000==0){
        st=FloatToStr(ts[s])+"\t"+DoubleToStr(energy[s])+"\t"+IntToStr(from[s])+"\t"+IntToStr(to[s])+"\t"+IntToStr(cumdx)+"\t"+IntToStr(MCsteps[s])+"\t"+FloatToStr(des[s])+"\n";
        FileWrite(f,st.c_str(),st.length());
      }
  }
  FileClose(f);
  FileClose(f2);
}
