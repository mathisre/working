//---------------------------------------------------------------------------
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string>

#pragma hdrstop

#include "systemclass.h"
#include "stringutils.h"
#include "fileutils.h"

// Remember to remove later
#include <iostream>
using namespace std;

//---------------------------------------------------------------------------

string IntToHexStr(unsigned int a,int dig)
{
	string s="";
	for(int i=0;i<(dig<<2);i+=4) s=string(1,b32char[(a>>i)%16])+s;
	return s;
}

string hexstr(unsigned char *a,int l)
{
  string s="";
  for(int i=0;i<l;i++)
   {
	 s=s+string(1,b32char[a[i]>>4])+string(1,b32char[a[i]&0x0F]);
   }
  return s;
}

//---------------------------------------------------------------------------

ESystem::ESystem(int size, int dim, double nu,int rs)
{
 int i,j,k;
 double r;

 idum=2354;iyy=0;iff=0;
 ir=new int[100];

 L=size;
 D=dim;
 tracepart=false;
 hasSPE=false;
 hasdg=false;
 traceL=0;
 traceN=0;
 tracepos=0;
 switch (dim)
 {
   case 1:N=L;break;
   case 2:N=L*L;break;
   case 3:N=L*L*L;break;
 }
 occnum=nu;
 NP=(int) (N*nu);
 dis=new double[N];
 indist=new double[N];
 indistDrag=new double[N];
 disgrid=new double[D*N];
 n=new char[N];
 posperm=new int[N];
 density=new double[N];
 partnum=new int[N];
 partpos=new int[NP];
 sitesMoved = new int[N];


 siteen=new double[N];
 spe=new double[N];
 Rseed=rs;
 ran2(rs);
 U=1.0;
 for(i=0;i<N;i++)
  {
   dis[i]=2*ran2(0)-1.0;
   //printf("%i\t%.20e\n",i,dis[i]);
   n[i]=0;
   indist[i]=0.0;
   siteen[i]=0.0;
   spe[i]=0.0;
   density[i]=0.0;
   posperm[i]=i;
   partnum[i]=-1;
   sitesMoved[i]=0;
  }
 reconfig(); //initalize the occupied sites


 //calculate the inverse distances
 // Blir dette noen gang brukt?
 if(D==1) {for(i=1;i<L;i++) indist[i]=1.0/((double) i);}
 else if(D==2)
  {
   for(j=0;j<L;j++)
	{
	 for(i=0;i<L;i++)
	  {
       r = i*i; r += j*j;
       if(r > 0.0) indist[j*L+i] = 1/sqrt(r);
	  }
	}
  }
 else if(D==3)
  {
   for(k=0;k<L;k++)
	{
   for(j=0;j<L;j++)
	{
	 for(i=0;i<L;i++)
	  {
	   r=i*i;r+=j*j;r+=k*k;
	   if(r>0.0) indist[(k*L+j)*L+i]=1/sqrt(r);
	  }
	}
	}
  }
}

ESystem::~ESystem()
{
 delete[] indist;
 delete[] disgrid;
 delete[] dis;
 delete[] siteen;
 delete[] spe;
 delete[] n;
 delete[] ir;
 delete[] posperm;
 delete[] density;
 delete[] partnum;
 delete[] partpos;
}

//---------------------------------------------------------------------------

double ESystem::ran2(int iseed)
{
//       based on ran2 from Numerical Recipes page 197 (fortran)
//       meant to emulate the unix function rand
//       positive iseed reinitializes, zero iseed - same string
// ir - class-array int[100], ma,ia,ic,rm - constants
// idum,iff,iyy - class-var (iff=0 before first run)
double rand;
int j;
  if((iseed>0) || (iff==0))
   {iff=1;
    idum=-iseed-3;
    idum=(ic-idum)%ma;       //-> idum < ma
    for(j=0;j<97;j++)  {
      idum=(ia*idum+ic)%ma;      //-> idum < ma
      ir[j]=idum;
     }
    idum=(ia*idum+ic)%ma;         //-> idum < ma
    iyy=idum;                     //-> iyy <ma
   }
  j=(97*iyy)/ma;                  //-> i < 97
  iyy=ir[j];
  rand=iyy*rm;
  idum=(ia*idum+ic)%ma;
  ir[j]=idum;
  if(rand==0.0) rand=1e-14; //double !!
  return rand;
}

//---------------------------------------------------------------------------

void ESystem::getbinstate(unsigned char *state,int len)
{int i,j;     //N should be multiple of 8!
 for(i=0;i<len;i++) state[i]=0;
 i=0;j=0;
 while((i<len) && (j<N))
  {
   if(n[j]==1) state[i]=state[i]+(1<<(j%8));
   j++;
   if((j%8)==0) i++;
  }
};

int ESystem::statecorr(unsigned char *s1,unsigned char *s2,int len,int *pos) //return the number of diffent bits in s1&s2
{
 int c=0;
 int i,j,p;
 unsigned char x_or;
 for(i=0;i<len;i++)
  {
	  x_or=s1[i]^s2[i];
	  p=i<<3;
	  for(j=0;j<8;j++)
	   if(((x_or>>j)&1)==1)
		 {c++;
		  if(pos!=NULL) pos[p+j]++;
		 }
  }
 return c;
}

string ESystem::getstatestr()
{
 int i,dig;
 string s="";
 for(i=0;i<N;i+=5)  //N has to be multiple of 5
  {
   dig=(n[i]<<4)+(n[i+1]<<3)+(n[i+2]<<2)+(n[i+3]<<1)+n[i+4];
   s=s+b32char[dig];
  }
 return s;
};

string ESystem::statestring()
{
 int i,j;
 // unsigned char c;
 unsigned int t,id=0;
 string s="";
 //generate an "id" for a state
 t=n[0];
 for(i=1;i<N;i++)
  {
   j=i%32;
   if(j==0)
	{
     id=id^t;
     t=0;
    }
   t=t+(n[i]<<j);
  }
 id=id^t;
 s=IntToHexStr(id)+"@"+IntToHexStr(Rseed);
 return s;
};

string ESystem::systemID(bool istate)
{
	string tmp=IntToHexStr(Rseed)+"_"+IntToHexStr((int) (U*10),2);
	if(istate) tmp=tmp+"_"+IntToHexStr(InitRS);
	return tmp;
}

//---------------------------------------------------------------------------

void ESystem::resample(int rs)
{
 int i;
 ran2(rs);
 Rseed=rs;

 for(i=0;i<N;i++) dis[i]=2*ran2(0)-1.0;
 //reconfig();
};

void ESystem::reconfigfromfile(string statefilename)
{
  int i,bsize,ne,ns; //ne=nr electrons ns=nr sites
  FILE *f;
  char c;
  char *buffer;

  if (statefilename=="") return;
  f=FileOpen(statefilename,fmOpenRead);
  if(f!=NULL)
    {
      bsize = FileSize(f);
      buffer=new char[bsize+1];
      FileRead(f,buffer,bsize);
      FileClose(f);
      buffer[bsize]=0;
      ne=0;
      ns=0;
      for(i=0;i<bsize;i++)
	{
	  c=buffer[i];
	  if (c=='0')
	    {
	      if (ns>N) 
		{printf("Esystem::reconfigfromfile error, ns>N ns: %i N: %i \n",ns,N);return;}
	      n[ns]=0;
	      partnum[ns]=-1;
	      ns+=1;

	    }
	  else if (c=='1')
	    {
	      if (ns>N) 
		{printf("Esystem::reconfigfromfile error, ns>N ns: %i N: %i \n",ns,N);return;}
	      if (ne>NP) 
		{printf("Esystem::reconfigfromfile error, ne>NP ne: %i NP: %i \n",ne,NP);return;}     
	      n[ns]=1;
	      partnum[ns]=ne;
	      partpos[ne]=ns;
	      ns+=1;
	      ne+=1;
	    } 
	}
      if (ns!=N) 
	{printf("Esystem::reconfigfromfile error, ns!=N ns: %i N: %i \n",ns,N);return;}
      if (ne!=NP) 
	{printf("Esystem::reconfigfromfile error, ne!=NP ne: %i NP: %i \n",ne,NP);return;}     
      delete[] buffer;
    }
  else {printf("Esystem::reconfigfromfile error, Error in opening file\n");return;}
}

void ESystem::readSPEfromfile(string spesfilename)
{
  int i,bsize,ne,ns; //ne=nr electrons ns=nr sites
  FILE *f;
  char c;
  char  g[80];

  char *buffer;

  if (spesfilename=="") return;
  //  std::ifstream fin(spesfilename);

    f=fopen(spesfilename.c_str(),"r");
  if(f!=NULL)
    {

      for (i=0;i<N;i++)
	{
	    fscanf(f, "%s", &g);
	    spe[i] = atof(g);
	
	    //  printf("%.15f\n", spe[i]);
	}
      
    }
  else {printf("Esystem::SPEfromfile error, Error in opening file\n");return;}
  fclose(f);

  hasSPE = true;

}

void ESystem::reconfig(int type)
{//initalize the occupied sites
// Blir kjort 2 ganger, forst med type 0 ogsaa med type 1

 int i,j,k;
 int p,pmin,jcount,x,y;
 double emin,en;
 int *parray;
 NP = (int) (occnum*N);
 delete[] partpos;
 partpos = new int[NP];
 for(i=0;i<N;i++) {n[i] = 0;partnum[i] = -1;} // reset!
 InitRS = type;



 if(type==0)  //disorder ground state
  {
   parray=new int[N];
   for(i=0; i<N; i++) parray[i]=i;
   sortdis(parray); // random order
   for(i=0;i<NP;i++)
    {
       p = parray[i];
       n[p] = 1;
       partnum[p] = i;
       partpos[i] = p;
	}
   delete[] parray;
  }
 else if(type==1)  //random
  {
   i=0;
   InitRS=ABS(idum);
   while(i<NP)
   {
     k = (int) (N*ran2(0));
     if(n[k] == 0)
     {
         n[k] = 1;
         partnum[k] = i;
         partpos[i] = k;
         i++;
     }
	}
  }
 else if(type==2) //ordered state - only for half filling, D=2
  {
   if((NP==N/2) && (D==2))
	{
	 k=0;
	 for(i=0;i<L;i++)
	  for(j=0;j<L;j++)
	   {if(((i+j)%2==0)) {p=i+L*j;n[p]=1;partnum[p]=k;
						  partpos[k]=p;k++;}}
	}
  }
 else if(type==3)
  for(i=0;i<NP;i++)
   {n[i]=1;partnum[i]=i;partpos[i]=i; }//simple
 else if(type==4) //one by one   , 2D
  {
	parray=new int[N];
	for(i=0;i<N;i++) parray[i]=i;
	sortdis(parray);
	p=parray[0];
	n[p]=1;   //gs for one particle
	partnum[p]=0;
	partpos[0]=p;
	i=1;
	while(i<NP)
	 {
	  pmin=-1;
	  //jmin=-1;
	  emin=1e10;
	  jcount=0;
	  for(j=0;j<N;j++)
	   {
		p=parray[j];
		if(n[p]==0)
		 {
		  x=p%L;
		  y=p/L;
		  n[p]=1;
		  en=getsiteen(x,y,0);
		  if(en<emin)
		   {pmin=p;emin=en;}
		  n[p]=0;
		  jcount++;
		 }
		if(jcount>250) j=N;
	   }
	  n[pmin]=1;
	  partnum[pmin]=i;partpos[i]=pmin;
	  i++;
     }
    delete[] parray;


  }
}

//---------------------------------------------------------------------------

void ESystem::setgriddis(double maxdist)
{
 hasdg=false;
 if(maxdist>0.25) maxdist=0.25;
 else if(maxdist<0) maxdist=0;
 if(maxdist>1e-6) hasdg=true;
 if(hasdg)
  {
	for(int i=0;i<D*N;i++) disgrid[i]=maxdist*(2*ran2(0)-1.0);
  }
}

//---------------------------------------------------------------------------

 void ESystem::pshuffle()
 {//reshuffle the posperm permutation array
  int x,y,i,p;
  for(i=0;i<N;i++)
   {
	x=(int) (N*ran2(0)+0.5);x=x%N;
	y=(int) (N*ran2(0)+0.5);y=y%N;
	p=posperm[x];
	posperm[x]=posperm[y];
	posperm[y]=p;
   }
 }


 int ESystem::getppos(int i,int &x,int &y,int &z)
 {

  int p=posperm[i];
  y=0;z=0;
  if(D==1) x=p;
  else if(D==2) {x=p%L;y=p/L;}
  else if(D==3) {x=p%L;y=p/L;z=y/L;y=y%L;}
  return p;
 }

//---------------------------------------------------------------------------

double ESystem::getdis(int x,int y,int z)
{
 if((D==1) || (y<0)) return dis[x];
 else if(D==2) return dis[x+L*y];
 else if(D==3) return dis[x+L*(y+L*z)];
 return 0.0;
}

char ESystem::getocc(int x,int y,int z)
{
 if((D==1) || (y<0)) return n[x];
 else if(D==2) return n[x+L*y];
 return n[x+L*(y+L*z)];
}

//return the index for indist wrt the two given points
int inline ESystem::getdistidx(int x,int y,int z,int a,int b,int c)
{
  int dist=0;
  getlatdist(x,y,z,a,b,c);
  if(D>2) dist=ABS(c);
  if(D>1) dist=ABS(b)+L*dist;
  return  ABS(a)+L*dist;;
}

//calculate distance of (a,b,c)-(i,j,k) into (a,b,c), note |a|=<L/2, dito for b,c
void inline ESystem::getlatdist(int x,int y,int z,int &a,int &b,int &c)
{
 if(y<0)
   {
	if(D>1) {y=x/L;x=x%L;b=a/L;a=a%L;}
	if(D>2) {z=y/L;y=y%L;c=b/L;b=b%L;}
   }
  a=a-x;if(ABS(a+a)>=L) {if(a>0) a=a-L;else a=a+L;}
  if(D>1) {b=b-y;if(ABS(b+b)>=L) {if(b>0) b=b-L;else b=b+L;};};
  if(D>2) {c=c-z;if(ABS(c+c)>=L) {if(c>0) c=c-L;else c=c+L;};};
}

//to make this function as sufficient as possible coordinates and index can be given
//q is calculated if hasdg  or calcq=true
double inline ESystem::getindist(int x,int y,int z,int p,int dx,int dy,int dz,int &q,bool calcq)
{

  double rin;
  int a,b,c;
  if((dx==0) && (dy==0) && (dz==0)) return -1;
  if(D==1) rin=indist[ABS(dx)];
  else if(D==2) rin=indist[ABS(dy)*L+ABS(dx)];
  else if(D==3) rin=indist[(ABS(dz)*L+ABS(dy))*L+ABS(dx)];
  q=-1; //not calculated

  if(!calcq && !hasdg) return rin;

  //calculate all the missing coordinates & indices
  if(x<0) //index must be used
	 { x=p;
	   if(D>1) {y=p/L;x=p%L;}
	   if (D>2) {z=y/L;y=y%L;}
	 }
  else if(p<0) //coordinates must be used
	 {if(D==1) p=x;
	  else if(D==2) p=x+L*y;
	  else if(D==3) p=x+L*(y+L*z);
	 }
  //calculate the "destination" lattice point
  a=x+dx;PBC(a,L);q=a;
  if(D>1) {b=y+dy;PBC(b,L);q=a+L*b;}
  if(D>2) {c=z+dz;PBC(c,L);q=a+L*(b+L*c);}
  //now be have p=(x,y,z); q=(a,b,c) [PBC adjusted] and grid distances (dx,dy,dz) for "q-p"
  if(hasdg)   //additional code to handle distorted lattices
   {
	double rx,ry,rz;
	if(D==1) rin=ABS(1/(dx+disgrid[q]-disgrid[p]));
	else if(D==2)
	 {
	  if((ABS(dx)<5) && (ABS(dy)<5)) //sqrt calculation
	   {rx=dx+disgrid[q+q]-disgrid[p+p];
		ry=dy+disgrid[q+q+1]-disgrid[p+p+1];
		rin=1/sqrt(rx*rx+ry*ry);
	   }
	  else
	   {
		 rx=disgrid[q+q]-disgrid[p+p];
		 ry=disgrid[q+q+1]-disgrid[p+p+1];
		 rin=rin*(1-rin*rin*(dx*rx+dy*ry));
	   }
	 }
	else if(D==3)
	 {
	  if((ABS(dx)<5) && (ABS(dy)<5) && (ABS(dz)<5)) //sqrt calculation
	   {rx=dx+disgrid[q+q]-disgrid[p+p];
		ry=dy+disgrid[q+q+1]-disgrid[p+p+1];
		rz=dz+disgrid[q+q+2]-disgrid[p+p+2];
		rin=1/sqrt(rx*rx+ry*ry+rz*rz);}
	  else
	   {
		 rx=disgrid[q+q]-disgrid[p+p];
		 ry=disgrid[q+q+1]-disgrid[p+p+1];
		 rz=disgrid[q+q+2]-disgrid[p+p+2];
		 rin=rin*(1-rin*rin*(dx*rx+dy*ry+dz*rz));
	   }
	 }
   }
  return rin;
}

//---------------------------------------------------------------------------

void ESystem::sortdis(int * pa)
{
 //sort the dis array ascending by changing the
 //permutation array pa
  qsort(dis,pa,0,N);
}

void ESystem::qsort(double *val,int *pa, int beg, int end)
{//sort the dis array ascending by changing the
 //permutation array pa
 //qsort(dis,parray,0,N);
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
//---------------------------------------------------------------------------

void ESystem::copyntodensity()
{
 int i;
 for (i = 0; i < N; i++) {
   density[i]=n[i];
 }
};

void ESystem::transformdensityton()
{
 int *pa=new int[N];
 int i,A;
 for(i=0;i<N;i++) {pa[i]=i;n[i]=0;}
 qsort(density,pa,0,N-1);
 A=(int) (occnum*N);
 for(i=0;i<A;i++) n[pa[N-i-1]]=1;
 delete[] pa;
};
//---------------------------------------------------------------------------


double ESystem::calcenergy()
{
 double E=0.0;
 for(int i=0;i<N;i++) E+=getsiteen(i,-1,0);
 return E;
};

double ESystem::calcenergyfromSPE()
{
 double E=0.0;
 for(int i=0;i<N;i++) {
   siteen[i] = (U*dis[i]*(n[i]+ occnum) + spe[i]*(n[i]-occnum))/2;
   E+=getsiteen(i,-1,0,false,false);
 }
 return E;
};


double ESystem::getenergy()
{
 double E=0.0;
 for(int i=0;i<N;i++) E+=siteen[i];
 return E;
};


double ESystem::addenergyFromSecondLayer(ESystem *es2)
{
 double E=0.0;
 for(int i=0;i<N;i++) E+=addsiteenFromSecondLayer(i,-1,0, es2);
 return E;
}


//---------------------------------------------------------------------------

double ESystem::getsiteen(int x, int y, int z,bool calc,bool getspe)
{
 double rin,sed=0.0,sec=0.0;
 int p,q,Ni,m1,m2,m3;
 int off,off2,l1,l2,l3;
 if(y<0) //array index given in x
  {
   p=x;
   if(D>1) {x=p%L;y=p/L;}
   if(D>2) {z=y/L;y=y%L;}
  }
 else
  {if(D==1) p=x;
	  else if(D==2) p=x+L*y;
	  else if(D==3) p=x+L*(y+L*z);
  }
 if(!calc) {if(getspe)
			 {if(hasSPE) return spe[p];else return calcSPE(p);}
			return siteen[p];}
 if(D==1)
  {
   for(m1=1;m1<=rM;m1++)
	{
		rin=getindist(x,0,0,p,m1,0,0,l1);
		if(rin>=rmaxi)
		 {if(l1<0)
		   {l1=x+m1;PBC(l1,L);
			l2=x-m1;PBC(l2,L);
			sec+=rin*(n[l1]+n[l2]-occnum-occnum);}
		  else sec+=rin*(n[l1]-occnum);
		 }
		if(hasdg)
		  {
			 rin=getindist(x,0,0,p,-m1,0,0,l2);  //if hasdg, l2 is valid
			 if(rin>=rmaxi) sec+=rin*(n[l2]-occnum);
		  }
	}
  }
 else if(D==2)
  { for(m2=-rM;m2<=rM;m2++)
		 {if(!hasdg) {off=ABS(m2)*L;l2=y+m2;PBC(l2,L);l2=l2*L;}
		  for(m1=-rM;m1<=rM;m1++)
		   {if(hasdg) rin=getindist(x,y,0,p,m1,m2,0,q,true);
			else {rin=indist[off+ABS(m1)];l1=x+m1;PBC(l1,L);q=l2+l1;}
			if(rin>=rmaxi) sec+=rin*(n[q]-occnum);
		   }
		 }
  }
 else if(D==3)
  { for(m1=-rM;m1<=rM;m1++)
		   {if(!hasdg) {off2=ABS(m3)*L;l3=z+m3;PBC(l3,L);l3=l3*L;}
			for(m2=-rM;m2<=rM;m2++)
			 {if(!hasdg) {off=(off2+ABS(m2))*L;l2=y+m2;PBC(l2,L);l2=(l3+l2)*L;}
			  for(m3=-rM;m3<=rM;m3++)
			   {if(hasdg) rin=getindist(x,y,z,p,m1,m2,m3,q,true);
				else {rin=indist[off+ABS(m1)];l1=x+m1;PBC(l1,L);q=l2+l1;}
				if(rin>=rmaxi) sec+=rin*(n[q]-occnum);
			   }
			 }
		   }
  }
 Ni=n[p];
 if(Ni!=0) sed=U*Ni*dis[p];
 siteen[p]=sed+0.5*sec*(Ni-occnum);
 if(hasSPE) {spe[p]=calcSPE(p);if(getspe) return spe[p];}
 if(getspe) return calcSPE(p);
 return siteen[p];
}

//---------------------------------------------------------------------------


double ESystem::addsiteenFromSecondLayer(int x, int y, int z, ESystem *es2)
{
 double rin,sec=0.0;
 int p,q,Ni,m1,m2;
 int off,l1,l2;
 if(y<0) //array index given in x
   {
     p=x;
     if(D>1) {x=p%L;y=p/L;}
     if(D>2) {z=y/L;y=y%L;}
   }
 else
   {if(D==1) p=x;
     else if(D==2) p=x+L*y;
     else if(D==3) p=x+L*(y+L*z);
   }
 if(D==1)
   {
     printf("ERROR!!! addsiteenFromSecondLayer: Not for 1D\n"); 
   }
 else if(D==2)
   {for(m2=-rM;m2<=rM;m2++)
       {off=ABS(m2)*L;l2=y+m2;PBC(l2,L);l2=l2*L;
	 for(m1=-rM;m1<=rM;m1++)
	   {
	     rin=indistDrag[off+ABS(m1)];l1=x+m1;PBC(l1,L);q=l2+l1;
	     if(rin>=rmaxi) sec+=rin*(es2->getocci(q)-es2->nu());
	   }
       }
   }
 else if(D==3)
   {
     printf("ERROR!!! addsiteenFromSecondLayer: Not for 3D\n"); 
   }
 Ni=n[p];
 siteen[p] += 0.5*sec*(Ni-occnum);
 if(hasSPE) {spe[p]=calcSPE(p);}
 return siteen[p];
}

//---------------------------------------------------------------------------



double ESystem::calcSPE(int p)
{
 return ((2*siteen[p] - U*dis[p]*((double)n[p] + occnum)) / ((double)n[p] - occnum));
}

void ESystem::updateSPE()
{
 for(int i=0;i<N;i++) spe[i]=calcSPE(i);
};

void ESystem::updateSPEi(int i)
{
  spe[i]=calcSPE(i);
}

//---------------------------------------------------------------------------

double ESystem::getMinE(int &p,bool usespe,bool occc,char occs)
{
 double en,m=1e100;
 int j=-1;
 for(int i=0;i<N;i++)
  {if(!occc || (n[i]==occs))
	{if(!usespe) en=siteen[i];
	 else if(hasSPE) en=spe[i];
	 else en=calcSPE(i);
	 if(en<m) {m=en;j=i;}
	}
  };
 p=j;
 return m;
};

double ESystem::getMaxE(int &p,bool usespe,bool occc,char occs)
{
 double en,m=-1e100;
 int j=-1;
 for(int i=0;i<N;i++)
  {if(!occc || (n[i]==occs))
	{if(!usespe) en=siteen[i];
	 else if(hasSPE) en=spe[i];
	 else en=calcSPE(i);
	 if(en>m) {m=en;j=i;}
	}
  };
 p=j;
 return m;
};

//---------------------------------------------------------------------------
//this needs to be the best optimized procedure...
double ESystem::flipsite(int x,int y,int z,int p)
{ // We are passing in (x,y,z,n) for each site
  double rin,en,sed=0.0,sec=0.0,dE=0.0,de;
  int NSite,Nn,m1,m2,m3,q;
  int off,off2,l1,l2,l3;
  if(x < 0) //then p must be valid
  {
    x = p;
    if(D > 1) {x = p%L; y = p/L;}
    if(D > 2) {z = y/L; y = y%L;}
  }
  else if(p < 0)
  {
      if(D == 1) p = x;
      else if(D == 2) p = x + L*y;
      else if(D == 3) p = x + L*(y + L*z);
  }
  NSite = n[p];
  Nn = 0;
  if(NSite == 0) Nn = 1;
  NSite = Nn-NSite; // 1 hvis n[p] = 0, -1 hvis n[p] = 1
  n[p] = Nn;

//  cout << "rM: " << rM << endl;
  if(D==2)
  {
      for(m2 = -rM; m2 <= rM; m2++) // loopen går fra -49 til +49 (over alle partikler!
         {
          if(!hasdg) { // hasgd = 0
              off = ABS(m2)*L;
              l2 = y+m2;
              PBC(l2,L); // shifts l2 by L if  outside grid
              l2 = l2*L;
          }
          for(m1 = -rM; m1 <= rM; m1++)
          {
              if(hasdg) rin = getindist(x,y,0,p,m1,m2,0,q,true);
              else {
                  rin = indist[off + ABS(m1)];
                  l1 = x + m1;
                  PBC(l1,L);
                  q = l2 + l1;
              }
            if(rin >= rmaxi)
             {
                en = rin*(n[q]-occnum);
                sec += en;
                de = 0.5*en*NSite;
                dE += de;
                siteen[q] += de;
                if(hasSPE) spe[q] = calcSPE(q);
			 }
		   }
		 }
  }

  if(Nn!=0) sed = U * Nn * dis[p];
  en = sed + 0.5 * sec * (Nn-occnum);
  dE += en-siteen[p];
  siteen[p] = en;
  if(hasSPE) spe[p] = calcSPE(p);
  return dE;
}


//if(D==1)
//{
// for(m1=-rM;m1<=rM;m1++)
//  {   if(hasdg) rin=getindist(x,0,0,p,m1,0,0,q,true);
//      else {rin=indist[ABS(m1)];q=x+m1;PBC(q,L);}
//      if(rin>=rmaxi)
//       {en=rin*(n[q]-occnum);
//        sec+=en;
//        siteen[q]+=0.5*en*Nd;
//        if(hasSPE) spe[q]=calcSPE(q);
//       }
//  }
//}
//else
//else if(D==3)
//{
//   for(m3=-rM;m3<=rM;m3++)
//         {if(!hasdg) {off2=ABS(m3)*L;l3=z+m3;PBC(l3,L);l3=l3*L;}
//          for(m2=-rM;m2<=rM;m2++)
//           {if(!hasdg) {off=(off2+ABS(m2))*L;l2=y+m2;PBC(l2,L);l2=(l3+l2)*L;}
//            for(m1=-rM;m1<=rM;m1++)
//             {if(hasdg) rin=getindist(x,y,z,p,m1,m2,m3,q,true);
//              else {rin=indist[off+ABS(m1)];l1=x+m1;PBC(l1,L);q=l2+l1;}
//              if(rin>=rmaxi)
//               {en=rin*(n[q]-occnum);
//                sec+=en;
//                siteen[q]+=0.5*en*Nd;
//                if(hasSPE) spe[q]=calcSPE(q);
//               }
//             }
//           }
//         }
//}

void  ESystem::setindistDrag(double d)
{
  int i;
  indistDrag[0] = 1/d/d;
  for(i=1;i<N;i++) 
    {
      indistDrag[i] = 1/sqrt(1/(indist[i]*indist[i])+d*d);
    }
};


//2D only
double ESystem::flipsiteDrag(int x,int y,int z,int p,ESystem *es2)
{
  double rin,en,sed=0.0,sec=0.0,dE=0.0,de;
  int Nd,Nn,m1,m2,m3,q;
  int off,off2,l1,l2,l3;
  if(x<0) //then p must be valid
   {
	x=p;
	if(D>1) {x=p%L;y=p/L;}
	if(D>2) {z=y/L;y=y%L;}
   }
  else if(p<0)
   {
	  if(D==1) p=x;
	  else if(D==2) p=x+L*y;
	  else if(D==3) p=x+L*(y+L*z);
   }
  Nd=n[p];
  Nn=0;if(Nd==0) Nn=1;
  Nd=Nn-Nd;
  n[p]=Nn;

  if(D==1)
  {
    printf("ERROR!!! flipsiteDrag: Not for 1D\n"); 
  }
  else if(D==2)
    {for(m2=-rM;m2<=rM;m2++)
	{if(!hasdg) {off=ABS(m2)*L;l2=y+m2;PBC(l2,L);l2=l2*L;}
	  for(m1=-rM;m1<=rM;m1++)
	    {
	      rin=indist[off+ABS(m1)];l1=x+m1;PBC(l1,L);q=l2+l1;
	      if(rin>=rmaxi)
		{
		  en=rin*(n[q]-occnum);
		  sec+=en;
		  de=0.5*en*Nd;dE+=de;
		  siteen[q]+=de;
		  if(hasSPE) updateSPEi(q);
		}
	      rin = indistDrag[off+ABS(m1)];
	      if(rin>=rmaxi)
		{
		  en=rin*(es2->getocci(q)-es2->nu());
		  sec+=en;
		  de=0.5*en*Nd;dE+=de;
		  es2->addToSiteen(q,de);
		  if(hasSPE) es2->updateSPEi(q);
		}
	    }
	}
    }
  else if(D==3)
    {
      printf("ERROR!!! flipsiteDrag: Not for 3D\n"); 
    }
  if(Nn!=0) sed=U*Nn*dis[p];
  en=sed+0.5*sec*(Nn-occnum);
  dE+=en-siteen[p];
  siteen[p]=en;
  if(hasSPE) spe[p]=calcSPE(p);
  return dE;
};



//calculate the total energy change, if n[(i,j,k)] <-> n[(a,b,c)]
double ESystem::hoppEdiff(int i,int j,int k,int a,int b,int c)
{
  int p,q,m;
  double inp,inq,rin;
  if(j<0) //means that array index is used in i & a, and not coordinates
   {
	p=i;q=a;
	if(D>1) {i=p%L;j=p/L;a=q%L;b=q/L;}
	if(D>2) {k=j/L;j=j%L;c=b/L;b=b%L;}
   }
  else
   {p=i+L*(j+L*k);q=a+L*(b+L*c);}
  if(n[p]==n[q]) return 0.0;

  getlatdist(i,j,k,a,b,c);
  rin=getindist(i,j,k,p,a,b,c,m);
  //rin=indist[getdistidx(i,j,k,a,b,c)];
  if(rin<=rmaxi) {rin=0;/* printf("rin<rmaxi %f %f\n",rin,rmaxi);*/}  //last term in energy difference only relevant if site q&p interact

  m=n[q]-n[p];
  if(hasSPE)
   {
	 if(m==1) return (spe[p]-spe[q]-rin); //no double multiplication
	 return (spe[q]-spe[p]-rin);
   }
  inp=1/((double)n[p]-occnum);
  inq=1/((double)n[q]-occnum);

  return (2*m*(siteen[p]*inp-siteen[q]*inq)+
		  U*m*(dis[q]*(occnum+n[q])*inq-dis[p]*(occnum+n[p])*inp)-
		  rin); //m*m omitted in the last term  10 multipl & 2 div
};


//calculate the total energy change, if n[(i,j,k)] <-> n[(a,b,c)]
double ESystem::hoppEdiffij(int i,int j)
{
  int p,q,m,a,b,c,k;
  double inp,inq,rin;
    p = i; q = j;
    if(D>1) {i=p%L; j=p/L; a=q%L; b=q/L;}
    if(D>2) {k=j/L; j=j%L; c=b/L; b=b%L;}
  if(n[p] == n[q]) return 0.0;

  getlatdist(i,j,k,a,b,c);
  rin = getindist(i,j,k,p,a,b,c,m);
  //rin=indist[getdistidx(i,j,k,a,b,c)];
  if(rin <= rmaxi) {rin=0;/* printf("rin<rmaxi %f %f\n",rin,rmaxi);*/}  //last term in energy difference only relevant if site q&p interact

  m = n[q] - n[p]; // -1
  if(hasSPE) // single particle energy
   {
     if(m == 1) return (spe[p] - spe[q] - rin); //no double multiplication
     return (spe[q] - spe[p] - rin); // spe maa oppdateres i hvert steg
   }
  inp = 1/((double)n[p] - occnum); //ocnum er fill ratio
  inq = 1/((double)n[q] - occnum);

  return (2*m*(siteen[p]*inp - siteen[q]*inq) +
          U*m*(dis[q]*(occnum + n[q])*inq - dis[p]*(occnum + n[p])*inp)- rin);
  //m*m omitted in the last term  10 multipl & 2 div
}


//calculate the total energy change, if i->j and k->l
double ESystem::hoppEdiffijkl(int i,int j,int k,int l,  int dx, int dy, int dx2, int dy2)
{
  int p,m,a,b,c,z,xi,xj,xk,xl,yi,yj,yk,yl;
  double rinij,rinjl,rinjk,rinik,rinil,rinkl;
  
  //dx,dy,dx2,dy2 skal alle veare mindre enn L/2, saa ingen problemer med Coulomb
  rinij=indist[ABS(dy)*L+ABS(dx)];
  rinkl=indist[ABS(dy2)*L+ABS(dx2)];

  xi=i%L;yi=i/L;
  xj=j%L;yj=j/L;
  xk=k%L;yk=k/L;
  xl=l%L;yl=l/L;
  
  if(n[i]==0 || n[j]==1 || n[k]==0 || n[l]==1) {printf("FEIL!!! prover aa gjore umulige flytt!! hoppEdiffijkl");// printf("%i %i %i %i \n",i,j,k,l);
    // printf("%i %i %i %i \n",n[i],n[j],n[k],n[l]);
 return 0.0;}

  a = xk; b= yk;
  getlatdist(xi,yi,z,a,b,c);
  rinik=getindist(xi,yi,z,p,a,b,c,m);

  a = xl; b= yl;
  getlatdist(xi,yi,z,a,b,c);
  rinil=getindist(xi,yi,z,p,a,b,c,m);

  getlatdist(xj,yj,z,xk,yk,c);
  rinjk=getindist(xj,yj,z,p,xk,yk,c,m);

  getlatdist(xj,yj,z,xl,yl,c);
  rinjl=getindist(xj,yj,z,p,xl,yl,c,m);

  if(rinij<=rmaxi) {rinij=0;} 
  if(rinik<=rmaxi) {rinik=0;} 
  if(rinil<=rmaxi) {rinil=0;} 
  if(rinjk<=rmaxi) {rinjk=0;} 
  if(rinjl<=rmaxi) {rinjl=0;} 
  if(rinkl<=rmaxi) {rinkl=0;} 

  if(hasSPE)
   {
     return (spe[j]-spe[i]+spe[l]-spe[k] -rinij-rinkl+rinjl-rinjk-rinil+rinik); 

   }

  printf("FEIL!!! Programmet skal aldri komme hit! hoppEdiffijkl \n");
  return 10000;
};





//calculate the total energy change, if n[(i,j,k)] <-> n[(a,b,c)]
double ESystem::hoppEdiffxyz(int i,int j,int k,int a,int b,int c)
{
  int p,q,m;
  double inp,inq,rin;
  p=i+L*(j+L*k);q=a+L*(b+L*c);
  if(n[p]==n[q]) return 0.0;

  getlatdist(i,j,k,a,b,c);
  rin=getindist(i,j,k,p,a,b,c,m);
  //rin=indist[getdistidx(i,j,k,a,b,c)];
  if(rin<=rmaxi) {rin=0;/* printf("rin<rmaxi %f %f\n",rin,rmaxi);*/}  //last term in energy difference only relevant if site q&p interact

  m=n[q]-n[p];
  if(hasSPE)
   {
	 if(m==1) return (spe[p]-spe[q]-rin); //no double multiplication
	 return (spe[q]-spe[p]-rin);
   }
  inp=1/((double)n[p]-occnum);
  inq=1/((double)n[q]-occnum);

  return (2*m*(siteen[p]*inp-siteen[q]*inq)+
		  U*m*(dis[q]*(occnum+n[q])*inq-dis[p]*(occnum+n[p])*inp)-
		  rin); //m*m omitted in the last term  10 multipl & 2 div
}


//do the actual hop, with recalculation of all site energies
//all the virtual mechanisms are replaced
#include <iostream>
double ESystem::hopp(int i,int j,int k,int a,int b,int c) // (i,-1,0,j,0,0)
{
  int p,q,m;
  double dE = 0.0;

  if(j < 0) //means that array index is used in i & a, and not coordinates
   {
    p = i; q = a;
    if(D>1) {i=p%L; j=p/L; a=q%L; b=q/L;}
    if(D>2) {k=j/L; j=j%L; c=b/L; b=b%L;}
   }
  else {p = i+L*(j+L*k); q = a+L*(b+L*c);}

  if(n[p] == n[q]) { return dE;} // This should never occur, because we have made sure i is occupied and j is vacant

  //dE=getenergy(); //hoppEdiff(i,j,k,a,b,c);

  dE = flipsite(i,j,k,p); //oppdateres spe og occunumber
  dE += flipsite(a,b,c,q);

  if(tracepart) // false
  {
    m = partnum[p];
    partnum[p] = partnum[q];
    partnum[q] = m;
    if(m>= 0) partpos[m]=q;
    m = partnum[p];
    if(m >= 0) partpos[m] = p;
  }
  sitesMoved[p] = 1;
  sitesMoved[q] = 1;

  return dE; //getenergy()-dE;
}



double ESystem::hoppDrag(int i,int j,int k,int a,int b,int c,ESystem *es2)
{
 int p,q,m;
 double dE=0.0;
  if(j<0) //means that array index is used in i & a, and not coordinates
   {
	p=i;q=a;
	if(D>1) {i=p%L;j=p/L;a=q%L;b=q/L;}
	if(D>2) {k=j/L;j=j%L;c=b/L;b=b%L;}
   }
  else
   {p=i+L*(j+L*k);q=a+L*(b+L*c);}
  if(n[p]==n[q]) return dE;
  //dE=getenergy(); //hoppEdiff(i,j,k,a,b,c);
  dE=flipsiteDrag(i,j,k,p,es2);
  dE+=flipsiteDrag(a,b,c,q,es2);
  if(tracepart)
  {
   m=partnum[p];
   partnum[p]=partnum[q];
   partnum[q]=m;
   if(m>=0) partpos[m]=q;
   m=partnum[p];
   if(m>=0) partpos[m]=p;
  }
  sitesMoved[p] = 1;
  sitesMoved[q] = 1;

  return dE; //getenergy()-dE;
}


//---------------------------------------------------------------------------

double ESystem::foursiteXdiff(int p,int q,int i,int j) //get energy difference for 4-exchange; siteonly array indices
{
	 return 0.0;
}

double ESystem::foursiteX(int p,int q,int i,int j)
{
	 return 0.0;
}

//---------------------------------------------------------------------------

void ESystem::starttrace(int timesteps,int tracenum,int skipsteps)
{int i,j,p;

 //init the trace recorder
 traceL=timesteps/skipsteps+1;
 traceskip=skipsteps;
 traceN=tracenum;
 tracepart=true;
 motiontrace=new int[traceN*traceL];
 tracepos=0;


 //label all particles
 for(i=0;i<NP;i++) partpos[i]=i;
 p=0;
 pshuffle();
 for(i=0;i<N;i++)
  {
   j=posperm[i];
   partnum[j]=-1;
   if(n[j]==1)
	{
	 partnum[j]=p;
	 partpos[p]=j;
	 p++;
	}
  }
}

int ESystem::recordpos(int ts)
{
 int i;
 if(((ts%traceskip)!=0) || (tracepos==traceL)) return tracepos;
 for(i=0;i<traceN;i++)
  motiontrace[i*traceL+tracepos]=partpos[i];
 tracepos++;
 return tracepos;
}


void ESystem::stoptrace()
{
  traceL=0;
  traceN=0;
  tracepos=0;
  tracepart=false;
  delete[] motiontrace;
}

//---------------------------------------------------------------------------

double ESystem::NNBoltzstep()
{
  int i,x,x1,y,y1,z;
  double d,dl,dr,dt,db,a,DE;
  bool pb;
  DE=0.0;
  if(D==2)
   {
	 pshuffle();
	 for(i=0;i<N;i++)
	  {
	   getppos(i,x,y,z);
	   if(n[y*L+x]==1)
		{
		 //hopp to left
		 x1=x-1;if(x1<0) {x1+=L;pb=true;} else pb=false;
		 if(pb && !xpbc) dl=0; //no probability to jump left
		 else dl=exp(-hoppEdiff(x,y,0,x1,y,0)/T-E);

		 x1=x+1;if(x1>=L) {x1-=L;pb=true;} else pb=false;
		 if(pb && !xpbc) dr=0; //no probability to jump right
		 else dr=exp(-hoppEdiff(x,y,0,x1,y,0)/T+E);

		 y1=y-1;if(y1<0) y1+=L;
		 dt=exp(-hoppEdiff(x,y,0,x,y1,0)/T);

		 y1=y+1;if(y1>=L) y1-=L;
		 db=exp(-hoppEdiff(x,y,0,x,y1,0)/T);

		 d=dl+dr+dt+db+1.0;
		 a=d*ran2(0);
		 x1=x;y1=y;
		 if(a<=dl) {x1=x-1;if(x1<0) x1+=L;}
		 else if(a<=dl+dr) {x1=x+1;if(x1>=L) x1-=L;}
		 else if(a<=d-db-1) {y1=y-1;if(y1<0) y1+=L;}
		 else if(a<=d-1) {y1=y+1;if(y1>=L) y1-=L;}
		 if((x1!=x) || (y1!=y))
		  {
		   DE=DE+hopp(x,y,0,x1,y1,0);
		  }
		}
	 }
  }
  return DE;
};

double ESystem::densityLE(double dt, int steps)
{
 int x,y,m1,m2,l1,l2,off,ns,p;
 double rin;
 double *dd=new double[N];
 copyntodensity();

 for(ns=0;ns<steps;ns++)
  {
   if(D==2)
	{
	 for(y=0;y<L;y++)
		for(x=0;x<L;x++)
		 {
		  p=y*L+x;
		  dd[p]=-U*dis[p];
		  for(m2=-rM;m2<=rM;m2++)
		   {
			off=ABS(m2);
			l2=y+m2;if(l2>=L) l2-=L;else if(l2<0) l2+=L;
			for(m1=-rM;m1<=rM;m1++)
			 {
				rin=indist[ABS(m1)+off*L];
				if(rin>=rmaxi)
				 {l1=x+m1;if(l1>=L) l1-=L;else if(l1<0) l1+=L;
				  dd[p]-=rin*(density[l1+l2*L]-occnum);
				 }
			 }
		   }
		 }
	}
   for(p=0;p<N;p++) density[p]=density[p]+dt*dd[p];
  }

 transformdensityton();
 delete[] dd;
 return calcenergy();
};



//---------------------------------------------------------------------------

void ESystem::fillit(int *usedres,unsigned char *fill,int pos,bool vert)
{//cluster[pos]>0 & fill[pos]==0 assumed
 //pbc only vert if vert=false, horiz else
 int x,y,p;
 fill[pos]=1;
 if(D==2)
  {
   x=pos%L;y=pos/L;
   if((y<L-1) || !vert)
	  {p=L*pbc(y+1)+x;if((fill[p]==0) && (usedres[pos+pos+1]>0)) fillit(usedres,fill,p,vert);}
   if((y>0) || !vert)
	  {p=L*pbc(y-1)+x;if((fill[p]==0) && (usedres[p+p+1]>0)) fillit(usedres,fill,p,vert);}
   if((x<L-1) || vert)
	  {p=L*y+pbc(x+1);if((fill[p]==0) && (usedres[pos+pos]>0)) fillit(usedres,fill,p,vert);}
   if((x>0) || vert)
	  {p=L*y+pbc(x-1);if((fill[p]==0) && (usedres[p+p]>0)) fillit(usedres,fill,p,vert);}
  }
};

bool ESystem::clusterconnected(unsigned char *cluster,int *usedres)
{
 //find out if cluster is critical
 //pass a marker line/plane through the sytem to left to right and top to bottom (2D)
 //--> misses reverse tracks --> depreciated
 //new version: recursive "float-fill"  using only neigbors connected by a resistor


 if(D!=2) return true;
 unsigned char *fill=new unsigned char[N];
 int i,j;
 bool foundone=false;

 //top to bottom first
 for(i=0;i<N;i++) fill[i]=0;
 for(i=0;i<L;i++)
  {if((cluster[i]>0) && (fill[i]==0)) fillit(usedres,fill,i,true);}
 for(i=0;i<L;i++) {if(fill[i+N-L]==1) foundone=true;}

 if(!foundone)
  {//left to right
   for(i=0;i<N;i++) fill[i]=0;
   for(i=0;i<L;i++)
	{j=i*L;if((cluster[j]>0) && (fill[j]==0)) fillit(usedres,fill,j,false);}
   for(i=1;i<=L;i++) {if(fill[i*L-1]==1) foundone=true;}
  }
 delete[] fill;

 return foundone;
};


double ESystem::calcresist(double *resists,bool vert)
{
//calc the resistance of a given (connected) random network
return 0.0;
};

double ESystem::addnextresist(unsigned char *cluster,int *usedresist,double *resval,int &num)    //must have length N
{
  //find the next lowest resistor
  double *resist;
  double min,dE;
  int x,y,i,j,k,p1,p2;

  dE=0;
  resist=new double[D*N];
  for(i=0;i<D*N;i++) resist[i]=1e20;

  if(D==2)
	{
	 for(y=0;y<L;y++)
	  for(x=0;x<L;x++)
	  {
	   p1=x+L*y;
	   i=pbc(x+1);j=y;
	   p2=L*j+i;
	   if((n[p1]!=n[p2]) && (usedresist[p1+p1]==0))
		{
		 resist[p1+p1]=exp(hoppEdiff(x,y,0,i,j,0));
		}
	   i=x;j=pbc(y+1);
	   p2=L*j+i;
	   if((n[p1]!=n[p2]) && (usedresist[p1+p1+1]==0))
		{
		 resist[p1+p1+1]=exp(hoppEdiff(x,y,0,i,j,0));
		}
	  }
	 }
	 min=1e19;k=-1;
	 for(i=0;i<D*N;i++) {if(resist[i]<min) {min=resist[i];k=i;}}
	 if(k>-1)
	  {usedresist[k]=num; //>=1
	   resval[k]=min;
	   num++;
	   if(D==2)
		{
		 p1=k/2;
		 x=p1%L;y=p1/L;
		 i=x;j=y;
		 if(k%2==0) i=pbc(x+1);
		 else j=pbc(y+1);
		 p2=L*j+i;
		 dE=hopp(x,y,0,i,j,0);
		 cluster[p1]=1;
		 cluster[p2]=1;
		}
	  }
  //if k=-1: num will not be increased -> error
  delete[] resist;
  return dE;
};

void ESystem::resitrace(double *resi,double *minres,int* pred,int pos,bool vert)
{
  int x,y,p;
  double r,m=minres[pos];
  r=1e20;
 if(D==2)
  {
   x=pos%L;y=pos/L;
   if((y<L-1) || !vert)
	{p=L*pbc(y+1)+x;r=resi[pos+pos+1];
		if(r<1e10)
		{
		  if(m+r<minres[p])
		   {
			 minres[p]=m+r;
			 pred[p]=pos;
			 resitrace(resi,minres,pred,p,vert);
		   }
		}
	}
   if((y>0) || !vert)
	{p=L*pbc(y-1)+x;r=resi[p+p+1];
		if(r<1e10)
		{
		  if(m+r<minres[p])
		   {
			 minres[p]=m+r;
			 pred[p]=pos;
			 resitrace(resi,minres,pred,p,vert);
		   }
		}
	}
   if((x<L-1) || vert)
	{p=L*y+pbc(x+1);r=resi[pos+pos];
		if(r<1e10)
		{
		  if(m+r<minres[p])
		   {
			 minres[p]=m+r;
			 pred[p]=pos;
			 resitrace(resi,minres,pred,p,vert);
		   }
		}
	}
   if((x>0) || vert)
	{p=L*y+pbc(x-1);r=resi[p+p];
		if(r<1e10)
		{
		  if(m+r<minres[p])
		   {
			 minres[p]=m+r;
			 pred[p]=pos;
			 resitrace(resi,minres,pred,p,vert);
		   }
		}
	}

  }
};

//find the minimal resistance path through the system for a given el-config
bool ESystem::findminresistpath(double *resists,int *path)
{
 bool blocked=true;
 int i,j,x,y,p1,p2;
 double min;
 double *minresist;
 int *predec;
 for(i=0;i<D*N;i++) resists[i]=1e20;

 calcenergy();

 if(D==2)
	{
	 for(y=0;y<L;y++)
	  for(x=0;x<L;x++)
	  {
	   p1=x+L*y;
	   i=pbc(x+1);j=y;
	   p2=L*j+i;
	   // if(n[p1]!=n[p2]) resists[p1+p1]=exp(vhopp(x,y,0,i,j,0));
	   resists[p1+p1]=exp(absval(getsiteen(i,j,0,false)-getsiteen(x,y,0,false)));
	   i=x;j=pbc(y+1);
	   p2=L*j+i;
	   // if(n[p1]!=n[p2]) resists[p1+p1+1]=exp(vhopp(x,y,0,i,j,0));
	   resists[p1+p1+1]=exp(absval(getsiteen(i,j,0,false)-getsiteen(x,y,0,false)));
	  }
	 }
  //now all resistors are calculated, find the path with lowest resistance
  //through the system next:
  minresist=new double[N];
  predec=new int[N];
  for(i=0;i<N;i++) {minresist[i]=1e20;predec[i]=-1;path[i]=-1;}

  //from left to right:
  for(i=0;i<L;i++) minresist[i*L]=0;
  for(i=0;i<L;i++) resitrace(resists,minresist,predec,i*L,false);
  min=1e10;
  for(i=0;i<L;i++)
   {
	p2=(i+1)*L-1;
	if(minresist[p2]<min) {p1=p2;min=minresist[p2];blocked=false;}
   }
  if(!blocked)
   {
	path[0]=p1;
	i=1;
	while((p1=predec[p1])>=0)
	 {
	  path[i]=p1;
	  i++;
	 }
   }
  /* only horizontally
  for(i=0;i<N;i++) {minresist[i]=1e20;predec[i]=-1;}
  //from top to bottom:
  for(i=0;i<L;i++) minresist[i]=0;
  for(i=0;i<L;i++) resitrace(resists,minresist,predec,i,true);
  minv=1e10;
  for(i=L;i>=1;i--)
   {
	p2=N-i;
	if(minresist[p2]<minv) {p1=p2;minv=minresist[p2];blocked=false;}
   }

  if((minv<min) && !blocked)
   {
	for(i=0;i<N;i++) path[i]=-1; //overwrite the path
	path[0]=p1;
	i=1;
	while((p1=predec[p1])>=0)
	 {
	  path[i]=p1;
	  i++;
	 }
   }
  */
  delete[] minresist;
  delete[] predec;
  return blocked;
}

void ESystem::nstofile(string filename)
{
  int  i;
  FILE *f;
  string st;

  f = FileCreate(filename);
  for (i=0;i<N;i++){
      st = IntToStr(n[i]);
      if ((i+1)%L == 0) st += "\n";
      else st += "\t";

      FileWrite(f,st.c_str(),st.length());
    }
  FileClose(f);

}

void ESystem::spestofile(string filename)
{
  int  i;
  FILE *f;
  string st;

  f=FileCreate(filename);
  for (i=0;i<N;i++){
      st = FloatToStr(spe[i]);

      if ((i+1)%L == 0) st += "\n";
      else st += "\t";
      FileWrite(f, st.c_str(), st.length());
    }
  FileClose(f);

}

void ESystem::movedtofile(string filename)
{
  int  i;
  FILE *f;
  string st;

  f=FileCreate(filename);
  for (i=0;i<N;i++){
      if(sitesMoved[i] == 1){
        st = IntToStr(i)+"\n";
        FileWrite(f,st.c_str(),st.length());
      }
  }
  FileClose(f);

}

