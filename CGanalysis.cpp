//---------------------------------------------------------------------------

#include <time.h>
#pragma hdrstop

#include "fileutils.h"
#include "stringutils.h"
#include "CGanalysis.h"
#include "systemclass.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)

//---------------------------------------------------------------------------
unsigned int CGcontrol::getsecs()
{
	 return clock()/CLOCKS_PER_SEC;
}
//---------------------------------------------------------------------------

void CGcontrol::setES(ESystem *e)
{
 es=e;
 //access some basic data from the system for internal use
 n=es->getnptr();
 dis=es->getdisptr();
 se=es->getseptr();
 D=es->dim();
 L=es->length();
 N=es->size();
 nu=es->nu(); // fill from params
 U=es->getU();
}
//---------------------------------------------------------------------------

void CGcontrol::addeecorr(double *ca,int cl,double avE)
{
  int i,j,x,y,n,c,Y;
  double norm,t;
  double *tmp;
  c=cl/2;
  n=0;
  if(es==NULL) return;
  if(D==2)
   {
	   tmp=new double[cl*cl];
	   for(i=0;i<cl*cl;i++) tmp[i]=0;

	   for(j=0;j<L;j++)
		 for(i=0;i<L;i++)
		 {
		  norm=se[j*L+i];       //if(ABS(norm)>1e-3) norm=1/norm;
		  for(y=0;y<cl;y++)
		   {Y=pbc(j+y-c)*L;
			 for(x=0;x<cl;x++)
			  {t=se[Y+pbc(i+x-c)]-norm;
			   tmp[y*cl+x]+=t*t;
			  }
		   }
		  n++;
		 }
	   for(i=0;i<cl*cl;i++) ca[i]=ca[i]+tmp[i]/n;
	   delete[] tmp;
   }
}

//---------------------------------------------------------------------------


double CGcontrol::findemin(int hopptyp,int rhopps,int &im, FMcallback callback)
{//lower energy using different methods given in hopptyp
 // # of random hops in rhopps
 // im is the algorithm type and return the maximal number of minimization steps needed
 // callback is am optional callback function to control the minimization process
  bool nnhopp,allhops,rhopp,fshopp,speopt;
  int state,stepst,i,j,Xc,XN,rM,alg;
  int p,q,x,z,y,x1,y1,z1,x2,y2,z2,xs,ys,zs;
  int *PP;
  unsigned int t,ts;
  double d,dm,oE,E,spe1,spe0;

  if(hopptyp==0) hopptyp=5; //NN+all

  nnhopp=((hopptyp&1)==1); //NN hops should always be set, could only be false for testing
  if(!nnhopp) printf("# Warning: no NN hopping in findemin performed.\n");

  rhopp=((hopptyp&2)==2);
  allhops=((hopptyp&4)==4);
  fshopp=((hopptyp&8)==8);
  speopt=((hopptyp&16)==16);
  alg=im; //the type of algorithm

  ts=getsecs();
  if(callback!=NULL) (*callback)(-2,0,0.0,0,ts,NULL);

  if(es==NULL) return 0.0;
  E=es->calcenergy();
  PP=es->getposperm();

  t=getsecs()-ts;
  if(callback!=NULL) (*callback)(-1,0,E,0,t,NULL);

  im=-1;   //returns the number of steps needed
  state=0; //current state of the search: 0: NN needs to be checked, 1: all NN are done, 2:all rand hops are done

  rM=es->getrM();

  XN=0; //total pair exchanges

  if((alg==11) && (callback!=NULL))   //save system energy, XN, spe histogram if alg=11
		  (*callback)(-11,XN,E,0,t,(void *) es);

  for(i=0;i<MAXFINDMIN;i++)
	 {
	  oE=E;
	  es->pshuffle();
	  Xc=0;  //pair exchanges in this iteration

	  stepst=10; //means: no exchange done

	  state=0; //nn
	  if(nnhopp)
	  {
		for(p=0;p<N;p++)
		 {
		   q=es->getppos(p,x,y,z);
		   dm=0.0;
		 
		   x1=x-1;if(x1<0) x1+=L;
		   d=es->hoppEdiff(x,y,z,x1,y,z);
		   dm=d; xs=x1;ys=y;zs=z;

		   x1=x+1;if(x1>=L) x1-=L;
		   d=es->hoppEdiff(x,y,z,x1,y,z);
		   if(d<dm) {dm=d;xs=x1;ys=y;zs=z;}

		   if(D>1)
			{
			 y1=y-1;if(y1<0) y1+=L;
			 d=es->hoppEdiff(x,y,z,x,y1,z);
			 if(d<dm) {dm=d;xs=x;ys=y1;zs=z;}

			 y1=y+1;if(y1>=L) y1-=L;
			 d=es->hoppEdiff(x,y,z,x,y1,z);
			 if(d<dm) {dm=d;xs=x;ys=y1;zs=z;}
			}
		   if(D>2)
			{
			 z1=z-1;if(z1<0) z1+=L;
			 d=es->hoppEdiff(x,y,z,x,y,z1);
			 if(d<dm) {dm=d;xs=x;ys=y;zs=z1;}

			 z1=z+1;if(z1>=L) z1-=L;
			 d=es->hoppEdiff(x,y,z,x,y,z1);
			 if(d<dm) {dm=d;xs=x;ys=y;zs=z1;}
			}
			if(dm<-1e-16) {
			  dm=es->hopp(x,y,z,xs,ys,zs);
			  E+=dm;
			  Xc++;
			  if(state<stepst) stepst=state; //remember that at least one NN hop was performed
			}
		   }
		  if((Xc==0) && (alg==11) && (callback!=NULL)) (*callback)(-12,XN,E,stepst,t,(void *) es);  //means we are in a local minimum
		 }
		 if(Xc==0) state=1;  //no nn

		 //  next: 4-site hops
		 if(fshopp && (state==1))
		  {
		 /*  for(p=0;p<N;p++)
			{
			 q=es->getppos(p,x,y,z);
			 dm=0.0;
			 if(D==2) //check only the "NN" aggregate
			 {
			   x1=es->pbc(x+1);y1=es->pbc(y+1);
			   d=es->foursiteXdiff(q,x1+y*L,x+y1*L,x1+y1*L);
			 }
			else if(D==3) //here 3 NN 4-site aggregates
			 {

			 }
			if(dm<-1e-16) {
			  //dm=es->hopp(x,y,z,xs,ys,zs);
			  E+=dm;
			  Xc++;
			  if(state<stepst) stepst=state; //remember that at least one 4site hop was performed
			}
		   }*/
		  }

		 if(Xc==0) state=2;  //no 4-site hop

		 //do random hops if no NN
		 if((state==2) && rhopp)
		  {
		   for(p=0;p<N;p++)
			{
			 q=es->getppos(p,x,y,z);
			 dm=0.0;
			 for(j=0;j<rhopps;j++)
			  {
			   x2=((int) (N*es->ran2(0)))%N;
			   d=es->hoppEdiff(q,-1,0,x2,0,0);
			   if(d<dm) {dm=d;xs=x2%L;ys=x2/L;zs=ys/L;ys=ys%L;}
			  }
			 if(dm<-1e-16) {
			  dm=es->hopp(x,y,z,xs,ys,zs);
			  E+=dm;
			  Xc++;
			  if(state<stepst) stepst=state; //remember that at least one rand. hop was performed
			}
		   }
		  }
		 if(Xc==0) state=3; //no NN or rhops

		
		 if(Xc==0) state=10; //none of the "simple" exchanges was used

		//end simpler routines

	   //stepst contains now the lowest ranked successful simple method
	   if(speopt && (stepst==10)) //now all simpler methods failed, use special SPE opimization
		{  es->setSPE(true); //switch on the internal (faster) spe calculation (which is slower for the above)
		   spe1=es->getMaxE(p,true,true,1); //get max occ spe
		   spe0=es->getMinE(q,true,true,0); //get min unocc spe
		   while(spe1>spe0)
			{  // d=es->hoppEdiff(p,-1,0,q,-1,0);
			   dm=es->hopp(p,-1,0,q,0,0);
			   E+=dm;
			   Xc++;
			   stepst=8;
			   spe1=es->getMaxE(p,true,true,1); //get max occ spe
			   spe0=es->getMinE(q,true,true,0); //get min unocc spe
			}
		  for(p=0;p<N;p++)
		   { dm=0.0;
			 q=es->getppos(p,x,y,z);
			 if(n[q]==1)
			  {
			   spe0=spe1-es->getsiteen(q,-1,0,false,true);   //>=0
			   if(spe0<=1.0)
				{
				  if(spe0>1e-8) j=(int) (1/spe0);
				  else j=rM;
				  if(j>rM) j=rM;
				  ys=0;zs=0;
				  if(D==1)
				   {for(x2=x-j;x2<=x+j;x2++)
					 {x1=pbc(x2);
					  if(n[x1]==0) {d=es->hoppEdiff(x,0,0,x1,0,0);
									if(d<dm) {dm=d;xs=x1;}}
					 }
				   }
				  else if(D==2)
				   {
					for(y2=y-j;y2<=y+j;y2++)
					 {y1=pbc(y2);
					  for(x2=x-j;x2<=x+j;x2++)
					   {x1=pbc(x2);
						if(n[x1+L*y1]==0) {d=es->hoppEdiff(x,y,0,x1,y1,0);
										   if(d<dm) {dm=d;xs=x1;ys=y1;}}
					   }
					 }
				   }
				  else if(D==3)
				   {
					for(z2=z-j;z2<=z+j;z2++)
					 {z1=pbc(z2);
					  for(y2=y-j;y2<=y+j;y2++)
					   {y1=pbc(y2);
						for(x2=x-j;x2<=x+j;x2++)
						 {x1=pbc(x2);
						  if(n[x1+L*(y1+L*z1)]==0) {d=es->hoppEdiff(x,y,z,x1,y1,z1);
													if(d<dm) {dm=d;xs=x1;ys=y1;zs=z1;}}
						 }
					   }
					 }
				   }
				  if(dm<0) {dm=es->hopp(x,y,z,xs,ys,zs);E+=dm;Xc++;stepst=8;}
				}
			  }
		   }
		}

	   if(allhops && (stepst==10)) //now everything failed: check all (can this be after speopt???)
		{ es->setSPE(true); //switch on the internal (faster) spe calculation
		  for(x=0;x<N;x++)
		   { dm=0.0;
			 p=PP[x];
			 for(q=0;q<N;q++)
			  {
				d=es->hoppEdiff(p,-1,0,q,0,0);
				if(d<dm) {dm=d;xs=q;}
			  }
			 if(dm<-1E-16)
			  {
			   dm=es->hopp(p,-1,0,xs,0,0);
			   E+=dm;
			   Xc++;
			   stepst=9; //some exchange occured
			  }
		   }
		}
	   XN+=Xc;
	   t=getsecs()-ts;
	   if(callback!=NULL)
		{(*callback)(i,Xc,E,stepst,t,NULL);
		 if(alg==11)   //save system energy, XN, spe histogram if alg=11
		  (*callback)(-11,XN,E,stepst,t,(void *) es);
		}
	   if(stepst==10) {im=i;i=MAXFINDMIN;}

   }
 es->setSPE(false); //switch the SPE off
 return E;
};

//---------------------------------------------------------------------------

//calculate the dipol moment
bool CGcontrol::dipol(double &Px,double &Py,double &Pz)
{
 int x,y,z,p;
 double hs=0.5*(L-1);
 Px=0.0;
 Py=0.0;
 Pz=0.0;
 if(es==NULL) return false;

 if(D==1)
  {
	for(x=0;x<L;x++) Px=Px+(n[x]-nu)*(x-hs);
  }
 else if(D==2)
  {
   for(y=0;y<L;y++)
	for(x=0;x<L;x++)
	 {p=L*y+x;
	  Px=Px+(n[p]-nu)*(x-hs);
	  Py=Py+(n[p]-nu)*(y-hs);
	 }
  }
  else if(D==3)
  {
   for(z=0;z<L;z++)
	for(y=0;y<L;y++)
	 for(x=0;x<L;x++)
	 {p=L*(L*z+y)+x;
	  Px=Px+(n[p]-nu)*(x-hs);
	  Py=Py+(n[p]-nu)*(y-hs);
	  Pz=Pz+(n[p]-nu)*(z-hs);
	 }
  }
 Px=Px/N;
 Py=Py/N;
 Pz=Pz/N;
 return true;
}

//---------------------------------------------------------------------------
//----------------------    FFT     ---------------------------------
//---------------------------------------------------------------------------


void CGcontrol::dfour1(double *data, unsigned long nn, int isign)
{
        unsigned long n,mmax,m,j,istep,i;
        double wtemp,wr,wpr,wpi,wi,theta;
        double tempr,tempi;

        n=nn << 1;
        j=1;
        for (i=1;i<n;i+=2) {
                if (j > i) {
                        SWAP(data[j],data[i]);
                        SWAP(data[j+1],data[i+1]);
                }
                m=n >> 1;
                while (m >= 2 && j > m) {
						j -= m;
                        m >>= 1;
                }
				j += m;
        }
        mmax=2;
        while (n > mmax) {
                istep=mmax << 1;
                theta=isign*(TWOPI/mmax);
                wtemp=sin(0.5*theta);
                wpr = -2.0*wtemp*wtemp;
                wpi=sin(theta);
                wr=1.0;
                wi=0.0;
                for (m=1;m<mmax;m+=2) {
                        for (i=m;i<=n;i+=istep) {
                                j=i+mmax;
                                tempr=wr*data[j]-wi*data[j+1];
                                tempi=wr*data[j+1]+wi*data[j];
                                data[j]=data[i]-tempr;
                                data[j+1]=data[i+1]-tempi;
                                data[i] += tempr;
                                data[i+1] += tempi;
                        }
                        wr=(wtemp=wr)*wpr-wi*wpi+wr;
                        wi=wi*wpr+wtemp*wpi+wi;
                }
                mmax=istep;
        }
}

//---------------------------------------------------------------------------

void CGcontrol::drealft(double *data, unsigned long n, int isign)
{
        unsigned long i,i1,i2,i3,i4,np3;
        double c1=0.5,c2,h1r,h1i,h2r,h2i;
        double wr,wi,wpr,wpi,wtemp,theta;

        theta=PI/(double) (n>>1);
        if (isign == 1) {
                c2 = -0.5;
                dfour1(data,n>>1,1);
        } else {
                c2=0.5;
				theta = -theta;
        }
        wtemp=sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi=sin(theta);
        wr=1.0+wpr;
        wi=wpi;
        np3=n+3;
        for (i=2;i<=(n>>2);i++) {
				i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
                h1r=c1*(data[i1]+data[i3]);
                h1i=c1*(data[i2]-data[i4]);
                h2r = -c2*(data[i2]+data[i4]);
                h2i=c2*(data[i1]-data[i3]);
                data[i1]=h1r+wr*h2r-wi*h2i;
                data[i2]=h1i+wr*h2i+wi*h2r;
                data[i3]=h1r-wr*h2r+wi*h2i;
                data[i4] = -h1i+wr*h2i+wi*h2r;
                wr=(wtemp=wr)*wpr-wi*wpi+wr;
                wi=wi*wpr+wtemp*wpi+wi;
        }
        if (isign == 1) {
                data[1] = (h1r=data[1])+data[2];
                data[2] = h1r-data[2];
		} else {
				data[1]=c1*((h1r=data[1])+data[2]);
				data[2]=c1*(h1r-data[2]);
				dfour1(data,n>>1,-1);
		}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

double CGcontrol::CGfindemin(int hopptyp,int &im)
{//lower energy by nn + rand hopps

  time_t starttime, currtime;
  double timeelapsed;

  bool nnhopp,nnnhopp,rhopp;
  int samec,x,x1,y,y1,i,j,xs,ys;
  double d,dm,oE,E;

  nnhopp=(hopptyp&1)==1;
  rhopp=(hopptyp&2)==2;
  nnnhopp=(hopptyp&4)==4;

  starttime =getsecs();

  printf("CGfindemin: starting calcbothenergies.\n");
  E = calcbothenergies();
  currtime = getsecs();
  timeelapsed = (currtime- starttime);
  printf("CGfindemin: calcbothenergies took %le s. \n", timeelapsed);

  //  printf("CGfindemin: starting calcenergy.\n");
  //  E=es->calcenergy();
  //  currtime = time(NULL);
  //  timeelapsed = difftime(currtime, starttime);
  //  printf("CGfindemin: calcenergy took %le s, starting calcaddenergy. \n", timeelapsed);
  //  es->calcaddenergy();
  //  currtime = time(NULL);
  //  timeelapsed = difftime(currtime, starttime);
  //  printf("CGfindemin: calcadd energy took %le s,starting minimum search.\n", timeelapsed);
 

  im=1;
  samec=0;
  for(i=0;i<MAXFINDMIN;i++)
    {
      oE=E;
	  for(x=0;x<L;x++)
	for(y=0;y<L;y++)
	  {
	    dm=0.0;
	    // NN hopping
	    if(nnhopp && (samec==0)) //if samec>0 : NN has not effect, since it's a local minimum
	      {
		x1=x-1;if(x1<0) x1+=L;
		d=es->hoppEdiff(x,y,0,x1,y,0);
		dm=d; xs=x1;ys=y;

		x1=x+1;if(x1>=L) x1-=L;
		d=es->hoppEdiff(x,y,0,x1,y,0);
		if(d<dm) {dm=d;xs=x1;ys=y;}

		y1=y-1;if(y1<0) y1+=L;
		d=es->hoppEdiff(x,y,0,x,y1,0);
		if(d<dm) {dm=d;xs=x;ys=y1;}

		y1=y+1;if(y1>=L) y1-=L;
		d=es->hoppEdiff(x,y,0,x,y1,0);
		if(d<dm) {dm=d;xs=x;ys=y1;}
	      }

	    if(nnnhopp && (samec==0)) //NNN hopping: 4 more vhopps
	      {
		x1=x+1;if(x1>=L) x1-=L;
		y1=y+1;if(y1>=L) y1-=L;
		d=es->hoppEdiff(x,y,0,x1,y1,0);
		if(d<dm) {dm=d;xs=x1;ys=y1;}

		x1=x-1;if(x1<0) x1+=L;
		y1=y+1;if(y1>=L) y1-=L;
		d=es->hoppEdiff(x,y,0,x1,y1,0);
		if(d<dm) {dm=d;xs=x1;ys=y1;}

		x1=x-1;if(x1<0) x1+=L;
		y1=y-1;if(y1<0) y1+=L;
		d=es->hoppEdiff(x,y,0,x1,y1,0);
		if(d<dm) {dm=d;xs=x1;ys=y1;}

		x1=x+1;if(x1>=L) x1-=L;
		y1=y-1;if(y1<0) y1+=L;
		d=es->hoppEdiff(x,y,0,x1,y1,0);
		if(d<dm) {dm=d;xs=x1;ys=y1;}
	      }

	    //random hopping
	    if(rhopp)
	      {
		for(j=0;j<4;j++)
		  {
			x1=x+(int) (L*es->ran2(0));if(x1>=L) x1-=L;
			y1=y+(int) (L*es->ran2(0));if(y1>=L) y1-=L;
			d=es->hoppEdiff(x,y,0,x1,y1,0);
			if(d<dm) {dm=d;xs=x1;ys=y1;}
		  }
	      }
	    
	    if(dm<0.0)
		  {
		   es->hopp(x,y,0,xs,ys,0);
		   E+=dm;
	      }
	  }
      if(oE==E) samec++;
      else samec=0;
	  if (samec==5) {im=i;i=MAXFINDMIN;}
	  if (i==0||i==MAXFINDMIN)
	{
	  currtime = getsecs();
	  timeelapsed = (currtime- starttime);
	  printf("CGfindmin: Step %d complete, %le s elapsed\n", im, timeelapsed);
	}    
    }
  return E;
};

//EoMK--CGfindemin

//MK--------------------------------------------------------

double CGcontrol::CGrandomizedfindemin(int hopptyp,int &im)
{//lower energy by nn + rand hopps

  time_t starttime, currtime;
  double timeelapsed;

  bool nnhopp,nnnhopp,rhopp;
  int samec,x,x1,y,y1,i,j, k,xs,ys;
  double d,dm,oE,E;

  nnhopp=(hopptyp&1)==1;
  rhopp=(hopptyp&2)==2;
  nnnhopp=(hopptyp&4)==4;

  starttime = getsecs();

  printf("CGfindemin: starting calcbothenergies.\n");
  E = calcbothenergies();
  currtime = getsecs();
  timeelapsed = (currtime - starttime);
  printf("CGfindemin: calcbothenergies took %le s. \n", timeelapsed);

  //  printf("CGfindemin: starting calcenergy.\n");
  //  E=es->calcenergy();
  //  currtime = time(NULL);
  //  timeelapsed = difftime(currtime, starttime);
  //  printf("CGfindemin: calcenergy took %le s, starting calcaddenergy. \n", timeelapsed);
  //  es->calcaddenergy();
  //  currtime = time(NULL);
  //  timeelapsed = difftime(currtime, starttime);
  //  printf("CGfindemin: calcadd energy took %le s,starting minimum search.\n", timeelapsed);
 

  im=1;
  samec=0;
  for(i=0;i<MAXFINDMIN;i++)
    {
      oE=E;
	  for(k=0;k<L*L;k++)
	{
	  x=(int) (L*es->ran2(0));if(x>=L) x-=L;
	  y=(int) (L*es->ran2(0));if(y>=L) y-=L;
	  dm=0.0;
	  // NN hopping
	  if(nnhopp && (samec==0)) //if samec>0 : NN has not effect, since it's a local minimum
	    {
		  x1=x-1;if(x1<0) x1+=L;
		  d=es->hoppEdiff(x,y,0,x1,y,0);
		  dm=d; xs=x1;ys=y;
	      
		  x1=x+1;if(x1>=L) x1-=L;
	      d=es->hoppEdiff(x,y,0,x1,y,0);
		  if(d<dm) {dm=d;xs=x1;ys=y;}
	      
		  y1=y-1;if(y1<0) y1+=L;
	      d=es->hoppEdiff(x,y,0,x,y1,0);
		  if(d<dm) {dm=d;xs=x;ys=y1;}
	      
		  y1=y+1;if(y1>=L) y1-=L;
	      d=es->hoppEdiff(x,y,0,x,y1,0);
		  if(d<dm) {dm=d;xs=x;ys=y1;}
	    }
	  
	  if(nnnhopp && (samec==0)) //NNN hopping: 4 more vhopps
	    {
		  x1=x+1;if(x1>=L) x1-=L;
		  y1=y+1;if(y1>=L) y1-=L;
	      d=es->hoppEdiff(x,y,0,x1,y1,0);
	      if(d<dm) {dm=d;xs=x1;ys=y1;}
	      
		  x1=x-1;if(x1<0) x1+=L;
		  y1=y+1;if(y1>=L) y1-=L;
	      d=es->hoppEdiff(x,y,0,x1,y1,0);
	      if(d<dm) {dm=d;xs=x1;ys=y1;}
	      
		  x1=x-1;if(x1<0) x1+=L;
		  y1=y-1;if(y1<0) y1+=L;
	      d=es->hoppEdiff(x,y,0,x1,y1,0);
	      if(d<dm) {dm=d;xs=x1;ys=y1;}
	      
		  x1=x+1;if(x1>=L) x1-=L;
		  y1=y-1;if(y1<0) y1+=L;
	      d=es->hoppEdiff(x,y,0,x1,y1,0);
	      if(d<dm) {dm=d;xs=x1;ys=y1;}
	    }
	  
	  //random hopping
	  if(rhopp)
	    {
	      for(j=0;j<4;j++)
		{
		  x1=x+(int) (L*es->ran2(0));if(x1>=L) x1-=L;
		  y1=y+(int) (L*es->ran2(0));if(y1>=L) y1-=L;
		  d=es->hoppEdiff(x,y,0,x1,y1,0);
		  if(d<dm) {dm=d;xs=x1;ys=y1;}
		}
	    }
	  
	  if(dm<0.0)
	    {
		  es->hopp(x,y,0,xs,ys,0);
	      E+=dm;
	    }
	}
      if(oE==E) samec++;
      else samec=0;
	  if (samec==5) {im=i;i=MAXFINDMIN;}
	  if (i==0||i==MAXFINDMIN)
	{
	  currtime = getsecs();
	  timeelapsed = (currtime- starttime);
	  printf("CGfindmin: Step %d complete, %le s elapsed\n", im, timeelapsed);
	}    
    }
  return E;
};
//EoMK randomized findmin


//MK--------------------------------------------------------

double CGcontrol::findStepInCG()
{
  //Assuming that Coulomb gap is established
  //Find and perform a step that lowers system energy
  //return energy change from step

  int CGrmax, rM;
  double Emax, Eloc, d, dm, CGrmaxi;
  int x,y,z,l1,l2,l3,m1,m2,m3,off, off2,xs,ys,zs;

  rM = es->getrM();
  Emax = getMAXoccen(x);

  if(D==1)
    {
	  for(x=0;x<L;x++)
	{
	  dm=0.0;
	  if (es->getocc(x,0,0)==1)
	    {
		  Eloc = getsiteadden(x,0,0,false);
	      CGrmaxi = Emax-Eloc; //big energy difference requires near neighbour to bridge...
	      if (ABS(CGrmaxi) < 1) //no shorter than nearest neighbor interaction
		{
		  if (CGrmaxi != 0) CGrmax=(int)((double)1/CGrmaxi);
		  else CGrmax = rM;
		  if (CGrmax > rM) CGrmax = rM;
		  
		  for(m1=-CGrmax;m1<=CGrmax;m1++)
		    {
		      l1=x+m1; if(l1>=L) l1-=L; else if (l1<0) l1+=L;
		      if (es->getocc(l1,0,0)==0)
			{
			  d=es->hoppEdiff(x,0,0,l1,0,0);
			  if(d<dm) {dm=d;xs=l1;}
			}
		    }
		  if (dm<0)
		    {
			  es->hopp(x,0,0,xs,0,0);
			  //printf("Saving hop... Energy change: %le\n",dm);
		      return dm;
		    }
		}
	    }
	}
    }
  else if (D==2)
    {
	  for(x=0;x<L;x++)
	{
	  for(y=0;y<L;y++)
	    {
	      dm=0.0;
	      if (es->getocc(x,y,0)==1)
		{
		  Eloc = getsiteadden(x,y,0,false);
		  CGrmaxi = Emax-Eloc; //big energy difference requires near neighbour to bridge...
		  if (ABS(CGrmaxi) < 1) //no shorter than nearest neighbor interaction
		    {
			  if (CGrmaxi != 0) CGrmax=(int)((double)1/CGrmaxi);
		      else CGrmax = rM;
		      if (CGrmax > rM) CGrmax = rM;
		      
		      for(m2=-CGrmax;m2<=CGrmax;m2++)
			{
			  off=ABS(m2);
			  l2=y+m2;if(l2>=L) l2-=L;else if(l2<0) l2+=L;
			  for(m1=-CGrmax;m1<=CGrmax;m1++)
			    {
			      l1=x+m1;if(l1>=L) l1-=L;else if(l1<0) l1+=L;
			      if (es->getocc(l1,l2,0)==0)
				{
				  d=es->hoppEdiff(x,y,0,l1,l2,0);
				  if(d<dm) {dm=d;xs=l1;ys=l2;}
				}
			    }
			}    
		      if (dm<0)
			{
			  es->hopp(x,y,0,xs,ys,0);
			  //printf("Saving hop... Energy change: %le\n",dm);
			   return dm;
			}
		    }
		}
	    }
	}
    }
  else if(D==3)
    {
	  for(x=0;x<L;x++)
	{
	  for(y=0;y<L;y++)
	    {
		  for(z=0;z<L;z++)
		{
		  dm=0.0;
		  if (es->getocc(x,y,z)==1)
		    {
			  Eloc = getsiteadden(x,y,z,false);
		      CGrmaxi = Emax-Eloc; //big energy difference requires near neighbour to bridge...
		      if (ABS(CGrmaxi) < 1) //no shorter than nearest neighbor interaction
			{
			  if (CGrmaxi != 0) CGrmax=(int)((double)1/CGrmaxi);
			  else CGrmax = rM;
			  if (CGrmax > rM) CGrmax = rM;

		      	  for(m1=-CGrmax;m1<=CGrmax;m1++)
			    {
			      off=ABS(m1);
			      l1=x+m1;if(l1>=L) l1-=L;else if(l1<0) l1+=L;
			      for(m2=-CGrmax;m2<=CGrmax;m2++)
				{
				  off2=ABS(m2);
				  l2=y+m2;if(l2>=L) l2-=L;else if(l2<0) l2+=L;
				  for(m3=-CGrmax;m3<=CGrmax;m3++)
				    {
				      l3=z+m3;if(l3>=L) l3-=L;else if(l3<0) l3+=L;
				      if (es->getocc(l1,l2,l3)==0)
					{
					  d=es->hoppEdiff(x,y,z,l1,l2,l3);
					  if(d<dm) {dm=d;xs=l1;ys=l2;zs=l3;}
					}
				    }
				}
			    }
			  if (dm<0)
			    {
				  es->hopp(x,y,z,xs,ys,zs);
			      //printf("Saving hop... Energy change: %le\n",dm);
			      return dm;
			    }
			}
		    }
		}
	    }
	}
    }
  return 0;
}
//EoMK findStepinCG


//MK--------------------------------------------------------------------------------------

double CGcontrol::checkAllSPjumps()
{
  //Assuming local minimum has been found, test all possible jumps
  //Perform jump that lowers energy
  //return energy of succesful jump
  int rM;
  double d, dm;
  int x,y,z,m1,m2,m3,off, off2,xs,ys,zs;

  rM = es->getrM();

  if(D==1)
    {
	  for(x=0;x<L;x++)
	{
	  dm=0.0;
	  if (es->getocc(x,0,0)==1)
	    {
		  for(m1=0;m1<L;m1++)
		{
		  if (es->getocc(m1,0,0)==0)
		    {
		      d=es->hoppEdiff(x,0,0,m1,0,0);
			  if(d<dm) {dm=d;xs=m1;}
		    }
		}
	      if (dm<0)
		{
		  es->hopp(x,0,0,xs,0,0);
		  //printf("Saving hop... Energy change: %le\n",dm);
		  return dm;
		}
	    }
	}
    }  
  else if (D==2)
    {
	  for(x=0;x<L;x++)
	{
	  for(y=0;y<L;y++)
	    {
	      dm=0.0;
	      if (es->getocc(x,y,0)==1)
		{
		  for(m2=0;m2<L;m2++)
		    {
			  for(m1=0;m1<L;m1++)
			{
			  if (es->getocc(m1,m2,0)==0)
			    {
				  d=es->hoppEdiff(x,y,0,m1,m2,0);
				  if(d<dm) {dm=d;xs=m1;ys=m2;}
				}
			}
		    }    
		  if (dm<0)
		    {
			  es->hopp(x,y,0,xs,ys,0);
			  //printf("Saving hop... Energy change: %le\n",dm);
		      return dm;
		    }
		}
	    }
	}
    }
  else if(D==3)
    {
	  for(x=0;x<L;x++)
	{
	  for(y=0;y<L;y++)
		{
		  for(z=0;z<L;z++)
		{
		  dm=0.0;
		  if (es->getocc(x,y,0)==1)
		    {
			  for(m1=0;m1<L;m1++)
			{
			  for(m2=0;m2<L;m2++)
			    {
				  for(m3=0;m3<L;m3++)
				{
				  if (es->getocc(m1,m2,m3)==0)
				    {
					  d=es->hoppEdiff(x,y,z,m1,m2,m3);
					  if(d<dm) {dm=d;xs=m1;ys=m2;zs=m3;}
				    }
				}
			    }
			}
		      if (dm<0)
			{
			  es->hopp(x,y,z,xs,ys,zs);
			  //printf("Saving hop... Energy change: %le\n",dm);
			  return dm;
			}
		    }
		}
	    }
	}
    }
  return 0;
}
//EoMK checkAllSPjumps


//MK----------------------------------------------------------------------------------

double CGcontrol::optimizeCG()
{
  //requires updated siteadden!!!
  //Assumes CG to have formed, but not being perfect
  //Return energy gained by optimization

  time_t starttime, currtime;
  double timeelapsed;

  int changes1, changes2,i,L,off,off2,D;
  double d,dm,oE,pE,E, dmtot1, dmtot2;

  int maxOptSteps = 1000000;

  starttime = getsecs();

  oE=0;
  pE=0;

  changes1=0;
  dmtot1=0;

  printf("optimizeCG started \n");
  for(i=0;i<maxOptSteps;i++)
    {
	  dm = findStepInCG();
      dmtot1+=dm;
      if (dm==0) {changes1=i;i=maxOptSteps;}
      //printf("Energy changed by %le in step %d\n",dm,i);
    }

  currtime = getsecs();
  timeelapsed = currtime- starttime;
  printf("optimizeCG: 1. step finished after %d steps in %le s.\n",changes1,timeelapsed);
  printf("Energy gain from 1. step optimization:  %le\n", dmtot1);

  dmtot2=0;
  changes2=0;

  for(i=0;i<maxOptSteps;i++)
    {
	  dm = checkAllSPjumps();
      dmtot2+=dm;
      if (dm==0) {changes2=i;i=maxOptSteps;}
      //printf("Energy changed by %le in step %d\n",dm,i);
    }
  currtime = getsecs();
  timeelapsed = currtime- starttime;
  printf("optimizeCG: 2. step finished after %d steps in %le s.\n",changes2,timeelapsed);
  printf("Energy gain from 2. step optimization:  %le\n", dmtot2);
  return dmtot1+dmtot2;
}
//EoMK 

//MK--------------------------------------------------------------------------------------

void CGcontrol::makeHistogram( int nrBoxes, int *histogram, int type, double highestE, double lowestE) //type=1: occ, 2: uocc, 3:all
{
  int x,y,z,box;
  double e, histBoxSize;


  histBoxSize = (highestE-lowestE)/nrBoxes;

  for (x=0;x<nrBoxes;x++) histogram[x]=0;

  if (D == 1)
    {
	  for (x=0;x<L;x++)
	{
	  if ( ((es->getocc(x,0,0) == 1)&&type==1) ||
	       ((es->getocc(x,0,0) == 0)&&type==2) ||
	       type==3)
	    {
		  e=getsiteadden(x,0,0,false);
	      box = (int)((e-lowestE)/histBoxSize+0.5);
	      if (box>nrBoxes-1) box=nrBoxes-1; //Error for Emax=>nrBoxes+1;
	      if (box<0) box=0;
	      histogram[box]++;
	    }
	}
    }
  else if (D == 2)
    {
	  for (y=0;y<L;y++)
	{
	  for (x=0;x<L;x++)
	    {
	      if ( ((es->getocc(x,y,0) == 1)&&type==1) ||
		   ((es->getocc(x,y,0) == 0)&&type==2) ||
		   type==3)
		{
		  e=getsiteadden(x,y,0,false);
		  box = (int)((e-lowestE)/histBoxSize+0.5);
		  if (box>nrBoxes-1) box=nrBoxes-1; //Error for Emax=>nrBoxes+1;
		  if (box<0) box=0;
		  histogram[box]++;
		}
	    }
	}
    }
  else if (D == 2)
    {
	  for (z=0;z<L;z++)
	{
	  for (y=0;y<L;y++)
	    {
		  for (x=0;x<L;x++)
		{
		  if ( ((es->getocc(x,0,0) == 1)&&type==1) ||
			   ((es->getocc(x,0,0) == 0)&&type==2) ||
			   type==3)
			{
			  e=getsiteadden(x,y,z,false);
			  box = (int)((e-lowestE)/histBoxSize+0.5);
			  if (box>nrBoxes-1) box=nrBoxes-1; //Error for Emax=>nrBoxes+1;
			  if (box<0) box=0;
			  histogram[box]++;
			}
		}
		}
	}
	}
}
//--------------------------------------------------------------------------------------

double CGcontrol::checkMicrocrystallinity(string filename)
{
  int *correlation1;

  int totalABS1;
  double av1;

  int x,y,z,xm,ym,zm,xn,yn,zn,l1,l2,l3,m1,m2,m3, ni, nj;

  xm=L/2;
  xn=xm;
  if ((xn+xm)==L) xn--;
  ym=xm;
  yn=ym;
  if ((yn+ym)==L) yn--;
  zm=xm;
  zn=zm;
  if ((zn+zm)==L) zn--;

  correlation1 = new int[N];

  for (x=0;x<N;x++) correlation1[x]=0;

  if(D==1)
	{
	  for(x=0;x<L;x++)
	{
	  ni = es->getocc(x,0,0);
	  for(m1=-xm;m1<=xn;m1++)
		{
		  l1=x+m1; if(l1>=L) l1-=L; else if (l1<0) l1+=L;
		  nj = es->getocc(l1,0,0);
		  if (ni==nj)
		{
		  correlation1[xm+m1]+=1;
		}
		  else {correlation1[xm+m1]-=1;}
		}
	}
	}
  else if (D==2)
	{
	  for(y=0;y<L;y++)
	{
	  for(x=0;x<L;x++)
		{
		  ni= es->getocc(x,y,0);
		  for(m2=-ym;m2<=yn;m2++)
		{
		  l2=y+m2;if(l2>=L) l2-=L;else if(l2<0) l2+=L;
		  for(m1=-xm;m1<=xn;m1++)
			{
			  l1=x+m1;if(l1>=L) l1-=L;else if(l1<0) l1+=L;
			  nj=es->getocc(l1,l2,0);
			  if (ni==nj) {
			correlation1[xm+m1+L*(ym+m2)]+=1;
			  }
			  else {correlation1[xm+m1+L*(ym+m2)]-=1;}
			}
		}
		}
	}
	}
  else if(D==3)
	{
	  for(x=0;x<L;x++)
	{
	  for(y=0;y<L;y++)
		{
		  for(z=0;z<L;z++)
		{
		  ni=es->getocc(x,y,z);
		  for(m1=-xm;m1<=xn;m1++)
			{
			  l1=x+m1;if(l1>=L) l1-=L;else if(l1<0) l1+=L;
			  for(m2=-ym;m2<=yn;m2++)
			{
			  l2=y+m2;if(l2>=L) l2-=L;else if(l2<0) l2+=L;
			  for(m3=-zm;m3<=zn;m3++)
				{
				  l3=z+m3;if(l3>=L) l3-=L;else if(l3<0) l3+=L;
				  nj = es->getocc(l1,l2,l3);
				  if (ni==nj)
				{
				  correlation1[xm+m1+L*(ym+m2)+L*L*(zm+m3)]+=1;
				}
				  else {correlation1[xm+m1+L*(ym+m2)+L*L*(zm+m3)]-=1;}
				}
			}
			}
		}
		}
	}
	}
  //printf("Correlation1\n");
  //  printCorrelation(D,L,N,correlation1);
  saveCorrelation(correlation1,filename);

  for (x=0;x<N;x++)
	{
	  totalABS1 += ABS(correlation1[x]);
	}
  av1=totalABS1/N;
  printf("Correlation1: %le\n",av1);

  return av1;
}
//--------------------------------------------------------------------------------------

void CGcontrol::saveHistogram(int *histogram, int nrBoxes, string filename)
{
  string output;
  int x;
  FILE *f;
  f=FileCreate(filename);
  output = "U: "+FloatToStr(U)+"\n";
  FileWrite(f,output.c_str(),output.length());
  for (x=0;x<nrBoxes;x++)
	{
      output = IntToStr(histogram[x])+"\t";
      FileWrite(f,output.c_str(),output.length());
    }
  
  FileClose(f);
}
//--------------------------------------------------------------------------------------

void CGcontrol::printCorrelation(int *correlation)
{
  int x, y, z;
  int highestABS, space;
  
  for (x=0;x<N;x++){if (ABS(correlation[x])>highestABS) highestABS = ABS(correlation[x]);} 
  
  space = IntToStr(highestABS).length()+2;
  
  if (D==1)
    {
      for (x=0;x<L;x++)
	{
	  printf("%*d",space,correlation[x]);
	}
      printf("\n");
    }
  else if (D==2)
    {
      for (y=0;y<L;y++)
	{
	  for (x=0;x<L;x++)
	    {
	      printf("%*d",space,correlation[x+L*y]);
	    }
	  printf("\n");
	}
    }
  else if (D==3)
    {
      for (z=0;z<L;z++)
	{
	  printf("level %d\n",z);
	  for (y=0;y<L;y++)
	    {
	      for (x=0;x<L;x++)
		{
		  printf("%*d",space,correlation[x+L*y+L*L*z]);
		}
	      printf("\n");
	    }
	  printf("\n");
	}
    }
}

//--------------------------------------------------------------------------------------


void CGcontrol::saveCorrelation(int *correlation, string filename)
{
  string output;
  int  x, y;
  FILE *f;
  f=FileCreate(filename);
  //  output = "U: ";
  //FileWrite(f,output.c_str(),output.length());
  //output = FloatToStr(p->disU)+"\n";
  //FileWrite(f,output.c_str(),output.length());
  output="";

  if (D==1)
    {
      for (x=0;x<L;x++)
	{
	  output+=(IntToStr(correlation[x])+"\t");

	}
      output+="\n";
      FileWrite(f,output.c_str(),output.length());
    }
  else if (D==2)
    { 
      for (y=0;y<L;y++)
	{
	  for (x=0;x<L;x++)
	    {
	      output+=(IntToStr(correlation[x+L*y])+"\t");
	    }
	  output+="\n";
	  FileWrite(f,output.c_str(),output.length());
	  output="";
	}
    }
  else if (D==3)
    {

    }
  
  FileClose(f);
}
//--------------------------------------------------------------------------------------

int histogram(double *data,int L,int *bins,int N,double dmin, double dmax)
{ int i,j,res;
  double x,delta;

  for(i=0;i<N;i++) bins[i]=0;
  delta=(dmax-dmin)/N;
  res=0;
  for(i=0;i<L;i++)
   {
	 x=(data[i]-dmin)/delta;
	 j=(int) x;
	 if((j>=0) && (j<N)) {res++;bins[j]++;}
   }
  return res;
}
