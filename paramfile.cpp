//---------------------------------------------------------------------------
// paramfile.cpp
//---------------------------------------------------------------------------
#include <stdio.h>
#include <math.h>
#include <string>

#pragma hdrstop

#include "fileutils.h"
#include "stringutils.h"

#include "paramfile.h"

//---------------------------------------------------------------------------

int paramfilereader::openread(string fn)
{
 int a,b,vs,i,bsize;
 FILE *f;
 char c;
 char *buffer;
 string str;
 if(fn!="") filename=fn;
 if(filename=="") return -1;
 f=FileOpen(filename,fmOpenRead);
 if(f!=NULL)
     {
      done();
      bsize = FileSize(f);
      buffer=new char[bsize+1];
      FileRead(f,buffer,bsize);
      FileClose(f);
      buffer[bsize]=0;
      vs=0;
      for(i=0;i<bsize;i++) {c=buffer[i];if((c==0x0D) || (c==0x0A)) buffer[i]=0;if(c=='=') vs++;}
      a=0;
      vlist=new string[vs+vs];

      while(a<bsize)
       {
        if(buffer[a]!=0)
         {str=&buffer[a];
          while((buffer[a]!=0) && (a<bsize)) a++;
          if((b=str.find('#'))!=string::npos) str.erase(b,str.length()-b); //remove comments
          if((b=str.find('='))!=string::npos) //paramter split
           {
            vlist[vnum+vnum]=TrimStr(str.substr(0,b));
            vlist[vnum+vnum+1]=TrimStr(str.substr(b+1,str.length()-b-1));
            vnum++;
           }
         }
        else a++;
       }
      delete[] buffer;
     }
 else return -1;


 return vnum;
}

int paramfilereader::cmdlineopenread(int n,char *cmdl[])
{
 if(n<2) return -2;
 if(FileExists(cmdl[1])) return openread(cmdl[1]);

 int i,a,b,vs;
 string str;

 vs=0;
 for(i=1;i<n;i++) {str=cmdl[i];if(str.find('=')!=string::npos) vs++;}
 if(vs>0)
 {
  done();
  vlist=new string[vs+vs];
  for(i=1;i<n;i++)
   {
    str=cmdl[i];
    if((b=str.find('='))!=string::npos)
    {
     vlist[vnum+vnum]=TrimStr(str.substr(0,b));
     vlist[vnum+vnum+1]=TrimStr(str.substr(b+1,str.length()-b-1));
     vnum++;
    }
   }
 }
 else return 0;
 return vnum;
}


string paramfilereader::getstring(string vname)
{
 for(int i=0;i<vnum+vnum;i+=2)
 {
    if(vname==vlist[i]) return vlist[i+1];
 }
 return "";
}


int paramfilereader::getint(string vname,int def)
{
 string s=getstring(vname);
 if(s!="") return StrToInt(s);
 return def;
}

//returns the number of read integers
int paramfilereader::getintarray(string vname,int *data,int maxl,char sep)
{
 int n,p;
 string s=getstring(vname);
 if (s=="") return 0;
 n=0;
 while(((p=s.find(sep,0))!=string::npos) && (n<maxl))
 {
    data[n]=StrToInt(s.substr(0,p));
    s.erase(0,p+1);
    n++;
 }
 if((s!="") && (n<maxl)) {data[n]=StrToInt(s);n++;}
 return n;
}

double paramfilereader::getdouble(string vname,double def)
{
 string s=getstring(vname);
 if(s!="") return StrToFloat(s);
 return def;
}


