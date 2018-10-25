//---------------------------------------------------------------------------
// paramfile
//---------------------------------------------------------------------------
#if !defined(paramfile_H)
  #define paramfile_H
  //---------------------------------------------------------------------------
  #include <stdio.h>
  #include <math.h>
  #include <string>

  using namespace std;
//---------------------------------------------------------------------------
class paramfilereader {

  public:
   paramfilereader(string fn="") {filename=fn;vnum=0;vlist=NULL;}
   ~paramfilereader() {done();}

   int openread(string fn="");

   int cmdlineopenread(int n,char *cmdl[]);

   string getstring(string vname);
   int getint(string vname,int def=0);
   int getintarray(string vname,int *data,int maxl,char sep=',');
   double getdouble(string vname,double def=0.0);


   void done() {if(vlist!=NULL) {delete[] vlist;};vlist=NULL;vnum=0;}

  private:
   string filename;
   string *vlist;
   int vnum;
};



#endif // paramfile_H
