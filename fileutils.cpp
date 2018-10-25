#include <stdio.h>
#include <string>
#include "fileutils.h"

using namespace std;

FILE* FileOpen(string FileName, int Mode)
{//long r=-1;
 FILE *f;
 switch(Mode) {

  case fmOpenRead:f=fopen(FileName.c_str(),"rb");break;
  case fmOpenWrite:f=fopen(FileName.c_str(),"r+b");break;
  case fmOpenReadWrite:f=fopen(FileName.c_str(),"r+b");break;
  case fmOpenAppend:f=fopen(FileName.c_str(),"a+b");break;
  case fmCreate:f=fopen(FileName.c_str(),"wb");break;
 }
 //if(f!=NULL) r= f;
 //printf("Debug: fopen: f=%d, r=%d\n",(int) f,r);
 return f;
}

void FileClose(FILE* Handle)
{
 if(Handle!=NULL) fclose( Handle);
}

FILE* FileCreate(string FileName)
{
 return FileOpen(FileName,fmCreate);

}

bool FileExists(string FileName)
{
 FILE* f=FileOpen(FileName,fmOpenRead);
 if(f!=NULL)
  {FileClose(f);
   return true;
  }
 return false;
} 

int FileRead(FILE* Handle, void *Buffer, int Count)
{
 return fread(Buffer,Count,1,Handle);
}

int FileWrite(FILE* Handle, const void *Buffer, int Count)
{
 return fwrite(Buffer,Count,1, Handle);
}
 
int FileSeek(FILE* Handle, int Offset, int Origin)
{
 if(Handle!=NULL)
  {
   fseek(Handle,Offset,Origin);
   return ftell( Handle);
  }
 return -1;
}

int FileSize(FILE* Handle)
{int s=-1;
 if(Handle!=NULL)
   {
    s = FileSeek(Handle,0,2);
    FileSeek(Handle,0,0);
   }
 return s;
}
void FileFlush(FILE* Handle)
{
 if(Handle!=NULL)
   {
    fflush( Handle);
   }
}

void DeleteFile(string FileName)
{
  remove(FileName.c_str());
}
