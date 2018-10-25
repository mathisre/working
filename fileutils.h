#ifndef fileutilsH
#define fileutilsH

#include <string>

using namespace std;

#define fmOpenRead 0x0
#define fmOpenWrite 0x1
#define fmOpenReadWrite 0x2
#define fmOpenAppend 0x3
#define fmCreate 0xFF

#ifdef __cplusplus          

extern "C" {

#endif
             
FILE* FileOpen(string FileName, int Mode);

void FileClose(FILE* Handle);

FILE* FileCreate(string FileName);

bool FileExists(string FileName);

int FileRead(FILE* Handle, void *Buffer, int Count);

int FileWrite(FILE* Handle, const void *Buffer, int Count);
 
int FileSeek(FILE* Handle, int Offset, int Origin);
   
int FileSize(FILE* Handle);

void FileFlush(FILE* Handle);

void DeleteFile(string FileName);

#ifdef __cplusplus          

}

#endif

/* Origin
0	Der Dateizeiger wird Offset Bytes nach dem Dateianfang positioniert. 
1	Der Dateizeiger wird Offset Bytes nach der aktuellen Position positioniert. 
2	Der Dateizeiger wird Offset Bytes vor dem Dateiende positioniert.

*/
//-- end unit ----------------------------------------------------------------
#endif
