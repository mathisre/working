#include <stdlib.h>
#include <stdio.h>
#include <string>
#include "stringutils.h"

using namespace std;


int StrToInt(string s)
{
 return atoi(s.c_str());
}

string IntToStr(int value)
{char s[16];
 snprintf(s,16,"%d",value);
 return string(s);
}

string FloatToStr(double value)
{char s[32];
 snprintf(s,32,"%le",value);
 return string(s);
}

string DoubleToStr(double value)
{char s[32];
 snprintf(s,32,"%.8le",value);
 return string(s);
}

double StrToFloat(string s)
{
 return atof(s.c_str());	
}

string TrimStr(string s)
{int a,b,l;
 l=s.length();
 a=0;while((a<l) && (s[a]==' ')) a++;	
 b=l-1;while((b>=0) && (s[b]==' ')) b--;
 if(a<=b) return s.substr(a,b-a+1);
 return "";
}

