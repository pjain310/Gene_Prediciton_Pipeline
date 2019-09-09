#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "GBase.h"
#ifdef __WIN32__
  #include <windows.h>
#endif  
static char msg[4069];
//************************* Debug helpers **************************
// Assert failed routine
void GAssert(const char* expression, const char* filename, unsigned int lineno){
  sprintf(msg,"%s(%d): GASSERT(%s) failed.\n",filename,lineno,expression);
  fprintf(stderr,"%s",msg);
  }
// Error routine (prints error message and exits!)
void GError(const char* format,...){
  #ifdef __WIN32__	
    va_list arguments;
    va_start(arguments,format);
    vsprintf(msg,format,arguments);
    va_end(arguments);
    OutputDebugString(msg);
    fprintf(stderr,"%s",msg); // if a console is available
    MessageBox(NULL,msg,NULL,MB_OK|MB_ICONEXCLAMATION|MB_APPLMODAL);
    DebugBreak();
  #else
    va_list arguments;
    va_start(arguments,format);
    vfprintf(stderr,format,arguments);
    va_end(arguments);
    //abort();
  #endif
    exit(1);
  }
// Warning routine (just print message without exiting)
void GMessage(const char* format,...){
  va_list arguments;
  va_start(arguments,format);
  vsprintf(msg,format,arguments);
  va_end(arguments);
  #ifdef __WIN32__
    OutputDebugString(msg);
  #endif
  fprintf(stderr,"%s",msg);fflush(stderr);
  }

/*************** Memory management routines *****************/
// Allocate memory
bool GMalloc(pointer* ptr,unsigned long size){
  //GASSERT(ptr);
  if (size!=0) *ptr=malloc(size);
  //GMessage("GMalloc here! %d\n", size);
  return *ptr!=NULL;
  }
  
// Allocate cleaned memory (0 filled)
bool GCalloc(pointer* ptr,unsigned long size){
  GASSERT(ptr);
  *ptr=calloc(size,1);
  return *ptr!=NULL;
  }
// Resize memory
bool GRealloc(pointer* ptr,unsigned long size){
  GASSERT(ptr);
  if (size==0) {
    GFree(ptr);
    return true;
    }
  void *p=realloc(*ptr,size);
  if (p) {
      *ptr=p;
      return true;
      }
  return false;
  }
// Free memory, resets ptr to NULL afterward
void GFree(pointer* ptr){
  GASSERT(ptr);
  if (*ptr) free(*ptr);
  *ptr=NULL;
  }

char* Gstrdup(const char* str) {
  char *copy;
  GMALLOC(copy, strlen(str)+1);
  strcpy(copy,str);
  return copy;
  }

char* Gsubstr(const char* str, char* from, char* to) {
 //extract (and allocate) a substring, including boundaries (from/to)
 int len=strlen(str);
 if (to==NULL) to=(char *)(str+len-1); //extract tail from 'from' char, including it
 if (!(from>=str && to>from && to<str+len)) return NULL;
 int newlen=to-from+1;
 char* subs;
 GMALLOC(subs, newlen);
 memcpy(subs, str, newlen-1);
 subs[newlen]='\0';
 return subs;
 }

void* Gmemscan(void *mem, unsigned int len,
                   void *part, unsigned int partlen) {
char* p;
unsigned int restlen=len-partlen+1;
void* oldp=mem;
while ( (p=(char*)memchr(oldp, ((char*)part)[0], restlen))!=NULL) {
  //located first char, try to match the rest:
  p++;
  if (memcmp(p, &((char*)part)[1], partlen-1)==0) return p-1;
  //no string match, prepare next iteration
  restlen-=(p-(char*)oldp);
  oldp=p;
  }//while
return NULL;
}

//rindex function is missing on some platforms ?
char* rstrchr(char* str, char ch) {  /* returns a pointer to the rightmost
  occurence of ch in str  */
 char *p;
 if (str==NULL) return NULL;
 p=str+strlen(str)-1;
 while (p>=str) {
    if (*p==ch) return p;
    p--;
    }
 return NULL;
 }

/* DOS/UNIX safer (?) fgets : to read a text line from a binary file and
  update the file position accordingly */
char* fgetline(char *buf, int n, FILE *stream, long& f_pos) {
  //reads a char at a time until \n and/or \r are encountered
  int i=0;
  while (i<n-1) {
    int c=getc(stream);
    if (c=='\n' || c=='\r') {
      buf[i]='\0';
      if (c=='\r') {
        if ((c=getc(stream))!='\n')
                 ungetc(c,stream);
            else f_pos++;
        }
      f_pos++;
      return buf;
      }
    if (c==EOF) {
       if (i==0) return NULL;
          else {
             buf[i]='\0';
             return buf;
             }
       }
      else f_pos++;
    buf[i]=(char)c;
    i++;
    }
  buf[n-1]='\0';
  return buf;
  }


//strchr but with a set of chars instead of only one
char* strchrs(char* s, const char* chrs) {
  if (s==NULL || chrs==NULL || *chrs=='\0' || *s=='\0')
         return NULL;
  unsigned int l=strlen(s);
  unsigned int r=strcspn(s, chrs);
  if (r==l) return NULL;
  return (s+r);
}

char* rstrfind(char* str, char* substr) { 
/* like rindex() for a string */
 int l,i;
 if (str==NULL || *str=='\0') return NULL;
 if (substr==NULL || *substr=='\0') return NULL;
 char* p=str+strlen(str)-strlen(substr); 
   //rightmost position that could match
 l=strlen(substr);  
 while (p>=str) {
    for (i=0; i<l && *(p+i) == *(substr+i); i++);
    if (i==l) return p; //found!
    p--;
    }
 return NULL;
}


char* rstrstr(char* rstart, char *lend, char* substr) {  /*like strstr, but starts searching
 from right end, going up to lend and returns a pointer to the last (right) 
 matching character in str */
 char *p;
 int l,i;
 l=strlen(substr);
 p=rstart-l+1;
 while (p>=lend) {
    for (i=0;i<l;i++) if (*(p+i) != *(substr+i)) break;
    if (i==l) return p+l-1;
    p--;
    }
 return NULL;
 }


//hash function used for strings in GHash
int strhash(const char* str){
  register int h=0;
  register int g;
  while (*str) {
    h=(h<<4)+*str++;
    g=h&0xF0000000;
    if(g) h^=g>>24;
    h&=0x0fffffff;
    }
  GASSERT(h<=0x0fffffff);
  return h;
  }

