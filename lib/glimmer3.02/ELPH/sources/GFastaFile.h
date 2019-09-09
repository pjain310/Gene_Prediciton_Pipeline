#ifndef GFILE_H
#define GFILE_H
#include "GBase.h"
#include <stdio.h>

#define CAPINC 64

class FastaSeq {  /* fasta record storage */
   public:
     int  id_len; /* allocated size of the sequence name string*/
     char    *id; /* id only, up to first space */
     char *descr; /* any comment on the defline, after the first space */
     int   d_len; /* allocated size of the description */
     int descrlen; /* real length of the description */
     //-------actual sequence :
     int s_len; /* allocated length of the sequence string */
     int   len; /* the actual string length of seq */
     char* seq; /* the sequence buffer itself */
     //----
     FastaSeq() {
       GMALLOC(id, CAPINC);
       id_len=CAPINC;
       id[0]='\0';
       GMALLOC(descr, CAPINC);
       descr[0]='\0';
       d_len=CAPINC;
       GMALLOC(seq, CAPINC<<1);
       seq[0]='\0';
       len=0;
       descrlen=0;
       s_len=CAPINC<<1;
       }
     ~FastaSeq() {
       GFREE(id);id_len=0;
       GFREE(descr);d_len=0;descrlen=0;
       GFREE(seq);s_len=0;len=0;
       }
     int getDescrLen() { return descrlen; }
     const char* getDescr()  { return (const char*) descr; }
     const char* getSeqName()  { return (const char*) id; }
     int getSeqNameLen()  { return s_len; }
}; 


typedef int charFunc(char c, int pos, FastaSeq* fseq); //char processing function
/* passes: 
   current sequence character (c) 
   position of the given character within the sequence (pos)
   FastaSeq pointer (useful for retrieving sequence defline info) (fseq)
  the return value is not used yet 
*/


//(for reading/writing variable length records, etc.)
enum fileMode {
 fmRead,
 fmWrite
 };

class GFastaFile {
  char* fname;
  FILE* fh;
  long fpos;//current file position
  fileMode fmode;

  int record_read_pos; //the input stream offset of the current record to be read
  int read_pos;        //the input stream offset of the current byte to be read           

 protected:
  void bad_fastafmt() {
      GError("Error parsing file '%s'. Not a Fasta file?\n", fname);
      }
  void check_eof(int c) {
      if (c == EOF) bad_fastafmt();
      }
 public:
  GFastaFile(const char* filename, fileMode filemode=fmRead) {
     fh=NULL;
     read_pos=0;
     record_read_pos=0;
     fmode=filemode;
     const char *mode=(filemode==fmRead) ? "rb" : "wb";
     if (filename == NULL || filename[0]=='\0') {
           fh = (filemode == fmRead) ? stdin : stdout;
           fname=NULL;
           }
         else {
           if ((fh = fopen(filename, mode)) == NULL)
               GError("Cannot open file '%s'!", filename);
           fname=Gstrdup(filename);
           }
       /*    
       GCALLOC(curseqid, CAPINC);
       curseqidlen=CAPINC;
       GCALLOC(curdescr, CAPINC);
       curdescrlen=CAPINC;*/
     }

   //attach a GFastaFile object to an already open handle  
   GFastaFile(FILE* fhandle, fileMode filemode=fmRead, const char* filename=NULL) {
     fh=fhandle;
     read_pos=ftell(fh);
     fmode=filemode;
     record_read_pos=read_pos;
     if (filename == NULL || filename[0]=='\0') {
           fname=NULL;
           }
         else
           fname=Gstrdup(filename);           
     }
     
     
   void reset() {
    if (fh!=NULL && fh!=stdout && fh!=stdin) {
       fseek(fh,0L, SEEK_SET);
       read_pos=0;
       record_read_pos=0;
       }
     else GError("Cannot use GFastaFile::reset() on stdin, stdout or NULL handles.\n");
    }
    
   void seek(int pos) {
    if (fh!=NULL && fh!=stdout && fh!=stdin)
      fseek(fh, pos, SEEK_SET);
     else GError("Cannot use GFastaFile::reset() on stdin, stdout or NULL handles.\n");
    }
  ~GFastaFile() {
    if (fh!=NULL && fh!=stdout && fh!=stdin) fclose(fh);
    fh=NULL;
    GFREE(fname);
    /*GFREE(curseqid);
    GFREE(curdescr);*/
    }

   int getReadPos() { return read_pos; } /* returns current read position in the 
              input stream (can be used within callback) */
   int ReadSeqPos() {return record_read_pos; } /* returns the input stream offset of the last fasta 
                                                record processed by getFastaSeq*/
   
   //reads the Fasta sequence header 
   /* the first character must be '>' for each call
    seq must be a pointer to a initialized FastaSeq structure
    if seq is NULL, the sequence is not actually read,
     but just skipped and the file pointer set accordingly, while
     the returned "pointer" will not a valid one, but just NULL (or not NULL if end of file was encountered)
   if callbackFn is NULL, the sequence is read entirely in memory in a FastaSeq.seq field
      otherwise only the defline is parsed into FastaSeq::id and FastaSeq::descr but actual 
       sequence letters is are passed one by one to the callback function 
      and the actual sequence is never stored in memory (unless the callback does it)
   */
   FastaSeq *getFastaSeq(bool& is_last, FastaSeq* seq, charFunc* callbackFn = NULL ) {
      int c, len;
      int* buflen; 
      char** buf;
      int before;
      record_read_pos=read_pos;
      c = getc(fh); read_pos++;
      if (c==EOF) return NULL;
      if (c != '>')
            bad_fastafmt();
      len = 0; //chars accumulated so far
      seq->seq[0]='\0';
      seq->descr[0]='\0';
      seq->id[0]='\0';
      seq->descrlen=0;
      seq->len=0;
      // -------- read the defline first
      if (seq==NULL) { /* navigate only! don't read/parse anything but the record delimiter*/
          before=1;
          while ((c = getc(fh)) != EOF && c != '\n') read_pos++; /* skip defline */
          check_eof(c); /* it's wrong to have eof here! */
          read_pos++; //to account for the '\n' read
          /*----- read the sequence now: */
          before=1; /* "newline before" flag */
          while ((c = getc(fh)) != EOF && c != '>') {
                read_pos++;
                before = (c=='\n')?1:0;
                }
          } /* fasta fmt navigation only, no seq storage */
       else { // sequence storage: 
          buflen=&seq->id_len;
          buf=&seq->id;
          before=1;
          while ((c = getc(fh)) != EOF && c != '\n') {
              read_pos++; 
              if (len >= *buflen-1) {
                      GREALLOC(*buf, *buflen + CAPINC);
                      *buflen+=CAPINC;
                      }
              if (before && (c<=32)) {
                 /* space encountered => seq_name finished */
                 before=0;
                 (*buf)[len]='\0';
                 buf=&seq->descr;
                 buflen=&seq->d_len;
                 len=0;
                 if (c!=1)  /* special case, nrdb concatenation */
                   continue; // skip this space
                 }
              (*buf)[len]=c;
              len++;
              }
          (*buf)[len]='\0'; /* terminate the comment string */
          if (buf==&seq->descr)
              seq->descrlen=len;
          check_eof(c); /* it's wrong to have eof here */
          read_pos++; // to account for the last end of line read
          /*----- read the actual sequence now: */
          len=0;
          before=1; //newline before indicator          
          if (callbackFn==NULL) { //load the whole sequence in FastaSeq 
             while ((c = getc(fh)) != EOF && c != '>') {
                read_pos++;
                //if (isspace(c) || c<31) 
                if (c<=32) {
                       before = (c=='\n')?1:0;
                       continue; /* skip spaces */
                       }
                if (len >= seq->s_len-1) {
                      GREALLOC(seq->seq, seq->s_len + CAPINC);
                      seq->s_len+=CAPINC;
                      }
                seq->seq[len] = c;
                before=0;
                len++;
                }
             seq->seq[len] = '\0';
             seq->len=len;
             } /* sequence storage */
          else { //use the callback for each letter, do not store the whole sequence in FastaSeq
             while ((c = getc(fh)) != EOF && c != '>') {
                read_pos++;
                //if (isspace(c) || c<31) 
                  if (c<=32) {
                       before = (c=='\n')?1:0;
                       continue; /* skip spaces */
                       }
                (*callbackFn)(c, len, seq); //call the user function for each letter
                before=0;       
                len++;
                }
             seq->len=len;
             } /* callback sequence reading (no storage)*/
        } /* sequence parsing */     
      if (c=='>') {
         if (!before) bad_fastafmt(); /* '>' must only be at start of line,
                                       never within the sequence ! */
         is_last=false; /* FALSE - not the last one */
         ungetc(c, fh); 
         }
        else is_last=true; /* TRUE - eof() here */
      return ((seq==NULL) ? (FastaSeq*)fh : seq); //alwayws return non NULL here!
  } //getFastaSeq



   //simplified call to ignore the is_last flag
  FastaSeq *getFastaSeq(FastaSeq* seq, charFunc* callbackFn = NULL) {
     bool b;
     return getFastaSeq(b, seq, callbackFn);
     }

   //only for writing
   void putFastaSeq(FastaSeq *fa, const int linelen=60) {
      writeFasta(fh, fa->id, fa->descr, fa->seq, linelen);      
      }

   static void writeFasta(FILE *fh, char* seqid, char* descr, char* seq, const int linelen=60) {
      char *s;
      int i, ilen;
      s = (seqid == NULL) ? (char*)"ANONYMOUS" : seqid;
      if (*s != '>') putc('>', fh);
      fwrite(s, 1, strlen(s), fh);
      i=(descr==NULL)? 0 : strlen(descr);
      if (i>0) {
        putc(' ',fh);
        fwrite(descr, 1, i, fh);
        }
      ilen = linelen;
      if (ilen>0) {
          int len=strlen(seq);

          for (i=0, s = seq; i < len; i++, s++, ilen++) {
                if (ilen == linelen) {
                     putc('\n', fh);
                     ilen = 0;
                     }
                putc(*s, fh);
          }
        putc('\n', fh);  
         }        
       else { //no line length limit
         fprintf(fh, "\n%s\n", seq);
         }
      fflush(fh);      
      }
};


#endif
