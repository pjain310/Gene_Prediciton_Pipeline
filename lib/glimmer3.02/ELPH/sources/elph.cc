#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "GArgs.h"
#include "GString.h"
#include "GFastaFile.h"
#include "motif.h"

#define usage "Motif finder program ELPH (Estimated Locations of Pattern Hits)\n\
Usage:\n\
  elph <multi-fasta_file> [options] OR\n\
  elph <multi-fasta_file-1> <multi-fasta_file-2> [-t <matrix>]\n\
  Options:\n\
  -h             : print this help\n\
  -b             : just brief results; don't print motif\n\
  -d             : find significance deterministically\n\
  -v             : find motif deterministically\n\
  -o <out_file>  : write output in <out_file> instead of stdout\n\
  -m <motif>     : use the given pattern <motif> to compute its best fit matrix\n\
                   to the data\n\
  -a             : use aac residues; default = nucleotides residues\n\
  -s <seed>      : sets the seed for the random generation\n\
  -p n           : n = no of iterations before deciding local maximum; plateau\n\
                   period variable\n\
  -x             : print maximum positions within sequences\n\
  -g             : find significance of motif\n\
  -t <matrix>    : test if there is significant difference between the two \n\
                   input files for a given motif matrix; <matrix> is the file\n\
                   containing the motif matrix\n\
  -l             : compute Least Likely Consensus (LLC) for given motif \n\
  -r             : in conjunction with -m option: motif is not necessarily in\n\
                   the closest edit distance from input motif\n\
  -e             : only when an additional file is used to test the significance\n\
                   of the motif: find only the motifs that exactly match the\n\
                   input pattern (-m or -t options) \n\
  -n [0..5]      : degree of Markov chain used to generate the random file \n\
                   used to test the significance of the motif\n\
                   default = 2\n\
  LEN=n          : n = length of motif\n\
  ITERNO=n       : n = no of iterations to compute the global maximum;\n\
                   default = 10\n\
  MAXLOOP=n      : n = no of iterations to compute the local maximum;\n\
                   default = 500\n\
  SGFNO=n        : n = no of iterations to compute significance of motif;\n\
                   default = 1000\n\
 "

/* I am eliminating the following option from help file because it has no effect because I am never shuffling the file
  -z             : when using the makov chain file generation, don't exhaust\n\
                   the residues\n\
*/


// global variables: 
int ITER_NO=10;
int MAX_LOOP=500;
int printmax=0;
int SignifNo=1000;
int runsignif=0;
int defLLC=0;
int closest=1;
int exact=0;
int mdeg=2;

GString pattern, infile, testfile, matrixfile;
FILE *outf;
int motiflen=0;
int seed;
int inlocmax=20;
int gdet=0;
int mdet=0;
int print=1;
int zparam=0;

seqType Process_Options(GArgs *args);

int main(int argc, char *argv[])
{
  seqType t;

  GArgs args(argc, argv, "ho:abreglvdzxt:p:s:m:n:LEN=ITERNO=MAXLOOP=SGFNO=");
  
  // == Process arguments.

  int e;
  if ((e=args.isError())>0) {
    fprintf(stderr, "Error at argument %d, %s\n",e, argv[e]);
    exit(1);
  }
  t=Process_Options(&args);
  // == Done argument processing.

  Motif *M;
  if(!testfile.is_empty()) { // if testfile is defined then only compute significance between the two files

    M = new Motif(infile,outf,t,matrixfile,pattern,motiflen,ITER_NO,MAX_LOOP,inlocmax,mdet);    

    if(!pattern.is_empty()) M->twofilesignif(gdet,testfile,SignifNo,print,pattern,mdeg,exact);
    else M->twofilesignif(gdet,testfile,SignifNo,print,NULL,mdeg,exact);

    if(print) M->printMotif(printmax,SignifNo,print);

  }
  else {

    if(!pattern.is_empty()) { // motif was defined => find the best possible matches to the
                            // given motif

      M = new Motif(infile,outf,t,pattern);
      if(defLLC) { 
	double llc=M->computeLLC(pattern,print);
	fprintf(outf,"LLC = %f\n",llc);
      }
      
    }
    else {
      /* if(mdet) searchformotifLLC(); */

      M = new Motif(infile,outf,t,motiflen);
      closest=0;
    }

    double globAlignProb;

    globAlignProb=M->findMotif(ITER_NO,MAX_LOOP,inlocmax,1,mdet);

    double info;
    /*info=M->InfoPar(globAlignProb);
      fprintf(outf,"MAP for motif: %.3f InfoPar=%.3f\n\n",globAlignProb,info);
      M->printMotif();*/
    
    // optimizing
    fprintf(stderr,"Optimizing...\n");
    globAlignProb=M->optimize(globAlignProb,info,closest);
    fprintf(outf,"\n\n**********************\n\nMotif after optimizing\n");
    fprintf(outf,"MAP for motif: %.3f InfoPar=%.3f\n\n",globAlignProb,M->InfoPar(globAlignProb));
    
    if(runsignif) {
      M->runforsignif(SignifNo,print,gdet,pattern,zparam,mdeg);
    }

    M->printMotif(printmax,SignifNo,print);
  }


}

seqType Process_Options(GArgs* args)
{
    
  if (args->startNonOpt()) { //parse the non-options arguments 
                          //(usually filenames)
        infile=args->nextNonOpt();
  }

  if (infile.is_empty() || args->getOpt('h')!=NULL) 
    GError("%s",usage); // the empty test is optional you can ignore it if you accept stdin input

  testfile=args->nextNonOpt();
  
  GString outfile=args->getOpt('o');
  if (!outfile.is_empty()) {
    outf=fopen(outfile, "w");
    if (outf==NULL)
      GError("Cannot open file %s for writing!\n",outfile.chars());
  }
  else outf = stdout;

  matrixfile=args->getOpt('t');

  GString param;
  
  pattern=args->getOpt('m');
  if(pattern.is_empty()) {
    param=args->getOpt("LEN");
    if(param.is_empty() && matrixfile==NULL) {
      fprintf(stderr,"You must introduce motif length : "); fflush(stderr);
      scanf("%d",&motiflen);
    }
    else motiflen=param.asInt();
  }

  param=args->getOpt("ITERNO");
  if(!param.is_empty()) ITER_NO=param.asInt();

  param=args->getOpt("SGFNO");
  if(!param.is_empty()) SignifNo=param.asInt();

  param=args->getOpt("MAXLOOP");
  if(!param.is_empty()) MAX_LOOP=param.asInt();

  param=args->getOpt('s');
  if(!param.is_empty()) seed=param.asInt();
  else seed=1;
  srand(seed);

  param=args->getOpt('p');
  if(!param.is_empty()) inlocmax=param.asInt();

  param=args->getOpt('n');
  if(!param.is_empty()) mdeg=param.asInt();
  if(mdeg<0 || mdeg>5) GError("Invalid markov degree: value should be between 0 and 5!");

  if(args->getOpt('e')!=NULL) {exact=1;gdet=1;}

  if(args->getOpt('x')!=NULL) printmax=1;

  if(args->getOpt('g')!=NULL) runsignif=1;

  if(args->getOpt('r')!=NULL) closest=0;

  if(args->getOpt('d')!=NULL) gdet=1;

  if(args->getOpt('z')!=NULL) zparam=1;

  if(args->getOpt('v')!=NULL) mdet=1;

  if(args->getOpt('b')!=NULL) print=0;

  if(args->getOpt('l')!=NULL) defLLC=1;

  seqType t;
  if(args->getOpt('a')!=NULL) t=aac; else t=nucl;
  
  return(t);

}

/*void searchformotifLLC() {

  Motif *M;

  double llcmax=-HUGE_VAL;
  GString seed;	
  for(int i1=0;i1<4;i1++)
    for(int i2=0;i2<4;i2++)
      for(int i3=0;i3<4;i3++)
	for(int i4=0;i4<4;i4++) {
	  pattern="";
	  pattern.appendfmt("%c%c%c%c",inttobase(i1,t),inttobase(i2,t),inttobase(i3,t),inttobase(i4,t));
	  fprintf(stderr,"%s\n",(const char*)pattern);
	  M=new Motif(infile,t,pattern);
	  double llc=M->computeLLC(0);
	  if(llc>llcmax) {
	    llcmax=llc;
	    seed=pattern;
	  }
	  delete M;
	  printf("\n");
	}
  fprintf(stderr,"length 4: seed=%s\n",(const char*)seed);
  printf("length 4: seed=%s\n",(const char*)seed);
  seed="ttgg";
  int k=5;
  while(1) {
    GString tseed;
    llcmax=-HUGE_VAL;
    for(int i=0;i<4;i++) {
      pattern.format("%c%s",inttobase(i,t),(const char*)seed);
      fprintf(stderr,"%s\n",(const char*)pattern);
      M=new Motif(infile,t,pattern);
      double llc=M->computeLLC(0);
      if(llc>llcmax) {
	llcmax=llc;
	tseed=pattern;
      }
      delete M;
      printf("\n");
      pattern.format("%s%c",(const char*)seed,inttobase(i,t));
      fprintf(stderr,"%s\n",(const char*)pattern);
      M=new Motif(infile,t,pattern);
      llc=M->computeLLC(0);
      if(llc>llcmax) {
	llcmax=llc;
	tseed=pattern;
      }
      delete M;
      printf("\n");
    }
    seed=tseed;
    fprintf(stderr,"length %d: seed=%s\n",k,(const char*)seed);
    printf("length %d: seed=%s\n",k,(const char*)seed);
    k++;
  }
}*/

