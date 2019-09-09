//  file motif.h

#ifndef MOTIFH
#define MOTIFH

#include <ctype.h>
#include <string.h>
#include "GFastaFile.h"

//#define MAXINT 0xFFFFFFFF

template <class T>
void copyarray(T *dest, T *source, int len)
{
  for(int i=0;i<len;i++) dest[i]=source[i];
}
  
enum seqType {
  nucl,
  aac
};

struct CumWeight {
  double cw;   // cumulative weight
  int pos;  // position in string
};

struct FreqPos {
  double freq;
  int pos;
};

class Background {
 public:
  FILE *outf;
  double eps;
  int sumc; // sum of background residues counts 
  int maxlen; // length of longest sequence in the data
  double *d; // dirichlet counts
  double D; // sum of dirichlet counts;
  seqType SType;
  int *c; // background residues counts 
  int **p; // background first-order probabilities
  int ***p2; // background second-order probabilities
  int ****p3; // background third-order probabilities
  int *****p4; // background 4th-order probabilities
  int ******p5; // background 5th-order probabilities
  int noseq; //no. of sequences
  int size; // no of residues in the alphabet
  Background(GFastaFile *F,FILE *outfile,const seqType t,double e=0.0000001);
  Background &initData(GFastaFile *F,const seqType t);
  double freq(int n);
  void printFreq();
  ~Background();
};


struct ProbMotif {
  int *pos; // positions where there might be a motif;
  int no; // no. of positions in a sequence of a probable motif
  int dist; // distance from given motif
};

class Motif {
  GFastaFile *F;
  Background *B;
  ProbMotif *PM;
  double ***Q;
  double **q;
  double sumq;
  double sumloglength;
  int **motiffreq; // keeps motifs probability
  FreqPos *motifprob;
  double sumlogweights;
  int defmotif;
  int **MCount; // includes background info
  int sumMCount0; // this is the sum of background frequencies
  double **MFreq;
  int len; // motif len
  void initMCount();
  int check(int pos, int i);
  void updateMCount(int pos,char *seq, seqType t, int dir);
  void updateBackgr(char *seq, int seqlen,seqType t, int dir);
  void updateMFreq(int no);
  void MatchSeqPatt(char *seq,const char *pattern, int i, int seqlen);
  int MatchSeqPattExact(char *seq, const char *pattern, int i, int seqlen,int ignore);
  void compWeights(double *A,CumWeight *Acum,char * fseq, int seqlen, int l);
  void compWeights(double *A,CumWeight *Acum,char * fseq, int seqlen);
  void complogweights(double *A, double sumA, int elemno);
  void createtwofile(const char *testfile,char **fkseq, int *posfreqlen);
  void createrandfile(char **fkseq, int *posfreqlen,int zparam,int mdeg);
  double scoreF(char *motif);
  int maxWeight(char *fseq,int seqlen);
  int sample(CumWeight *Acum, int Alen);
  char randbase(int pb,int pb2,int pb3,int pb4,int pb5,int *res_c, int yres_c,int zparam);
  int nearoptim(CumWeight *Acum, int Alen, double perc);
  void GibbsSample(int SignfRunNo, char **fkseq,int print);
  void testsignif(int det,char **fkseq, int *posfreqlen,int print,int exact);
  int computemaxPMno();
  void computeQq();
  void computeQq(char **fkseq);
  double computeTR(char *seq, int pos);
  double LLC(int print);
  double computeBR(char *seq, int pos);
  void removeoneword();
  void removeoneword(char **fkseq);
 public:
  int *Align;
  double findMDet(int print);
  void print_consensus(char * consens);
  double optimize(double prevMax, double previnfo, int wpattern);
  double findMotif(int iter_no, int max_loop,int inlocmax,int init, int det);
  void initRandAlign();
  void predUpdate();
  double AlignProb();
  double InfoPar(double F);
  void computeMotif(int *MaxAlign);
  void printMotif(int printmax,int SignifNo,int print=1);
  Motif(const char *infile,FILE *outf,seqType t, int motiflen);
  Motif(const char *infile,FILE *outf,seqType t,const char *pattern);
  Motif(const char *infile,FILE *outf,seqType t, const char *matrixfile, const char *pattern, int motiflen, int iter_no, int max_loop,int inlocmax, int det);
  void runforsignif(int SignfRunNo, int print,int det, const char *pattern,int zparam,int mdeg);
  void detrunforsignif(int print, const char *pattern,int zparam,int mdeg);
  void twofilesignif(int det,const char *testfile, int SignfRunNo,int print, const char *pattern,int mdeg,int exact);
  int getNoseq();
  void printFreq();
  double computeLLC(const char *pattern,int print);
  ~Motif();
};

int basetoint(char c,seqType t);
char inttobase(int b,seqType t);
int strdist(char *s1,char *s2, int len);
int comp(const void *a, const void *b);
int compf(const void *a, const void *b);
double test_student(FILE *outf,double mean, double tmean, double sd, int n);
double Normal(double Z);
double test_wilcoxon(FILE *outf,double W, int n, int na, int nb);
void apply_test_wilcoxon_pair(FILE *outf,int diffno, FreqPos *diff, int print=1);
double WilcoxPair(double Winput, unsigned long int N);
void apply_test_student(FILE *outf,int maxno,FreqPos *maxf);
void computeranks(double *rank,int maxno,FreqPos *maxf);
int makeargv(char *s, char *delimiters, char ***argvp);

#endif
