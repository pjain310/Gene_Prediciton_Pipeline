//--------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include "motif.h"
#include "GFastaFile.h"
//--------------------------------------------------------------------------

// t-distribution table
// lines are df: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 inf
// columns are: 0.40 0.25 0.10 0.05 0.025 0.01 0.005 0.0005
//
static  double t_value[31][8]=
  { {0.324920, 1.000000, 3.077684, 6.313752, 12.70620, 31.82052, 63.65674, 636.6192} ,
    {0.288675, 0.816497, 1.885618, 2.919986,  4.30265,  6.96456,  9.92484,  31.5991},
    {0.276671, 0.764892, 1.637744, 2.353363,  3.18245,  4.54070,  5.84091,  12.9240},
    {0.270722, 0.740697, 1.533206, 2.131847,  2.77645,  3.74695,  4.60409,   8.6103},
    {0.267181, 0.726687, 1.475884, 2.015048,  2.57058,  3.36493,  4.03214,   6.8688},
    {0.264835, 0.717558, 1.439756, 1.943180,  2.44691,  3.14267,  3.70743,   5.9588},
    {0.263167, 0.711142, 1.414924, 1.894579,  2.36462,  2.99795,  3.49948,   5.4079},
    {0.261921, 0.706387, 1.396815, 1.859548,  2.30600,  2.89646,  3.35539,   5.0413},
    {0.260955, 0.702722, 1.383029, 1.833113,  2.26216,  2.82144,  3.24984,   4.7809},
    {0.260185, 0.699812, 1.372184, 1.812461,  2.22814,  2.76377,  3.16927,   4.5869},
    {0.259556, 0.697445, 1.363430, 1.795885,  2.20099,  2.71808,  3.10581,   4.4370},
    {0.259033, 0.695483, 1.356217, 1.782288,  2.17881,  2.68100,  3.05454,   4.3178},
    {0.258591, 0.693829, 1.350171, 1.770933,  2.16037,  2.65031,  3.01228,   4.2208},
    {0.258213, 0.692417, 1.345030, 1.761310,  2.14479,  2.62449,  2.97684,   4.1405},
    {0.257885, 0.691197, 1.340606, 1.753050,  2.13145,  2.60248,  2.94671,   4.0728},
    {0.257599, 0.690132, 1.336757, 1.745884,  2.11991,  2.58349,  2.92078,   4.0150},
    {0.257347, 0.689195, 1.333379, 1.739607,  2.10982,  2.56693,  2.89823,   3.9651},
    {0.257123, 0.688364, 1.330391, 1.734064,  2.10092,  2.55238,  2.87844,   3.9216},
    {0.256923, 0.687621, 1.327728, 1.729133,  2.09302,  2.53948,  2.86093,   3.8834},
    {0.256743, 0.686954, 1.325341, 1.724718,  2.08596,  2.52798,  2.84534,   3.8495},
    {0.256580, 0.686352, 1.323188, 1.720743,  2.07961,  2.51765,  2.83136,   3.8193},
    {0.256432, 0.685805, 1.321237, 1.717144,  2.07387,  2.50832,  2.81876,   3.7921},
    {0.256297, 0.685306, 1.319460, 1.713872,  2.06866,  2.49987,  2.80734,   3.7676},
    {0.256173, 0.684850, 1.317836, 1.710882,  2.06390,  2.49216,  2.79694,   3.7454},
    {0.256060, 0.684430, 1.316345, 1.708141,  2.05954,  2.48511,  2.78744,   3.7251},
    {0.255955, 0.684043, 1.314972, 1.705618,  2.05553,  2.47863,  2.77871,   3.7066},
    {0.255858, 0.683685, 1.313703, 1.703288,  2.05183,  2.47266,  2.77068,   3.6896},
    {0.255768, 0.683353, 1.312527, 1.701131,  2.04841,  2.46714,  2.76326,   3.6739},
    {0.255684, 0.683044, 1.311434, 1.699127,  2.04523,  2.46202,  2.75639,   3.6594},
    {0.255605, 0.682756, 1.310415, 1.697261,  2.04227,  2.45726,  2.75000,   3.6460},
    {0.253347, 0.674490, 1.281552, 1.644854,  1.95996,  2.32635,  2.57583,   3.2905} };


//
// Background::Background
// Background constructor for a given fasta file and sequence type (aac or nt)
//
Background::Background(GFastaFile *F,FILE *outfile,const seqType t,double e)
{
  outf=outfile;
  eps=e;
  initData(F,t);
}

//
// Background::freq
// computes the frequence of a residue
//
double Background::freq(int n)
{
  return (c[n]+d[n])/(sumc+D);
}

//
// Background::printFreq
// prints the background frequencies
//
void Background::printFreq()
{
  fprintf(outf,"Background probability model:\n");
  for(int n=0;n<size;n++) 
    fprintf(outf,"   %c   ",inttobase(n,SType));
  fprintf(outf,"\n");
  for(int n=0;n<size;n++) 
    fprintf(outf,"  %5.2f",freq(n));
  fprintf(outf,"\n");
  //  printf("sumc=%d D=%f\n",sumc,D);
  fprintf(outf,"\nBackground counts:\n");
  for(int n=0;n<size;n++) 
    fprintf(outf,"%c: %d\n",inttobase(n,SType),c[n]);
  fprintf(outf,"\n\n");
}

//
// Background::initData
// computes all Background frequencies
//
Background &Background::initData(GFastaFile *F, const seqType t)
{
  SType=t;
  size=(SType==nucl) ? 4 : 22;
  c= new int[size];
  d= new double[size];
  p= new int*[size];
  p2 = new int**[size];
  p3 = new int***[size];
  p4 = new int****[size];
  p5 = new int*****[size];
  noseq=0;
  for(int i1=0;i1<size;i1++) {
    c[i1]=0;
    p[i1] = new int[size];
    p2[i1]=new int*[size];
    p3[i1]=new int**[size];
    p4[i1]=new int***[size];
    p5[i1]=new int****[size];
    for(int i2=0;i2<size;i2++) {
      p[i1][i2]=0;
      p2[i1][i2]=new int[size];
      p3[i1][i2]=new int*[size];
      p4[i1][i2]=new int**[size];
      p5[i1][i2]=new int***[size];
      for(int i3=0;i3<size;i3++) {
	p2[i1][i2][i3]=0;
	p3[i1][i2][i3]=new int[size];
	p4[i1][i2][i3]=new int*[size];
	p5[i1][i2][i3]=new int**[size];
	for(int i4=0;i4<size;i4++) {
	  p3[i1][i2][i3][i4]=0;
	  p4[i1][i2][i3][i4]=new int[size];
	  p5[i1][i2][i3][i4]=new int*[size];
	  for(int i5=0;i5<size;i5++) {
	    p4[i1][i2][i3][i4][i5]=0;
	    p5[i1][i2][i3][i4][i5]=new int[size];
	    for(int i6=0;i6<size;i6++) {
	      p5[i1][i2][i3][i4][i5][i6]=0;
	    }
	  }
	}
      }
    }
  }
  sumc=0;
  F->reset();
  FastaSeq fseq;
  maxlen=0;
  while(F->getFastaSeq(&fseq)) {

    //fprintf(stderr,"Read seq of len %d: %s\n",fseq.len,fseq.id);fflush(stdout);

    if(fseq.len>maxlen) maxlen=fseq.len;
    noseq++;
    int pb=basetoint(fseq.seq[0],t);
    if(pb!=-1) c[pb]++;
    int b2=-1;
    int b3=-1;
    int b4=-1;
    int b5=-1;

    for(int i=1;i<fseq.len;i++) {
      int b=basetoint(fseq.seq[i],t);
      if(b!=-1) {
	c[b]++;

	if(pb!=-1) {
	  p[b][pb]++;

	  if(b2!=-1) {
	    p2[b][pb][b2]++;

	    if(b3!=-1) {
	      p3[b][pb][b2][b3]++;

	      if(b4!=-1){
		p4[b][pb][b2][b3][b4]++;

		if(b5 !=-1) {
		  p5[b][pb][b2][b3][b4][b5]++;

		}
	      }
	    }
	  }
	}
      }
      b5=b4;
      b4=b3;
      b3=b2;
      b2=pb;
      pb=b;
    }
  }

  D=0;sumc=0;
  for(int i=0;i<size;i++) {
    d[i]=eps*(c[i]+1);
    D+=d[i];
    sumc+=c[i];
  }

  return *this;
}

// 
// Background::~Background
// Background destructor
//
Background::~Background() {
  delete [] c;
  delete [] d;
  for(int i=0;i<size;i++) {
    delete [] p[i];
    for( int j=0;j<size;j++) {
      delete [] p2[i][j];
      for(int k=0;k<size;k++) {
	delete [] p3[i][j][k];
	for(int u=0;u<size;u++) {
	  delete [] p4[i][j][k][u];
	  for(int v=0;v<size;v++) {
	    delete [] p5[i][j][k][u][v];
	  }
	  delete [] p5[i][j][k][u];
	}
	delete [] p4[i][j][k];
	delete [] p5[i][j][k];
      }
      delete [] p3[i][j];
      delete [] p4[i][j];
      delete [] p5[i][j];
    }
    delete [] p2[i];
    delete [] p3[i];
    delete [] p4[i];
    delete [] p5[i];
  }
  delete [] p;
  delete [] p2;
  delete [] p3;
  delete [] p4;
  delete [] p5;
}

//
// Motif::Motif
// motif constructor for a given file, sequence type and motif length
//
Motif::Motif(const char *infile, FILE *outf,seqType t, int motiflen) 
{
  F = new GFastaFile(infile);
  B = new Background(F,outf,t);

  defmotif=0;
  len=motiflen;
  sumloglength=0;
  sumlogweights=0;

  Align=new int[B->noseq];
  MCount = new int*[len+1];                             // MCount=new int[len+1][size];
  MFreq = new double*[len+1];                             //MFreq=new double[len+1][size]; 

  for(int i=0;i<len+1;i++) { 
    MCount[i] = new int[B->size];
    MFreq[i] = new double[B->size];     
  }

  motiffreq=NULL;
  motifprob=NULL;

  q = new double *[len];
  Q = new double **[len];

  for(int i=0;i<len;i++) {
    q[i] = new double[B->size];
    Q[i] = new double *[B->size];
    for(int j=0;j<B->size;j++) Q[i][j] = new double[B->size];
  }

  F->reset();

  FastaSeq fseq;
  while(F->getFastaSeq(&fseq)) {

    if(fseq.len<len) GError("Motif length is greater then input sequence %s\n",fseq.id);
    sumloglength+=log(fseq.len-len+1);

  }

}

//
// Motif::randbase
// randomly generates a residue; if pb is specified then a 1order Markov chain is
// used to generate the next resifue; if pb is -1 then a 0order Mrkov chain is used
// if yres_c is 1 then the generation effect is equivalent to shuffling the 
// residues from the file on which the distributions are computed
//
char Motif::randbase(int pb,int pb2,int pb3,int pb4, int pb5,int *res_c,int yres_c,int zparam) 
{
  int csumc=0;
  int i;
  
  while(!csumc) {
    if(pb<0) {
      for(i=0;i<B->size;i++) if(res_c[i]) csumc+=B->c[i];
    }
    else
      if(pb2<0)  {  
	for(i=0;i<B->size;i++) if(res_c[i]) csumc+=B->p[i][pb];
	if(!csumc) pb=-1;
      }
      else 
	if(pb3<0) { 
	  for(i=0;i<B->size;i++) if(res_c[i]) csumc+=B->p2[i][pb][pb2];
	  if(!csumc) pb2=-1;
	}
	else 
	  if(pb4<0) { 
	    for(i=0;i<B->size;i++) if(res_c[i]) csumc+=B->p3[i][pb][pb2][pb3];
	      if(!csumc) pb3=-1;
	  }
	  else 
	    if(pb5<0) { 
	      for(i=0;i<B->size;i++) if(res_c[i]) csumc+=B->p4[i][pb][pb2][pb3][pb4];
	      if(!csumc) pb4=-1;
	    }
	    else {  
	      for(i=0;i<B->size;i++) if(res_c[i]) csumc+=B->p5[i][pb][pb2][pb3][pb4][pb5];
	      if(!csumc) pb5=-1;
	    }
  }
  //fprintf(stderr,"1.1.1 csumc=%d pb=%d\n",csumc,pb);
  
  int r = 1 + rand() % csumc;
  

  i=-1;
  csumc=0;
  while(r>csumc) {

    i++;
    if(!zparam || yres_c) while(res_c[i]==0) i++;
    if(pb<0) csumc+=B->c[i];
    else if(pb2<0) csumc+=B->p[i][pb];
    else if(pb3<0) csumc+=B->p2[i][pb][pb2];
    else if(pb4<0) csumc+=B->p3[i][pb][pb2][pb3];
    else if(pb5<0) csumc+=B->p4[i][pb][pb2][pb3][pb4];
    else csumc+=B->p5[i][pb][pb2][pb3][pb4][pb5];
  }


  if(yres_c) res_c[i]--;
  char ret=inttobase(i,B->SType);

  return(ret);
}

//
// Motif::createtwofile
// creates test sequences by appending testfile to the input file;
// at most the number of sequences in input file are appended from
// the testfile
//
void Motif::createtwofile(const char *testfile,char **fkseq, int *posfreqlen)
{
  int testlen;


  GFastaFile *T;
  T = new GFastaFile(testfile);

  F->reset();
  FastaSeq fseq1,fseq2;
  
  int i=0;
  while(F->getFastaSeq(&fseq1)) {
    
    testlen=fseq1.len;
    posfreqlen[i]=fseq1.len;

    if(T->getFastaSeq(&fseq2)) {
      testlen+=fseq2.len;
    }
    
    motiffreq[i] = new int[testlen];
    fkseq[i] = new char[testlen+1];
    
    for(int j=0; j<testlen; j++) {
      motiffreq[i][j]=0;
    }
    
    strcpy(fkseq[i],fseq1.seq);
    if(testlen>fseq1.len) strcat(fkseq[i],fseq2.seq);
    
    updateBackgr(fkseq[i]+fseq1.len,testlen-fseq1.len,B->SType,1);
    
    i++;
    

  }

  delete T;
}

// 
// Motif::createrandfile
// create test sequences by appending sequences of the same length as the
// ones in input file; the test sequences are randomly generated according
// to a 1order Markov chain that respects the distribution of the residues in
// the input file
//
void Motif::createrandfile(char **fkseq, int *posfreqlen,int zparam,int mdeg)
{
  int *res_c;
  int i;

  res_c= new int[B->size];
  for(i=0;i<B->size;i++) res_c[i]=B->c[i];
  
  F->reset();
  FastaSeq fseq;
  i=0;
  while(F->getFastaSeq(&fseq)) {

    motiffreq[i] = new int[2*fseq.len];
    fkseq[i] = new char[2*fseq.len+1];
    strcpy(fkseq[i],fseq.seq);

    fkseq[i][fseq.len]=randbase(-1,-1,-1,-1,-1,res_c,0,zparam);
    //fkseq[i][fseq.len]=randbase(-1,-1,-1,-1,-1,res_c,1,zparam); // shuffling sequences only

    int pb=basetoint(fkseq[i][fseq.len],B->SType);
    int pb2=-1;
    int pb3=-1;
    int pb4=-1;
    int pb5=-1;
    motiffreq[i][0]=0;
    motiffreq[i][fseq.len]=0;

    for(int j=1; j<fseq.len; j++) {

      switch (mdeg) {
      case 0: fkseq[i][j+fseq.len]=randbase(-1,-1,-1,-1,-1,res_c,0,zparam); break;
      case 1: fkseq[i][j+fseq.len]=randbase(pb,-1,-1,-1,-1,res_c,0,zparam); break; // first oder markov
      case 2: fkseq[i][j+fseq.len]=randbase(pb,pb2,-1,-1,-1,res_c,0,zparam); break; // second order markov
      case 3: fkseq[i][j+fseq.len]=randbase(pb,pb2,pb3,-1,-1,res_c,0,zparam); break; // third order markov
      case 4: fkseq[i][j+fseq.len]=randbase(pb,pb2,pb3,pb4,-1,res_c,0,zparam); break; // forth order markov
      case 5: fkseq[i][j+fseq.len]=randbase(pb,pb2,pb3,pb4,pb5,res_c,0,zparam); break; // 5th order order markov
      }
      //fkseq[i][j+fseq.len]=randbase(-1,-1,-1,-1,-1,res_c,1,zparam); // shuffling sequences; 0th order marko

      pb5=pb4;
      pb4=pb3;
      pb3=pb2;
      pb2=pb;
      pb=basetoint(fkseq[i][j+fseq.len],B->SType);
      motiffreq[i][j]=0;
      motiffreq[i][j+fseq.len]=0;
    }

    fkseq[i][2*fseq.len]='\0';

    posfreqlen[i]=fseq.len;

    updateBackgr(fkseq[i]+fseq.len,fseq.len,B->SType,1);

    i++;
  }

  delete [] res_c;
}



//
// Motif::runforsignif
// computes the significance of the motif by appending
// randomly generated sequences of the same length to the
// sequences in the input file;
// if det is specified then a deterministic computation of the 
// significance is used; otherwise the significance is computed
// by running Gibbs sampling and computing the frequencies of
// choosing the motif in the input file
//
void Motif::runforsignif(int SignfRunNo,int print, int det, const char *pattern,int zparam,int mdeg)  // give option at run time to specify these parameters
{
  char **fkseq;
  int *posfreqlen;

  posfreqlen= new int[B->noseq];

  fkseq = new char*[B->noseq];
  motiffreq = new int*[B->noseq];

  createrandfile(fkseq,posfreqlen,zparam,mdeg);
    
  if(det) { 
    char consensus[len+1];

    if(pattern==NULL || pattern[0]=='0') {
      PM=new ProbMotif[B->noseq];
      defmotif=1;

    }
    else strcpy(consensus,pattern);
    
    for(int i=0;i<B->noseq;i++) {    
      MatchSeqPatt(fkseq[i],consensus,i,strlen(fkseq[i]));
    }
    
    testsignif(1,fkseq,posfreqlen,print,0);
    
    if(pattern==NULL || pattern[0]=='0') defmotif=0;

  }
  else {

    // run Gibbs sampling SignfRunNo times and compute frequencies of choosing motifs
    GibbsSample(SignfRunNo,fkseq,print);
  
    // compute motif ranking and significance
    testsignif(0,fkseq,posfreqlen,print,0);

  }

  for(int i=0;i<B->noseq;i++) delete [] fkseq[i];
  delete [] fkseq;
    
  delete [] posfreqlen;

}

//
// apply_test_wilcoxon_pair
// applies the wilcoxon pair test on a vector of differences
//
void apply_test_wilcoxon_pair(FILE *outf,int diffno, FreqPos *diff,int print) 
{
  double *rank;

  GMALLOC(rank,diffno*sizeof(double));
  computeranks(rank,diffno,diff);

  //if(print) fprintf(outf,"\n\nWilcoxon Pairs Signed-Ranks Test:\n\n");
  fprintf(outf,"\n\nWilcoxon Pairs Signed-Ranks Test:\n\n");

  double mean=0;
  for(int i=0;i<diffno;i++) {
    mean+=rank[i];
  }
  mean/=diffno;

  //if(print) fprintf(outf,"Mean of ranks in sequences: %f\n",mean);
  fprintf(outf,"Mean of ranks in sequences: %f\n",mean);

  // Determine test parameters
  double Wp = 0;
  double Wm = 0;
  for(int i=0;i<diffno;i++) {
    if(rank[i]>0) Wp += rank[i];
    else Wm += fabs(rank[i]);
  }

  if(Wp<=Wm || diffno<=1)  { 
    //if(print) fprintf(outf,"The motif has no significant difference from the shoufled sequences.\n"); else fprintf(outf,"P-value: 1\n"); 
    fprintf(outf,"The motif has no significant difference from the shoufled sequences.\n");
  }
  else {

    double Wmax = (Wp > Wm) ? Wp : Wm;

    // Determine level of significance
    double p = 0;
    if(diffno <= 16) {   // Exact
      p = WilcoxPair(Wmax, diffno);
      fprintf(outf,"Z=0 N=%d\n",diffno);
    }
    else  {  // Normal approximation
      
      double Z = (Wmax - 0.5 - (double)diffno*(diffno+1)/4)/sqrt((double)diffno*(diffno+1)*(2*diffno+1)/24);
      p = Normal(Z);
      
      //if(print) fprintf(outf,"Z=%f N=%d\n",Z,diffno);
      fprintf(outf,"Z=%f N=%d\n",Z,diffno);

    }

    p=p/2;

    fprintf(outf,"P-value: %e\n",p);
  }

  GFREE(rank);

}


double WilcoxPair(double Winp, unsigned long int N)
{
  unsigned long int W, MaxW, No, C;
  unsigned long int Srank, i, j;
  double p;

  MaxW = N*(N+1)/2;
  if(Winp < MaxW/2)Winp = MaxW - Winp;
  W = (unsigned long int) Winp;   
  if(W != Winp)++W;  
  
  No = (unsigned long int) pow(2, N); 
  
  C = 0;
  for(i=0; i < No; ++i) {
    Srank = 0;
    for(j=0; j < N; ++j) {
      if((i >> j) & 1)Srank += j + 1;  
    }
    if(Srank >= W)++C;
  }

  p = 2*((double)C) / ((double)No);

  return p;
}

// 
// apply_test_student
// applies the student test on a vector of differences
//
void apply_test_student(FILE *outf,int maxno,FreqPos *maxf)
{
  //double rank[maxno];
  //computeranks(rank,maxno,maxf);
  
  fprintf(outf,"\n\nStudent t Test:\n\n");

  double mean=0;
  for(int i=0;i<maxno;i++) {
    //fprintf(stderr,"%d: %f\n",i,maxf[i].freq);
    mean+=maxf[i].freq;   // here I replaced rank[i] by maxf[i].freq
  }
  mean/=maxno;

  fprintf(outf,"Mean of differences in sequences: %f\n",mean);
  if(mean<=0 || maxno<=1) fprintf(outf,"The motif has no significant difference from the shoufled sequences.\n");
  else {
    double sd=0;
    for(int i=0;i<maxno;i++) {
      sd+=pow(maxf[i].freq-mean,2); // here I replaced rank[i] by maxf[i].freq
    }
    sd/=(maxno-1);
    sd=sqrt(sd);

    fprintf(outf,"SD of differences in sequences: %f\nN=%d\n",sd,maxno);
 
    test_student(outf,mean,0,sd,maxno);
  }
} 

//
// computeranks
// computes the ranks on a vector of values
//
void computeranks(double *rank,int maxno,FreqPos *maxf)
{
  qsort(maxf,maxno,sizeof(FreqPos),compf);

  double previous=0;
  int start=0;
  for(int i=0;i<maxno;i++) {
    if(fabs(maxf[i].freq)==previous) {
      double mean_rank=(double)(start+i+2)/2;
      for(int j=start;j<=i;j++) if(maxf[j].freq) rank[j] = maxf[j].freq>0 ? mean_rank : -mean_rank;
      else rank[j]=0;
    }
    else {
      rank[i]= maxf[i].freq>0 ? i+1 : -(i+1);
      previous=fabs(maxf[i].freq);
      start=i;
    }
  }

  /* printf("Ranks (val)  :\n");
     for(int i=0;i<maxno; i++) printf("%f (%d)  ",rank[i],maxf[i].freq);
     printf("\n"); */
}

//
// test_wilcoxon
// the actual wilcoxon test
//
double test_wilcoxon(FILE *outf,double W, int n, int n1, int n2)
{
  double pairs=n*(n+1)/2;

  W = W > pairs/2 ? pairs - W : W;

  int min_n = n1 < n2 ? n1 : n2;
  double C = (W >= pairs/2) ? -0.5 : 0.5;

  double Z=(W+C-min_n*(n+1)/2)/sqrt(n1*n2*(n+1)/12);
  Z=fabs(Z);

  double p = Normal(Z);

  fprintf(outf,"Motif significance approximation under Wilcoxon test: %f ",p/2);
  if(min_n<=10) fprintf(outf," - this significance value is unreliable.");
  fprintf(outf,"\n");

  return(p/2);
}


double Normal(double Z) 
{
  double b[5]={0.319381530, -0.356563782, 1.781477937, -1.821255978, 1.330274429};
  double pi = 3.14159265358979323846;
  double t = 1/(1+0.2316419*Z);

  double powt = t;
  double sum = 0;

  for(int i=0; i < 5; i++) {
    sum += b[i]*powt;
    powt *= t;
  };

  double P = 2*sum*exp(-Z*Z/2.0)/(sqrt(2*pi));
  
  return(P);
}

// 
// test_student
// the actual student test
//
double test_student(FILE *outf,double mean, double tmean, double sd, int n)
{
  double t=(mean-tmean)*sqrt((double)n)/sd;
  
  int line=30;

  if(n-1<=30) line=n-2;

  int col=7;

  while(col>=0 && t<t_value[line][col]) col--;

  switch(col) {
  case -1: 
    fprintf(outf,"P-value is bigger than 0.4 : the motif has no significant difference from the shoufled sequences.\n");
    return(1);
  case 0: 
    fprintf(outf,"0.25 < P-value <= 0.40 : the motif has no significant difference from the shoufled sequences at significance level 0.25\n");
    return(0.4);
  case 1: 
    fprintf(outf,"0.10 < P-value <= 0.25 : the motif has no significant difference from the shoufled sequences at significance level 0.1\n"); 
    return(0.25);
  case 2: 
    fprintf(outf,"0.05 < P-value <= 0.1 : the motif is statistically significant at significance level 0.1\n");
    return(0.1);
  case 3:
    fprintf(outf,"0.025 < P-value <= 0.05 : the motif is statistically significant at significance level 0.05\n");
    return(0.05);
  case 4:
    fprintf(outf,"0.01 < P-value <= 0.025 : the motif is statistically significant at significance level 0.025\n");
    return(0.025);
  case 5:
    fprintf(outf,"0.005 < P-value <= 0.01 : the motif is statistically significant at significance level 0.01\n");
    return(0.01);
  case 6:
    fprintf(outf,"0.0005 < P-value <= 0.005 : the motif is statistically significant at significance level 0.005\n");
    return(0.005);
  case 7: 
    fprintf(outf,"P-value <= 0.0005 : the motif is highly statistically significant.\n"); 
    return(0.0005);
  }
  
  return(1);
}

//
// compf
// sorting function on FreqPos
//
int compf(const void *a, const void *b)
{
  if(fabs(((FreqPos *)a)->freq)>fabs(((FreqPos *)b)->freq)) return(1);
  else if(fabs(((FreqPos *)a)->freq)==fabs(((FreqPos *)b)->freq)) return(0);
  else return(-1);

  // else return(0); using just this might make the system crash on Solaris
}

//
// Motif::nearoptim
// this function selects only those positions which are not less then perc from the maximum value
//
int Motif::nearoptim(CumWeight *Acum, int Alen, double perc)
{
  if(Alen==1) return(Alen);

  int retlen=Alen;
  double threshold=Acum[Alen-1].cw-Acum[Alen-2].cw-(Acum[Alen-1].cw-Acum[Alen-2].cw-Acum[0].cw)*perc;

  int i=Alen-2;
  while(i>=0) {
    double val=Acum[i].cw;
    if(i) val-=Acum[i-1].cw;
    if(val < threshold) {
      retlen=Alen-1-i;
      for(int j=0;j<retlen;j++) {
	Acum[j].cw=Acum[i+1+j].cw;
	Acum[j].pos=Acum[i+1+j].pos;
      }
      return(retlen);
    }
    i--;
  }

  return(retlen);
}

//   
// Motif::InfoPar
// computes information per parameter for a given motif with a computed MAP scF
//
double Motif::InfoPar(double scF)
{
  double info=0;
  //double myF=0;
  int i=0;
  double  *A;
  CumWeight *Acum;
  
  GMALLOC(A,(B->maxlen)*sizeof(double));
  GMALLOC(Acum,(B->maxlen)*sizeof(CumWeight));

  F->reset();
  FastaSeq fseq;

  sumlogweights=0;

  while(F->getFastaSeq(&fseq)) {

      compWeights(A,Acum,fseq.seq,fseq.len);
      complogweights(A,Acum[fseq.len-len].cw,fseq.len-len+1);    
    
      //myF+=log(A[Align[i]]); this computes the product of weights

    i++;
  }
       
  // var 1: info=(scF-sumloglength-sumlogweights)/((B->size-1)*len);

  // var 2: info=(sumloglength+sumlogweights)/(B->noseq);

  // var 3:   info=scF/((B->noseq)*(B->size-1)*len)+(sumloglength+sumlogweights)/(B->noseq);

  //var 4:   info=(scF+sumloglength+sumlogweights)/((B->noseq)*(B->size-1)*len);

  //var 5:   

  info=(scF-sumloglength-sumlogweights)/((B->noseq)*(B->size-1)*len);

  info=info/log(2);

  //printf("F=%.2f info=%.2f\n",scF/((B->noseq)*(B->size-1)*len), info);fflush(stdout);


  GFREE(A);
  GFREE(Acum);

  return(info);
}

// 
// Motif::Motif
// Motif constructor when two files are given at input to test the significance
//
Motif::Motif(const char *infile,FILE *outf,seqType t, const char *matrixfile, const char *pattern, int motiflen, int iter_no, int max_loop,int inlocmax, int det) 
{
  F = new GFastaFile(infile);
  B = new Background(F,outf,t);

  motiffreq=NULL;
  motifprob=NULL;

  Align=new int[B->noseq];

  defmotif=0;
  len=motiflen;
  sumloglength=0;
  sumlogweights=0;

  F->reset();
  FastaSeq fseq;  

  motiflen=strlen(pattern);
  if(motiflen) len=motiflen;

  if(!len && strlen(matrixfile)) {

    FILE *M;
        
    M = fopen(matrixfile,"r");
    if(M == NULL) {
      fprintf (stderr, "ERROR:  Unable to open file %s.\n", matrixfile);
      exit (-1);
    }
    
    char line[5000];
    while( fgets(line,5000,M) != NULL) {

      if(!strcmp(line,"Motif counts:\n")) {

	for(int i=0;i<4;i++) {
	  fgets(line,5000,M);
	  
	  char **myargv;
	  char delim[]=" ";
	  int no=makeargv(line,delim,&myargv);

	  if(no<0) {
	    fprintf(stderr,"ERROR: Wrong matrix file format: %s\n",matrixfile);
	    exit(-1);
	  }
	  else {
	    if(!i) {
	      len=no-2;
	      MCount = new int*[len+1];                             // MCount=new int[len+1][size];
	      for(int j=0;j<len+1;j++) MCount[j] = new int[B->size];

	      MFreq = new double*[len+1];                             //MFreq=new double[len+1][size]; 
	      for(int j=0;j<len+1;j++) MFreq[j] = new double[B->size];     

	      initMCount();
	    }

	    for(int j=1;j<=len;j++) {

	      int tok=atoi(myargv[j]);
	      MCount[0][i]-=tok;
	      sumMCount0-=tok;
	      MCount[j][i]+=tok;
	    }
	  }
	}
	updateMFreq(B->noseq);
	break;
      }
    }

    fclose(M);

    q = new double *[len];
    Q = new double **[len];
      
    for(int i=0;i<len;i++) {
      q[i] = new double[B->size];
      Q[i] = new double *[B->size];
      for(int j=0;j<B->size;j++) Q[i][j] = new double[B->size];
    }



    int i=0;
    while(F->getFastaSeq(&fseq)) {

      if(fseq.len<len) GError("Motif length is greater then input sequence %s\n",fseq.id);
      sumloglength+=log(fseq.len-len+1);
      
      int pos=maxWeight(fseq.seq,fseq.len); 
      if(pos==-1) { 
	fprintf(stderr,"No motif position possible in sequence %s. Ignoring sequence...\n",fseq.id);
	continue;
      }
      Align[i]=pos;
      i++;
    } 

  }
  else {

    // init parameters;

    PM=new ProbMotif[B->noseq];
    for(int i=0;i<B->noseq;i++) PM[i].pos=NULL;

    q = new double *[len];
    Q = new double **[len];
      
    for(int i=0;i<len;i++) {
      q[i] = new double[B->size];
      Q[i] = new double *[B->size];
      for(int j=0;j<B->size;j++) Q[i][j] = new double[B->size];
    }


    MCount = new int*[len+1];                             // MCount=new int[len+1][size];
    for(int i=0;i<len+1;i++) MCount[i] = new int[B->size];

    MFreq = new double*[len+1];                             //MFreq=new double[len+1][size]; 
    for(int i=0;i<len+1;i++) MFreq[i] = new double[B->size];     


    if(motiflen) {

      initMCount();
      for(int j=1;j<len+1;j++) {
	int k=basetoint(pattern[j-1],B->SType);
	for(int i=0;i<4;i++) { 
	  MCount[j][i]=1;
	  MCount[0][i]-=1;
	}
	MCount[j][k]=B->noseq-3;
	MCount[0][k]-=B->noseq-4;
	sumMCount0-=B->noseq;
      }
      updateMFreq(B->noseq);

      int i=0;
      while(F->getFastaSeq(&fseq)) {

	if(fseq.len<len) GError("Motif length is greater then input sequence %s\n",fseq.id);
	sumloglength+=log(fseq.len-len+1);
      
	int pos=maxWeight(fseq.seq,fseq.len); 
	if(pos==-1) { 
	  fprintf(stderr,"No motif position possible in sequence %s. Ignoring sequence...\n",fseq.id); 
	  continue;
	}
	Align[i]=pos;
	i++;
      } 

    }
    else {

      while(F->getFastaSeq(&fseq)) {
	
	if(fseq.len<len) GError("Motif length is greater then input sequence %s\n",fseq.id);
	sumloglength+=log(fseq.len-len+1);

      }

      fprintf(stderr,"Compute motif of input file...\n");
      double globAlignProb;
    
      globAlignProb=findMotif(iter_no,max_loop,inlocmax,1,det);

      double info;
      // optimizing
      fprintf(stderr,"Optimizing...\n");
      globAlignProb=optimize(globAlignProb,info,0);
      fprintf(outf,"MAP for motif: %.3f InfoPar=%.3f\n\n",globAlignProb,InfoPar(globAlignProb));
    }
  }
}

//
// Motif::twofilesignif
// compute significance of motif when the test sequences are defined by joining
// the input file and the test file together; if det is set to 1 the significance
// is computed deterministically otherwise it's computed by running the Gibbs 
// sampler
//
void Motif::twofilesignif(int det,const char *testfile, int SignfRunNo,int print, const char *pattern,int mdeg,int exact) 
{
  char **fkseq;
  int *posfreqlen;

  fkseq = new char*[B->noseq];
  motiffreq = new int*[B->noseq];
  posfreqlen= new int[B->noseq];

  createtwofile(testfile,fkseq,posfreqlen);

  if(det) {

    char consensus[len+1];
    
    if(pattern==NULL || pattern[0]=='0') {
      print_consensus(consensus);
      PM=new ProbMotif[B->noseq];
      defmotif=1;
    }
    else strcpy(consensus,pattern);

    fprintf(B->outf,"Pattern: %s ",consensus);fflush(stdout);


    for(int i=0;i<B->noseq;i++) {    
      if(exact) {
	if(!MatchSeqPattExact(fkseq[i],consensus,i,strlen(fkseq[i]),1)) {
	  strcpy(fkseq[i],fkseq[-1+B->noseq]);
	  B->noseq=-1+B->noseq;
	  i--;
	}
      }
      else MatchSeqPatt(fkseq[i],consensus,i,strlen(fkseq[i]));
    }

    //fprintf(stderr,"noseq=%d",B->noseq);
    
    testsignif(1,fkseq,posfreqlen,print,exact);

    if(pattern==NULL || pattern[0]=='0') defmotif=0;

  }
  else {
    // run Gibbs sampling SignfRunNo times and compute frequencies of choosing motifs
    GibbsSample(SignfRunNo,fkseq,print);
  
    // compute motif ranking and significance
    testsignif(0,fkseq,posfreqlen,print,0);

  }

  for(int i=0;i<B->noseq;i++) delete [] fkseq[i];
  delete [] fkseq;
  
  delete [] posfreqlen;
}

//
// Motif::Motif
// motif constructor for a given file, sequence type and pattern
//
Motif::Motif(const char *infile,FILE *outf,seqType t, const char *pattern) 
{

  F = new GFastaFile(infile);
  B = new Background(F,outf,t);

  defmotif=1;

  Align=new int[B->noseq];
  PM=new ProbMotif[B->noseq];
  for(int i=0;i<B->noseq;i++) PM[i].pos=NULL;
  len=strlen(pattern);

  motiffreq=NULL;
  motifprob=NULL;

  sumloglength=0;
  sumlogweights=0;

  MCount = new int*[len+1];                             // MCount=new int[len+1][size];
  MFreq = new double*[len+1];                           // MFreq=new double[len+1][size]; 
  q = new double *[len];
  Q = new double **[len];

  for(int i=0;i<len+1;i++) {
    MCount[i] = new int[B->size];
    MFreq[i] = new double[B->size];     
  }
  for(int i=0; i<len;i++) {
    q[i] = new double[B->size];
    Q[i] = new double *[B->size];
    for(int j=0;j<B->size;j++) Q[i][j] = new double[B->size];
  }

  F->reset();
  FastaSeq fseq;
  int i=0;

  while(F->getFastaSeq(&fseq)) {

    if(fseq.len<len) GError("Motif length is greater then input sequence %s\n",fseq.id);

    MatchSeqPatt(fseq.seq,pattern,i,fseq.len); // this function computes PM[i]
    sumloglength+=log(fseq.len-len+1);
   
    i++;

  }

}

//
// Motif::findMDet
// find motif deterministically
//
double Motif::findMDet(int print)
{
  double optinfo=0;;
  double diff=1;
  int t=0;
  initMCount();

  while(diff>B->eps) {
    
    int i=0;
    F->reset();
    FastaSeq fseq;
    double maxinfo;
    while(F->getFastaSeq(&fseq)) {
      maxinfo=-HUGE_VAL;
      int maxpos=-1;
      if(t) updateMCount(Align[i],fseq.seq,B->SType,-1);
      int nelem;
      if(defmotif) nelem=PM[i].no;
      else nelem=fseq.len-len;
      for(int j=0;j<nelem;j++) {
	int pos;
	if(defmotif) pos=PM[i].pos[j];
	else pos=j;
	updateMCount(pos,fseq.seq,B->SType,1);
	int no=i+1;
	if(t) no=B->noseq;
	updateMFreq(no);
	double info=0;
	for(int ji=1;ji<=len;ji++) {
	  for(int ii=0;ii<B->size;ii++)
	    info+=MFreq[ji][ii]*log(MFreq[ji][ii]/MFreq[0][ii])/log(2);
	}
	if(info>maxinfo) { maxinfo=info; maxpos=pos;}
	updateMCount(pos,fseq.seq,B->SType,-1);
      }
      Align[i]=maxpos;
      updateMCount(Align[i],fseq.seq,B->SType,1);
      if(print) fprintf(stderr,"%d: %4.2f          \r",i+1,maxinfo);
    
      //printFreq();
      i++;
    }

    if(print) fprintf(stderr,"\n");
    //fprintf(stderr,"maxinfo=%f optinfo=%f diff=%f \n",maxinfo,optinfo,diff);

    diff=maxinfo-optinfo;
    optinfo=maxinfo;


    //double globAlignProb=AlignProb();
    //printf("MAP for motif: %.3f InfoPar=%.3f\n\n",globAlignProb,InfoPar(globAlignProb));
    //printFreq();
    t++;
  }

  double globAlignProb=AlignProb();

  /* fprintf(stderr,"Optimizing...\n");
     double info;
     globAlignProb=optimize(globAlignProb,info);
     if(print) {
     printf("MAP for motif: %.3f InfoPar=%.3f\n\n",globAlignProb,InfoPar(globAlignProb));
     printFreq();
     }*/

  return(globAlignProb);
}

//
// Motif::computeLLC
// compute the Least Likely Consensus score and the IC for a given pattern
// 
double Motif::computeLLC(const char *pattern, int print)
{
  double llc;

  int maxnoinseq=computemaxPMno();

  while(maxnoinseq>1) {
    //    fprintf(stderr,"Max no of elements in seqs: %d           \r",maxnoinseq);
    computeQq();
    removeoneword();
    maxnoinseq--;
  }
  computeQq();
  fprintf(stderr,"\n");
  llc=LLC(print);

  for(int i=0; i<B->noseq; i++) {
    Align[i]=PM[i].pos[0];
  }

  computeMotif(Align);

  if(print) {
    char consensus[len+1];
    print_consensus(consensus);
    fprintf(B->outf," Consensus: %s ",consensus);
  }
  
  // redo PM

  F->reset();
  FastaSeq fseq;
  int i=0;
  
  while(F->getFastaSeq(&fseq)) {
    
    if(fseq.len<len) GError("Motif length is greater then input sequence %s\n",fseq.id);

    MatchSeqPatt(fseq.seq,pattern,i,fseq.len); // this function computes PM[i]
   
    i++;

  }

  return(llc);

}

//
// Motif::print_consensus
// computes consensus for a given matrix score
//
void Motif::print_consensus(char *consens)
{
  for(int i=1;i<=len;i++) {
    int max=0;
    int pos=0;
    for(int j=0;j<B->size;j++) 
      if(MCount[i][j]>max) { 
	pos=j;
	max=MCount[i][j];
      }
    consens[i-1]=inttobase(pos,B->SType);
  }
  consens[len]='\0';
}

//
// Motif::detrunforsignif
// find motif deterministically, supposing that a pattern is given
//
void Motif::detrunforsignif(int print, const char *pattern,int zparam,int mdeg)
{
  char **fkseq;
  int *posfreqlen;

  posfreqlen= new int[B->noseq];

  fkseq = new char*[B->noseq];
  motiffreq = new int*[B->noseq];

  createrandfile(fkseq,posfreqlen,zparam,mdeg);

  for(int i=0;i<B->noseq;i++) {    
    MatchSeqPatt(fkseq[i],pattern,i,strlen(fkseq[i]));
  }

  testsignif(1,fkseq,posfreqlen,print,0);

  for(int i=0;i<B->noseq;i++) delete [] fkseq[i];
  delete [] fkseq;

  delete [] posfreqlen;
}

// 
// Motif::LLC
// computes LLC score
//
double Motif::LLC(int print)
{
  double IC=0;
  double I=0;
  double *sq;

  sq=new double[len];
  
  for(int i=0;i<len;i++) {
    sq[i]=0;
    for(int j=0;j<B->size;j++) 
      sq[i]+=q[i][j];
    for(int j=0;j<B->size;j++) {
      I+=q[i][j]/sq[i]*log(q[i][j]/sq[i]);
      IC+=q[i][j]/sq[i]*log(q[i][j]*(B->sumc+B->eps)/(sq[i]*(B->c[j]+B->eps)));
    }
  }

  I/=log(2);
  IC/=log(2);

  double BR=0;
  F->reset();
  FastaSeq fseq;
  int i=0;
  while(F->getFastaSeq(&fseq)) {
    BR+=computeBR(fseq.seq,PM[i].pos[0]);
    i++;
  }
  BR/=log(2)*B->noseq;


  //double LLCval=log(B->noseq)*(IC-BR)/log(2);
  double LLCval=I-BR;
  
  if(print) fprintf(B->outf,"IC = %f I= %f BR = %f ",IC,I,BR); 

  delete [] sq;
  
  return(LLCval);
  //return(IC);
}

//
// Motif::computeBR
// computes Background Rareness score needed by LLC
//
double Motif::computeBR(char *seq, int pos) 
{
  int pb=basetoint(seq[pos],B->SType);
  double BR=0;
  if(pb!=-1) BR=log((B->c[pb]+B->eps)/(B->sumc+B->eps));

  for(int i=1;i<len;i++) {
    int b=basetoint(seq[pos+i],B->SType);
    if(b!=-1 && pb != -1) BR+=log((B->p[b][pb]+B->eps)/(B->c[pb]+B->eps));
    pb=b;
  }

  return(BR);
}

//
// Motif::removeoneword
// from each sequence in the input file with more then one candidate for the motif 
// removes the word that satisifes less the consensus 
//
void Motif::removeoneword()
{
  F->reset();
  
  FastaSeq fseq;
  int i=0;
  while(F->getFastaSeq(&fseq)) {
    if(PM[i].no>1) {
      double minTR=HUGE_VAL;
      int minpos=-1;
      for(int j=0; j<PM[i].no;j++) {
	double TR=computeTR(fseq.seq,PM[i].pos[j]);
	if(TR<minTR) {
	  minTR=TR;
	  minpos=j;
	}
      }
      PM[i].pos[minpos]=PM[i].pos[PM[i].no-1];
      PM[i].no--;
    }
    i++;
  }
}

//
// Motif::removeoneword
// from each sequence in the input+test(random) file with more then one candidate for the motif 
// removes randomly on of the words that satisify less the consensus 
//
void Motif::removeoneword(char **fkseq)
{
  for(int i=0;i<B->noseq;i++) {

    if(PM[i].no>1) {

      int *minpos;
      minpos= new int[PM[i].no];
      int minposno=0;

      double minTR=HUGE_VAL;
      for(int j=0; j<PM[i].no;j++) {
	double TR=computeTR(fkseq[i],PM[i].pos[j]);
	if(TR<minTR) { // removes the first word that satisfies less the consensus
	  minTR=TR;
	  minpos[0]=j;
	  minposno=1;
	}
	else if(TR==minTR) {
	  minpos[minposno++]=j;
	}
      }

      int k =rand() % minposno;

      //fprintf(stderr,"k=%d minposno=%d\n",k,minposno);

      PM[i].pos[minpos[k]]=PM[i].pos[PM[i].no-1];
      PM[i].no--;

      delete [] minpos;

    }

  }
}

//
// Motif::computeTR
// computes the Transition score from the consensus score matrix of the motif;
// needed by the removeoneword procedure
//
double Motif::computeTR(char *seq, int pos) 
{
  int pb=basetoint(seq[pos],B->SType);
  double TR=0;
  if(pb!=-1) TR=log(q[0][pb]/sumq)-log((B->c[pb]+B->eps)/(B->sumc+B->eps)); 
  for(int i=1;i<len;i++) {
    int b=basetoint(seq[pos+i],B->SType);
    //    if(b!=-1 && pb!=-1) TR+=log(Q[i][b][pb]/q[i][pb])-log((B->p[b][pb]+B->eps)/(B->c[pb]+B->eps));
    if(b!=-1 && pb!=-1) TR+=log(Q[i][b][pb]/q[i-1][pb])-log((B->p[b][pb]+B->eps)/(B->c[pb]+B->eps));
    pb=b;
  }
  return(TR);
}

//
// Motif::computeQq
// computes the first order transition score matrix of the motif from the input file
//
void Motif::computeQq() 
{
  sumq=B->eps;
  for(int i=0;i<len;i++) {
    for(int j=0;j<B->size;j++) {
      q[i][j]=B->eps;
      for(int k=0;k<B->size;k++) {
	Q[i][j][k]=B->eps;
      }
    }
  }

  F->reset();

  FastaSeq fseq;
  int i=0;
  sumq=0;
  while(F->getFastaSeq(&fseq)) {

    for(int j=0;j<PM[i].no;j++) {

      int pb=basetoint(fseq.seq[PM[i].pos[j]],B->SType);
      if(pb!=-1) q[0][pb]++;
      sumq++;

      for(int k=1;k<len;k++) {

	int b=basetoint(fseq.seq[PM[i].pos[j]+k],B->SType);
	
	if(b!=-1) {
	  q[k][b]++;
	  if(pb!=-1) Q[k][b][pb]++;
	}
	pb=b;

      }

    }
    i++;
  }
}
	
//
// Motif::computeQq
// computes the first order transition score matrix of the motif from the input+test files
//
void Motif::computeQq(char **fkseq) 
{

  sumq=B->eps;
  for(int i=0;i<len;i++) {
    for(int j=0;j<B->size;j++) {
      q[i][j]=B->eps;
      for(int k=0;k<B->size;k++) {
	Q[i][j][k]=B->eps;
      }
    }
  }


  sumq=0;
  for(int i=0; i<B->noseq;i++) {

    for(int j=0;j<PM[i].no;j++) {

      int pb=basetoint(fkseq[i][PM[i].pos[j]],B->SType);
      if(pb!=-1) q[0][pb]++;
      sumq++;

      for(int k=1;k<len;k++) {

	int b=basetoint(fkseq[i][PM[i].pos[j]+k],B->SType);
	
	if(b!=-1) {
	  q[k][b]++;
	  if(pb!=-1) Q[k][b][pb]++;
	}
	pb=b;

      }

    }
  }


}

//
// Motif::computemaxPMno
// computes the maximum number of elements of PM in all sequences
//
int Motif::computemaxPMno()
{
  int max=0;
  for(int i=0;i<B->noseq;i++) {
    if(PM[i].no>max) {
      max=PM[i].no;
    }
  }
  return(max);
}

//
// Motif::computeMotif
// copies the argument array into the motif array and updates the 
// motif counts
//
void Motif::computeMotif(int *MaxAlign)  
{
  if(MaxAlign != Align) copyarray(Align,MaxAlign,B->noseq);

  initMCount();

  F->reset();
  FastaSeq fseq;
  int i=0;

  while(F->getFastaSeq(&fseq)) {
    updateMCount(Align[i],fseq.seq,B->SType,1);
    i++;
  }

  // compute residue frequencies now;
  updateMFreq(B->noseq);

}

//
// Motif::check
// checks if pos is between the PM[i].pos positions; returns 1 if yes, 0 otherwise
//
int Motif::check(int pos, int i)
{
  int found=0;
  int j=0;

  while(!found && j<PM[i].no) {
    if(pos==PM[i].pos[j]) found=1;
    j++;
  }
  
  return(found);
}

//
// Motif::optimize
// optimizes the MAP probability computed by the findmotif function; 
// the optimization can be done according to the presence of the pattern (wpattern=1)
// or not 
//
double Motif::optimize(double prevMax, double previnfo, int wpattern)
{
  double probMotif;
  int optimize,maxweightpos,prevalign;
  FastaSeq fseq;
  int *copyAlign=new int[B->noseq];

  do {

    copyarray(copyAlign,Align,B->noseq);

    optimize=0;

    F->reset();
    int i=0;

    while(F->getFastaSeq(&fseq)) {
      maxweightpos=maxWeight(fseq.seq,fseq.len);

      if(maxweightpos==-1) continue;

      if(wpattern && defmotif && !check(maxweightpos,i)) maxweightpos=Align[i]; // if maximum weight positions is not in the neighbouring vicinity of the motif as defined by PM then take as maximum the Align[i] position; this way you don't allow the motif to randomly go further from the initial consensus

      if(maxweightpos != Align[i] ) { 
	updateMCount(Align[i],fseq.seq,B->SType,-1);
	prevalign=Align[i];
	Align[i]=maxweightpos;
	updateMCount(Align[i],fseq.seq,B->SType,1);
	updateMFreq(B->noseq);
	probMotif=AlignProb();
	if(probMotif>prevMax) {
	  optimize=1;
	  prevMax=probMotif;
	}
	else {
	  updateMCount(Align[i],fseq.seq,B->SType,-1);
	  Align[i]=prevalign;
	  updateMCount(Align[i],fseq.seq,B->SType,1);
	  updateMFreq(B->noseq);
	}
      }
      i++;
    }


    /*if(optimize) {
      double info=InfoPar(probMotif);
      if(info>previnfo) {
	fprintf(stderr,"New align prob=%.2f Info/par=%.2f\n",prevMax,info);
	previnfo=info;
      }
      else {
	optimize=0;
	computeMotif(copyAlign);
      }
      }*/

    if(optimize) fprintf(stderr,"New align prob=%.2f\n",probMotif);

  }
  while(optimize);
  
  delete [] copyAlign;

  return(prevMax);
}

//
// Motif::maxWeight
// computes the position that achieves the maximum weight in the sequence
// according to the motif score matrix
//
int Motif::maxWeight(char *fseq,int seqlen)
{
  double Qx,Px;
  double weight,maxweight=0;
  int maxpos=-1;

  for(int i=0;i<seqlen-len+1;i++) {

    Qx=1;Px=1;

    for(int j=1;j<=len;j++) {
      int val=basetoint(fseq[i+j-1],B->SType);
      if(val==-1) { Qx=0;Px=1;break;}
      Qx*=MFreq[j][val];
      Px*=MFreq[0][val];
    }

    weight=Qx/Px;
    
    if(weight>maxweight) {
      maxpos=i;
      maxweight=weight;
    }
  }
  
  //  if(maxpos==-1) GError("No valid motif positions (without N's found for sequence %s\n",fseq);  // initial variant

  return(maxpos);
}

// 
// Motif::printMotif
// prints the motif
//
void Motif::printMotif(int printmax,int SignifNo,int print)
{
  fprintf(B->outf,"Motif found:\n\n");

  printFreq();
  
  F->reset();
  FastaSeq fseq;
  int i=0;
  int copied;

  char align[len+14];
  char motif[len+1];
  char *spaces;

  if(len>5) {
    spaces=new char[len-3];
    for(int j=0;j<len-4;j++) spaces[j]=' ';
    spaces[len-4]='\0';
  }
  else {
    spaces=new char[2];
    strcpy(spaces," ");
  }

  if(print) {
    fprintf(B->outf,"    Seq.no  Pos ***** Motif%s***** Prob    D Seq.Id\n",spaces);

    while(F->getFastaSeq(&fseq)) {
    
      strncpy(motif,fseq.seq+Align[i],len);
      motif[len]='\0';
      copied=0;

      for(int j=Align[i]-5;j<=Align[i]+len+5;j++) 
	if(j<0) align[j-Align[i]+5]=' ';
	else if(j>=fseq.len) align[j-Align[i]+7]=' ';
	else if(j<Align[i]) align[j-Align[i]+5]=tolower(fseq.seq[j]);
	else if (j>=Align[i]+len) align[j-Align[i]+7]=tolower(fseq.seq[j]);
	else { //align[j-Align[i]+5]=toupper(fseq.seq[j]); 
	  if(!copied) {
	    align[5]='\0';
	    strcat(align," ");
	    strcat(align,motif);
	    strcat(align," ");
	    copied=1;
	  }
	  align[j-Align[i]+6]=toupper(align[j-Align[i]+6]);
	}
      align[len+12]='\0';
      
    
      int dist = defmotif ? PM[i].dist : -1;
      fprintf(B->outf,"%10d %4d %s %6.4f %2d %s\n",i+1,Align[i]+1,align,scoreF(motif),dist,fseq.id);
      
      i++;
    }

    if(motifprob != NULL) {
      fprintf(B->outf,"\n\n******\n\nMost frequent sites matching motif:\n(sites appearing more then half a time are marked by *)\n\n");fflush(stdout);

      fprintf(B->outf,"     Seq.no  Pos ***** Motif%s***** Prob    D Seq.Id\n",spaces);fflush(stdout);
      F->reset();
      i=0;

      while(F->getFastaSeq(&fseq)) {

	strncpy(motif,fseq.seq+motifprob[i].pos,len);
	motif[len]='\0';
	copied=0;
      
	for(int j=motifprob[i].pos-5;j<=motifprob[i].pos+len+5;j++) 
	  if(j<0) align[j-motifprob[i].pos+5]=' ';
	  else if(j>=fseq.len) align[j-motifprob[i].pos+7]=' ';
	  else if(j<motifprob[i].pos) align[j-motifprob[i].pos+5]=tolower(fseq.seq[j]);
	  else if (j>=motifprob[i].pos+len) align[j-motifprob[i].pos+7]=tolower(fseq.seq[j]);
	  else { //align[j-motifprob[i].pos+5]=toupper(fseq.seq[j]); 
	    if(!copied) {
	      align[5]='\0';
	      strcat(align," ");
	      strcat(align,motif);
	      strcat(align," ");
	      copied=1;
	    }
	    align[j-motifprob[i].pos+6]=toupper(align[j-motifprob[i].pos+6]);
	  }
	align[len+12]='\0';
      
	int dist = defmotif ? PM[i].dist : -1;
	char c=' ';
	if(motifprob[i].freq >= (double)SignifNo/2) c='*';
	fprintf(B->outf,"%c%10d %4d %s %6.4f %2d %s\n",c,i+1,motifprob[i].pos+1,align,motifprob[i].freq/(double)SignifNo,dist,fseq.id);fflush(stdout);
	i++;
      }
    }
  }

  // print max positions within motif
  
  
  if(printmax) {
    fprintf(B->outf,"\n\n******\n\nPositions within sequences for which the maximum weight according to motif is reached:\n\n");

    fprintf(B->outf,"    Seq.no  Pos ***** Motif%s***** Prob    D Seq.Id\n",spaces);
    F->reset();
    i=0;
    int pos;
    
    while(F->getFastaSeq(&fseq)) {
    
      pos=maxWeight(fseq.seq,fseq.len); 
      if(pos==-1) continue;
      strncpy(motif,fseq.seq+pos,len); 
      motif[len]='\0';
      copied=0;
      
      for(int j=pos-5;j<=pos+len+5;j++) 
	if(j<0) align[j-pos+5]=' ';
	else if(j>=fseq.len) align[j-pos+7]=' ';
	else if(j<pos) align[j-pos+5]=tolower(fseq.seq[j]);
	else if(j>=pos+len) align[j-pos+7]=tolower(fseq.seq[j]);
	else { //align[j-Align[i]+5]=toupper(fseq.seq[j]);
	  if(!copied) {
	    align[5]='\0';
	    strcat(align," ");
	    strcat(align,motif); 
	    strcat(align," ");
	    copied=1;
	  }
	  align[j-pos+6]=toupper(align[j-pos+6]);
	}
      align[len+12]='\0';
      
    
      int dist = defmotif ? PM[i].dist : -1;
      fprintf(B->outf,"%10d %4d %s %6.4f %2d %s\n",i+1,pos+1,align,scoreF(motif),dist,fseq.id);
    
      i++;
    }
  }

  if(spaces != NULL) delete [] spaces;

}

//
// Motif::scoreF
// computes the probability of seeing that motif according to the probabilities in
// the motif score matrix
//
double Motif::scoreF(char *motif)
{
  double F=1;

  for(int i=1;i<=len;i++) {
    int val=basetoint(motif[i-1],B->SType);
    //assert(val!=-1);
    if(!MFreq[0][val]) GError("Background frequency q[0][%d] is 0!\n",val);
    //if(MFreq[i][val]) F+=log(MFreq[i][val]/MFreq[0][val]);
    if(MFreq[i][val]) F*=MFreq[i][val];
  }

  return(F);
}
  
//
// Motif::initMCount
// initializes all the frequencies defining the background and the motif
//
void Motif::initMCount()
{
  sumMCount0=0;
  for(int i=0;i<B->size;i++) {
    MCount[0][i]=B->c[i];
    sumMCount0+=MCount[0][i];
    for(int j=1;j<=len;j++) MCount[j][i]=0;
  }
}

// 
// Motif::initRandAlign
// computes an initial (random) alignment of the sequences in the input files
//
void Motif::initRandAlign()
{

  initMCount();

  F->reset();
  FastaSeq fseq;
  int i=0;

  while(F->getFastaSeq(&fseq)) {
  
    if(defmotif) { 
      int r = rand() % PM[i].no;
      Align[i]=PM[i].pos[r];
    }
    else Align[i]=rand() % (fseq.len-len);


    updateMCount(Align[i],fseq.seq,B->SType,1);
    i++;

  }

  // compute residue frequencies now;
  // updateMFreq();

}

//
// Motif::updateMFreq
// updates the probabilities of the motif according to the frequencies 
// computed so far
//
void Motif::updateMFreq(int no)
{
  
  for(int i=0;i<B->size;i++) {
    if(!sumMCount0) GError("Sum count of frequencies is 0!\n");
    MFreq[0][i]=((double)MCount[0][i]+B->d[i])/(sumMCount0+B->D);
    for(int j=1;j<=len;j++) {
      if(B->noseq<=1) GError("Two few sequencies in data.\n");
      //MFreq[j][i]=((double)MCount[j][i])/(B->noseq-1);
      MFreq[j][i]=(MCount[j][i]+B->d[i])/(no+B->D);
    }
  }

}

// 
// Motif::printFreq
// prints all background and motif frequencies
//
void Motif::printFreq()
{
  B->printFreq();
  //printf("summcount=%d\n",sumMCount0);
  
  fprintf(B->outf,"Motif probability model:\n");
  fprintf(B->outf,"Pos: ");
  for(int j=1;j<=len;j++) fprintf(B->outf," %4d ",j); 
  fprintf(B->outf,"\n");
  for(int i=0;i<B->size;i++) {
    fprintf(B->outf,"%c    ",inttobase(i,B->SType));
    for(int j=1;j<=len;j++)  fprintf(B->outf," %4.2f ",MFreq[j][i]);      
    fprintf(B->outf,"\n");
  }
  for(int j=0;j<=len;j++) fprintf(B->outf,"------");
  fprintf(B->outf,"\nInfo ");
  for(int j=1;j<=len;j++) {
    double info=0;
    for(int i=0;i<B->size;i++)
      info+=MFreq[j][i]*log(MFreq[j][i]/MFreq[0][i])/log(2);
    fprintf(B->outf," %4.2f ",info);
  }

  fprintf(B->outf,"\n\nMotif counts:\n");
  for(int i=0;i<B->size;i++) {
    fprintf(B->outf,"%c: ",inttobase(i,B->SType));
    for(int j=1;j<=len;j++) fprintf(B->outf,"%15d ",MCount[j][i]);
    fprintf(B->outf,"\n");
  }
  fprintf(B->outf,"\n\n");
}

//
// Motif::~Motif
// Motif destructor
//
Motif::~Motif()
{ 
  delete [] Align;
  for(int i=0;i<len+1;i++) {
    delete [] MCount[i];
    delete [] MFreq[i];
  }
  delete [] MCount;
  delete [] MFreq;
  if(PM != NULL) {
    for(int i=0;i<B->noseq;i++) {
      GFREE(PM[i].pos);
    }
    delete [] PM;
    }
  delete F;
  if(motiffreq != NULL) {
    for(int i=0;i<B->noseq;i++) delete [] motiffreq[i];
    delete [] motiffreq;
  }
  if(motifprob != NULL) delete [] motifprob;
  if(q != NULL) {
    for(int i=0;i<len;i++) {
      delete [] q[i];
      for(int j=0;j<B->size;j++) delete [] Q[i][j];
      delete [] Q[i];
    }
    delete [] q;
    delete [] Q;
  }
  delete B;
}

//
// Motif::MatchSeqPatt
// computes (updates) the PM vector for a given pattern and sequence
//
void Motif::MatchSeqPatt(char *seq, const char *pattern, int i, int seqlen) 
{
  int mindist=len+1;
  int d;
  int allocsize=10;
  int no=0;
  char *copy;



  copy=new char[len+1];
  strcpy(copy,pattern);

  for(int j=0;j<seqlen-len+1;j++) {


    d=strdist(seq+j,copy,len);

    if(d<=mindist) {

      if(d==mindist) {
	if(no % allocsize == 0) {
	  GREALLOC(PM[i].pos,(no+allocsize)*sizeof(int));
	}
	PM[i].pos[no]=j;
	no++;
      }
      else {

	if(PM[i].pos!=NULL) GFREE(PM[i].pos);

	GMALLOC(PM[i].pos,allocsize*sizeof(int));
	PM[i].pos[0]=j;
	no=1;
	mindist=d;
      }
    }
  }

  delete []copy;

  PM[i].no=no;
  PM[i].dist=mindist;

  /*fprintf(stderr,"PM[%d].pos: ",i);
    for(int j=0;j<no;j++) fprintf(stderr,"%d ",PM[i].pos[j]);
    fprintf(stderr,"\n");exit(0);*/
}    


//
// Motif::MatchSeqPattExact
// computes (updates) the PM vector for a given pattern and sequence; finds only exact matches
// if ignore is 1 then positions at the end before middle of sequence needs to be ignored
int Motif::MatchSeqPattExact(char *seq, const char *pattern, int i, int seqlen,int ignore) 
{ 
  int d;
  int allocsize=10;
  int no=0;
  char *copy;



  copy=new char[len+1];
  strcpy(copy,pattern);

  if(ignore) {
    for(int j=0;j<seqlen/2-len+1;j++) {

      d=strdist(seq+j,copy,len);

      if(d==0) {
	if(no % allocsize == 0) {
	  GREALLOC(PM[i].pos,(no+allocsize)*sizeof(int));
	}
	PM[i].pos[no]=j;
	no++;
      }
    }
    for(int j=seqlen/2;j<seqlen-len+1;j++) {

      d=strdist(seq+j,copy,len);

      if(d==0) {
	if(no % allocsize == 0) {
	  GREALLOC(PM[i].pos,(no+allocsize)*sizeof(int));
	}
	PM[i].pos[no]=j;
	no++;
      }
    }
  }
  else 
    for(int j=0;j<seqlen-len+1;j++) {

      
      d=strdist(seq+j,copy,len);

      if(d==0) {
	if(no % allocsize == 0) {
	  GREALLOC(PM[i].pos,(no+allocsize)*sizeof(int));
	}
	PM[i].pos[no]=j;
	no++;
      }
    }


  delete []copy;

  PM[i].no=no;
  PM[i].dist=0;

  //  if(no!=0) fprintf(stderr,"i=%d pos=%d seq=\n%s\n",i,PM[i].pos[0],seq);

  /*printf("PM[%d].pos: ",i);
    for(int j=0;j<no;j++) printf("%d ",PM[i].pos[j]);
    printf("\n");*/

  if(no==0) return(0);
  else return(1);

}    

//
// Motif::compWeights
// computes the weights of all positions in a sequence
//
void Motif::compWeights(double *A, CumWeight *Acum, char * fseq, int seqlen)
{
  double Qx,Px;

  for(int i=0;i<seqlen-len+1;i++) {

    Qx=1;Px=1;

    for(int j=1;j<=len;j++) {
      int val=basetoint(fseq[i+j-1],B->SType);
      if(val==-1) { Qx=0;Px=1; break;}
      Qx*=MFreq[j][val];
      Px*=MFreq[0][val];
    }
    if(!Px) GError("Px is 0!");

    A[i]=Qx/Px;

    Acum[i].cw = A[i];
    Acum[i].pos = i;

  }

  qsort(Acum, seqlen-len+1, sizeof(CumWeight), comp);

  //printf("A[%d]=%f",Acum[0].pos,Acum[0].cw);

  for(int i=1;i<seqlen-len+1;i++) {

    //printf("A[%d]=%f",Acum[i].pos,Acum[i].cw);

    Acum[i].cw+=Acum[i-1].cw;
  }

}

//
// Motif::complogweights
// computes the sum of the log weights needed by the InfoPar function
//
void Motif::complogweights(double *A, double sumA, int elemno) 
{
  double norm;
  
  for(int i=0;i<elemno;i++) {
    norm=A[i]/sumA;
    if(norm!=0) sumlogweights+=norm*log(norm);
  }
}
 
//
// Motif::compWeights
// computes the weights of all positions in a sequence when a pattern motif is defined
//
void Motif::compWeights(double *A, CumWeight *Acum, char * fseq, int seqlen, int l)
{
  double Qx,Px;

  for(int i=0;i<PM[l].no;i++) {

    Qx=1;Px=1;

    for(int j=1;j<=len;j++) {
      int val=basetoint(fseq[PM[l].pos[i]+j-1],B->SType);
      if(val==-1) { Qx=0; Px=1; break;}
      Qx*=MFreq[j][val];
      Px*=MFreq[0][val];
    }
    A[i]=Qx/Px;
    if(!Px) GError("Px is 0!\n");

    Acum[i].cw = A[i];
    Acum[i].pos = i;
  }

  qsort(Acum, PM[l].no, sizeof(CumWeight), comp);

  //printf("A[%d]=%f",Acum[0].pos,Acum[0].cw);

  for(int i=1;i<PM[l].no;i++) {

    //printf("A[%d]=%f",Acum[i].pos,Acum[i].cw);

    Acum[i].cw+=Acum[i-1].cw;
  }

  

}

//
// comp 
// sort function for the CumWeight type
//
int comp(const void *a, const void *b)
{ 
  if(((CumWeight *)a)->cw < ((CumWeight *)b)->cw) return(1);
  else if (((CumWeight *)a)->cw==((CumWeight *)b)->cw) return(0);
  else return(-1); // return(0) here can cause the system to crash on Solaris

}  

//
// Motif::predUpdate
// the predictive update step in the Gibbs sampler
//
void Motif::predUpdate()
{
  F->reset();

  FastaSeq fseq;
  int i=0;

  double  *A;
  CumWeight *Acum;

  GMALLOC(A,(B->maxlen)*sizeof(double));
  GMALLOC(Acum,(B->maxlen)*sizeof(CumWeight));

  while(F->getFastaSeq(&fseq)) {

    updateMCount(Align[i],fseq.seq,B->SType,-1);
    updateBackgr(fseq.seq,fseq.len,B->SType,-1);
    updateMFreq(B->noseq-1);
    
    if(defmotif) compWeights(A,Acum,fseq.seq,fseq.len,i);
    else compWeights(A,Acum,fseq.seq,fseq.len);

    int Alen = defmotif ? PM[i].no : fseq.len-len+1;

    int r=sample(Acum,Alen);

    Align[i]= defmotif ? PM[i].pos[r] : r;

    updateMCount(Align[i],fseq.seq,B->SType,1);
    updateBackgr(fseq.seq,fseq.len,B->SType,1);
    
    i++;
  }

  GFREE(A);
  GFREE(Acum);

}

// 
// Motif::getNoseq
// a function that returns the no of sequences in the input file
//
int Motif::getNoseq()
{
  return(B->noseq);
}

//
// Motif::AlignProb
// computes the MAP for the current alignment
//
double Motif::AlignProb()
{
  double F=0;
  
  updateMFreq(B->noseq);

  for(int i=1;i<=len;i++) 
    for(int j=0;j<B->size;j++) {
      if(!MFreq[0][j]) GError("Background frequency q[0][%d] is 0!\n",j);
      if(MFreq[i][j]) F+=MCount[i][j]*log(MFreq[i][j]/MFreq[0][j]);
    }

  return(F);

}
    
//
// Motif::sample
// (weighted samples) randomly selects a postion according to a given set of weights deistribution
//
int Motif::sample(CumWeight *Acum, int Alen)
{
  int r = rand();
  
  double x = r;
  x = x/RAND_MAX * Acum[Alen-1].cw;

  //printf(" r=%d rmax=%d Acum[%d]=%f x=%f\n",r,RAND_MAX,Alen-1,Acum[Alen-1].cw,x);fflush(stdout);
  
  int i=0;
  while(x>Acum[i].cw) i++;

  if(i>=Alen) GError("Random number exceeds range!\n");
   

  return(Acum[i].pos);
}

//
// Motif::updateBackgr
// updates the Background frequencies by removing or adding one sequence
//
void Motif::updateBackgr(char *fseq,int seqlen,seqType t, int dir)
{
  for(int i=0;i<seqlen;i++) {
    int val=basetoint(fseq[i],t);
    if(val!=-1) MCount[0][val]+=dir;
  }
  sumMCount0+=dir*seqlen;
}
  
//
// Motif::updateMCount
// updates the background and motif frequencies by removing or adding a candidate 
// motif to the background
//
void Motif::updateMCount(int pos,char *fseq,seqType t, int dir)
{
  for(int j=1;j<=len;j++) {
    int val=basetoint(fseq[pos+j-1],t);
    if(val!=-1) {
      MCount[j][val]+=dir;
      MCount[0][val]-=dir;
      
      sumMCount0-=dir;
    }
  }
}

//
// Motif::testsignif
// called by runforsignif or twofilesignif to actually compute the 
// statistical significance after the test sequences have been created; 
// can be deterministic (when det=1) or not
//
void Motif::testsignif(int det,char **fkseq, int *posfreqlen,int print,int exact)
{
  FreqPos diff[B->noseq];
  unsigned long int diffno=0;
  FreqPos *pos = new FreqPos[B->noseq];


  if(det) {

    int maxnoinseq=computemaxPMno();
    
    if(!exact) {
      while(maxnoinseq>1) {
	//fprintf(stderr,"Test: Max no of elements in seqs: %d           \r",maxnoinseq);
	computeQq(fkseq);
	removeoneword(fkseq);
	maxnoinseq--;
      }
      fprintf(stderr,"\n");
    }
    computeQq(fkseq);

    double TRposall=0;
    double TRnegall=0;

    for(int i=0;i<B->noseq;i++) {
      double TRpos=-HUGE_VAL;
      double TRneg=-HUGE_VAL;
      int fkseqlen=strlen(fkseq[i]);
      if(exact) {
	int j=0;
	TRpos=0;
	diff[diffno].pos=0;
	while(j<PM[i].no && PM[i].pos[j]<posfreqlen[i]-len) {
	  TRpos+=computeTR(fkseq[i],PM[i].pos[j]);
	  diff[diffno].pos=PM[i].pos[j];
	  j++;
	}
	while(j<PM[i].no && PM[i].pos[j]<posfreqlen[i]) { 
	  j++;
	}
	TRneg=0;
	while(j<PM[i].no && PM[i].pos[j]<fkseqlen-len) {
	  TRneg+=computeTR(fkseq[i],PM[i].pos[j]);
	  j++;
	}
      }
      else {
	for(int j=0;j<posfreqlen[i]-len;j++) {
	  double TR=computeTR(fkseq[i],j);
	  if(TR>TRpos) { 
	    TRpos=TR; 
	    //diff[i].pos=j; 
	    diff[diffno].pos=j; 
	  }
	}
	for(int j=posfreqlen[i];j<fkseqlen-len;j++) {
	  double TR=computeTR(fkseq[i],j);
	  if(TR>TRneg) { TRneg=TR; }
	}
      }
      if(exact) {
	TRposall+=TRpos;
	TRnegall+=TRneg;
      }
      
      //if(TRpos!=TRneg) {
	diff[diffno++].freq=TRpos-TRneg;
	//}
	
	//fprintf(stderr,"%d: TRpos=%f TRneg=%f diff=%f\n",i,TRpos, TRneg,diff[i].freq);
    }

    if(exact) { 
      fprintf(B->outf,"TR: %f %f\n",TRposall,TRnegall);
      if(TRposall==0 && TRnegall==0) { fprintf(stderr,"Motif not found!\n");exit(0);}
    }

    if(!print) fprintf(B->outf," ");
    //    diffno=B->noseq;
  }
  else { 
    FreqPos neg[B->noseq];

    for(int i=0; i<B->noseq; i++) {
      pos[i].freq=0;
      neg[i].freq=0;
      pos[i].pos=0;
      neg[i].pos=0;
      
      int motiffreqlen=strlen(fkseq[i]);

      for(int j=0; j<posfreqlen[i]; j++) {
	if(motiffreq[i][j]>pos[i].freq) {
	  pos[i].freq=motiffreq[i][j];
	  pos[i].pos=j;
	}
      }

      for(int j=posfreqlen[i]; j<motiffreqlen; j++) {
	if(motiffreq[i][j]>neg[i].freq) {
	  neg[i].freq=motiffreq[i][j];
	  neg[i].pos=j;
	}
      }

      if(pos[i].freq != neg[i].freq) {
	diff[diffno].freq=pos[i].freq-neg[i].freq;
	diff[diffno++].pos=pos[i].pos;
      }

    }
  }
  if(print) fprintf(B->outf,"**********************\n\nSignificance of motif:\n");

  //if(print) apply_test_student(B->outf,diffno, diff);
  apply_test_student(B->outf,diffno, diff);
  apply_test_wilcoxon_pair(B->outf,diffno,diff,print);

  if(!det) motifprob=pos;
  if(print) fprintf(B->outf,"\n**********************\n\n\n");
}

//
// Motif::GibbsSample
// runs Gibbs Sampler in a limited way - just sampling according to the pre-computed motif
// ; needed in computing the significance of a motif
//
void Motif::GibbsSample(int SignfRunNo, char **fkseq,int print)
{
  double  *A;
  CumWeight *Acum;
  
  GMALLOC(A,(B->maxlen)*sizeof(double));
  GMALLOC(Acum,(B->maxlen)*sizeof(CumWeight));

  for(int j=0;j<SignfRunNo;j++) {  

    if(print) fprintf(stderr, "Finding significance of motif: %d/%d                  \r",j+1,SignfRunNo);
    
    for(int i=0;i<B->noseq; i++) {
      
      // this variant makes continous updates to motifs

      //updateMCount(Align[i],fkseq[i],B->SType,-1);
      //int fkseqlen=strlen(fkseq[i]);
      //updateBackgr(fkseq[i],fkseqlen,B->SType,-1);
      //updateMFreq(B->noseq-1);
    
      //compWeights(A,Acum,fkseq[i],fkseqlen);
      //int Alen=fkseqlen-len+1;
      //int Alen=nearoptim(Acum,fkseqlen-len+1,nearoptperc);  

      //int r=sample(Acum,Alen);

      //Align[i]= r;
      //motiffreq[i][r]++;
    
      //updateMCount(Align[i],fkseq[i],B->SType,1);
      //updateBackgr(fkseq[i],fkseqlen,B->SType,1);

      
      // this variant is update free

      int fkseqlen=strlen(fkseq[i]);
      
      compWeights(A,Acum,fkseq[i],fkseqlen/2);
      int Alen=fkseqlen/2-len+1;
      int r=sample(Acum,Alen);
      motiffreq[i][r]++;
      compWeights(A,Acum,fkseq[i]+fkseqlen/2,fkseqlen/2);
      r=sample(Acum,Alen);
      motiffreq[i][r+fkseqlen/2]++;
    }

  }
  fprintf(stderr,"\n");
  GFREE(A);
  GFREE(Acum);
}

//
// basetoint
// transform a residue to its integer equivalent
//
int basetoint(char c,seqType t) 
{
  int j;

  if(t==nucl) {
    switch(c) {
    case 'a':
    case 'A': return(0);
    case 'c':
    case 'C': return(1);
    case 'g':
    case 'G': return(2);
    case 't': 
    case 'u': 
    case 'U':
    case 'T': return(3);
    case 'm':
    case 'M': return(0);
      j=(int) (rand()/(RAND_MAX+1.0));
      if(j) return(0);
      else return(1);
    case 'r':
    case 'R': return(0);
      j=(int) (rand()/(RAND_MAX+1.0));
      if(j) return(0);
      else return(2);
    case 'w':
    case 'W': return(0);
      j=(int) (rand()/(RAND_MAX+1.0));
      if(j) return(0);
      else return(3);
    case 's':
    case 'S':return(1);
      j=(int) (rand()/(RAND_MAX+1.0));
      if(j) return(1);
      else return(2);
    case 'y':
    case 'Y':return(1);
      j=(int) (rand()/(RAND_MAX+1.0));
      if(j) return(1);
      else return(3);
    case 'k':
    case 'K':return(2);
      j=(int) (rand()/(RAND_MAX+1.0));
      if(j) return(2);
      else return(3);
    case 'v':
    case 'V':return(0);
      j=(int) (2.0*rand()/(RAND_MAX+1.0));
      switch(j) {
      case 0: return(0);
      case 1: return(1);
      case 2: return(2);
      }
    case 'h':
    case 'H':return(0);
      j=(int) (2.0*rand()/(RAND_MAX+1.0));
      switch(j) {
      case 0: return(0);
      case 1: return(1);
      case 2: return(3);
      }
    case 'd':
    case 'D':return(0);
      j=(int) (2.0*rand()/(RAND_MAX+1.0));
      switch(j) {
      case 0: return(0);
      case 1: return(2);
      case 2: return(3);
      }
    case 'b':
    case 'B':return(1);
      j=(int) (2.0*rand()/(RAND_MAX+1.0));
      switch(j) {
      case 0: return(1);
      case 1: return(2);
      case 2: return(3);
      }
    case 'x':
    case 'X':
    case 'n':
    case 'N':return(-1);
      j=(int) (3.0*rand()/(RAND_MAX+1.0));
      switch(j) {
      case 0: return(0);
      case 1: return(1);
      case 2: return(2);
      case 3: return(3);
      }
    default: GError("Unrecognized character: %c\n",c);
    }
  }

  if(t==aac) {
    switch(c) {
    case 'a': case 'A': return(0);
    case 'b': case 'B': return(1);
    case 'c': case 'C': return(2);
    case 'd': case 'D': return(3);
    case 'e': case 'E': return(4);
    case 'f': case 'F': return(5);
    case 'g': case 'G': return(6);
    case 'h': case 'H': return(7);
    case 'i': case 'I': return(8);
    case 'k': case 'K': return(9);
    case 'l': case 'L': return(10);
    case 'm': case 'M': return(11);
    case 'n': case 'N': return(12);
    case 'p': case 'P': return(13);
    case 'q': case 'Q': return(14);
    case 'r': case 'R': return(15);
    case 's': case 'S': return(16);
    case 't': case 'T': return(17);
    case 'v': case 'V': return(18);
    case 'w': case 'W': return(19);
    case 'y': case 'Y': return(20);
    case 'z': case 'Z': return(21);
    case  'x': case 'X': return(-1);
      j=(int) (21.0*rand()/(RAND_MAX+1.0));
      return(j);
    default: GError("Unrecognized character: %c\n",c);
    }
  }

  return(-1);
}

// 
// inttobase
// transforms an integer to its equivalent residue
//
char inttobase(int b,seqType t)
{
  if(t==nucl) {
    switch(b) {
    case 0: return('a');
    case 1: return('c');
    case 2: return('g');
    case 3: return('t');
    default: return('n');
    }
  }
  if(t==aac) {
    switch(b) {
      case 0: return('A');
      case 1: return('B');
      case 2: return('C');
      case 3: return('D');
      case 4: return('E');
      case 5: return('F');
      case 6: return('G');
      case 7: return('H');
      case 8: return('I');
      case 9: return('K');
      case 10: return('L');
      case 11: return('M');
      case 12: return('N');
      case 13: return('P');
      case 14: return('Q');
      case 15: return('R');
      case 16: return('S');
      case 17: return('T');
      case 18: return('V');
      case 19: return('W');
      case 20: return('Y');
      case 21: return('Z');
    default: return('X');
    }
  }
  return(-1);
}

//
// strdist
// computes the edit distance between two strings of equal length
//
int strdist(char *s1,char *s2, int len)
{
  int d=0;
  for(int i=0;i<len;i++) 
    d += (toupper(s1[i]) == toupper(s2[i])) ? 0 : 1;
  return d;
} 

//
// makeargv
// tokenize an input string according to the specified delimiters
//
int makeargv(char *s, char *delimiters, char ***argvp)
{
  char *t, *snew;
  int numtokens,i;

  snew=s+strspn(s,delimiters);
  if((t=(char *)calloc(strlen(snew)+1,sizeof(char)))==NULL)
    {
      *argvp=NULL;
      numtokens=-1;
    }
  else
    {
      strcpy(t,snew);
      if(strtok(t,delimiters)==NULL)
	numtokens=0;
      else
	for(numtokens=1;strtok(NULL,delimiters)!=NULL;numtokens++) ;
      if((*argvp=(char **)calloc(numtokens+1,sizeof(char *)))==NULL)
	{
	  free(t);
	  numtokens=-1;
	}
      else
	{
	  if(numtokens>0)
	    {
	      strcpy(t,snew);
	      **argvp=strtok(t,delimiters);
	      for(i=1;i<numtokens+1;i++)
		*((*argvp)+i)=strtok(NULL,delimiters);
	    }
	  else
	    {
	      **argvp=NULL;
	      free(t);
	    }
	}
    }

  return numtokens;
}

//
// Motif::findMotif
// search for motif in a limited space (when defmotif is 1) or in all space
// if det is set to 1 then the motif is computed deterministically; otherwise
// the Gibbs sampling method is used
//
double Motif::findMotif(int iter_no, int max_loop, int inlocmax, int init, int det)
{
  double globMaxProb=0; 
  int *MaxAlign;

  if(det) {
    globMaxProb=findMDet(0);
    MaxAlign=Align;
  }
  else {
    
    MaxAlign = new int[B->noseq];  
    int *locAlign;
 
    locAlign=new int[B->noseq];

    int i,j;
    
    for(i=1; i<=iter_no; i++) {

      //printf("**** %d ****\n",i);

      if(init) initRandAlign();

      double locMaxProb=0;
      double prevProb=0;
      int not_loc_max=inlocmax;
      double locMaxInfo=0;

      j=0;

      
      while(not_loc_max && j<max_loop) {

	predUpdate();

	double alignProb=AlignProb(); 
	double infopar=InfoPar(alignProb);

	//fprintf(stderr,"[%d]%d: alignProb=%.2f Info/param=%.2f    \r",i,j,alignProb,infopar);
			
	if(alignProb>locMaxProb) {
	  fprintf(stderr,"[%d]%d: alignProb=%.2f Info/param=%.2f diff=%.2f    \r",i,j,alignProb,infopar,alignProb-locMaxProb);
	  //fprintf(stderr,"[%d]%d: alignProb=%.2f diff=%.2f        \r",i,j,alignProb,alignProb-locMaxProb);
	  locMaxProb=alignProb;
	  locMaxInfo=infopar;
	  copyarray(locAlign,Align,B->noseq);
	  not_loc_max=inlocmax;
	}
	else /*if(alignProb == prevProb)*/ not_loc_max--;
	
	prevProb=alignProb;

	j++;
      }

      fprintf(stderr,"\n");
      //fprintf(stderr,"[%d]: alignProb=%.2f Info/param=%.2f                        \n",i,locMaxProb,locMaxInfo);
      //fprintf(stderr,"[%d]: alignProb=%.2f                         \n",i,locMaxProb);

      //printf("Iter no: %d maxprob=%f\n",i,globMaxProb);fflush(stdout);

      if(locMaxProb == globMaxProb) 
	break;
      else if(locMaxProb > globMaxProb) {
	globMaxProb = locMaxProb;
	copyarray(MaxAlign,locAlign,B->noseq);
      }
      
      
    }
  }

  computeMotif(MaxAlign);
  return(globMaxProb);

}
