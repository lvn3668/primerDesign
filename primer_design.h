#include<iostream.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<stdlib.h>
#include<ctype.h>
#define min_percent 10
#define max_percent 50
#define gap_percent 20
#define mismatches 0
#define mismatch_score -1
#define match_score 2
#define complement(x) (x) == 'A' ? (x)='T' : (x)=='T' ? (x) = 'A' : (x)=='G'? (x)='C' : (x) =='C' ? (x)='G'

long int len;
char* seq;
long int minlength, maxlength,gaplength;
float cutoff_GC_min, cutoff_GC_max,cutoff_MT_min,cutoff_MT_max;
int max_length_of_primer,min_length_of_primer;
char* getstr(int,int,char*,char*);
int pal(char*, int);
void find_probes(int,int,int);
char* reverse(char*);
int check_complementality(char*,char*);
