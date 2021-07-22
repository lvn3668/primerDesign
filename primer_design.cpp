// Author: Lalitha Viswanathan
// Org: Tata Consultancy Services
//read DNA sequence
//search for percentage GC content
//serach for MT compuation
// Find putative primer targets based on GC Content,  melting temperature and max length of primer / dist from target seq 
#include "primer_design.h"
#include "pal.cpp"
#include "readDNA.cpp"
// Reads in a file with GC Min max cutoffs, MT min max cutoffs 
int main(int argc, char* argv[])
{
		if(argc!=11)
		{
		printf("<File contaiing DNA sequence in which probe is required> <min length of primer>" + 
			"<max length of primer> <distance from target sequence> <min cutoff GC content> <max cutoff GC content>"+ 
			"<min cutoff melting temperature> <max cutoff melting temperature> <from--pos> <to-pos>\n");
		exit(1);
		}
		char* fname=argv[1];
		min_length_of_primer=atoi(argv[2]);
		max_length_of_primer=atoi(argv[3]);
		int dist=atoi(argv[4]);
		cutoff_GC_min=atoi(argv[5]);
		cutoff_GC_max=atoi(argv[6]);
		cutoff_MT_min=atoi(argv[7]);
		cutoff_MT_max=atoi(argv[8]);
		int from=atoi(argv[9]);
		int to=atoi(argv[10]);
		readDNA(fname);
		// Check if required primer length is greater than length of sequence
		if((min_length_of_primer >len)||(max_length_of_primer > len))
		{
				printf("length of primer required is greater than length of sequence\n");
				exit(1);	
		}
		if((from-dist <0)||(from <0))
		{
				printf("cannot find probe at %d distance from target sequence\n",from-dist);
				exit(1);
		}
		if(to+dist > len) //takes care of to being greater than length as well
		{
				printf("cannot find probe at %d distance from target sequence. %d is greater than length of sequence %d\n",to+dist,len);
				exit(1);
		}
		if((to < 0)|| (from<0))
		{
				printf("invalid entry\n");
				exit(1);
		}
		// Cutoff GC min max 
		if((cutoff_GC_min < 0) || (cutoff_MT_min < 0))
		{
				printf("cutoff values cannot be less than zero\n");
				exit(1);
		}
		find_probes(from,to,dist);
	printf("after finding probes!!!\n");
	printf("Enter 2 sequences for checking for complementality\n");
	char fwd_primer[100],rev_primer[100];
	scanf("%s %s",&fwd_primer,&rev_primer);
	if(check_complementality(fwd_primer,rev_primer))
	{
			printf("This pair cannot be used as a valid primer pair\n");
	}
	else
	{
			//valid primer pair
			printf("This can be used as a valid primer pair\n");
	}
}
/// pal.cpp

char* getstr(int start,int end,char* string,char* target)
{
		int q=0,counter;
		for(counter=start;counter<=end;counter++)
		{
		    string[q]=target[counter];
	  	    q++;		
		}
		string[q]='\0';
		return string;
}	
// Find probes
void find_probes(int from, int to,int dist)
{
		float no_AT, no_GC,no_ATrev,no_GCrev;
		float MT, MTrev;
		int loc;
		FILE* fp=fopen("fwd_revprimers.txt","w");
for(int length_of_primer=min_length_of_primer;length_of_primer<=max_length_of_primer;length_of_primer++)
		{
		fprintf(fp,"\nFORWARD PRIMERS\n");
				loc=from-length_of_primer-dist;
				//a is at maximum distance from "from"; i
				//ncreasing "a" reduces distance from  "from"
				//for a given length of primer, searh for primers at varying distances from target location
				//printf("forward primer finding\n");
					for(int a=loc;a+length_of_primer<from;a++)
					{
				//printf("a=%d from=%d a+length_of_primer=%d length_of_primer=%d\n",a,from,a+length_of_primer,length_of_primer);
				//get the forward primer
				
						char* sequence=(char*)malloc(sizeof(char)*(length_of_primer+2));
						sequence=getstr(loc,loc+length_of_primer,sequence,seq);
						no_AT=0; no_GC=0;
						for(int counter=0;counter<strlen(sequence);counter++)
						{
								switch(sequence[counter])
								{
									case 'A':
									case 'T':	   	
											no_AT++; break;
									case 'G':
									case 'C':
											no_GC++;break;		
								}
						}
		  	    //compute MT
				MT=4*(no_GC)+2*(no_AT);
				//MT is greater than cutoff and GC content is greater than cutoff and there are no inverted repeats, then its a forward primer
				//printf("GC percentage is %f MT is %f %d\n",(no_GC/length_of_primer)*100,MT,pal(sequence,length_of_primer));
				float GC_percent=(no_GC/length_of_primer)*100;
				if( (GC_percent>=cutoff_GC_min) && (GC_percent<=cutoff_GC_max) && (MT>=cutoff_MT_min)&&(MT<cutoff_MT_max))
				{
						if(!pal(sequence,length_of_primer))
						{
								fprintf(fp,"\n%s length= %d GC%=%5.3f Melt. Temp.=%5.3f dist. from target=%d",sequence,length_of_primer,(no_GC/length_of_primer)*100,MT,from-a-length_of_primer);
						}
				}
					}
						fprintf(fp,"\n\nREVERSE PRIMERS\n");
						//search for reverse primer
						int loc1=to+length_of_primer+dist;
						for(int p=loc1;p>to+length_of_primer;p--)
						{
								//get complementary sequence
								char* sequence1=(char *)malloc(sizeof(char)*(length_of_primer+3));
								sequence1=getstr(p-length_of_primer,p,sequence1,seq);
								//complement the sequence
								for(int c=0;c<strlen(sequence1);c++)
								{
								if(sequence1[c]=='A')
									sequence1[c]='T';
								else if(sequence1[c]=='T')
							 	    sequence1[c]='A';
								else if(sequence1[c]=='G')
									sequence1[c]='C';
								else
									sequence1[c]='G';	
								}		
								no_ATrev=0;no_GCrev=0;
								for(int b=0;b<strlen(sequence1);b++)
								{
										switch(sequence1[b])
										{
												case 'A':
												case 'T':
													no_ATrev++;
													break;
												case 'G':
												case 'C':
													no_GCrev++;	
													break;
										}
								}
							MTrev=4*(no_GCrev)+2*(no_ATrev);
				//reverse n print the complementary sequence;
							if((((float)(no_GCrev/length_of_primer)*100)>=cutoff_GC_min) &&(((float)(no_GCrev/length_of_primer)*100)<=cutoff_GC_max)&& (MTrev>=cutoff_MT_min)&&(MTrev<=cutoff_MT_max))
						{
								if(!pal(sequence1,length_of_primer))
									fprintf(fp,"\n%s length=%d GC%=%5.3f Melt. Temp.=%5.3f Dist. from target=%d",reverse(sequence1),length_of_primer,(no_GCrev/length_of_primer)*100,MTrev,p-to-length_of_primer);
						}
						
						free(sequence1);
					}
		}
}
//reverse sequence
char* reverse(char* pointertosequence)
{
		char* temp=(char *)malloc(sizeof(char)*strlen(pointertosequence));
		int t=0,n=0;
		for(t=strlen(pointertosequence)-1,n=0;t>=0;t--)
		{
				temp[n++]=pointertosequence[t];
		}
		temp[n]='\0';
	return temp;
}
// check complementality
int check_complementality(char* sequence1, char* sequence2)
{
//get the complement of the string b
		// if string b is smaller than string a ie. reverse primer is smaller than forward
		// primer, then can compare only strlen(reverse primer) characters
		int cutoff=85,length_count;
		char* temp_str=(char*)malloc(sizeof(char)*strlen(sequence2));
		int temp;
		for(temp=0;temp<strlen(sequence2);temp++)
				{
						switch(sequence2[temp])
						{
								case 'A':
										temp_str[temp]='T';
										break;
								case 'T':
										temp_str[temp]='A';
										break;
								case 'G':
										temp_str[temp]='C';
										break;
								case 'C':
										temp_str[temp]='G';
										break;
						}
				}
		temp_str[temp]='\0';
		//check for complementality
		// a is forward primer
		// temp_str is complement of reverse primer b
		// reverse primer is stored in 3' to 5' format w/o the complementality done
		// temp_str is complementary seq 
		// since it is also in 3' to 5' format, the first strlen(a) characters of the sequences need
		// to be checked for MATCHING
		// i.e. if they match, then complement of temp_str(reverse primer) is complementary
		// to forward primer
		// else if there r enuff mismatches then, complement of temp_str i.e. the reverse primer 
		// is not complement of forward primer
if(strlen(sequence1)<strlen(sequence2))
		length_count=strlen(sequence1);
else
		length_count=strlen(sequence2);
int mat=0, mismat=0;
printf("%s %s %d \n",sequence1,temp_str,length_count);
for(int counter=0;counter<length_count;counter++)
{
		if(sequence1[counter]==temp_str[counter])
				mat++;
		else
				mismat++;
}
// if more matches are there, then complementality is higher 
// else lower complementality
		if(((mat*2-(mismat*1))/length_count)*100>cutoff)
		{
			return 1;
		}
		return 0;
}
