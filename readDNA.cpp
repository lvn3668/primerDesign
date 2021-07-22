// Author Lalitha Viswanathan
// Org Tata Consultancy Services
//assuming that the file contains only 1 sequence for which probes are to be designed
void readDNA(char *filename)
{
	  FILE* fp;
	  char str[500],ch;
	  long int size=1,i;
	  if(!(fp=fopen(filename,"r")))
	  {
		  printf("Filename %s does not exist\n",filename);
		  exit(1);
	  }
	  fseek(fp,0L,SEEK_END);
	  long int temp=ftell(fp);
	  seq=(char *)(malloc(sizeof(char)*temp));
	  fseek(fp,0L,SEEK_SET);
	 //read the first line containing ">" 
	  fgets(str,500,fp);
	  strcpy(seq,"");
	  while(!feof(fp))
	  {	
			  	 ch=fgetc(fp);
   				 if((ch=='A')||(ch=='T')||(ch=='U')||(ch=='C')||(ch=='G'))
				 seq[size++]=ch;
	  }
	  len=size;
	  //for(int a=0;a<len;a++)
	//	printf("%c",seq[a]);
	  seq[len]='\0';
	//  printf("%c %c %c is the sequence",seq[len-2],seq[len-1],seq[len]);
	  len=size;
	//exit(1);
}
