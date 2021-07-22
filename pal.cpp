
int pal(char* temp,int length_of_primer)
{
		long int nomisinv;
	float scoreinv;
		char* str1,*str2;
		int flag=0;
		minlength=(min_percent/100)*length_of_primer;
		maxlength=(max_percent/100)*length_of_primer;
		gaplength=(gap_percent/100)*length_of_primer;
		
		for(long int a=0;a<strlen(temp)-minlength;a++)
		    for(long int length=minlength;length<=maxlength;length++)
		        {
					 //get 1st string
					 str1=(char *)malloc(sizeof(char)*(length+1));
					 str1=getstr(a,a+length-1,str1,temp);
					 float cutoff=(length*match_score)-(mismatches*mismatch_score);
						for(long int gapl=0;gapl<=gaplength;gapl++)
						{
							long int b=a+length+gapl;
							//get second string
					        str2=(char *)malloc(sizeof(char)*(length+1));
						    str2=getstr(b,b+length-1,str2,temp);
							scoreinv=0,nomisinv=0;
							//search for inverted repeats
							for(long int q=0,r=length-1;q<length,r>=0;q++,r--)
							{
									 if(str1[q]==str2[r])
									 {}
									 else
									 {nomisinv++;}
							}
							scoreinv=(length*match_score)-(nomisinv*mismatch_score);	
							//if score is greater than or equal to cutoff,then string is palindrome 
							if(scoreinv>cutoff)
									{
									flag=1;		
									break;	
									}
						free(str2);
						}
						if(flag)
								return 1;
					free(str1);
				}
		return flag;
}

