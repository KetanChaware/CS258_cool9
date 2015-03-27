#include<stdio.h>
#include<omp.h>
#include<string.h>
#define threads_num 4
#define var_num 8

struct variable
{
    char* name;
    double values[3];
};

//Struct for storing first two max looping variable values
struct ip_var
{
    char* name;
    double value;
};


int main(int argc,char* argv[])
{
    int i,j,threads_obtained;
    struct variable input[var_num];
    
    char* X = argv[2];

	omp_set_num_threads(threads_num);
    
/*****************************************************************************************
Section for Reading INPUT FILE
*****************************************************************************************/                                                          

    FILE* input_file = fopen(argv[1],"r");

    for(i=0;i<var_num;i++)
    {
        fscanf(input_file,"%ms",&input[i].name);
        for(j=0;j<3;j++)
        {
            fscanf(input_file,"%lf",&input[i].values[j]);
        }
    }

    
///////////////////////////////////////////////////////////////////////////////////////////    
    
    
/******************************************************************************************
Section for finding first two maximum looping variables and storing them in struct ip_var
******************************************************************************************/

    int max1=0,max2=1; //Variables to store Indices for first two max looping variables  

    for(i=1 ; i<var_num; i++)
    {
        if((input[i].values[1] - input[i].values[0])/input[i].values[2]  >  (input[max1].values[1] - input[max1].values[0])/input[max1].values[2])
        {
            max2 = max1;
            max1 = i;
        }
    }

    int loop1 = (input[max1].values[1] - input[max1].values[0])/input[max1].values[2];
    int loop2 = (input[max2].values[1] - input[max2].values[0])/input[max2].values[2] ? (input[max2].values[1] - input[max2].values[0])/input[max2].values[2]:1;

    struct variable temp;
    temp = input[0];
    input[0] = input[max1];
    input[max1] = temp;
    temp = input[1];
    input[1] = input[max2];
    input[max2] = temp;

    struct ip_var ip_file[(loop1+1)*(loop2+1)][2];

    int k=0;
    double x,y;
    
    for(x=input[0].values[0]; x<=input[0].values[1]; x += input[0].values[2])
    {
        for(y=input[1].values[0]; y<=input[1].values[1]; y += input[1].values[2])
        {
            ip_file[k][0].name = input[0].name;
            ip_file[k][0].value = x;
            ip_file[k][1].name = input[1].name;
            ip_file[k++][1].value = y;
        }
    }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	

/****************************************************************************************************************
Parallel Section which creates number of threads and then each thread calls main code of micrOMEGAs.
Each thread has its own file for input and output and loops modify input file and output is appended
in output file everytime.
*****************************************************************************************************************/
    
    #pragma omp parallel
    {
        int pf;
        double pi[var_num-1];
        threads_obtained = omp_get_num_threads();
        int id = omp_get_thread_num();
        FILE* infiles;
        FILE* outfiles;
        char fname[16];
        char foname[16];

        for(pf=id;  pf<k; pf+=threads_obtained)
        {
            for(pi[1] = input[2].values[0]; pi[1] <= input[2].values[1]; pi[1] += input[2].values[2])
            {
                for(pi[2] = input[3].values[0]; pi[2] <= input[3].values[1]; pi[2] += input[3].values[2])
                {
                    for(pi[3] = input[4].values[0]; pi[3] <= input[4].values[1]; pi[3] += input[4].values[2])
                    {
                        for(pi[4] = input[5].values[0]; pi[4] <= input[5].values[1]; pi[4] += input[5].values[2])
                        {
                            for(pi[5] = input[6].values[0]; pi[5] <= input[6].values[1]; pi[5] += input[6].values[2])
                            {
                            	for(pi[6] = input[7].values[0]; pi[6] <= input[7].values[1]; pi[6] += input[7].values[2] )
                            	{
                            		sprintf(fname,"iFiles%d.txt",id);
                                	sprintf(foname,"oFiles%d.txt",id);
                                	outfiles = fopen(foname,"a");
                                	
                                	fprintf(outfiles,"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,",ip_file[pf][0].value,ip_file[pf][1].value,pi[1],pi[2],pi[3],pi[4],pi[5],pi[6]);
                                	fclose(outfiles);
                                	
                                	infiles = fopen(fname,"w");

                                	fprintf(infiles,"%s  %lf \n",input[0].name,ip_file[pf][0].value);
                                	fprintf(infiles,"%s  %lf \n",input[1].name,ip_file[pf][1].value);
                                	fprintf(infiles,"%s  %lf \n",input[2].name,pi[1]);
                                	fprintf(infiles,"%s  %lf \n",input[3].name,pi[2]);
                                	fprintf(infiles,"%s  %lf \n",input[4].name,pi[3]);
                                	fprintf(infiles,"%s  %lf \n",input[5].name,pi[4]);
                                	fprintf(infiles,"%s  %lf \n",input[6].name,pi[5]);
                                	fprintf(infiles,"%s  %lf \n",input[7].name,pi[6]);                                

                                	fclose(infiles);

                                	char s[128];

                                	sprintf(s,"./main %s %s>> %s",fname,X,foname);

                                	system(s);
                            	}
                            }
                        }
                    }
                }
            }
        }
    }
    
//////////////////////////////////////////////////////////////////////////////////////////////////////



/*****************************************************************************************************
Section for combining all output files corresponding to indivisual thread into a single ".csv" file
which can be opened in Gnumeric or Excel or any such program.
*****************************************************************************************************/

    FILE* output = fopen("output_final.csv","w");

    fprintf(output,"%s,%s,%s,%s,%s,%s,%s,%s,ERROR\n",input[0].name,input[1].name,input[2].name,input[3].name,input[4].name,input[5].name,input[6].name,input[7].name);

    fclose(output);
    char cmd[1024]="";
    char temp1[16];

    for(i=0; i<threads_obtained; i++)
    {
        sprintf(temp1,"oFiles%d.txt ",i);
        strcat(cmd,temp1);
    }

    char command[1024];

    sprintf(command,"cat %s >> output_final.csv",cmd);

    system(command);
    
    for(i=0; i<threads_obtained; i++)
    {
        sprintf(temp1,"iFiles%d.txt ",i);
        strcat(cmd,temp1);
    }
    
    sprintf(command,"rm -f %s",cmd);
    system(command);
    
// End of Program
    
////////////////////////////////////////////////////////////////////////////////////////////////////

}
