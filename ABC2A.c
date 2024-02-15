/**************************
**Program:	ABC2A.c
**Description:	Reconstruct the all-atom structure from 3-bead CG conformation(s)
**Author:	YZ Shi & Hao Wu
**Date:		2023.2
**Update_version: v1.0(2023.3) for multiple fragment
**                v2.0(2023.8) bond check
**                v3.0(2023.10) clash remove
**Compile:	gcc -Wall -o3 ABC2A.c -o ABC2A -lm
**Run:		./ABC2A &pdb_name(e.g., 1zih_CG) &fragment_path
***************************/
#include<stdio.h>
#include<math.h>
#include<dirent.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<time.h>
#include<string.h>
#include<stdlib.h>
#define SIZE 14           //Max_size of str 
#define Max_atom 10000    //Max_num of atoms
#define Max_frag_num 12
FILE *out_AA,*out_rmsd;
char	frag_path[60];
//char	frag_path[60]="fragment/";
int	N_frag[5]={0,6,6,6,6};  //[base_type:XAUGC-01234] 
//int	N_frag[5]={0,1,1,1,1};
char	atom_name[Max_atom][SIZE],res_name[Max_atom][SIZE],chain[Max_atom][SIZE];	
float	x_frag[5][Max_frag_num][5][50];  //[base_type:AUGC-1234][frag_no][xyz-123][N]
int	N_frag_atom[5][Max_frag_num];
char	atom_frag[5][Max_frag_num][50][SIZE];
char	res_frag[5][Max_frag_num][50][SIZE];
float	x_frag_CG[5][Max_frag_num][6][5];   //[base_type:AUGC-1234][frag_no][CPCNPC-012345][xyz-123]
char 	*str1="ATOM",*str2="C4'",*str3="P",*str4="A",*str5="G",*str6="C",*str7="U",*str8="N1",*str9="N9";
int	NN[Max_atom/3][50],NN_num[Max_atom/3];  //Neighbour of each nucleotide
int main(int argc, char *argv[])
{
   clock_t start,end;
   start=clock();
   if (argc<3) {
      printf("You must run like: ./ABC2A $pdb_name $fragment_path (e.g., ./ABC2A 1zih_CG fragment/)\n");}
   sprintf(frag_path,"%s",argv[2]);
   char  pdb_file[20],AA_PDB[30]; 
   sprintf(pdb_file,"%s.pdb",argv[1]);   printf("pdb_file: %s\n",pdb_file);
   int	Read_conf_PDB();
   void	Read_frag_AUGC(),Assemble();
   float	x[5][Max_atom];
   int N0=Read_conf_PDB(pdb_file,x);  //Read CG pdb file
   Read_frag_AUGC(); 	      //Read fragment
   
   sprintf(AA_PDB,"%s_AA.pdb",argv[1]);
   out_AA = fopen(AA_PDB, "w+");  
   char rmsd_file[20];
   strcpy(rmsd_file,argv[1]);
   strcat(rmsd_file,"_rmsd.txt");
   out_rmsd=fopen(rmsd_file,"w+");   
   Assemble(N0,x);
   end = clock();
   float spendtime = (float)(end-start)/(CLOCKS_PER_SEC);
   printf("Finished with the spendtime of %f s.\n",spendtime);
   fprintf(out_rmsd,"#Time(s):	%.5f",spendtime);
   fclose(out_AA); fclose(out_rmsd);
   return 0;
}
/*********************************************************************************/
void Parse_aline(char str[100],char str_temp[15][SIZE])
{
   const char sep[2]=" ";
   char *token;
   int i=1;
   token = strtok(str,sep);
   while(token != NULL) {
      strcpy(str_temp[i],token);
      i++;
      token = strtok(NULL,sep);
   } //separate the string str into str_temp by space 
}
float Cal_distance(float x0[5],float x1[5])
{
   float dist=0.0;
   for (int i=1;i<=3;i++) {dist += pow(x0[i]-x1[i],2);}
   return sqrt(dist);
}
/*************************Get Neighbour**************************/
void Get_neighbour(char atom[Max_atom][SIZE],float x0[5][Max_atom],int N)
{
   int k;
   float cutoff=12.5;
   for (int i=N-1;i>3;i--) { 
       if (strcmp(atom[i],"C4'")==0) { k=1; 
          for (int j=2;j<=i-3;j++) {
              if (strcmp(atom[j],"C4'")==0) { 
                 float dist=0.0;
                 for (int l=1;l<=3;l++) {dist += pow(x0[l][i]-x0[l][j],2);}
                 if (sqrt(dist)<=cutoff) { 
                    NN[(i+3)/3][k]=(j+1)/3; //printf("i = %d --> (%d) %d %f\n",(i+1)/3,k,NN[(i+3)/3][k],sqrt(dist));
                    k++;}}
            } 
          NN_num[(i+3)/3]=k-1;
       } }
}
/*****************************Read_Native_PDB_file**************************/
int Read_conf_PDB(char file_name[20],float x[5][Max_atom])
{    
   FILE *inpdb;     
   inpdb=fopen(file_name,"r+");   
   char  atom0[SIZE],atom_name0[Max_atom][SIZE],res_name0[Max_atom][SIZE],chain0[Max_atom][SIZE];
   float x0[Max_atom],y0[Max_atom],z0[Max_atom]; 
  // int  res_id0[Max_atom];
   char aline[100],str_line[15][SIZE];   
   int i=1;
   while(!feof(inpdb)) {
       fgets(aline,100,inpdb);
       strncpy(atom0,aline,4);  atom0[4]='\0';
       if (strcmp(atom0,str1)!=0)  	continue;
       Parse_aline(aline,str_line);
       if (strcmp(atom0,str_line[1])!=0) {
          printf("Warning! The Parse_aline would be not right!\n"); exit(0);}
       if (strcmp(str_line[4],str4)!=0&&strcmp(str_line[4],str5)!=0
         &&strcmp(str_line[4],str6)!=0&&strcmp(str_line[4],str7)!=0) continue;  
       if (i==1&&strcmp(str_line[3],str3)!=0) continue;   //avoid that the first is not P   
       strcpy(atom_name0[i],str_line[3]); strcpy(res_name0[i],str_line[4]);
       strcpy(chain0[i],str_line[5]);     //res_id0[i]=atoi(str_line[6]);
       x0[i]=atof(str_line[7]);y0[i]=atof(str_line[8]);z0[i]=atof(str_line[9]);
       if (i>3&&strcmp(chain0[i],chain0[i-1])!=0) {
          printf("Note: This PDB include multi_chains(first is rebuilt)\n"); break;}
       i++;
    }
   fclose(inpdb);   
   int	Num_atom=i-1;   
   for (int j=1;j<=Num_atom;j++) {
       strcpy(atom_name[j],atom_name0[j]);         
       strcpy(res_name[j],res_name0[j]);    
       strcpy(chain[j],chain0[j]);    //res_id[i]=res_id0[i]; 
       x[1][j]=x0[j]; x[2][j]=y0[j]; x[3][j]=z0[j];   //Used for rmsd calculation
       //printf("%d	%s	%s	%f\n",j,atom_name[j],res_name[j],x[1][j]);
    }
   
   if (strcmp(atom_name[Num_atom],"N1")==0||strcmp(atom_name[Num_atom],"N9")==0) {
       float vectorCP[5];
       for (int j=1;j<=3;j++) {vectorCP[j]=x[j][Num_atom-2]-x[j][Num_atom-4];}
       Num_atom += 1;
       strcpy(atom_name[Num_atom],"P");         
       strcpy(res_name[Num_atom],res_name[Num_atom-2]);        
       for (int j=1;j<=3;j++) {x[j][Num_atom]=x[j][Num_atom-2]+vectorCP[j];}
       //printf("%s %s %f %f %f\n",atom_name[Num_atom],res_name[Num_atom],x[1][Num_atom],x[2][Num_atom],x[3][Num_atom]); 
   }
   else if (strcmp(atom_name[Num_atom],"C4'")==0) {Num_atom -= 1;}
   else {}
   
   Get_neighbour(atom_name,x,Num_atom);
   printf("Enjoy: End of conf_PDB-Read;\nThe length of RNA [N_all_atom: %d length: %dnt)\n",Num_atom,Num_atom/3);
   return Num_atom;
}
/***********************************************************************/
void Read_one_frag(char file[100],int base_type,int frag_No)   //base_type: A-1; U-2; G-3; C-4
{
   FILE *frag;   
   int	i=1;
   char	aline[100],str_line[15][SIZE]; 
   frag = fopen(file,"r+");
   while(!feof(frag)) {
      fgets(aline,100,frag);
      Parse_aline(aline,str_line);
      if (strcmp(str_line[1],"END")==0)  	break;
      if (strcmp(str_line[1],str1)!=0)  	continue;
      strcpy(atom_frag[base_type][frag_No][i],str_line[3]); 
      strcpy(res_frag[base_type][frag_No][i],str_line[4]);
      x_frag[base_type][frag_No][1][i]=atof(str_line[7]);
      x_frag[base_type][frag_No][2][i]=atof(str_line[8]);
      x_frag[base_type][frag_No][3][i]=atof(str_line[9]);
      if (strcmp(str_line[3],str3)==0) {  //atom is P
         if (i<10) { //the first P
            for (int k=1;k<=3;k++) {x_frag_CG[base_type][frag_No][1][k]=x_frag[base_type][frag_No][k][i];}}
         else {      //the last P
            for (int k=1;k<=3;k++) {x_frag_CG[base_type][frag_No][4][k]=x_frag[base_type][frag_No][k][i];}}
        }
      if (strcmp(str_line[3],str2)==0) {  //atom is C4'
         if (i<5) {  //the first C4'
            for (int k=1;k<=3;k++) {x_frag_CG[base_type][frag_No][0][k]=x_frag[base_type][frag_No][k][i];}}
         else if (i<20) {   //the middle C4'
            for (int k=1;k<=3;k++) {x_frag_CG[base_type][frag_No][2][k]=x_frag[base_type][frag_No][k][i];}}
         else {
            for (int k=1;k<=3;k++) {x_frag_CG[base_type][frag_No][5][k]=x_frag[base_type][frag_No][k][i];}}
        }
      if ((base_type==2||base_type==4)&&strcmp(str_line[3],str8)==0) {  //atom is N1
         for (int k=1;k<=3;k++) {x_frag_CG[base_type][frag_No][3][k]=x_frag[base_type][frag_No][k][i];}
        }
      if ((base_type==1||base_type==3)&&strcmp(str_line[3],str9)==0) {  //atom is N9
         for (int k=1;k<=3;k++) {x_frag_CG[base_type][frag_No][3][k]=x_frag[base_type][frag_No][k][i];}
        }
    //printf("%d	%s	%s	%f\n",i,atom_frag[base_type][frag_No][i],res_frag[base_type][frag_No][i],x_frag[base_type][frag_No][1][i]);
      i++;
    }
   //for (int j=1;j<=5;j++) {printf("%d %d %d %f\n",j,base_type,frag_No,x_frag_CG[base_type][frag_No][j][1]);}
   N_frag_atom[base_type][frag_No]=i-1; //printf("base_type: %s %d %d\n",file,frag_No,N_frag_atom[base_type][frag_No]);
   fclose(frag);
}
void Read_frag_AUGC()
{
   char filename[100];
   char base[5]="XAUGC";
   for (int i=1;i<=4;i++) {
       for (int j=1;j<=N_frag[i];j++) {
           sprintf(filename,"%s%c/%c%d.pdb",frag_path,base[i],base[i],j);
           Read_one_frag(filename,i,j);}
   }
   printf("Finish reading AUGC_fragment files\n");
}
/*****************Calculate the RMSD & Replace the CG by AA******************/
/*************************Calculate eigenvalue***************************/
void Calc_eigenvalue(float Tao[5][5],float eigenvalue[5])
{
   float 	max = 10000., sigma0 = 0;
   int	n = 0, mm = 0;  
   float	lambda = -100., unit = 0.001; 
   do { mm++;   
      float sigma = (Tao[1][1]-lambda)*(Tao[2][2]-lambda)*(Tao[3][3]-lambda) + Tao[1][2]*Tao[2][3]*Tao[3][1] + Tao[1][3]*Tao[2][1]*Tao[3][2] - Tao[1][3]*Tao[3][1]*(Tao[2][2]-lambda) - Tao[1][2]*Tao[2][1]*(Tao[3][3]-lambda) - Tao[2][3]*Tao[3][2]*(Tao[1][1]-lambda);                         
      if (sigma * sigma0 < 0) { //if two neighbour sigma with opposite sign, then the solution between them--Zero Theory 
         n++; eigenvalue[n] = lambda - unit/2.;   //the n-th eigenvalue   
         if (eigenvalue[n]<=0) {eigenvalue[n]=0.000001;}  //To avoid NAN in the next. Should be carefull!!!         
        }        
      sigma0 = sigma;  lambda += unit; 
   }  while(!(n==3||lambda > max) && mm<=1090000);  //1030000   
   if (n<3 || mm==1090000) {eigenvalue[1]=0.008;eigenvalue[2]=0.18;eigenvalue[3]=0.6;}
}
float Rotate(float target[5][5],float ref[5][5],float Rot[5][5])  //RMSD
{
  int    i,j,k; 
  float  w[4]={0.05,0.05,0.05,0.01}; //w=0.05;
  float  xnew_CG[5][5],xmiddle[5],ymiddle[5],xx_CG[5][5],yy_CG[5][5];
  float  eigenvalue[5],dis;
  float  R[5][5],R_T[5][5],Tao[5][5],A[5][5],B[5][5];
/****************initialize position & matrix******************/
  for (i=1;i<=4;i++)   {
      xmiddle[i]=0.0; ymiddle[i]=0.0;
      for (j=1;j<=4;j++)  {
          R[i][j]=0.0;  R_T[i][j]=0.0; Tao[i][j]=0.0; 
          A[i][j]=0.0;  B[i][j]=0.0;   Rot[i][j]=0.0;
          xnew_CG[i][j]=0.;} 
  } 
/**************************Translation************************/
  for (i=1;i<=4;i++) {
      for (j=1;j<=3;j++) {
          xmiddle[j] += target[j][i]; ymiddle[j] += ref[j][i]; }} 
  for (i=1;i<=3;i++) {xmiddle[i] /= 4.0;ymiddle[i] /= 4.0;} //center 
  for (i=1; i<=4; i++) {   
      for (j=1; j<=3; j++) {xx_CG[j][i] = target[j][i]-xmiddle[j]; yy_CG[j][i] = ref[j][i]-ymiddle[j]; } 
  }  //put the two center at (0,0,0)  
/*******************Calaulate eigenvalue********************/
  for (i=1; i<=3; i++) {
      for (j=1; j<=3; j++) {            
          for (k=1; k<=4; k++) {R[i][j] += yy_CG[i][k]*xx_CG[j][k]*w[k-1];}}}
  for (i=1; i<=3; i++) {
      for (j=1; j<=3; j++) {R_T[i][j] = R[j][i];} }   //transpose   
  for (i=1; i<=3; i++) {
      for (j=1; j<=3; j++) {
          for (k=1; k<=3; k++) {Tao[i][j] += R_T[i][k]*R[k][j];} } } 
  dis=0.0;
  Calc_eigenvalue(Tao,eigenvalue);
/*******************Calaulate rotation matrix****************/
   for (i=1; i<=3; i++){
       A[1][i] = 1.;   
       A[3][i] = (Tao[2][1]+(Tao[2][2]-eigenvalue[i])*(eigenvalue[i]-Tao[1][1])/Tao[1][2])/(Tao[1][3]*(Tao[2][2]-eigenvalue[i])/Tao[1][2]-Tao[2][3]);    
       A[2][i] = (eigenvalue[i]-Tao[1][1]-Tao[1][3]*A[3][i])/Tao[1][2];            
   }       
   for (i=1; i<=3; i++)  {        
       dis = sqrt(A[1][i]*A[1][i] + A[2][i]*A[2][i] + A[3][i]*A[3][i]);
       A[1][i] /= dis;    A[2][i] /= dis;   A[3][i] /= dis;
   }    
   for (i=1; i<=3; i++)  { //Calculation of matrix B   
       for (j=1; j<=3; j++) {
           for(k=1; k<=3; k++) { B[j][i] += R[j][k]*A[k][i]/sqrt(eigenvalue[i]); }  }
   }
   for (i=1; i<=3; i++) {
       for (j=1; j<=3; j++) { 
           for (k=1; k<=3; k++) { Rot[i][j] += B[i][k]*A[j][k];} }
   }
/*******************rotat the conformation to match****************/
   float RMSD=0.0,dis_r,xnew_CG_i;
   for (k=1; k<=4; k++) {dis_r = 0.;
       for (i=1; i<=3; i++) {xnew_CG_i=0.; 
           for (j=1; j<=3; j++) {xnew_CG_i += Rot[i][j]*xx_CG[j][k];}  
           xnew_CG[i][k]=xnew_CG_i;
           dis_r += (xnew_CG[i][k]-yy_CG[i][k])*(xnew_CG[i][k]-yy_CG[i][k]);}
       RMSD += dis_r;    
    }    
    RMSD = sqrt(RMSD/4.0);    
    fprintf(out_rmsd,"%.5f\n",RMSD);  fflush(out_rmsd);  
    return RMSD;   
}
void Replace_AA(float x_AA[5][50],float xx_CG[5][5],float yy_CG[5][5],float Rot[5][5],float xnew[5][50],int N_AA)
{   
    float	xmiddle[5]={0.0},ymiddle[5]={0.0},xnew_i=0.0;
    int	i,j,k;
    for (j=1;j<=3;j++) {
        for (i=1;i<=4;i++) {
            xmiddle[j] += xx_CG[j][i]; 
            ymiddle[j] += yy_CG[j][i];}
        xmiddle[j] /= 4.0; //center 
        ymiddle[j] /= 4.0;}
    for (k=1;k<=N_AA;k++) {
        for (i=1;i<=3;i++) {x_AA[i][k] -= xmiddle[i]; xnew[i][k]=0.;}}
    for (k=1;k<=N_AA;k++) {   
        for (i=1;i<=3;i++) {xnew_i=0.0;  
            for (j=1;j<=3;j++) {xnew_i += Rot[i][j]*x_AA[j][k];}
            xnew[i][k]=xnew_i+ymiddle[i];} }  
}
/***********************************************************************/
int Rebuild_one_nt(int index,float ref_x[5][5],int base,float Rot[5][5],float final_rmsd[Max_atom/20]) 
{
   float target_x[5][5],rmsd=100.0,U_temp[Max_frag_num][5][5];
   int   num_best_frag=1;
   for (int i=1;i<=4;i++)  {for (int j=1;j<=4;j++)  {Rot[i][j]=0.0;} }
   //memset(Rot,0.0,sizeof(Rot));
   for (int k=1;k<=N_frag[base];k++)  {       
       for (int i=1;i<=3;i++) {
           for (int j=1;j<=4;j++) {
               target_x[i][j]=x_frag_CG[base][k][j][i];} }
       fprintf(out_rmsd,"%d	%d	",index,k); fflush(out_rmsd);       
       float rmsd_temp = Rotate(target_x,ref_x,U_temp[k]);       
       if (rmsd_temp<rmsd) {rmsd=rmsd_temp; num_best_frag=k;}
    }
   final_rmsd[index]=rmsd;
   for (int i=1;i<=3;i++) {
       for (int j=1;j<=3;j++) {Rot[i][j]=U_temp[num_best_frag][i][j];}}
   return num_best_frag;
}
/*****Check clash between atoms in new and old fragments & remove it*******/
void Update_coor(float x0[Max_atom/3][50][5],int nu1,int nu1_a,int nu2,int nu2_a,int n[Max_atom/3])
{
   float move_vector[5],d_mv,x1[5],x2[5],rot_theta;
   float dN1,dN1_temp,tran_vector[5],d_tv,axis[5],d_axis,d_tv_axis;
   d_mv=0.0;dN1=0.0;dN1_temp=0.0;
   for (int i=1;i<=3;i++) {
       move_vector[i]=x0[nu1][nu1_a][i]-x0[nu2][nu2_a][i];
       d_mv += pow(move_vector[i],2);
       dN1 += pow(x0[nu1][nu1_a][i]-x0[nu1][13][i],2);}
   dN1=sqrt(dN1);
   //printf("nu1: %d (%d) %d (%d) %f %f\n",nu1,nu1_a,nu2,nu2_a,sqrt(d_mv),dN1);
   for (int i=1;i<=3;i++) {
       x1[i]=x0[nu1][nu1_a][i]+move_vector[i]/sqrt(d_mv)*6.0; //translate 3A along move_vector
       x2[i]=x1[i]-x0[nu1][13][i];
       dN1_temp += pow(x2[i]-x0[nu1][13][i],2);}
   d_tv=0.0;d_axis=0.0;d_tv_axis=0.0;
   for (int i=1;i<=3;i++) {
       x2[i] /= sqrt(dN1_temp)*dN1;
       tran_vector[i]=x2[i]-x0[nu1][nu1_a][i];  
       d_tv += pow(tran_vector[i],2);
       axis[i]=x0[nu1][13][i]-x0[nu1][6][i];  //N1/N9-C4'
       d_axis += pow(axis[i],2);
       d_tv_axis += tran_vector[i]*axis[i];} 
   for (int i=1;i<=3;i++) {axis[i] /= sqrt(d_axis);}
   rot_theta=acos(d_tv_axis/(sqrt(d_tv)*sqrt(d_axis)));
   if (rot_theta>3.14/2) {rot_theta -= 3.14/2;}   
   //rot_theta = 3.14/4.;
   //printf("theta: %f\n",rot_theta);
   
   float rot_matrix[5][5];
   rot_matrix[1][1]=axis[1]*axis[1]*(1-cos(rot_theta))+cos(rot_theta);
   rot_matrix[1][2]=axis[1]*axis[2]*(1-cos(rot_theta))-axis[3]*sin(rot_theta);
   rot_matrix[1][3]=axis[1]*axis[3]*(1-cos(rot_theta))+axis[2]*sin(rot_theta);
   rot_matrix[2][1]=axis[2]*axis[1]*(1-cos(rot_theta))+axis[3]*sin(rot_theta);
   rot_matrix[2][2]=axis[2]*axis[2]*(1-cos(rot_theta))+cos(rot_theta);
   rot_matrix[2][3]=axis[2]*axis[3]*(1-cos(rot_theta))-axis[1]*sin(rot_theta);
   rot_matrix[3][1]=axis[3]*axis[1]*(1-cos(rot_theta))-axis[2]*sin(rot_theta);
   rot_matrix[3][2]=axis[3]*axis[2]*(1-cos(rot_theta))+axis[1]*sin(rot_theta);
   rot_matrix[3][3]=axis[3]*axis[3]*(1-cos(rot_theta))+cos(rot_theta);   
   for (int k=14;k<=n[nu1];k++) {
       float x3[5];
       for (int i=1;i<=3;i++) {x3[i]=x0[nu1][k][i]-x0[nu1][13][i];}
       x0[nu1][k][1]=rot_matrix[1][1]*x3[1]+rot_matrix[1][2]*x3[2]+rot_matrix[1][3]*x3[3]+x0[nu1][13][1];
       x0[nu1][k][2]=rot_matrix[2][1]*x3[1]+rot_matrix[2][2]*x3[2]+rot_matrix[2][3]*x3[3]+x0[nu1][13][2];
       x0[nu1][k][3]=rot_matrix[3][1]*x3[1]+rot_matrix[3][2]*x3[2]+rot_matrix[3][3]*x3[3]+x0[nu1][13][3];
    } 
}
void Remove_clash(float x0[Max_atom/3][50][5],int N,int n[Max_atom/3],float rmsd_k[Max_atom/3]) 
{
   float cutoff=1.5;
   for (int i=1;i<=NN_num[N];i++) {   
       int NN_n=NN[N][i]; //printf("%d %d %d\n",N,i,NN_n);
       for (int j=13;j<=n[NN_n];j++) {
           for (int k=13;k<=n[N];k++) { //just consider the base atoms
               float dist=0.0; 
               for (int l=1;l<=3;l++) {dist += pow(x0[N][k][l]-x0[NN_n][j][l],2);}   
               if (sqrt(dist)<=cutoff) {
                  if (rmsd_k[N]>=rmsd_k[NN_n]) {Update_coor(x0,N,k,NN_n,j,n);}
                  else {Update_coor(x0,NN_n,j,N,k,n);} }
              } } }
}
void Out_AA_PDB(int res_id,int N,int atom_id[50],char atom_type[50][SIZE],char res_type[SIZE],float AA_x[50][5])
{
    for (int i=1;i<=N;i++) {
        fprintf(out_AA,"%s  %5d  %-4s %2s %c%4d    %8.3f%8.3f%8.3f  %s  %s\n","ATOM",atom_id[i],atom_type[i],res_type,'A',res_id,AA_x[i][1],AA_x[i][2],AA_x[i][3],"1.00","0.00");  
        fflush(out_AA); }
}
void Assemble(int N0,float x[5][Max_atom])
{
   int	i,j,k,nn=0,base,nucl_AA,num_frag,N_AA;
   float	AA_x[Max_atom/3][50][5],final_rmsd[Max_atom/3];
   char	AA_atom[Max_atom/3][50][SIZE],AA_res[Max_atom/3][50][SIZE];
   int	AA_N[Max_atom/3],atom_id[Max_atom/20][50];
   for (i=1;i<N0;i++) {
       if (fmod(i,3)==0) {
          if (strcmp(atom_name[i],str8)!=0 && strcmp(atom_name[i],str9)!=0) {
              printf("Warning: the given CG pdb file is not standard(i.e., P,C,N1/N9) %d %s\n",i/3,atom_name[i]);break;}
          float y_CG[5][5],x_CG[5][5],Rot[5][5]={{0.0}},x_AA[5][50];
          for (j=1;j<=3;j++) {   //coor of four CG beads 
              y_CG[j][1]=x[j][i-2];  //P
              y_CG[j][2]=x[j][i-1];  //C4'
              y_CG[j][3]=x[j][i];    //N
              y_CG[j][4]=x[j][i+1];  //P
             } 
          if (strcmp(res_name[i],"A")==0)      {base=1;} 
          else if (strcmp(res_name[i],"U")==0) {base=2;}
          else if (strcmp(res_name[i],"G")==0) {base=3;}
          else if (strcmp(res_name[i],"C")==0) {base=4;}
          else {printf("Error: There is a usual nucleotide\n"); break;}        
          num_frag = Rebuild_one_nt(i/3,y_CG,base,Rot,final_rmsd); 
          N_AA=N_frag_atom[base][num_frag]; 
          for (j=1;j<=3;j++) {
              for (k=1;k<=4;k++)    {x_CG[j][k]=x_frag_CG[base][num_frag][k][j];}
              for (k=1;k<=N_AA;k++) {x_AA[j][k]=x_frag[base][num_frag][j][k];}  }          
          float xnew[5][50]; //coor of rebuilt AA         
          Replace_AA(x_AA,x_CG,y_CG,Rot,xnew,N_AA);
          nucl_AA=i/3;
          for (k=2;k<=N_AA-3;k++) {  //remove the first C4' & last P, PO1, PO2
              atom_id[nucl_AA][k-1]=nn+k-1;     //the atom order in reconstructed all-atom file
              strcpy(AA_atom[nucl_AA][k-1],atom_frag[base][num_frag][k]); //all-atom name
              strcpy(AA_res[nucl_AA][k-1],res_frag[base][num_frag][k]); //res name (e.g., AUGC)     
              for (j=1;j<=3;j++) {AA_x[nucl_AA][k-1][j]=xnew[j][k];}                      
             }           
          nn += N_AA-4;   //fragement includes other C4',P,PO1,PO2 beyond this nucleotide (i.e., N_AA-4)   
          AA_N[nucl_AA]=N_AA-4;  //Num of atoms in i/3-th nucleotide 
          //Check bond connection & repair
          float v0[5],dPO3; 
          int O3_id = 10;                      
          if (nucl_AA>1) { 
             for (k=1;k<=AA_N[nucl_AA-1];k++) {//Get the last O3'
                 if (strcmp(AA_atom[nucl_AA-1][k],"O3'")==0) {O3_id=k; break;}}  
             dPO3=0.0;
             for (j=1;j<=3;j++) {v0[j]=0.0;
                 v0[j] = AA_x[nucl_AA-1][O3_id][j]-AA_x[nucl_AA][1][j];  //bond O3(i-1)-P(i)
                 dPO3 += pow(v0[j],2);}
             dPO3 = sqrt(dPO3); 
             if (dPO3>1.8 || dPO3<1.4) {
                for (j=1;j<=3;j++) {v0[j] /= dPO3; //unit vector
                    AA_x[nucl_AA-1][O3_id][j] -= v0[j]*(dPO3-1.6)/2.;  //change O3 coor
                    for (int l=1;l<=3;l++) {  //Change coor of P, PO1, PO2
                        AA_x[nucl_AA][l][j] += v0[j]*(dPO3-1.6)/2.;} } 
                }
           }                 
      }
   }
  for (i=1;i<N0;i++) {
      if (fmod(i,3)==0) {
         Remove_clash(AA_x,i/3,AA_N,final_rmsd); }}
  for (i=1;i<N0;i++) {
       if (fmod(i,3)==0) {
          nucl_AA=i/3;
          Out_AA_PDB(nucl_AA,AA_N[nucl_AA],atom_id[nucl_AA],AA_atom[nucl_AA],res_name[i],AA_x[nucl_AA]); }}
}


