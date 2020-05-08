#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>





struct particle_data
{
    float Pos[3];
    float Vel[3];
    float Mass;
    int Type;

    float Tho,U,Temp,Ne;
}p[N];//Basic information of particles  

int *Id;

float Distance(particle_data A,particle_data B)
{   
    return sqrt(pow((A.Pos[0]-B.Pos[0]),2)+pow((A.Pos[1]-B.Pos[1]),2)+pow((A.Pos[2]-B.Pos[2]),2));
}//calculate the distance between two particles//计算两个粒子之间距离的函数

ShellSort(float Arr[],int length)
{
    int increasement = length;
    int i,j,k;
    do{
       increasement=increasement/3+1;
       for(i=0;i<increasement;i++){
                                   for(j=i+increasement;j<length;j+=increasement)
                                      {
                                       if(Arr[j]<Arr[j-increasement]&&fabs(Arr[j]-Arr[j-increasement])>1e-6)
                                         {
                                          float t=Arr[j];
                                          for(k=j-increasement;k>=0&&t<Arr[k]&&fabs(Arr[k]-t)>1e-6;k-=increasement)
                                              {
                                               Arr[k+increasement]=Arr[k];
                                              }
                                           Arr[k+increasement]=t;   
                                         }                                                
                                      }
                                  }
      }while(increasement>1);
}//希尔排序

particle_data *BH(particle_data *p,int z,char *str)//p为进行操作的particle_data结构体的头指针，x为判定所用临界数密度，z为粒子数
{   
    FILE *f2;//后面将筛选出来的粒子输入进文件分f2中
    particle_data q[z],*r;//q数组用来记录筛选出的粒子，r用来作函数的返回值    
    float a[z-1];//用来记录粒子间距离，因为自身与自身距离无意义，故数组长度为z-1
    float b=0.5;//以某粒子为中心，在该半径区域内进行粒子密度计算
    float m=p[0].Mass;
    int i,j,k,l;//用来循环
    int n=0;//用来标记筛选出粒子个数
    int d=100;  
                    
    for(i=0;i<z;i++){   k=0;
                        for(j=0;j<z;j++){
                                         
                                         if(j!=i){
                                                 a[k]=Distance(p[i],p[j]);
                                                 k++;//将i，j两粒子之间距离依次算出，并存放进浮点数组a
                                                 }
                                        }
                        ShellSort(a,z-1);//对浮点数组a进行从小到大的希尔排序
                        l=0;//每次循环前，l,k需要赋值为0
                        for(j=0;j<z-1;j++) {
                                            if(((a[j]<b))) 
                                              {l++;}
                                           }
                                         //判断该半径内有多少粒子，用l来记录粒子数（因为上一步骤中已经跳过了对i，i粒子即同一粒子计算距离，所以循环数为z-1）
                                         
                        if(l!=0&&l>d){ for(j=0;j<l;j++) { if(((j+2)*3*m/(4*pi)/pow(a[j],3))>(7.98553*10^(-5)/(a[j])^2))
                                                                            {q[n]=p[i];n++;break;}
                                                   }
                                     }/*当l不为0即半径b内有粒子时，判断(j+2)*m/a[j]^3是否大于等于临界数密度，此处取以b为半边长的正方体内，与取球体相近，n为记录该，j是因为作为数组的a[j]元素，实际上是整个数组的j+1个元素，再算上中心处粒子自身，所以j+2
                                         质量密度判据数值为根据史瓦西半径反推而得，前面常数单位为Msun/m，a[j]单位为具体模拟所采用的单位，例如kpc,m单位为Msun
                                      */
                        /*if(l>x){q[n]=p[i];n++;}*/
                    }//上述操作对所有粒子遍历，一旦满足，则作为q[n]记录下来
    
    for(i=0;i<n;i++) printf("x%d=%f,y%d=%f,z%d=%f\n",i+1,q[i].Pos[0],i+1,q[i].Pos[1],i+1,q[i].Pos[2]);//print the coordinates of target particles//这一步是在终端中直接看出筛选出来多少粒子，无需打开文件，可以删去
    f2=fopen(str,"w");
    for(i=0;i<n;i++) {fprintf(f2,"%.2f ",q[i].Pos[0]);fprintf(f2,"%.2f ",q[i].Pos[1]);fprintf(f2,"%.2f\n",q[i].Pos[2]);}
    fclose(f2);//将筛选出粒子坐标输入文件
    r = (particle_data *)malloc(n* sizeof(struct particle_data)); //为r分配空间
    for(i=0;i<n;i++) r[i] = q[i];  //赋值    
    return r;//函数返回r,r是长度为n恰好装下所有筛选出粒子的结构提数组头指针
          
}
int main(int argc, char *argv[]){
    //arg
    char input_fname[200],output_fname[200],name1[200],name2[200];
    char c;
    int type,snapshot_number,files;
    int N;
    FILE *f1;//f1用来存放所有粒子数据的文件，如不从外界读取初始条件，则要保存给定的所有粒子位置信息
    
    int i,j,k;//用来循环
    if(argc==1) 
    {
      printf("Too few arguments to function int main()");return 0;
    }
   
    
    else 
    {
      
                          for(k=0;k<argc-1;k++){
                                           sprintf(input_fname,"%s",argv[k+1]);
                                           sprintf(output_fname,"%s_target_particles",argv[k+1]);
                                           N=0;//用来标记粒子数;
                                           FILE *fr=fopen(input_fname,"r");//此文件指针用来计数有多少粒子
                                           if(fr==NULL)
                                             {
                                              printf("No file available\n"); return 0;
                                             }
                                           while(c=fgetc(fr)!=EOF) {if(c=='\n') N++;}
                                                                                N++;//最后一行粒子也算进去;
                                                                                N=N-14;//snapshots_reader输出的文本文件前面14行为一些注释内容，所以减去;
      
                                           FILE *fp=fopen(input_fname,"r");
                                           

                                           if(fp==NULL)
                                             {
                                              printf("No file available\n"); return 0;
                                             }
                                            
                                           for (i=0;i<N;i++)
                                             {                 
                                                               fscanf(fp,"%f",&Id[i]);
                                                               fscanf(fp,"%f",&p[i].Type);
                                              for(j=0;j<3;j++) fscanf(fp,"%f",&p[i].Pos[j]);
                                              for(j=0;j<3;j++) fscanf(fp,"%f",&p[i].Vel[j]);
                                                               fscanf(fp,"%f",&p[i].Mass);
                                                               fscanf(fp,"%f",&p[i].Rho);
                                                               fscanf(fp,"%f",&p[i].U);
                                                               fscanf(fp,"%f",&p[i].Temp);
                                                               fscanf(fp,"%f",&P[i].Ne);
                                             }
                                            /*for(int i = 1; i <= NumPart; i++){
    fprintf(pcd, "%d %d %f %f %f %f %f %f %f %10e %f %f %f\n",\
    Id[i], P[i].Type, P[i].Pos[0],P[i].Pos[1],P[i].Pos[2], P[i].Vel[0], P[i].Vel[1], P[i].Vel[2], P[i].Mass,\
    P[i].Rho, P[i].U, P[i].Temp, P[i].Ne);
  }*/
                                            fclose(fp);
                                            BH(p,N,output_fname);
                                            printf("%d\n\n\n\n",argc);

                                          }
                         
    }
    
    return 0; 
   
}
