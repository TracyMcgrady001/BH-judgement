#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#define m 1//mass of single particles//粒子质量
#define N 10000//number of particles//例子总数
#define D 5000//critical mass density//质量密度判据



struct particle_data
{
    float Pos[3];
    float Vel[3];
    float Mass;
    int Type;

    float Tho,U,Temp,Ne;
}p[N];//Basic information of particles  



float Distance(particle_data A,particle_data B)
{   
    return sqrt(pow((A.Pos[0]-B.Pos[0]),2)+pow((A.Pos[1]-B.Pos[1]),2)+pow((A.Pos[2]-B.Pos[2]),2));
}//calculate the distance between two particles//计算两个粒子之间距离的函数

void ShellSort(float Arr[],int length)
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

particle_data *BH(particle_data *p,float x,int z)//p为进行操作的particle_data结构体的头指针，x为判定所用临界数密度，z为粒子数
{   
    FILE *f2;//后面将筛选出来的粒子输入进文件分f2中
    particle_data q[z],*r;//q数组用来记录筛选出的粒子，r用来作函数的返回值    
    float a[z-1];//用来记录粒子间距离，因为自身与自身距离无意义，故数组长度为z-1
    float b=0.0001;//以某粒子为中心，在该半径区域内进行粒子密度计算
    int i,j,k,l;//用来循环
    int n=0;//用来标记筛选出粒子个数
      
                    
    for(i=0;i<z;i++){   k=0;
                        for(j=0;j<z;j++){
                                         
                                         if(j!=i){
                                                 a[k]=Distance(p[i],p[j]);
                                                 k++;//将i，j两粒子之间距离依次算出，并存放进浮点数组a
                                                 }
                                        }
                        ShellSort(a,z-1);//对浮点数组a进行从小到大的希尔排序
                        l=0;//每次循环前，l,k需要赋值为0
                        for(j=0;j<z-1;j++) {if(((a[j]<b)||(fabs(a[j]-b)<1e-6))) 
                                              {l++;}
                                           }//判断该半径内有多少粒子，用l来记录粒子数（因为上一步骤中已经跳过了对i，i粒子即同一粒子计算距离，所以循环数为z-1）
                                           
                        if(l!=0){ for(j=0;j<l;j++) { if((((j+2)*m/pow(a[j],3)>x)||(fabs((j+2)*m/pow(a[j],3)-x)<1e-6))) 
                                                        {q[n]=p[i];n++;}
                                                   }
                                }//当l不为0即半径b内有粒子时，判断(j+2)*m/a[j]^3是否大于等于临界数密度，此处取以b为半边长的正方体内，与取球体相近，n为记录该，j是因为作为数组的a[j]元素，实际上是整个数组的j+1个元素，再算上中心处粒子自身，所以j+2

                    }//上述操作对所有粒子遍历，一旦满足，则作为q[n]记录下来
    
    for(i=0;i<n;i++) printf("x%d=%f,y%d=%f,z%d=%f\n",i+1,q[i].Pos[0],i+1,q[i].Pos[1],i+1,q[i].Pos[2]);//print the coordinates of target particles//这一步是在终端中直接看出筛选出来多少粒子，无需打开文件，可以删去
    f2=fopen("/home/tracy/Workreport/objectparticles.txt","w");
    for(i=0;i<n;i++) {fprintf(f2,"%.2f ",q[i].Pos[0]);fprintf(f2,"%.2f ",q[i].Pos[1]);fprintf(f2,"%.2f\n",q[i].Pos[2]);}
    fclose(f2);//将筛选出粒子坐标输入文件
    r = (particle_data *)malloc(n* sizeof(struct particle_data)); //为r分配空间
    for(i=0;i<n;i++) r[i] = q[i];  //赋值    
    return r;//函数返回r,r是长度为n恰好装下所有筛选出粒子的结构提数组头指针
          
}
int main(){
    FILE *f1;//f1用来存放所有粒子数据的文件，如不从外界读取初始条件，则要保存给定的所有粒子位置信息
    
    int i,j;//用来循环
    
    FILE *fp=fopen("./random1_028_posi.txt","r");
    if(fp==NULL)
    {
     printf("No coordinates available\n"); return 0;
    }
    for (i=0;i<N;i++){
        for(j=0;j<3;j++) fscanf(fp,"%f",&p[i].Pos[j]);
    }
    
    fclose(fp);//input coordinates of particles from a certain text//粒子位置信息从文件中读取，文件指针为fp
    
    /*srand((unsigned)time(NULL));
    for(i=0;i<N;i++){
                      p[i].Pos[0]=rand()%101,p[i].Pos[1]=rand()%101,p[i].Pos[2]=rand()%101;
                    }   *///the particles are distributed evenly in a cubic which is 10 units in length//此处注释掉的是用于测试的随机位置分布粒子
    BH(p,D,N);//p，D，N作为参数传入函数
    printf("\n\n\n\n");//为了美观

    /*for(i=0;i<N;i++) printf("x%d=%f,y%d=%f,z%d=%f\n",i+1,p[i].Pos[0],i+1,p[i].Pos[1],i+1,p[i].Pos[2]);*///输出所有粒子信息，但是太长，注释掉

    f1=fopen("/home/tracy/Workreport/allparticles.txt","w");
    for(i=0;i<N;i++) {fprintf(f1,"%.2f ",p[i].Pos[0]);fprintf(f1,"%.2f ",p[i].Pos[1]);fprintf(f1,"%.2f\n",p[i].Pos[2]);}
    fclose(f1);//将所有粒子位置信息输入文件f1中
    return 0; 
   
}
