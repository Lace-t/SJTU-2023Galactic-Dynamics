#include<stdio.h>
#include<math.h>
#include<memory.h>

#define Mb 409.0  //M_gal=2.235*10^7
#define Md 2856.0
#define Mh 1018.0
#define bb 0.23 //kpc
#define ad 4.22
#define bd 0.292
#define ah 2.562
#define L  200.0
#define gamma 2.0

#define vf 0.0102202 //将速度从 10 km/s -> kpc/Myr
#define step 1000//写入间隔步数

void bar(double &potb,double &fxb,double &fyb,double &fzb,double x,double y,double z){
    potb=-Mb/sqrt(pow(x,2.0)+pow(y,2.0)+pow(z,2.0)+pow(bb,2.0));
    double com=-Mb/pow(pow(x,2.0)+pow(y,2.0)+pow(z,2.0)+pow(bb,2.0),1.5);
    fxb=x*com;
    fyb=y*com;
    fzb=z*com;
}

void disk(double &potd,double &fxd,double &fyd,double &fzd,double x,double y,double z){
    potd=-Md/sqrt(pow(x,2.0)+pow(y,2.0)+pow(ad+sqrt(pow(z,2.0)+pow(bd,2.0)),2.0));
    double com=-Md/pow(pow(x,2.0)+pow(y,2.0)+pow(ad+sqrt(pow(z,2.0)+pow(bd,2.0)),2.0),1.5);
    fxd=x*com;
    fyd=y*com;
    fzd=z*com*(ad+sqrt(pow(bd,2.0)+pow(z,2.0)))/sqrt(pow(bd,2.0)+pow(z,2.0));
}

void halo(double &poth,double &fxh,double &fyh,double &fzh,double x,double y,double z){
    poth=Mh/ah*(1.0/(gamma-1.0)*log((1.0+pow(sqrt(pow(x,2.0)+pow(y,2.0)+pow(z,2.0))/ah,gamma-1.0))/(1.0+pow(L/ah,gamma-1.0)))-pow(L/ah,gamma-1.0)/(1.0+pow(L/ah,gamma-1)));
    double com=-Mh/ah*pow(sqrt(pow(x,2.0)+pow(y,2.0)+pow(z,2.0))/ah,gamma-1.0)/(pow(x,2.0)+pow(y,2.0)+pow(z,2.0))/(1.0+pow(sqrt(pow(x,2.0)+pow(y,2.0)+pow(z,2.0))/ah,gamma-1.0));
    fxh=com*x;
    fyh=com*y;
    fzh=com*z;
}

void integrator(double h,int n,double d[][8]){//Leapfrog integrator
    int i,j,k;
    FILE *fp;
    double potb,potd,poth,potf;
    double fxb,fxd,fxh,fxf;
    double fyb,fyd,fyh,fyf;
    double fzb,fzd,fzh,fzf;
    //首先计算初始势
    bar(potb,fxb,fyb,fzb,d[0][1],d[0][2],d[0][3]);
    disk(potd,fxd,fyd,fzd,d[0][1],d[0][2],d[0][3]);
    halo(poth,fxh,fyh,fzh,d[0][1],d[0][2],d[0][3]);
    d[0][7]=potb+potd+poth;
    //写入第一行数据
    fp=fopen("pot.txt","w");
    for(i=0;i<8;i++){
        fprintf(fp,"%lf ",d[0][i]);
    }
    fprintf(fp,"\n");
    fclose(fp);
    //开始演化
    for(j=0;j<n/10000;j++){
        for(i=1;i<10001;i++){
            d[i][0]=d[i-1][0]+h;
            d[i][1]=d[i-1][1]+d[i-1][4]*h*vf+0.5*pow(h*vf,2.0)*(fxb+fxd+fxh);
            d[i][2]=d[i-1][2]+d[i-1][5]*h*vf+0.5*pow(h*vf,2.0)*(fyb+fyd+fyh);
            d[i][3]=d[i-1][3]+d[i-1][6]*h*vf+0.5*pow(h*vf,2.0)*(fzb+fzd+fzh);
            d[i][4]=d[i-1][4]+0.5*h*vf*(fxb+fxd+fxh);
            d[i][5]=d[i-1][5]+0.5*h*vf*(fyb+fyd+fyh);
            d[i][6]=d[i-1][6]+0.5*h*vf*(fzb+fzd+fzh);

            bar(potb,fxb,fyb,fzb,d[i][1],d[i][2],d[i][3]);
            disk(potd,fxd,fyd,fzd,d[i][1],d[i][2],d[i][3]);
            halo(poth,fxh,fyh,fzh,d[i][1],d[i][2],d[i][3]);
            d[i][4]=d[i][4]+0.5*h*vf*(fxb+fxd+fxh);
            d[i][5]=d[i][5]+0.5*h*vf*(fyb+fyd+fyh);
            d[i][6]=d[i][6]+0.5*h*vf*(fzb+fzd+fzh);
            d[i][7]=potb+potd+poth;
            if(i%1000==999){
                printf("step=%d\n",i+10000*j+1);
            }
        }
        //每10000行写入一次数据
        fp=fopen("pot.txt","a+");
        for(i=step;i<10001;i=i+step){
            for(k=0;k<8;k++){
                fprintf(fp,"%lf ",d[i][k]);
            }
            fprintf(fp,"\n");
        }
        fclose(fp);
        //写完后将数组复原
        for(i=0;i<8;i++){
            d[0][i]=d[10000][i];
        }
    }
}

int main(){
    double h;
    int n;

    h=0.001;//时间间隔，以Myr为单位
    n=(int)1000/h;
    double d[10001][8];//t,x,y,z,vx,vy,vz,pot
    memset(d,0.0,sizeof(d));
    //初始化
    d[0][0]=0.0;d[0][1]=8.2;d[0][2]=0.0;d[0][3]=0.025;
    d[0][4]=-1.11;d[0][5]=25.224;d[0][6]=0.725;

    integrator(h,n,d);
    return 0;
}
