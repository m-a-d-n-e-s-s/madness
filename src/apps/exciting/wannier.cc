#include <iostream>
#include <math.h>

extern "C" void readinput_();
extern "C" void wann_init1_();
extern "C" void wann_unk_(int* n,double* vpl,double* vrc,double* val);

int main(int argn, char** argv) {
    //k-mesh division
    int nk[3]={4,4,4};
    //number of k-points
    int nkpt;
    //k-point in lattice coordinates
    double vkl[3];
    //vectors of a 3D-plot box
    double p1[3]={12.0,0.0,0.0};
    double p2[3]={0.0,12.0,0.0};
    double p3[3]={0.0,0.0,12.0};
    //cener point of the box
    double p0[3]={0.0,0.0,0.0};
    //number of r-mesh points
    int nrmesh[3]={50,50,50};
    //index of Wannier function
    int n=3;
    
    double vr[3],val[2],val0[2];
    
    readinput_();
    wann_init1_();
    
    nkpt=nk[0]*nk[1]*nk[2];
    
    for (int i1=0;i1<nrmesh[0];i1++) {
      for (int i2=0;i2<nrmesh[1];i2++) {
        for (int i3=0;i3<nrmesh[2];i3++) {
          //point inside 3D box
          for (int i=0;i<3;i++) vr[i]=(1.0*i1/nrmesh[0]-0.5)*p1[i] +
                                      (1.0*i2/nrmesh[1]-0.5)*p2[i] +
                                      (1.0*i3/nrmesh[2]-0.5)*p3[i] + p0[i];
          val[0]=val[1]=0.0;
          //compute \sum_{k}u_{nk}(r)
          for (int j1=0;j1<nk[0];j1++) {
            for (int j2=0;j2<nk[1];j2++) {
              for (int j3=0;j3<nk[2];j3++) {
                //k-point in lattice coordinates
                vkl[0]=1.0*j1/nk[0];
                vkl[1]=1.0*j2/nk[1];
                vkl[2]=1.0*j3/nk[2];
                //get value of u_{nk}(r)
                wann_unk_(&n,vkl,vr,val0);
                val[0]+=val0[0];
                val[1]+=val0[1];
              }
            }
          }
          val[0]=val[0]/nkpt;
          val[1]=val[1]/nkpt;
          std::cout<<vr[0]<<" "<<vr[1]<<" "<<vr[2]<<" "<<sqrt(val[0]*val[0]+val[1]*val[1])<<std::endl;
        }
      }
    }
}                
