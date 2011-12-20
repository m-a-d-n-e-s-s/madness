#include <madness_config.h>
#define ENABLE_CUBLAS 1
#include <tensor/tensor.h>
#include <tensor/cu_mtxmq.h>
#include <stdio.h>

extern "C" void device_synchronize(void **,unsigned int);
using namespace madness;

static double tttt;
#define STARTt_TIMER  tttt=wall_time();
#define ENDt_TIMER(msg) tttt=wall_time()-tttt;  printf("timer: %20.20s %8.10fs \n", msg, tttt)

int main(){

  STARTt_TIMER;
  sleep(2);
  ENDt_TIMER("timer init");
 
  STARTt_TIMER;
  lsk(1); 
  device_synchronize(0,0);
  ENDt_TIMER("1");
  
  STARTt_TIMER;
  lsk(1); 
  device_synchronize(0,0);
  ENDt_TIMER("1");
  
  STARTt_TIMER;
  lsk(2); 
  device_synchronize(0,0);
  ENDt_TIMER("2");
  
  STARTt_TIMER; 
  lsk(10);
  device_synchronize(0,0);
  ENDt_TIMER("10");
  
  STARTt_TIMER; 
  lsk(100);
  device_synchronize(0,0);
  ENDt_TIMER("100");
  
  STARTt_TIMER; 
  lsk(1000);
  device_synchronize(0,0);
  ENDt_TIMER("1000");

  STARTt_TIMER; 
  lsk(10000);
  device_synchronize(0,0);
  ENDt_TIMER("10000");

  STARTt_TIMER;
  lsk1(1);
  device_synchronize(0,0);
  ENDt_TIMER("access 1");

  STARTt_TIMER;
  lsk1(2);
  device_synchronize(0,0);
  ENDt_TIMER("access 2");

  STARTt_TIMER;
  lsk1(10);
  device_synchronize(0,0);
  ENDt_TIMER("access 10");

  STARTt_TIMER;
  lsk1(100);
  device_synchronize(0,0);
  ENDt_TIMER("access 100");

  STARTt_TIMER;
  lsk1(1000);
  device_synchronize(0,0);
  ENDt_TIMER("access 1000");

  STARTt_TIMER;
  lsk1(10000);
  device_synchronize(0,0);
  ENDt_TIMER("access 10000");
}
