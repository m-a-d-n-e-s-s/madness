#include <iostream>

extern "C" void readinput_();
extern "C" void gndstate_();

int main(int argc, char** argv)
{
  readinput_();
  gndstate_();
  return 0;
}
