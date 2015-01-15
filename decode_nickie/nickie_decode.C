#include <iostream>
#include "nickie_decode.h"
using namespace std;

int main(){

  int   len      = 5;
  int*  len_pt   = &len;
  float input[5] = {0., 0., 0., 0., 0.};
  int   steps    = 31;
  //float max      = .0025;
  float max      = .06;
  int   in       = 1;
  int   out      = 0;
  float step     = max * 2. / ( steps - 1. );
  for( int i = 0; i < steps; i++){
    input[in]   = -max + step * i;
    float fp[5] = {x_sp_cfp_(input, len_pt), t_sp_cfp_(input, len_pt), y_sp_cfp_(input, len_pt),
		   p_sp_cfp_(input, len_pt), l_sp_cfp_(input, len_pt)};
    cout << fp[out] << " " << input[in] << endl; 
  }

  return 0;

}
