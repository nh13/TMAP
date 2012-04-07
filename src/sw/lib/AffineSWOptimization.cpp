#include <cstring>
#include <sstream>
#include "solution1.h"
#include "solution2.h"
#include "solution3.h"
#include "solution4.h"
#include "solution5.h"
#include "solution6.h"
#include "solution7.h"
#include "solution8.h"
#include "solution9.h"
#include "solution10.h"
#include "AffineSWOptimization.h"

using namespace std;

AffineSWOptimization::AffineSWOptimization(int type) {
    myType = type;
    switch(type) {
      case 1:
        s = new Solution1();
        break;
      case 2:
        s = new Solution2();
        break;
      case 3:
        s = new Solution3();
        break;
      case 4:
        s = new Solution4();
        break;
      case 5:
        s = new Solution5();
        break;
      case 6:
        s = new Solution6();
        break;
      case 7:
        s = new Solution7();
        break;
      case 8:
        s = new Solution8();
        break;
      case 9:
        s = new Solution9();
        break;
      case 10:
        s = new Solution10();
        break;
      default:
        fprintf(stderr, "Error: unknown type: %d\n", type);
        exit(1);
    }
}

int AffineSWOptimization::process(string b, string a, int qsc, int qec,
                                     int mm, int mi, int o, int e, int dir,
                                     int *opt, int *te, int *qe, int *n_best) {
    return s->process(b, a, qsc, qec, mm, mi, o, e, dir, opt, te, qe, n_best);
}

AffineSWOptimization::~AffineSWOptimization()
{
  //s->~Solution();
  delete s;
}

