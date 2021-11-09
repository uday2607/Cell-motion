#include <algorithm>
#include <vector>
#include <iostream>
#include <eigen3/Eigen/Dense>
using namespace Eigen;
using namespace std;

class isEqual{
  public:
    vector<int>* ET;
    double E;
    isEqual(vector<int> *e, double a){ET = e; E = a;}
    void operator()(double i){static int it = 0;
                    if (i == E) ET->push_back(it); it++;}
};

class isGreater{
  public:
    vector<int>* ET;
    double E;
    isGreater(vector<int> *e, double a){ET = e; E = a;}
    void operator()(double i){static int it = 0;
                    if (i > E) ET->push_back(it); it++;}
};

int main(int argc,char **argv){
    ArrayXd P = ArrayXd::Random(20);
    vector<int> GT;
    for_each(&P(1),&P(10),isGreater(&GT, 0.0));
    cout<<P<<endl;
    for(int i=0;i<GT.size();++i)cout<<GT[i]<<" ";
    cout << endl;
    Map<ArrayXi> H(&GT[0], GT.size());;
    ArrayXd U = P(seq(1, 10))(H);
    cout<<U<<endl;
    return 0;
}
