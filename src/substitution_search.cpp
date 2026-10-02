// Permutation annealing. Space (23) is fixed; masks never enter input grams.
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>
using namespace std;
using Key=array<int,24>;
struct Gram {array<int,4> c; double w;};
int main(int argc,char**argv){
 if(argc!=8){cerr<<"lm.bin grams.txt rank.txt seed restarts steps output\n";return 2;}
 vector<double> lm(24*24*24*24);ifstream lf(argv[1],ios::binary);
 lf.read((char*)lm.data(),lm.size()*sizeof(double));if(!lf)return 3;
 vector<Gram> grams;ifstream gf(argv[2]);Gram g;
 while(gf>>g.c[0]>>g.c[1]>>g.c[2]>>g.c[3]>>g.w)grams.push_back(g);
 if(grams.empty())return 4;
 Key ranked;ifstream rf(argv[3]);for(auto&x:ranked)rf>>x;if(!rf||ranked[23]!=23)return 5;
 int restarts=stoi(argv[5]),steps=stoi(argv[6]);mt19937_64 rng(stoull(argv[4]));
 uniform_real_distribution<double> unif(0,1);uniform_int_distribution<int> sym(0,22);
 auto value=[&](int i,const Key& k){auto &v=grams[i].c;return lm[((k[v[0]]*24+k[v[1]])*24+k[v[2]])*24+k[v[3]]];};
 auto full=[&](const Key&k){double z=0;for(int i=0;i<(int)grams.size();i++)z+=grams[i].w*value(i,k);return z;};
 double mass=0;for(auto&v:grams)mass+=v.w;
 vector<int> affected[23][23];
 for(int a=0;a<23;a++)for(int b=a+1;b<23;b++)for(int i=0;i<(int)grams.size();i++){
  bool found=false;for(int c:grams[i].c)if(c==a||c==b)found=true;
  if(found)affected[a][b].push_back(i);
 }
 ofstream out(argv[7]);out<<setprecision(17);
 for(int restart=0;restart<restarts;restart++){
  Key key=ranked;if(restart)shuffle(key.begin(),key.begin()+23,rng);
  vector<double> cached(grams.size());for(int i=0;i<(int)grams.size();i++)cached[i]=value(i,key);
  double score=full(key),best=score;Key bestkey=key;
  auto delta=[&](int a,int b){double d=0;for(int i:affected[a][b])d+=grams[i].w*(value(i,key)-cached[i]);return d;};
  auto accept=[&](int a,int b,double d){score+=d;for(int i:affected[a][b])cached[i]=value(i,key);if(score>best){best=score;bestkey=key;}};
  for(int cycle=0;cycle<3;cycle++){
   key=bestkey;score=full(key);for(int i=0;i<(int)grams.size();i++)cached[i]=value(i,key);
   for(int step=0;step<steps;step++){
    int a=sym(rng),b=sym(rng);if(a==b)continue;if(a>b)swap(a,b);
    double temp=mass*0.03*pow(0.0001/0.03,double(step)/(steps-1));
    swap(key[a],key[b]);double d=delta(a,b);
    if(d>=0||unif(rng)<exp(d/temp))accept(a,b,d);else swap(key[a],key[b]);
   }
  }
  key=bestkey;score=full(key);for(int i=0;i<(int)grams.size();i++)cached[i]=value(i,key);
  for(int round=0;round<1000;round++){
   double dmax=1e-8;int aa=-1,bb=-1;
   for(int a=0;a<23;a++)for(int b=a+1;b<23;b++){
    swap(key[a],key[b]);double d=delta(a,b);swap(key[a],key[b]);
    if(d>dmax){dmax=d;aa=a;bb=b;}
   }
   if(aa<0)break;
   swap(key[aa],key[bb]);accept(aa,bb,dmax);
  }
  best=full(bestkey);out<<restart<<' '<<best/mass;
  for(int x:bestkey)out<<' '<<x;
  out<<'\n';out.flush();
  cerr<<"restart "<<restart<<" score "<<best/mass<<'\n';
 }
}
