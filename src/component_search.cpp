// Experiment014: same frozen011 search with a dynamic plaintext alphabet.
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
struct Gram{array<int,4> c;double n;};
int main(int argc,char**argv){
 if(argc!=11){cerr<<"lm.bin grams counts initial seed restarts steps cycles output alphabet_with_space\n";return 2;}
 int A=stoi(argv[10]),L=A-1;if(A<2||A>64)return 8;
 vector<double> lm(A*A*A*A);ifstream lf(argv[1],ios::binary);lf.read((char*)lm.data(),lm.size()*8);if(!lf)return 3;
 vector<Gram> grams;Gram g;ifstream gf(argv[2]);while(gf>>g.c[0]>>g.c[1]>>g.c[2]>>g.c[3]>>g.n)grams.push_back(g);
 ifstream cf(argv[3]);int G;cf>>G;vector<int> freq(G);for(int&n:freq)cf>>n;if(!cf||grams.empty())return 4;
 vector<int> initial(G+1);ifstream ini(argv[4]);for(int&x:initial)ini>>x;if(!ini||initial[G]!=L)return 5;
 int N=accumulate(freq.begin(),freq.end(),0);vector<double> nlog(N+1);for(int i=1;i<=N;i++)nlog[i]=i*log(double(i));
 double constant=0;for(int n:freq)constant+=nlog[n];double mass=0;for(auto&x:grams)mass+=x.n;
 vector<vector<int>> affected(G);for(int i=0;i<int(grams.size());i++){
  auto c=grams[i].c;sort(c.begin(),c.end());int last=-1;
  for(int a:c)if(a<G&&a!=last){affected[a].push_back(i);last=a;}
 }
 mt19937_64 rng(stoull(argv[5]));uniform_real_distribution<double> uniform(0,1);
 uniform_int_distribution<int> symbol(0,G-1),letter(0,L-1);discrete_distribution<int> weighted(freq.begin(),freq.end());
 int restarts=stoi(argv[6]),steps=stoi(argv[7]),cycles=stoi(argv[8]);ofstream out(argv[9]);out<<setprecision(17);
 auto val=[&](int i,const vector<int>& k){auto&c=grams[i].c;return lm[((k[c[0]]*A+k[c[1]])*A+k[c[2]])*A+k[c[3]]];};
 auto full=[&](const vector<int>& k){double s=constant;vector<int> totals(L,0);for(int c=0;c<G;c++)totals[k[c]]+=freq[c];
  for(int n:totals)s-=nlog[n];
  for(int i=0;i<int(grams.size());i++)s+=grams[i].n*val(i,k);
  return s;};
 vector<int> marks(grams.size(),0);int tick=0;
 for(int restart=0;restart<restarts;restart++){
  vector<int> key=initial;
  if(restart){
   vector<int> permutation(L); iota(permutation.begin(),permutation.end(),0);shuffle(permutation.begin(),permutation.end(),rng);
   for(int c=0;c<G;c++)key[c]=permutation[key[c]];
   // Alter group membership as well as names for half of the random starts.
   if(restart%2==0)for(int c=0;c<G;c++)if(uniform(rng)<.3)key[c]=letter(rng);
  }
  double current=full(key),best=current;vector<int> bestkey=key;vector<int> totals(L,0);
  vector<double> cache(grams.size());
  auto reset=[&](){fill(totals.begin(),totals.end(),0);for(int c=0;c<G;c++)totals[key[c]]+=freq[c];for(int i=0;i<int(grams.size());i++)cache[i]=val(i,key);current=full(key);};reset();
  vector<int> changed,old,impacted;vector<int> increments(L,0);
  auto proposal=[&](const vector<int>& which,const vector<int>& replacement,double temperature,bool force){
   changed=which;old.clear();impacted.clear();fill(increments.begin(),increments.end(),0);tick++;
   for(int j=0;j<int(changed.size());j++){
    int c=changed[j];old.push_back(key[c]);increments[key[c]]-=freq[c];increments[replacement[j]]+=freq[c];key[c]=replacement[j];
    for(int i:affected[c])if(marks[i]!=tick){marks[i]=tick;impacted.push_back(i);}
   }
   double delta=0;
   for(int l=0;l<L;l++)if(increments[l])delta+=nlog[totals[l]]-nlog[totals[l]+increments[l]];
   for(int i:impacted)delta+=grams[i].n*(val(i,key)-cache[i]);
   bool take=force||(temperature>0 ? delta>=0||uniform(rng)<exp(delta/temperature):delta>1e-8);
   if(take){
    current+=delta;for(int l=0;l<L;l++)totals[l]+=increments[l];for(int i:impacted)cache[i]=val(i,key);
    if(current>best){best=current;bestkey=key;}
   }else for(int j=0;j<int(changed.size());j++)key[changed[j]]=old[j];
   return take;
  };
  for(int cycle=0;cycle<cycles;cycle++){
   key=bestkey;reset();
   for(int step=0;step<steps;step++){
    double temperature=mass*.02*pow(.00005/.02,double(step)/(steps-1));double mode=uniform(rng);
    if(mode<.8){int c=uniform(rng)<.5?symbol(rng):weighted(rng),v=letter(rng);if(key[c]!=v)proposal({c},{v},temperature,false);}
    else if(mode<.95){int a=symbol(rng),b=symbol(rng);if(a!=b&&key[a]!=key[b])proposal({a,b},{key[b],key[a]},temperature,false);}
    else{
     int a=letter(rng),b=letter(rng);if(a==b)continue;vector<int> cs,vs;
     for(int c=0;c<G;c++)if(key[c]==a||key[c]==b){cs.push_back(c);vs.push_back(key[c]==a?b:a);}
     if(!cs.empty())proposal(cs,vs,temperature,false);
    }
   }
  }
  key=bestkey;reset();
  for(int pass=0;pass<100;pass++){
   bool improved=false;
   for(int c=0;c<G;c++)for(int v=0;v<L;v++)if(v!=key[c])improved=proposal({c},{v},0,false)||improved;
   if(!improved)break;
  }
  best=full(bestkey);out<<restart<<' '<<best/mass;for(int x:bestkey)out<<' '<<x;out<<'\n';out.flush();
  cerr<<"restart "<<restart<<" objective "<<best/mass<<'\n';
 }
}
