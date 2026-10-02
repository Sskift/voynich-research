// Experiment015 v1. Complete-word likelihood, fixed one/two-letter units.
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
constexpr int L=23,A=24,U=552,MG=64,MN=128;
struct Word {vector<int> c; int n;};
int main(int argc,char**argv){
 if(argc!=12){cerr<<"lm terms words counts initials seed restarts steps cycles polish output\n";return 2;}
 vector<double> lm(A*A*A*A),term((MG+1)*(MN+1));
 ifstream lf(argv[1],ios::binary),tf(argv[2],ios::binary);lf.read((char*)lm.data(),lm.size()*8);tf.read((char*)term.data(),term.size()*8);if(!lf||!tf)return 3;
 vector<Word> words;ifstream wf(argv[3]);int len,n;
 while(wf>>len>>n){if(len<1||len>MG||n<1)return 4;Word w;w.n=n;w.c.resize(len);for(int&c:w.c)wf>>c;if(!wf)return 4;words.push_back(w);}
 ifstream cf(argv[4]);int G;cf>>G;vector<int> freq(G);for(int&v:freq)cf>>v;if(!cf||words.empty())return 5;
 for(auto&w:words)for(int c:w.c)if(c<0||c>=G)return 5;
 int N=accumulate(freq.begin(),freq.end(),0);vector<double> nlog(N+1);for(int i=1;i<=N;i++)nlog[i]=i*log(double(i));double constant=0;for(int v:freq)constant+=nlog[v];
 vector<vector<int>> affected(G);for(int i=0;i<int(words.size());i++){auto cs=words[i].c;sort(cs.begin(),cs.end());cs.erase(unique(cs.begin(),cs.end()),cs.end());for(int c:cs)affected[c].push_back(i);}
 array<array<int,2>,U> unit;array<int,U> ulen;
 for(int u=0;u<U;u++){ulen[u]=u<L?1:2;unit[u]={u<L?u:(u-L)/L,u<L?-1:(u-L)%L};}
 auto uid=[&](int a,int b){return b<0?a:L+a*L+b;};
 auto val=[&](int i,const vector<int>&k){int a=L,b=L,c=L,m=0;double s=0;for(int x:words[i].c){int u=k[x];for(int j=0;j<ulen[u];j++){int d=unit[u][j];s+=lm[((a*A+b)*A+c)*A+d];a=b;b=c;c=d;m++;}}s+=lm[((a*A+b)*A+c)*A+L];return s+term[words[i].c.size()*(MN+1)+m];};
 auto full=[&](const vector<int>&k){double s=constant;vector<int> totals(U,0);for(int c=0;c<G;c++)totals[k[c]]+=freq[c];for(int v:totals)s-=nlog[v];for(int i=0;i<int(words.size());i++)s+=words[i].n*val(i,k);return s;};
 mt19937_64 rng(stoull(argv[6]));uniform_real_distribution<double> uniform(0,1);uniform_int_distribution<int> symbol(0,G-1),letter(0,L-1),anyunit(0,U-1);discrete_distribution<int> weighted(freq.begin(),freq.end());
 int restarts=stoi(argv[7]),steps=stoi(argv[8]),cycles=stoi(argv[9]),passes=stoi(argv[10]);ifstream ini(argv[5]);ofstream out(argv[11]);out<<setprecision(17);
 vector<int> marks(words.size(),0);int tick=0;
 for(int restart=0;restart<restarts;restart++){
  vector<int> key(G);for(int&v:key){ini>>v;if(v<0||v>=U)return 6;}if(!ini)return 6;
  double current=full(key),best=current;vector<int> bestkey=key,totals(U,0);vector<double> cache(words.size());
  auto reset=[&](){fill(totals.begin(),totals.end(),0);for(int c=0;c<G;c++)totals[key[c]]+=freq[c];for(int i=0;i<int(words.size());i++)cache[i]=val(i,key);current=full(key);};reset();
  vector<int> old,impacted,changedgroups;vector<int> increments(U,0);vector<double> newvals;
  auto proposal=[&](const vector<int>&which,const vector<int>&replacements,double temperature){
   old.clear();impacted.clear();changedgroups.clear();fill(increments.begin(),increments.end(),0);tick++;
   for(int j=0;j<int(which.size());j++){int c=which[j];old.push_back(key[c]);increments[key[c]]-=freq[c];increments[replacements[j]]+=freq[c];key[c]=replacements[j];for(int i:affected[c])if(marks[i]!=tick){marks[i]=tick;impacted.push_back(i);}}
   double delta=0;for(int u=0;u<U;u++)if(increments[u]){delta+=nlog[totals[u]]-nlog[totals[u]+increments[u]];changedgroups.push_back(u);}
   newvals.clear();for(int i:impacted){double v=val(i,key);newvals.push_back(v);delta+=words[i].n*(v-cache[i]);}
   bool take=temperature>0 ? delta>=0||uniform(rng)<exp(delta/temperature):delta>1e-8;
   if(take){current+=delta;for(int u:changedgroups)totals[u]+=increments[u];for(int j=0;j<int(impacted.size());j++)cache[impacted[j]]=newvals[j];if(current>best){best=current;bestkey=key;}}
   else for(int j=0;j<int(which.size());j++)key[which[j]]=old[j];
   return take;
  };
  for(int cycle=0;cycle<cycles;cycle++){
   key=bestkey;reset();
   for(int step=0;step<steps;step++){
    double t=N*.02*pow(.00005/.02,double(step)/(steps-1)),mode=uniform(rng);
    if(mode<.80){int c=uniform(rng)<.5?symbol(rng):weighted(rng),u=key[c],v;
     if(mode<.45){int a=unit[u][0],b=unit[u][1];if(b<0||uniform(rng)<.5)a=letter(rng);else b=letter(rng);v=uid(a,b);}
     else if(mode<.60){if(u<L)v=uniform(rng)<.5?uid(u,letter(rng)):uid(letter(rng),u);else v=unit[u][uniform(rng)<.5?0:1];}
     else v=anyunit(rng);
     if(v!=u)proposal({c},{v},t);
    }else if(mode<.95){int a=symbol(rng),b=symbol(rng);if(a!=b&&key[a]!=key[b])proposal({a,b},{key[b],key[a]},t);}
    else{int a=letter(rng),b=letter(rng);if(a==b)continue;vector<int> cs,vs;for(int c=0;c<G;c++){int x=unit[key[c]][0],y=unit[key[c]][1];auto swapname=[&](int z){return z==a?b:z==b?a:z;};int v=uid(swapname(x),swapname(y));if(v!=key[c]){cs.push_back(c);vs.push_back(v);}}if(!cs.empty())proposal(cs,vs,t);}
   }
  }
  key=bestkey;reset();int used=0;
  for(int pass=0;pass<passes;pass++){bool improved=false;used++;for(int c=0;c<G;c++)for(int u=0;u<U;u++)if(u!=key[c])improved=proposal({c},{u},0)||improved;if(!improved)break;}
  double checked=full(bestkey);if(abs(checked-best)>1e-5){cerr<<"CACHE MISMATCH "<<checked<<' '<<best<<'\n';return 7;}
  out<<restart<<' '<<checked/N;for(int v:bestkey)out<<' '<<v;out<<'\n';out.flush();cerr<<"restart "<<restart<<" objective "<<checked/N<<" polish_passes "<<used<<'\n';
 }
}
