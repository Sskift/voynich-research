// Experiment022 shared types; protocol must be frozen before real experiments.
// Public structural oracle, never a true key file.
#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace nu {
constexpr int A=23, T=6, R=3, S=R*T*A, H=24*24;
constexpr int WEIGHT[T]={20,8,8,8,4,4};
constexpr double NEG=-std::numeric_limits<double>::infinity();
inline void need(bool ok,const std::string&message){if(!ok)throw std::runtime_error(message);}

struct Inventory {
 int G=0, duplicate=-1;
 std::vector<int> role,multiplicity;
 std::vector<std::pair<int,int>> forbidden;
};
struct Edge {int from=0,to=0,first=-1,second=-1,token=-1;double constant=0;};
struct Row {int positions=0,tokens=0;std::vector<Edge> edges;};
struct Lattice {std::vector<Row> rows;};
struct Model {
 std::vector<double> logp;
 double at(int h,int c) const{return logp.at(h*24+c);}
};
struct Key {
 std::array<int,S> slots{};
 std::vector<int> letter,mass;
 std::array<int,A*A> blocked{};
 bool rebuild(const Inventory&in){
  letter.assign(in.G,-1);mass.assign(in.G,0);std::vector<int> counts(in.G,0);
  for(int p=0;p<S;p++){
   int g=slots[p],r=p/(T*A),a=p%A,t=(p/A)%T;
   if(g<0||g>=in.G||in.role[g]!=r)return false;
   if(letter[g]>=0&&letter[g]!=a)return false;
   letter[g]=a;mass[g]+=WEIGHT[t];counts[g]++;
  }
  for(int g=0;g<in.G;g++)if(counts[g]!=in.multiplicity[g]||letter[g]<0||mass[g]<=0)return false;
  blocked.fill(0);
  for(auto [p,s]:in.forbidden)blocked[letter[p]*A+letter[s]]+=mass[p]*mass[s];
  for(int n:blocked)if(n>=52*52||n<0)return false;
  return true;
 }
 double atom_log(int g)const{return std::log(double(mass[g])/52);}
 double rejection_log(int p,int s)const{
  return -std::log(double(52*52-blocked[letter[p]*A+letter[s]])/(52*52));
 }
};
struct DecodedRow {
 std::vector<int> edges,letters,atoms;
 double score=NEG;
 long long states=0,transitions=0;
};
struct Decoded {std::vector<DecodedRow> rows;double score=0;long long max_states=0;};

inline Inventory read_inventory(const std::string&path){
 std::ifstream f(path);need(bool(f),"Cannot open public inventory");
 Inventory in;int n=0;need(bool(f>>in.G>>in.duplicate),"Truncated inventory header");
 need(in.G==S-1&&in.duplicate>=0&&in.duplicate<137,"Invalid inventory shape");
 in.role.resize(in.G);in.multiplicity.resize(in.G);
 for(int g=0;g<in.G;g++){
  need(bool(f>>in.role[g]>>in.multiplicity[g]),"Truncated inventory atom");
  need(in.role[g]==(g<137?0:(g<275?1:2)),"Invalid canonical role order");
  need(in.multiplicity[g]==(g==in.duplicate?2:1),"Invalid atom multiplicity");
 }
 need(bool(f>>n)&&n>=0&&n<=138*138,"Invalid forbidden count");
 for(int i=0;i<n;i++){
  int p=-1,s=-1;need(bool(f>>p>>s),"Truncated forbidden pair");
  need(p>=0&&p<in.G&&s>=0&&s<in.G&&in.role[p]==1&&in.role[s]==2,"Invalid forbidden pair");
  in.forbidden.emplace_back(p,s);
 }
 std::set<std::pair<int,int>> unique(in.forbidden.begin(),in.forbidden.end());need(unique.size()==in.forbidden.size(),"Repeated forbidden pair");
 std::string extra;need(!(f>>extra),"Extra inventory content");
 return in;
}
inline Lattice read_lattice(const std::string&path,const Inventory&in){
 std::ifstream f(path);need(bool(f),"Cannot open observation lattice");
 int n=0;need(bool(f>>n)&&n>0&&n<=4096,"Invalid row count");Lattice out;
 for(int i=0;i<n;i++){
  Row row;int e=0;need(bool(f>>row.positions>>row.tokens>>e),"Truncated lattice row");
  need(row.positions>1&&row.positions<=100000&&row.tokens>0&&row.tokens<row.positions&&e>0&&e<=2000000,"Invalid lattice row");
  for(int j=0;j<e;j++){
   Edge edge;need(bool(f>>edge.from>>edge.to>>edge.first>>edge.second>>edge.token>>edge.constant),"Truncated lattice edge");
   need(edge.from>=0&&edge.to>edge.from&&edge.to<row.positions&&edge.first>=0&&edge.first<in.G,"Invalid lattice edge");
   need(edge.token>=0&&edge.token<row.tokens&&std::isfinite(edge.constant)&&edge.constant<=0,"Invalid edge metadata");
   need((edge.second==-1&&in.role[edge.first]==0)||(edge.second>=0&&edge.second<in.G&&in.role[edge.first]==1&&in.role[edge.second]==2),"Invalid U/PS roles");
   row.edges.push_back(edge);
  }
  out.rows.push_back(std::move(row));
 }
 std::string extra;need(!(f>>extra),"Extra lattice content");return out;
}
inline Model read_model(const std::string&path){
 Model m;m.logp.resize(24*24*24);std::ifstream f(path,std::ios::binary);
 f.read(reinterpret_cast<char*>(m.logp.data()),m.logp.size()*sizeof(double));need(bool(f),"Truncated model");
 for(double x:m.logp)need(std::isfinite(x)&&x<=0,"Invalid model log probability");
 need(f.peek()==std::char_traits<char>::eof(),"Extra model content");
 for(int h=0;h<H;h++){
  double total=0;for(int c=0;c<24;c++)total+=std::exp(m.at(h,c));
  need(std::abs(total-1)<1e-10,"Unnormalized model row");
 }
 return m;
}
inline Key read_key(const std::string&path,const Inventory&in){
 Key k;std::ifstream f(path);for(int &g:k.slots)f>>g;
 need(bool(f)&&k.rebuild(in),"Invalid supplied candidate slots");
 std::string extra;need(!(f>>extra),"Extra key content");return k;
}
inline double edge_channel(const Edge&e,const Key&k){
 double score=e.constant+k.atom_log(e.first);
 if(e.second>=0)score+=k.atom_log(e.second)+k.rejection_log(e.first,e.second);
 return score;
}

// Planned implementation in naibbe_unknown_decode.cpp, separate from optimizer.
Decoded decode(const Inventory&,const Lattice&,const Model&,const Key&,long long max_states=1000000);
std::vector<std::vector<int>> initial_parse(const Lattice&);
void decoder_self_check();
}
