// Experiment022 draft shared types; decoder/search are not implemented yet.
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
 Inventory in;int n;f>>in.G>>in.duplicate;need(in.G>0&&in.G<=S,"Invalid inventory G");
 in.role.resize(in.G);in.multiplicity.resize(in.G);
 for(int g=0;g<in.G;g++)f>>in.role[g]>>in.multiplicity[g];
 f>>n;need(n>=0&&n<=S*S,"Invalid forbidden count");
 for(int i=0;i<n;i++){int p,s;f>>p>>s;need(p>=0&&p<in.G&&s>=0&&s<in.G&&in.role[p]==1&&in.role[s]==2,"Invalid forbidden pair");in.forbidden.emplace_back(p,s);}
 need(bool(f),"Truncated inventory");
 std::set<std::pair<int,int>> unique(in.forbidden.begin(),in.forbidden.end());need(unique.size()==in.forbidden.size(),"Repeated forbidden pair");
 return in;
}
inline Lattice read_lattice(const std::string&path,const Inventory&in){
 std::ifstream f(path);need(bool(f),"Cannot open observation lattice");
 int n;f>>n;need(n>0&&n<100000,"Invalid row count");Lattice out;
 for(int i=0;i<n;i++){
  Row row;int e;f>>row.positions>>row.tokens>>e;
  need(row.positions>1&&row.tokens>0&&e>0,"Invalid lattice row");
  for(int j=0;j<e;j++){
   Edge edge;f>>edge.from>>edge.to>>edge.first>>edge.second>>edge.token>>edge.constant;
   need(edge.from>=0&&edge.to>edge.from&&edge.to<row.positions&&edge.first>=0&&edge.first<in.G,"Invalid lattice edge");
   need(edge.token>=0&&edge.token<row.tokens&&std::isfinite(edge.constant),"Invalid edge metadata");
   need((edge.second==-1&&in.role[edge.first]==0)||(edge.second>=0&&edge.second<in.G&&in.role[edge.first]==1&&in.role[edge.second]==2),"Invalid U/PS roles");
   row.edges.push_back(edge);
  }
  out.rows.push_back(std::move(row));
 }
 need(bool(f),"Truncated lattice");return out;
}
inline Model read_model(const std::string&path){
 Model m;m.logp.resize(24*24*24);std::ifstream f(path,std::ios::binary);
 f.read(reinterpret_cast<char*>(m.logp.data()),m.logp.size()*sizeof(double));need(bool(f),"Truncated model");
 for(double x:m.logp)need(std::isfinite(x)&&x<=0,"Invalid model log probability");
 return m;
}
inline Key read_key(const std::string&path,const Inventory&in){
 Key k;std::ifstream f(path);for(int &g:k.slots)f>>g;
 need(bool(f)&&k.rebuild(in),"Invalid supplied candidate slots");return k;
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
