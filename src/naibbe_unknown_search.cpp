// Experiment022: legal six-table key search; no true-key or plaintext input.
#include "naibbe_unknown.hpp"
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <map>
#include <random>
#include <sstream>
#include <fcntl.h>
#include <unistd.h>
#include <cerrno>
#include <cstring>

namespace nu {
namespace search {
using Clock=std::chrono::steady_clock;
using Paths=std::vector<std::vector<int>>;
constexpr double ABS_TOL=1e-7, REL_TOL=1e-12;
constexpr std::uint64_t RESTART_STRIDE=0x9e3779b97f4a7c15ULL;
constexpr const char* SCOPE="MAX_ATOM_PATH_APPROX_WITH_REPLACEMENT_NOT_MARGINAL";

double seconds(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}
double tolerance(double score){return std::max(ABS_TOL,REL_TOL*std::abs(score));}
void close(double actual,double expected,const std::string&where){
 need(std::isfinite(actual)&&std::isfinite(expected)&&std::abs(actual-expected)<=tolerance(expected),
      where+" score mismatch: "+std::to_string(actual)+" vs "+std::to_string(expected));
}

// Distribution helpers are explicit so the RNG stream does not depend on the
// implementation of std::uniform_int_distribution or std::shuffle.
struct RNG {
 std::mt19937_64 engine;
 explicit RNG(std::uint64_t seed):engine(seed){}
 std::uint64_t index(std::uint64_t n){
  need(n>0,"Empty RNG range");std::uint64_t threshold=(std::uint64_t(0)-n)%n,x;
  do{x=engine();}while(x<threshold);return x%n;
 }
 double unit(){return double(engine()>>11)*0x1.0p-53;}
 template<class V> void shuffle(V&v){for(std::size_t i=v.size();i>1;--i)std::swap(v[i-1],v[index(i)]);}
};

struct Components {
 double language=0,emission=0,rejection=0,constant=0;
 long long atoms=0,units=0,rows=0;
 double total()const{return language+emission+rejection+constant;}
};

// The direct reference deliberately walks every parsed edge and every emitted
// letter.  It does not use the gram cache, collapsed pair matrix or deltas.
Components direct(const Inventory&in,const Lattice&lat,const Model&lm,const Key&key,const Paths&paths){
 need(paths.size()==lat.rows.size(),"Parse row count mismatch");Components out;
 for(std::size_t r=0;r<paths.size();++r){
  const Row&row=lat.rows[r];int cursor=0,history=A*24+A;
  need(!paths[r].empty(),"Empty parsed row");
  auto emit=[&](int g){int a=key.letter.at(g);out.language+=lm.at(history,a);
   history=(history%24)*24+a;out.emission+=std::log(double(key.mass.at(g))/52);++out.atoms;};
  for(int id:paths[r]){
   need(id>=0&&id<int(row.edges.size()),"Parse edge outside lattice");const Edge&e=row.edges[id];
   need(e.from==cursor,"Disconnected parse");cursor=e.to;out.constant+=e.constant;++out.units;
   emit(e.first);
   if(e.second>=0){
    emit(e.second);int a=key.letter[e.first],b=key.letter[e.second],blocked=0;
    for(auto [p,s]:in.forbidden)if(key.letter[p]==a&&key.letter[s]==b)blocked+=key.mass[p]*key.mass[s];
    need(blocked>=0&&blocked<2704,"Invalid direct rejection normalizer");
    out.rejection-=std::log(double(2704-blocked)/2704);
   }
  }
  need(cursor==row.positions-1,"Parse does not reach row end");
  out.language+=lm.at(history,A);++out.rows;
 }
 need(std::isfinite(out.total()),"Nonfinite direct objective");return out;
}

Paths paths_of(const Decoded&decoded){Paths result;for(const auto&r:decoded.rows)result.push_back(r.edges);return result;}

Components verify_decoded(const Inventory&in,const Lattice&lat,const Model&lm,const Key&key,const Decoded&decoded){
 need(decoded.rows.size()==lat.rows.size(),"Incomplete decoded rows");
 for(std::size_t r=0;r<lat.rows.size();++r){
  std::vector<int> atoms,letters;
  for(int id:decoded.rows[r].edges){
   need(id>=0&&id<int(lat.rows[r].edges.size()),"Decoded edge index invalid");
   const Edge&e=lat.rows[r].edges[id];atoms.push_back(e.first);
   if(e.second>=0)atoms.push_back(e.second);
  }
  for(int g:atoms)letters.push_back(key.letter[g]);
  need(atoms==decoded.rows[r].atoms&&letters==decoded.rows[r].letters,"Decoded labels disagree with path");
  Lattice one;one.rows.push_back(lat.rows[r]);
  close(decoded.rows[r].score,direct(in,one,lm,key,Paths{decoded.rows[r].edges}).total(),"Decoded row");
 }
 Components out=direct(in,lat,lm,key,paths_of(decoded));close(decoded.score,out.total(),"Decoded total");return out;
}

struct Gram {std::array<int,3> atoms;long long count;};
struct Pair {int first,second;long long count;};
struct Change {
 Key key;Components components;std::array<long long,A*A> q{};
 std::vector<std::pair<int,double>> gram_values;
};

class FixedParse {
 const Inventory&in;const Lattice&lat;const Model&lm;const Paths&paths;
 std::vector<Gram> grams;std::vector<Pair> pairs;
 std::vector<std::vector<int>> affected_grams,affected_pairs;
 std::vector<double> values;
 std::vector<std::uint64_t> gram_marks,pair_marks;std::uint64_t tick=0;
 double gram_value(const Gram&g,const Key&k)const{
  auto letter=[&](int a){return a==in.G?A:k.letter[a];};
  return lm.at(letter(g.atoms[0])*24+letter(g.atoms[1]),letter(g.atoms[2]));
 }
public:
 Key key;Components components;std::vector<long long> atom_counts;
 std::array<long long,A*A> q{};
 FixedParse(const Inventory&i,const Lattice&l,const Model&m,const Paths&p,const Key&k)
 :in(i),lat(l),lm(m),paths(p),key(k),atom_counts(i.G,0){
  need(paths.size()==lat.rows.size(),"Fixed parse row count mismatch");
  std::map<std::array<int,3>,long long> gram_counts;
  std::map<std::pair<int,int>,long long> pair_counts;
  for(std::size_t r=0;r<paths.size();++r){
   int prev2=in.G,prev1=in.G,cursor=0;need(!paths[r].empty(),"Empty fixed parse");
   auto emit=[&](int g){++atom_counts[g];++gram_counts[{prev2,prev1,g}];prev2=prev1;prev1=g;++components.atoms;};
   for(int id:paths[r]){
    need(id>=0&&id<int(lat.rows[r].edges.size()),"Invalid fixed parse edge");
    const Edge&e=lat.rows[r].edges[id];need(e.from==cursor,"Disconnected fixed parse");cursor=e.to;
    components.constant+=e.constant;++components.units;emit(e.first);
    if(e.second>=0){emit(e.second);++pair_counts[{e.first,e.second}];}
   }
   need(cursor==lat.rows[r].positions-1,"Incomplete fixed parse");
   ++gram_counts[{prev2,prev1,in.G}];++components.rows;
  }
  affected_grams.resize(in.G);affected_pairs.resize(in.G);
  for(const auto&item:gram_counts){
   int id=int(grams.size());grams.push_back({item.first,item.second});
   std::set<int> unique(item.first.begin(),item.first.end());
   for(int g:unique)if(g!=in.G)affected_grams[g].push_back(id);
   double v=gram_value(grams.back(),key);values.push_back(v);components.language+=item.second*v;
  }
  for(const auto&item:pair_counts){
   int id=int(pairs.size()),p1=item.first.first,p2=item.first.second;pairs.push_back({p1,p2,item.second});
   affected_pairs[p1].push_back(id);affected_pairs[p2].push_back(id);
   q[key.letter[p1]*A+key.letter[p2]]+=item.second;
  }
  for(int g=0;g<in.G;++g)components.emission+=atom_counts[g]*key.atom_log(g);
  for(int cell=0;cell<A*A;++cell)if(q[cell])components.rejection-=q[cell]*std::log(double(2704-key.blocked[cell])/2704);
  gram_marks.assign(grams.size(),0);pair_marks.assign(pairs.size(),0);
  audit();
 }
 Change prepare(const Key&candidate){
  Change out;out.key=candidate;out.components=components;out.q=q;++tick;
  for(int g=0;g<in.G;++g){
   if(candidate.mass[g]!=key.mass[g])out.components.emission+=atom_counts[g]*std::log(double(candidate.mass[g])/key.mass[g]);
   if(candidate.letter[g]==key.letter[g])continue;
   for(int id:affected_grams[g])if(gram_marks[id]!=tick){
    gram_marks[id]=tick;double value=gram_value(grams[id],candidate);
    out.components.language+=grams[id].count*(value-values[id]);out.gram_values.emplace_back(id,value);
   }
   for(int id:affected_pairs[g])if(pair_marks[id]!=tick){
    pair_marks[id]=tick;const Pair&p=pairs[id];
    out.q[key.letter[p.first]*A+key.letter[p.second]]-=p.count;
    out.q[candidate.letter[p.first]*A+candidate.letter[p.second]]+=p.count;
   }
  }
  // A B change applies to ALL instances in its letter-pair cell, even when
  // none contains the moved surface atom. Q changes and B changes both matter.
  for(int cell=0;cell<A*A;++cell)if(out.q[cell]!=q[cell]||candidate.blocked[cell]!=key.blocked[cell]){
   need(out.q[cell]>=0,"Negative collapsed pair count");
   double old=-std::log(double(2704-key.blocked[cell])/2704);
   double value=-std::log(double(2704-candidate.blocked[cell])/2704);
   out.components.rejection+=out.q[cell]*value-q[cell]*old;
  }
  need(std::isfinite(out.components.total()),"Nonfinite incremental objective");return out;
 }
 void commit(Change&&change){
  key=std::move(change.key);components=change.components;q=change.q;
  for(auto [id,value]:change.gram_values)values[id]=value;
 }
 double audit(){
  Components exact=direct(in,lat,lm,key,paths);double error=std::abs(components.total()-exact.total());
  close(components.total(),exact.total(),"Fixed parse incremental/direct");
  close(components.language,exact.language,"Language component");
  close(components.emission,exact.emission,"Emission component");
  close(components.rejection,exact.rejection,"Rejection component");
  close(components.constant,exact.constant,"Constant component");
  need(components.atoms==exact.atoms&&components.units==exact.units&&components.rows==exact.rows,"Component counts differ");
  // Reset only after a passing comparison; never hide an out-of-tolerance drift.
  components=exact;return error;
 }
};

enum Move {SLOT_SWAP=0,ROLE_COLUMN=1,GLOBAL_COLUMN=2,DUPLICATE_BUNDLE=3};
const char* move_name(int m){const char*names[]={"slot_swap","role_column_swap","global_column_swap","duplicate_bundle"};return names[m];}

Key propose(const Inventory&in,const Key&current,RNG&rng,int&kind){
 Key candidate=current;int draw=int(rng.index(100));
 kind=draw<70?SLOT_SWAP:draw<85?ROLE_COLUMN:draw<95?GLOBAL_COLUMN:DUPLICATE_BUNDLE;
 if(kind==SLOT_SWAP){
  int role=int(rng.index(R)),a=int(rng.index(T*A)),b=int(rng.index(T*A-1));if(b>=a)++b;
  std::swap(candidate.slots[role*T*A+a],candidate.slots[role*T*A+b]);
 }else if(kind==ROLE_COLUMN||kind==GLOBAL_COLUMN){
  int role=kind==ROLE_COLUMN?int(rng.index(R)):0;
  int a=int(rng.index(A)),b=int(rng.index(A-1));if(b>=a)++b;
  for(int r=role;r<(kind==ROLE_COLUMN?role+1:R);++r)
   for(int t=0;t<T;++t)std::swap(candidate.slots[r*T*A+t*A+a],candidate.slots[r*T*A+t*A+b]);
 }else{
  int a=int(rng.index(A)),t1=int(rng.index(T)),t2=int(rng.index(T-1));if(t2>=t1)++t2;
  int x=t1*A+a,y=t2*A+a;
  for(int target:{x,y})if(candidate.slots[target]!=in.duplicate){
   int source=-1;
   for(int p=0;p<T*A;++p)if(candidate.slots[p]==in.duplicate&&p!=x&&p!=y){source=p;break;}
   need(source>=0,"Cannot move duplicate bundle");std::swap(candidate.slots[target],candidate.slots[source]);
  }
 }
 return candidate;
}

void validate_family(const Inventory&in){
 need(in.G>0&&int(in.role.size())==in.G&&int(in.multiplicity.size())==in.G,"Invalid family arrays");
 need(in.duplicate>=0&&in.duplicate<in.G&&in.role[in.duplicate]==0&&in.multiplicity[in.duplicate]==2,"Expected U duplicate");
 std::array<int,R> capacity{};
 for(int g=0;g<in.G;++g){
  need(in.role[g]>=0&&in.role[g]<R,"Invalid atom role");
  need(in.multiplicity[g]==(g==in.duplicate?2:1),"Unexpected atom multiplicity");capacity[in.role[g]]+=in.multiplicity[g];
 }
 for(int n:capacity)need(n==T*A,"Role slot capacity mismatch");
}

Key randomized(const Inventory&in,RNG&rng,long long&attempts){
 for(int attempt=0;attempt<10000;++attempt){
  ++attempts;Key k;k.slots.fill(-1);
  int a=int(rng.index(A)),t1=int(rng.index(T)),t2=int(rng.index(T-1));if(t2>=t1)++t2;
  k.slots[t1*A+a]=k.slots[t2*A+a]=in.duplicate;
  for(int r=0;r<R;++r){
   std::vector<int> atoms;for(int g=0;g<in.G;++g)if(in.role[g]==r&&g!=in.duplicate)atoms.push_back(g);
   rng.shuffle(atoms);std::size_t next=0;
   for(int p=r*T*A;p<(r+1)*T*A;++p)if(k.slots[p]<0){need(next<atoms.size(),"Random family overflow");k.slots[p]=atoms[next++];}
   need(next==atoms.size(),"Random family underflow");
  }
  if(k.rebuild(in))return k;
 }
 throw std::runtime_error("RANDOM_INITIALIZATION_BUDGET_EXCEEDED");
}

struct Statistics {
 long long proposals=0,accepted=0,illegal=0,no_op=0,direct_checks=0;
 std::array<long long,4> move_proposals{},move_accepted{},move_illegal{};
 double max_direct_error=0,elapsed_seconds=0;
};
struct Annealed {Key key;Components components;double start=0;Statistics stats;};

Annealed anneal(const Inventory&in,const Lattice&lat,const Model&lm,const Paths&paths,const Key&initial,RNG&rng,int steps){
 auto started=Clock::now();FixedParse cache(in,lat,lm,paths,initial);Annealed out;
 out.start=cache.components.total();out.key=initial;double best=out.start;
 const double scale=double(cache.components.atoms+cache.components.rows);
 for(int step=0;step<steps;++step){
  ++out.stats.proposals;int kind;Key candidate=propose(in,cache.key,rng,kind);++out.stats.move_proposals[kind];
  if(candidate.slots==cache.key.slots){++out.stats.no_op;}
  else if(!candidate.rebuild(in)){++out.stats.illegal;++out.stats.move_illegal[kind];}
  else{
   Change change=cache.prepare(candidate);double delta=change.components.total()-cache.components.total();
   double progress=steps==1?1.:double(step)/double(steps-1);
   double temperature=scale*.02*std::pow(.00005/.02,progress);
   bool accept=delta>=0;
   if(!accept){double u=rng.unit();accept=u==0||std::log(u)<delta/temperature;}
   if(accept){
    cache.commit(std::move(change));++out.stats.accepted;++out.stats.move_accepted[kind];
    double score=cache.components.total();
    if(score>best||(score==best&&cache.key.slots<out.key.slots)){best=score;out.key=cache.key;}
   }
  }
  if((step+1)%1000==0||step+1==steps){
   out.stats.max_direct_error=std::max(out.stats.max_direct_error,cache.audit());++out.stats.direct_checks;
  }
 }
 out.components=direct(in,lat,lm,out.key,paths);close(best,out.components.total(),"Best fixed parse");
 // Returning the starting key is allowed when all accepted moves were worse.
 need(out.components.total()+tolerance(out.start)>=out.start,"Annealing lost its initial best");
 out.stats.elapsed_seconds=seconds(started);return out;
}

struct RoundResult {
 int round=0;double fixed_start=0,fixed_best=0,decoded=0,elapsed_seconds=0;
 Key key;Components components;Statistics stats;long long max_states=0;
};
struct RestartResult {
 int restart=0,best_round=-1;std::uint64_t seed=0;std::string initialization;
 long long initialization_attempts=0;bool initialized=false,complete=false;
 Key initial_key,best_key;Components initial_components,best_components;
 double initial_objective=NEG,best_objective=NEG,elapsed_seconds=0;
 std::vector<RoundResult> rounds;
};
struct FitResult {
 std::uint64_t seed=0;int restarts_count=0,rounds_count=0,steps=0,selected_restart=-1,selected_round=-1;
 Key selected_key;Components selected_components;double selected_objective=NEG,elapsed_seconds=0;
 std::vector<RestartResult> restarts;
};

void fit(const Inventory&in,const Lattice&lat,const Model&lm,const Key&initial,FitResult&result){
 validate_family(in);auto started=Clock::now();
 for(int restart=0;restart<result.restarts_count;++restart){
  result.restarts.emplace_back();RestartResult&r=result.restarts.back();auto rstarted=Clock::now();
  r.restart=restart;r.seed=result.seed+RESTART_STRIDE*std::uint64_t(restart);RNG rng(r.seed);Key key=initial;
  if(restart==0)r.initialization="PYTHON_INITIAL_UNCHANGED";
  else if(restart%2){
   r.initialization="PYTHON_INITIAL_PLUS_64_VALID_NONIDENTICAL_RANDOM_MOVES";int valid=0;
   for(int attempt=0;attempt<6400&&valid<64;++attempt){
    ++r.initialization_attempts;int kind;Key next=propose(in,key,rng,kind);
    if(next.slots!=key.slots&&next.rebuild(in)){key=std::move(next);++valid;}
   }
   need(valid==64,"PERTURBED_INITIALIZATION_BUDGET_EXCEEDED");
  }else{r.initialization="UNIFORM_ROLE_SLOT_SHUFFLE_WITH_UNKNOWN_U_DUPLICATE_PAIR";key=randomized(in,rng,r.initialization_attempts);}
  r.initial_key=key;Decoded decoded=decode(in,lat,lm,key);r.initial_components=verify_decoded(in,lat,lm,key,decoded);
  r.initial_objective=r.initial_components.total();r.initialized=true;
  r.best_key=key;r.best_components=r.initial_components;r.best_objective=r.initial_objective;
  auto update_selected=[&](int round,const Key&k,const Components&components){
   double score=components.total();
   if(score>result.selected_objective||(score==result.selected_objective&&k.slots<result.selected_key.slots)){
    result.selected_objective=score;result.selected_key=k;result.selected_components=components;
    result.selected_restart=restart;result.selected_round=round;
   }
  };
  update_selected(-1,key,r.initial_components);
  for(int round=0;round<result.rounds_count;++round){
   auto round_started=Clock::now();Paths paths=paths_of(decoded);
   Annealed fitted=anneal(in,lat,lm,paths,key,rng,result.steps);
   Key candidate=fitted.key;Decoded next=decode(in,lat,lm,candidate);
   Components components=verify_decoded(in,lat,lm,candidate,next);
   need(components.total()+tolerance(fitted.components.total())>=fitted.components.total(),"Full decode below fixed-parse path");
   RoundResult rr;rr.round=round;rr.fixed_start=fitted.start;rr.fixed_best=fitted.components.total();
   rr.decoded=components.total();rr.key=candidate;rr.components=components;rr.stats=fitted.stats;
   rr.max_states=next.max_states;rr.elapsed_seconds=seconds(round_started);r.rounds.push_back(rr);
   if(components.total()>r.best_objective||(components.total()==r.best_objective&&candidate.slots<r.best_key.slots)){
    r.best_objective=components.total();r.best_components=components;r.best_key=candidate;r.best_round=round;
   }
   update_selected(round,candidate,components);
   // Every new fixed parse begins with the complete decode under the current
   // restart's best full-objective key, never a fragment/profile-selected key.
   if(candidate.slots==r.best_key.slots){key=candidate;decoded=std::move(next);}
   else{key=r.best_key;decoded=decode(in,lat,lm,key);verify_decoded(in,lat,lm,key,decoded);}
   std::cerr<<"restart "<<restart<<" round "<<round<<" training_objective "<<std::setprecision(17)<<rr.decoded<<"\n";
  }
  r.complete=true;r.elapsed_seconds=seconds(rstarted);
 }
 need(int(result.restarts.size())==result.restarts_count,"Missing restart");
 for(const auto&r:result.restarts)need(r.complete&&int(r.rounds.size())==result.rounds_count,"Incomplete restart");
 need(result.selected_restart>=0&&std::isfinite(result.selected_objective),"No selected training key");
 result.elapsed_seconds=seconds(started);
}

std::string quoted(const std::string&s){
 std::ostringstream out;out<<'"';
 for(unsigned char c:s){
  if(c=='"'||c=='\\')out<<'\\'<<c;
  else if(c=='\n')out<<"\\n";else if(c=='\r')out<<"\\r";else if(c=='\t')out<<"\\t";
  else if(c<32)out<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<int(c)<<std::dec;
  else out<<c;
 }
 return out.str()+'"';
}
template<class V>void integers(std::ostream&o,const V&v){o<<'[';bool first=true;for(auto n:v){if(!first)o<<',';first=false;o<<n;}o<<']';}
void components_json(std::ostream&o,const Components&c){
 o<<"{\"language\":"<<c.language<<",\"emission\":"<<c.emission<<",\"rejection\":"<<c.rejection
  <<",\"constant\":"<<c.constant<<",\"total\":"<<c.total()<<",\"atoms\":"<<c.atoms<<",\"units\":"<<c.units<<",\"rows\":"<<c.rows<<'}';
}
void statistics_json(std::ostream&o,const Statistics&s){
 o<<"\"proposals\":"<<s.proposals<<",\"accepted\":"<<s.accepted<<",\"illegal\":"<<s.illegal
  <<",\"no_op\":"<<s.no_op<<",\"direct_checks\":"<<s.direct_checks<<",\"max_direct_error\":"<<s.max_direct_error
  <<",\"anneal_elapsed_seconds\":"<<s.elapsed_seconds<<",\"moves\":{";
 for(int m=0;m<4;++m){if(m)o<<',';o<<quoted(move_name(m))<<":{\"proposals\":"<<s.move_proposals[m]
  <<",\"accepted\":"<<s.move_accepted[m]<<",\"illegal\":"<<s.move_illegal[m]<<'}';}o<<'}';
}

std::string fit_json(const FitResult&f,bool complete,const std::string&error=""){
 std::ostringstream o;o<<std::setprecision(17);
 o<<"{\"status\":"<<quoted(complete?"FIT_COMPLETE_TRAIN_ONLY":"FAILED")<<",\"complete\":"<<(complete?"true":"false")
  <<",\"score_scope\":"<<quoted(SCOPE)<<",\"seed\":"<<f.seed<<",\"restarts_count\":"<<f.restarts_count
  <<",\"rounds_count\":"<<f.rounds_count<<",\"steps_per_round\":"<<f.steps<<",\"completed_restarts\":";
 int completed=0;for(const auto&r:f.restarts)completed+=r.complete;o<<completed;
 if(!error.empty())o<<",\"error\":"<<quoted(error);
 o<<",\"selection_rule\":\"Highest complete training decode objective; exact ties lexicographic slots; no truth\""
  <<",\"rng_rule\":\"mt19937_64(seed + 0x9e3779b97f4a7c15*restart modulo 2^64); explicit rejection-sampled bounded integers; top53-bit uniforms; Fisher-Yates; acceptance draw only for negative delta\""
  <<",\"temperature_rule\":\"(fixed_parse_atoms+rows)*.02*(.00005/.02)^(step/(steps-1)); one-step progress=1\""
  <<",\"verification\":{\"interval_proposals\":1000,\"absolute_tolerance\":"<<ABS_TOL<<",\"relative_tolerance\":"<<REL_TOL
  <<",\"rule\":\"max(absolute_tolerance,relative_tolerance*abs(direct_score)); components also checked\"}"
  <<",\"selected_restart\":"<<f.selected_restart<<",\"selected_round\":"<<f.selected_round;
 if(f.selected_restart>=0){o<<",\"selected_slots\":";integers(o,f.selected_key.slots);
  o<<",\"selected_objective\":"<<f.selected_objective<<",\"selected_components\":";components_json(o,f.selected_components);}
 else o<<",\"selected_slots\":null,\"selected_objective\":null,\"selected_components\":null";
 o<<",\"restarts\":[";
 for(std::size_t i=0;i<f.restarts.size();++i){
  if(i)o<<',';
  const auto&r=f.restarts[i];o<<"{\"restart\":"<<r.restart<<",\"seed\":"<<r.seed
   <<",\"initialization\":"<<quoted(r.initialization)<<",\"initialization_attempts\":"<<r.initialization_attempts
   <<",\"complete\":"<<(r.complete?"true":"false");
  if(r.initialized){o<<",\"initial_slots\":";integers(o,r.initial_key.slots);o<<",\"initial_objective\":"<<r.initial_objective
   <<",\"initial_components\":";components_json(o,r.initial_components);o<<",\"best_round\":"<<r.best_round
   <<",\"best_objective\":"<<r.best_objective<<",\"best_slots\":";integers(o,r.best_key.slots);}
  o<<",\"rounds\":[";
  for(std::size_t j=0;j<r.rounds.size();++j){
   if(j)o<<',';
   const auto&x=r.rounds[j];o<<"{\"round\":"<<x.round<<",\"fixed_start_objective\":"<<x.fixed_start
    <<",\"fixed_best_objective\":"<<x.fixed_best<<",\"decoded_objective\":"<<x.decoded<<",\"components\":";
   components_json(o,x.components);o<<",\"slots\":";integers(o,x.key.slots);o<<',';statistics_json(o,x.stats);
   o<<",\"max_decode_states\":"<<x.max_states<<",\"elapsed_seconds\":"<<x.elapsed_seconds<<'}';
  }
  o<<"],\"elapsed_seconds\":"<<r.elapsed_seconds<<'}';
 }
 o<<"],\"elapsed_seconds\":"<<f.elapsed_seconds<<"}\n";return o.str();
}

std::string decode_json(const Key&key,const Decoded&decoded,const Components&components,double elapsed){
 std::ostringstream o;o<<std::setprecision(17)<<"{\"status\":\"DECODE_COMPLETE\",\"complete\":true,\"score_scope\":"<<quoted(SCOPE)<<",\"slots\":";
 integers(o,key.slots);o<<",\"score\":"<<decoded.score<<",\"components\":";components_json(o,components);
 o<<",\"max_states\":"<<decoded.max_states<<",\"rows\":[";
 for(std::size_t i=0;i<decoded.rows.size();++i){if(i)o<<',';const auto&r=decoded.rows[i];
  o<<"{\"edges\":";integers(o,r.edges);o<<",\"atoms\":";integers(o,r.atoms);o<<",\"letters\":";integers(o,r.letters);
  o<<",\"score\":"<<r.score<<",\"states\":"<<r.states<<",\"transitions\":"<<r.transitions<<'}';}
 o<<"],\"elapsed_seconds\":"<<elapsed<<"}\n";return o.str();
}

void exclusive_write(const std::string&path,const std::string&text){
 int fd=::open(path.c_str(),O_WRONLY|O_CREAT|O_EXCL,0644);
 need(fd>=0,"Exclusive output create failed: "+path+" ("+std::strerror(errno)+")");
 std::size_t written=0;
 while(written<text.size()){
  ssize_t n=::write(fd,text.data()+written,text.size()-written);
  if(n<0&&errno==EINTR)continue;
  if(n<=0){int e=errno;::close(fd);throw std::runtime_error("Output write failed: "+std::string(std::strerror(e)));}
  written+=std::size_t(n);
 }
 need(::close(fd)==0,"Output close failed");
}

struct Toy {Inventory inventory;Key key;Model model;Lattice lattice;Paths paths;};
Toy toy(){
 Toy x;Inventory&in=x.inventory;in.G=S-1;in.duplicate=0;in.role.resize(in.G);in.multiplicity.assign(in.G,1);in.multiplicity[0]=2;
 x.key.slots.fill(-1);x.key.slots[0]=x.key.slots[A]=0;int g=1;
 for(int p=0;p<T*A;++p)if(x.key.slots[p]<0)x.key.slots[p]=g++;
 for(int r=1;r<R;++r)for(int p=r*T*A;p<(r+1)*T*A;++p){x.key.slots[p]=g;in.role[g++]=r;}
 need(g==in.G,"Toy family construction");
 for(int k=0;k<9;++k)in.forbidden.emplace_back(x.key.slots[T*A+k],x.key.slots[2*T*A+k]);
 need(x.key.rebuild(in),"Illegal toy key");validate_family(in);
 x.model.logp.resize(24*24*24);
 for(int h=0;h<H;++h){double total=0;std::array<double,24>w{};
  for(int c=0;c<24;++c){w[c]=1+(h*7+c*11+(h/24)*c)%17;total+=w[c];}
  for(int c=0;c<24;++c)x.model.logp[h*24+c]=std::log(w[c]/total);
 }
 RNG rng(72211);
 for(int r=0;r<9;++r){Row row;row.positions=13+r;row.tokens=row.positions-1;std::vector<int> path;
  for(int pos=0;pos<row.positions-1;++pos){
   int u=x.key.slots[rng.index(T*A)],p,s;
   do{p=x.key.slots[T*A+rng.index(T*A)];s=x.key.slots[2*T*A+rng.index(T*A)];}
   while(std::find(in.forbidden.begin(),in.forbidden.end(),std::make_pair(p,s))!=in.forbidden.end());
   int start=int(row.edges.size());
   row.edges.push_back({pos,pos+1,u,-1,pos,pos+1==row.positions-1?0.:std::log(17./36)});
   row.edges.push_back({pos,pos+1,p,s,pos,std::log(19./36)});path.push_back(start+(pos%3==0));
   if(pos+2<row.positions)row.edges.push_back({pos,pos+2,p,s,pos,std::log(19./36)+std::log(.03)});
  }
  x.lattice.rows.push_back(row);x.paths.push_back(path);
 }
 return x;
}

void self_check(){
 auto started=Clock::now();Toy x=toy();RNG rng(9022022);FixedParse cache(x.inventory,x.lattice,x.model,x.paths,x.key);
 long long checked=0,illegal=0;double max_error=0;std::array<long long,4> moves{};std::set<std::tuple<int,int,int>> duplicate_places;
 for(int i=0;i<20000;++i){
  int kind;Key candidate=propose(x.inventory,cache.key,rng,kind);++moves[kind];
  if(!candidate.rebuild(x.inventory)){++illegal;continue;}
  Change change=cache.prepare(candidate);Components exact=direct(x.inventory,x.lattice,x.model,candidate,x.paths);
  double error=std::abs(change.components.total()-exact.total());max_error=std::max(max_error,error);
  close(change.components.total(),exact.total(),"Synthetic candidate total");
  close(change.components.language,exact.language,"Synthetic candidate LM");
  close(change.components.emission,exact.emission,"Synthetic candidate emission");
  close(change.components.rejection,exact.rejection,"Synthetic candidate rejection");++checked;
  if(i%5!=0)cache.commit(std::move(change)); // Exercise rejected updates as well.
  if(i%1000==0)cache.audit();
  int a=-1,t1=-1,t2=-1;
  for(int p=0;p<T*A;++p)if(cache.key.slots[p]==x.inventory.duplicate){a=p%A;if(t1<0)t1=p/A;else t2=p/A;}
  duplicate_places.emplace(a,t1,t2);
 }
 cache.audit();for(long long n:moves)need(n>0,"Untested move class");
 need(duplicate_places.size()>100,"Duplicate bundle did not explore letters/table pairs");

 // Changing an UNUSED forbidden-prefix atom must change the rejection score
 // of a DIFFERENT, observed P/S pair in the same collapsed letter cell.
 Inventory in=x.inventory;in.forbidden={{x.key.slots[T*A],x.key.slots[2*T*A]}};Key old=x.key;need(old.rebuild(in),"Toy old key");
 int p=x.key.slots[T*A+A],s=x.key.slots[2*T*A+A];Lattice lat;Row row;row.positions=2;row.tokens=1;
 row.edges.push_back({0,1,p,s,0,0.});lat.rows.push_back(row);Paths paths{{0}};FixedParse indirect(in,lat,x.model,paths,old);
 Key next=old;std::swap(next.slots[T*A],next.slots[T*A+1]);need(next.rebuild(in),"Toy next key");
 need(indirect.atom_counts[old.slots[T*A]]==0&&indirect.atom_counts[old.slots[T*A+1]]==0,"Indirect test atoms were used");
 Change change=indirect.prepare(next);Components exact=direct(in,lat,x.model,next,paths);
 close(change.components.total(),exact.total(),"Indirect Z update");
 need(std::abs(change.components.rejection-indirect.components.rejection)>1e-3,"Indirect Z effect missing");
 close(change.components.language,indirect.components.language,"Indirect LM unchanged");
 close(change.components.emission,indirect.components.emission,"Indirect emission unchanged");

 Inventory unsupported=x.inventory;unsupported.forbidden.clear();
 for(int t=0;t<T;++t)for(int u=0;u<T;++u)unsupported.forbidden.emplace_back(x.key.slots[T*A+t*A],x.key.slots[2*T*A+u*A]);
 Key bad=x.key;need(!bad.rebuild(unsupported),"Zero-support key accepted");
 Key split=x.key;std::swap(split.slots[0],split.slots[1]);need(!split.rebuild(x.inventory),"Split duplicate accepted");
 long long attempts=0;for(int i=0;i<20;++i){Key k=randomized(x.inventory,rng,attempts);need(k.rebuild(x.inventory),"Random key illegal");}

 RNG arng(772);Annealed a=anneal(x.inventory,x.lattice,x.model,x.paths,x.key,arng,2100);
 need(a.stats.proposals==2100&&a.stats.direct_checks==3,"Annealing audit interval mismatch");
 need(a.stats.accepted+a.stats.illegal+a.stats.no_op<=a.stats.proposals,"Invalid proposal bookkeeping");
 decoder_self_check();Decoded d=decode(x.inventory,x.lattice,x.model,a.key);verify_decoded(x.inventory,x.lattice,x.model,a.key,d);
 FitResult integrated;integrated.seed=78123;integrated.restarts_count=3;integrated.rounds_count=2;integrated.steps=200;
 fit(x.inventory,x.lattice,x.model,x.key,integrated);
 need(integrated.restarts.size()==3,"Synthetic fit restart count");
 double maximum=NEG;
 for(const auto&r:integrated.restarts){
  need(r.complete&&r.rounds.size()==2,"Synthetic fit incomplete");maximum=std::max(maximum,r.initial_objective);
  for(const auto&rr:r.rounds){maximum=std::max(maximum,rr.decoded);need(rr.stats.proposals==200,"Synthetic fit step count");}
 }
 close(integrated.selected_objective,maximum,"Selected complete training maximum");
 Decoded selected=decode(x.inventory,x.lattice,x.model,integrated.selected_key);
 close(verify_decoded(x.inventory,x.lattice,x.model,integrated.selected_key,selected).total(),integrated.selected_objective,"Selected key redecode");
 std::string serialized=fit_json(integrated,true);
 need(serialized.find("\"completed_restarts\":3")!=std::string::npos&&serialized.find("\"status\":\"FIT_COMPLETE_TRAIN_ONLY\"")!=std::string::npos,"Fit serializer completion");
 need(fit_json(integrated,false,"synthetic failure").find("\"status\":\"FAILED\"")!=std::string::npos,"Failure serializer status");
 std::cout<<std::setprecision(17)<<"{\"status\":\"PASS\",\"synthetic_candidate_checks\":"<<checked
  <<",\"synthetic_illegal_proposals\":"<<illegal<<",\"max_incremental_direct_error\":"<<max_error
  <<",\"duplicate_letter_table_configurations\":"<<duplicate_places.size()
  <<",\"indirect_Z_update_checked\":true,\"zero_support_rejection_checked\":true,\"split_duplicate_rejection_checked\":true"
  <<",\"anneal_proposals\":"<<a.stats.proposals<<",\"anneal_direct_checks\":"<<a.stats.direct_checks
  <<",\"integrated_restarts\":3,\"integrated_rounds\":2,\"selected_key_redecode_checked\":true"
  <<",\"elapsed_seconds\":"<<seconds(started)<<",\"real_inputs_read\":false}\n";
}
} // namespace search
} // namespace nu

int main(int argc,char**argv){
 using namespace nu;using namespace nu::search;std::string output;FitResult result;bool fitting=false;
 try{
  if(argc==2&&std::string(argv[1])=="self-check"){self_check();return 0;}
  need(argc>=2,"Usage: fit inventory lattice lm.bin initial-slots seed restarts rounds steps output.json | decode inventory lattice lm.bin slots output.json | self-check");
  std::string command=argv[1];
  need((command=="fit"&&argc==11)||(command=="decode"&&argc==7),"Invalid command or argument count");
  output=argv[argc-1];need(!std::filesystem::exists(output),"Output already exists: "+output);
  Inventory in=read_inventory(argv[2]);validate_family(in);Lattice lattice=read_lattice(argv[3],in);Model lm=read_model(argv[4]);Key key=read_key(argv[5],in);
  if(command=="decode"){
   auto started=Clock::now();Decoded decoded=decode(in,lattice,lm,key);Components components=verify_decoded(in,lattice,lm,key,decoded);
   exclusive_write(output,decode_json(key,decoded,components,seconds(started)));return 0;
  }
  auto number=[](const char*s){std::string value=s;std::size_t used=0;need(!value.empty()&&value.find_first_not_of("0123456789")==std::string::npos,"Invalid numeric argument");
   auto n=std::stoull(value,&used);need(used==value.size(),"Invalid numeric suffix");return n;};
  result.seed=number(argv[6]);auto nr=number(argv[7]),rounds=number(argv[8]),steps=number(argv[9]);
  need(nr>0&&nr<=1000&&rounds>0&&rounds<=1000&&steps>0&&steps<=100000000,"Invalid search budget");
  result.restarts_count=int(nr);result.rounds_count=int(rounds);result.steps=int(steps);fitting=true;
  fit(in,lattice,lm,key,result);exclusive_write(output,fit_json(result,true));return 0;
 }catch(const std::exception&e){
  std::cerr<<"FAILED: "<<e.what()<<'\n';
  if(!output.empty()&&!std::filesystem::exists(output)){
   try{exclusive_write(output,fitting?fit_json(result,false,e.what()):
    "{\"status\":\"FAILED\",\"complete\":false,\"score_scope\":"+quoted(SCOPE)+",\"error\":"+quoted(e.what())+"}\n");}
   catch(const std::exception&write_error){std::cerr<<"Failure record not written: "<<write_error.what()<<'\n';}
  }
  return 1;
 }
}
