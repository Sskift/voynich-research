// Exploratory post-022 diagnostic. No truth input and no iterative fitting.
// Reuse the frozen and independently checked objective/cache without editing it.
#define main frozen022_main
#include "naibbe_unknown_search.cpp"
#undef main

namespace nu { namespace neighborhood {
using namespace search;

template<class F> void each_neighbor(const Inventory& in,const Key& key,F visit){
 for(int r=0;r<R;r++)for(int a=0;a<T*A;a++)for(int b=a+1;b<T*A;b++){
  Key k=key;std::swap(k.slots[r*T*A+a],k.slots[r*T*A+b]);visit(k,0);
 }
 for(int r=0;r<R;r++)for(int a=0;a<A;a++)for(int b=a+1;b<A;b++){
  Key k=key;for(int t=0;t<T;t++)std::swap(k.slots[r*T*A+t*A+a],k.slots[r*T*A+t*A+b]);visit(k,1);
 }
 for(int a=0;a<A;a++)for(int b=a+1;b<A;b++){
  Key k=key;for(int r=0;r<R;r++)for(int t=0;t<T;t++)std::swap(k.slots[r*T*A+t*A+a],k.slots[r*T*A+t*A+b]);visit(k,2);
 }
 for(int a=0;a<A;a++)for(int t1=0;t1<T;t1++)for(int t2=t1+1;t2<T;t2++){
  Key k=key;int x=t1*A+a,y=t2*A+a;
  for(int target:{x,y})if(k.slots[target]!=in.duplicate){
   int source=-1;for(int p=0;p<T*A;p++)if(k.slots[p]==in.duplicate&&p!=x&&p!=y){source=p;break;}
   need(source>=0,"Duplicate relocation source missing");std::swap(k.slots[target],k.slots[source]);
  }
  visit(k,3);
 }
}

void check(){
 Toy x=toy();std::array<long long,4> count{};long long legal=0;
 Paths paths=paths_of(decode(x.inventory,x.lattice,x.model,x.key));
 FixedParse cache(x.inventory,x.lattice,x.model,paths,x.key);
 double max_error=0;
 each_neighbor(x.inventory,x.key,[&](Key k,int kind){
  ++count[kind];if(!k.rebuild(x.inventory))return;++legal;
  Change change=cache.prepare(k);double reference=direct(x.inventory,x.lattice,x.model,k,paths).total();
  max_error=std::max(max_error,std::abs(change.components.total()-reference));
  close(change.components.total(),reference,"Neighborhood candidate direct");
 });
 need(count==std::array<long long,4>{28359,759,253,345},"Neighborhood cardinality");
 std::cout<<std::setprecision(17)<<"{\"status\":\"PASS\",\"synthetic_neighbors\":29716,\"legal\":"<<legal
  <<",\"maximum_score_error\":"<<max_error<<",\"real_inputs_read\":false}\n";
}

void scan(const Inventory&in,const Lattice&lat,const Model&lm,const Key&key,const std::string&output){
 auto started=Clock::now();Decoded decoded=decode(in,lat,lm,key);Paths paths=paths_of(decoded);
 Components baseline=verify_decoded(in,lat,lm,key,decoded);FixedParse cache(in,lat,lm,paths,key);
 Key best=key;double best_value=baseline.total();int best_kind=-1;
 std::array<long long,4> counts{},legal{},positive{},same{};double error=0;
 each_neighbor(in,key,[&](Key candidate,int kind){
  ++counts[kind];if(candidate.slots==key.slots){++same[kind];return;}
  if(!candidate.rebuild(in))return;
  ++legal[kind];
  Change change=cache.prepare(candidate);double value=change.components.total();
  if(value>baseline.total()+1e-7)++positive[kind];
  if(value>best_value+1e-7||(std::abs(value-best_value)<=1e-7&&value>baseline.total()+1e-7&&candidate.slots<best.slots)){
   best=candidate;best_value=value;best_kind=kind;
  }
 });
 Components direct_best=direct(in,lat,lm,best,paths);error=std::abs(best_value-direct_best.total());
 close(best_value,direct_best.total(),"Neighborhood best direct");
 Decoded best_decoded=decode(in,lat,lm,best);Components complete=verify_decoded(in,lat,lm,best,best_decoded);
 need(complete.total()+1e-7>=direct_best.total(),"Full decode below feasible path");
 std::ostringstream out;out<<std::setprecision(17)<<"{\"status\":\"COMPLETE_EXPLORATORY_TRAINING_NEIGHBORHOOD\",\"truth_read\":false"
  <<",\"score_scope\":"<<quoted(SCOPE)<<",\"baseline_objective\":"<<baseline.total()
  <<",\"best_fixed_path_objective\":"<<direct_best.total()<<",\"best_full_objective\":"<<complete.total()
  <<",\"best_move_kind\":"<<best_kind<<",\"best_score_error\":"<<error<<",\"initial_slots\":";integers(out,key.slots);
 out<<",\"best_slots\":";integers(out,best.slots);out<<",\"neighborhoods\":[";
 for(int k=0;k<4;k++){
  if(k)out<<',';
  out<<"{\"kind\":"<<quoted(move_name(k))<<",\"enumerated\":"<<counts[k]
   <<",\"same_slots\":"<<same[k]<<",\"legal_changed\":"<<legal[k]<<",\"strict_improvements\":"<<positive[k]<<'}';
 }
 out<<"],\"elapsed_seconds\":"<<seconds(started)<<"}\n";exclusive_write(output,out.str());
}
}}

int main(int argc,char**argv){
 try{
  if(argc==2&&std::string(argv[1])=="self-check"){nu::neighborhood::check();return 0;}
  nu::need(argc==6,"Usage: inventory lattice model slots output.json | self-check");
  nu::need(!std::filesystem::exists(argv[5]),"Refusing output overwrite");
  auto in=nu::read_inventory(argv[1]);auto lat=nu::read_lattice(argv[2],in);
  auto model=nu::read_model(argv[3]);auto key=nu::read_key(argv[4],in);
  nu::neighborhood::scan(in,lat,model,key,argv[5]);return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
