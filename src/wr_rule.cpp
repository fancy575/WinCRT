// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <limits>
using namespace Rcpp;

/*
  General K-tier win-ratio comparison (faithful to your 2-tier semantics)

  Expected DataFrame columns:
    - id    : subject id (numeric, assumed integer-like)
    - z     : group/arm (constant within id)
    - y     : time
    - tier  : 0 = censoring, 1..K = event tiers (larger code => HIGHER priority)

  Pair rule for (i, j):
    Let mC = min(C_i, C_j), where C_i is the latest y with tier==0 (∞ if none).
    For r = K, K-1, ..., 1:
      ti = earliest y with tier==r for i (∞ if none)
      tj = earliest y with tier==r for j
      If min(ti, tj) > mC: continue (no decision at this tier)
      Else if min(ti, tj) == mC: tie (stop)
      Else:
        - if ti < tj and ti < mC → i loses (j's opponent wins)
        - if tj < ti and tj < mC → i wins
        - else (ti == tj < mC) → tie (stop)
    If no tier decides → tie.

  Output: NumericVector length 4*n
    [1..n]       : wins per id
    [n+1..2n]    : losses per id
    [2n+1..3n]   : ties per id
    [3n+1..4n]   : ties vs opposite group per id
*/

static inline long long id_key(double x){ return (long long)x; }

// [[Rcpp::export]]
NumericVector wr_hier(DataFrame df){
  // Required columns (new names)
  if (!df.containsElementNamed("id") ||
      !df.containsElementNamed("z")  ||
      !df.containsElementNamed("y")  ||
      !df.containsElementNamed("tier")) {
    stop("DataFrame must contain columns: id, z, y, tier");
  }

  NumericVector id_col  = df["id"];
  IntegerVector z_col   = df["z"];
  NumericVector y_col   = df["y"];
  IntegerVector tier    = df["tier"];

  const int m = df.nrow();
  if ((int)id_col.size()!=m || (int)z_col.size()!=m ||
      (int)y_col.size()!=m  || (int)tier.size()!=m) {
    stop("Columns id, z, y, tier must have the same length");
  }

  // Map ids to 0..n-1 (first appearance order) and find K = max tier>=1
  std::unordered_map<long long,int> id2idx;
  id2idx.reserve((size_t)(m*1.3));
  std::vector<long long> ids; ids.reserve(m);

  int K = 0;
  for(int r=0;r<m;++r){
    long long key = id_key(id_col[r]);
    if (id2idx.find(key)==id2idx.end()){
      int idx = (int)ids.size();
      id2idx.emplace(key, idx);
      ids.push_back(key);
    }
    if (tier[r] > K) K = tier[r];
  }

  const int n = (int)ids.size();
  NumericVector out(4*n, 0.0);
  if (n==0) return out;

  const double INF = std::numeric_limits<double>::infinity();

  // Per-id storage
  std::vector<int> Z(n,0);
  std::vector<char> Zset(n,0);
  std::vector<double> C(n, INF); // latest censor time (tier==0)
  std::vector< std::vector<double> > T(n, std::vector<double>(std::max(1,K)+1, INF)); // T[i][r], r=1..K

  // Build per-id summaries from long data
  for(int r=0;r<m;++r){
    int i = id2idx[id_key(id_col[r])];
    double t  = y_col[r];
    int tr    = tier[r];
    int z     = z_col[r];

    if(!Zset[i]) { Z[i]=z; Zset[i]=1; }
    else if (Z[i] != z) stop("Within-id z must be constant (id=%lld)", ids[i]);

    if (tr == 0){
      if (C[i]==INF || t > C[i]) C[i] = t;        // latest censoring
    } else if (tr >= 1){
      if (tr <= K && t < T[i][tr]) T[i][tr] = t;  // earliest tier time
    } else {
      stop("tier must be >= 0");
    }
  }

  // Pairwise counting
  std::vector<long long> W(n,0), L(n,0), Ties(n,0);
  std::vector<long long> Wop(n,0), Lop(n,0);

  for(int i=0;i<n;++i){
    for(int j=i+1;j<n;++j){
      const double mC = std::min(C[i], C[j]);
      int winner = 0; // +1: i wins; -1: i loses; 0: tie

      for(int r=K; r>=1; --r){
        const double ti = T[i][r];
        const double tj = T[j][r];
        const double mn = (ti < tj ? ti : tj);

        if (mn > mC){
          // No event at this tier before earliest censor → go to lower tier
          continue;
        } else if (mn == mC){
          // Event coincides with earliest censor → cannot descend → tie
          winner = 0; break;
        } else {
          // Some tier-r event occurs before earliest censor
          if (ti < tj && ti < mC){
            winner = -1; // i’s earlier (worse for i) ⇒ i loses
          } else if (tj < ti && tj < mC){
            winner = +1; // j’s earlier ⇒ i wins
          } else {
            // ti == tj < mC → tie at this tier; do NOT drop to lower
            winner = 0;
          }
          break;
        }
      }

      if (winner == 0){
        Ties[i]++; Ties[j]++;
      } else if (winner > 0){
        W[i]++; L[j]++;
        if (Z[i] != Z[j]) { Wop[i]++; Lop[j]++; }
      } else {
        W[j]++; L[i]++;
        if (Z[i] != Z[j]) { Wop[j]++; Lop[i]++; }
      }
    }
  }

  // Opposite-group counts for tie_opposite
  std::unordered_map<int,int> gcnt;
  gcnt.reserve(n*2);
  for(int i=0;i<n;++i) gcnt[ Z[i] ]++;

  // Assemble output
  for(int i=0;i<n;++i){
    out[i]         = (double)W[i];
    out[n + i]     = (double)L[i];
    out[2*n + i]   = (double)Ties[i];

    long long n_op = (long long)n - (long long)gcnt[ Z[i] ];
    long long tie_op = n_op - Wop[i] - Lop[i];
    if (tie_op < 0) tie_op = 0;
    out[3*n + i]   = (double)tie_op;
  }

  return out;
}
