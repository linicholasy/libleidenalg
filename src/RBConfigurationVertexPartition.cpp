#include "RBConfigurationVertexPartition.h"

#include <cmath>

namespace {

  // Deviation of a community population from the target.
  //
  // relative == 0 (default): p - target, the absolute form. The resulting
  // quadratic barrier scales as target^2, so it grows as the number of
  // communities falls and the target rises.
  //
  // relative != 0: p/target - 1, the scale-free form. The barrier is then
  // expressed in fractions of the target and is constant in the community
  // count.
  inline double pop_deviation(double p, double pop_target, int relative)
  {
    if (relative)
      return pop_target > 0.0 ? p / pop_target - 1.0 : 0.0;
    return p - pop_target;
  }

  // Deviation of the community count from the target.
  //
  // floor_only == 0 (default): |current - target|, penalising both directions.
  //
  // floor_only != 0: max(0, target - current), a one-way floor. The term is
  // exactly zero while the count is at or above the target, so it leaves the
  // early steps of the optimisation - which start from a very large number of
  // communities - undistorted.
  inline int community_count_deviation(int current_count, int target, int floor_only)
  {
    int deficit = target - current_count;
    if (floor_only)
      return deficit > 0 ? deficit : 0;
    return deficit > 0 ? deficit : -deficit;
  }

  // Scale applied to the summed population penalty.
  //
  // relative == 2: divide by the community count, so the penalty is a mean
  // rather than a sum.  Every vertex is assigned, so the mean community
  // population is exactly the target and the result is the squared
  // coefficient of variation of the community populations.
  //
  // That quantity is free of K, which the summed form is not.  Along a
  // near-equal trajectory any single merge makes one community twice the
  // target, a relative deviation of 1, and the summed form charges the same
  // for it at K = 2000 as at K = 10.  A lambda that controls the spread at
  // the end is then a wall against the first merge.  Dividing by K charges
  // lambda/K instead, so the early barrier falls with the community count
  // while the same relative error at the end still costs lambda/K_final.
  // `exponent` generalises that division to lambda/K^a.  a = 1 is the mean
  // form just described, and is the only count-neutral value.  a > 1 lowers
  // the early barrier further relative to the late one, so a lambda large
  // enough to hold the spread at the final K can still be passable at the
  // start; the price is that the meaning of lambda then depends on K.
  inline double pop_penalty_scale(int relative, size_t K, double exponent)
  {
    if (relative != 2 || K == 0)
      return 1.0;
    if (exponent == 1.0)                // the common case, and pow is not free
      return 1.0 / K;
    return 1.0 / pow((double) K, exponent);
  }

}

RBConfigurationVertexPartition::RBConfigurationVertexPartition(Graph* graph,
      vector<size_t> const& membership, double resolution_parameter) :
        LinearResolutionParameterVertexPartition(graph,
        membership, resolution_parameter),
        pop_lambda(0.0), pop_lambda2(0.0), pop_threshold(0.0),
        target_communities(0), community_count_lambda(0.0),
        eg_lambda(0.0), eg_lambda2(0.0), eg_target(0.0),
        cont_lambda(0.0), cont_lambda2(0.0),
        pop_relative(0), pop_count_exponent(1.0),
        community_count_floor_only(0)
{ }

RBConfigurationVertexPartition::RBConfigurationVertexPartition(Graph* graph,
      vector<size_t> const& membership) :
        LinearResolutionParameterVertexPartition(graph,
        membership),
        pop_lambda(0.0), pop_lambda2(0.0), pop_threshold(0.0),
        target_communities(0), community_count_lambda(0.0),
        eg_lambda(0.0), eg_lambda2(0.0), eg_target(0.0),
        cont_lambda(0.0), cont_lambda2(0.0),
        pop_relative(0), pop_count_exponent(1.0),
        community_count_floor_only(0)
{ }

RBConfigurationVertexPartition::RBConfigurationVertexPartition(Graph* graph,
      double resolution_parameter) :
        LinearResolutionParameterVertexPartition(graph, resolution_parameter),
        pop_lambda(0.0), pop_lambda2(0.0), pop_threshold(0.0),
        target_communities(0), community_count_lambda(0.0),
        eg_lambda(0.0), eg_lambda2(0.0), eg_target(0.0),
        cont_lambda(0.0), cont_lambda2(0.0),
        pop_relative(0), pop_count_exponent(1.0),
        community_count_floor_only(0)
{ }

RBConfigurationVertexPartition::RBConfigurationVertexPartition(Graph* graph) :
        LinearResolutionParameterVertexPartition(graph),
        pop_lambda(0.0), pop_lambda2(0.0), pop_threshold(0.0),
        target_communities(0), community_count_lambda(0.0),
        eg_lambda(0.0), eg_lambda2(0.0), eg_target(0.0),
        cont_lambda(0.0), cont_lambda2(0.0),
        pop_relative(0), pop_count_exponent(1.0),
        community_count_floor_only(0)
{ }

RBConfigurationVertexPartition::RBConfigurationVertexPartition(Graph* graph,
      vector<size_t> const& membership, double resolution_parameter,
      double pop_lambda, double pop_lambda2, double pop_threshold,
      int target_communities, double community_count_lambda,
      double eg_lambda, double eg_lambda2, double eg_target,
      double cont_lambda, double cont_lambda2,
      int pop_relative, int community_count_floor_only,
      double pop_count_exponent) :
        LinearResolutionParameterVertexPartition(graph,
        membership, resolution_parameter),
        pop_lambda(pop_lambda), pop_lambda2(pop_lambda2), pop_threshold(pop_threshold),
        target_communities(target_communities), community_count_lambda(community_count_lambda),
        eg_lambda(eg_lambda), eg_lambda2(eg_lambda2), eg_target(eg_target),
        cont_lambda(cont_lambda), cont_lambda2(cont_lambda2),
        pop_relative(pop_relative), pop_count_exponent(pop_count_exponent),
        community_count_floor_only(community_count_floor_only)
{ }

RBConfigurationVertexPartition::RBConfigurationVertexPartition(Graph* graph,
      double resolution_parameter,
      double pop_lambda, double pop_lambda2, double pop_threshold,
      int target_communities, double community_count_lambda,
      double eg_lambda, double eg_lambda2, double eg_target,
      double cont_lambda, double cont_lambda2,
      int pop_relative, int community_count_floor_only,
      double pop_count_exponent) :
        LinearResolutionParameterVertexPartition(graph, resolution_parameter),
        pop_lambda(pop_lambda), pop_lambda2(pop_lambda2), pop_threshold(pop_threshold),
        target_communities(target_communities), community_count_lambda(community_count_lambda),
        eg_lambda(eg_lambda), eg_lambda2(eg_lambda2), eg_target(eg_target),
        cont_lambda(cont_lambda), cont_lambda2(cont_lambda2),
        pop_relative(pop_relative), pop_count_exponent(pop_count_exponent),
        community_count_floor_only(community_count_floor_only)
{ }

RBConfigurationVertexPartition::~RBConfigurationVertexPartition()
{ }

RBConfigurationVertexPartition* RBConfigurationVertexPartition::create(Graph* graph)
{
  return new RBConfigurationVertexPartition(graph, this->resolution_parameter,
                                            this->pop_lambda, this->pop_lambda2, this->pop_threshold,
                                            this->target_communities, this->community_count_lambda,
                                            this->eg_lambda, this->eg_lambda2, this->eg_target,
                                            this->cont_lambda, this->cont_lambda2,
                                            this->pop_relative, this->community_count_floor_only,
                                            this->pop_count_exponent);
}

RBConfigurationVertexPartition* RBConfigurationVertexPartition::create(Graph* graph, vector<size_t> const& membership)
{
  return new RBConfigurationVertexPartition(graph, membership, this->resolution_parameter,
                                            this->pop_lambda, this->pop_lambda2, this->pop_threshold,
                                            this->target_communities, this->community_count_lambda,
                                            this->eg_lambda, this->eg_lambda2, this->eg_target,
                                            this->cont_lambda, this->cont_lambda2,
                                            this->pop_relative, this->community_count_floor_only,
                                            this->pop_count_exponent);
}

/*****************************************************************************
  Returns the difference in modularity if we move a node to a new community
*****************************************************************************/
double RBConfigurationVertexPartition::diff_move(size_t v, size_t new_comm)
{
  #ifdef DEBUG
    cerr << "double RBConfigurationVertexPartition::diff_move(" << v << ", " << new_comm << ")" << endl;
  #endif
  size_t old_comm = this->_membership[v];
  double diff = 0.0;
  double total_weight = this->graph->total_weight()*(2.0 - this->graph->is_directed());
  if (total_weight == 0.0)
    return 0.0;
  if (new_comm != old_comm)
  {
    #ifdef DEBUG
      cerr << "\t" << "old_comm: " << old_comm << endl;
    #endif
    double w_to_old = this->weight_to_comm(v, old_comm);
    #ifdef DEBUG
      cerr << "\t" << "w_to_old: " << w_to_old << endl;
    #endif
    double w_from_old = this->weight_from_comm(v, old_comm);
    #ifdef DEBUG
      cerr << "\t" << "w_from_old: " << w_from_old << endl;
    #endif
    double w_to_new = this->weight_to_comm(v, new_comm);
    #ifdef DEBUG
      cerr << "\t" << "w_to_new: " << w_to_new << endl;
    #endif
    double w_from_new = this->weight_from_comm(v, new_comm);
    #ifdef DEBUG
      cerr << "\t" << "w_from_new: " << w_from_new << endl;
    #endif
    double k_out = this->graph->strength(v, IGRAPH_OUT);
    #ifdef DEBUG
      cerr << "\t" << "k_out: " << k_out << endl;
    #endif
    double k_in = this->graph->strength(v, IGRAPH_IN);
    #ifdef DEBUG
      cerr << "\t" << "k_in: " << k_in << endl;
    #endif
    double self_weight = this->graph->node_self_weight(v);
    #ifdef DEBUG
      cerr << "\t" << "self_weight: " << self_weight << endl;
    #endif
    double K_out_old = this->total_weight_from_comm(old_comm);
    #ifdef DEBUG
      cerr << "\t" << "K_out_old: " << K_out_old << endl;
    #endif
    double K_in_old = this->total_weight_to_comm(old_comm);
    #ifdef DEBUG
      cerr << "\t" << "K_in_old: " << K_in_old << endl;
    #endif
    double K_out_new = this->total_weight_from_comm(new_comm) + k_out;
    #ifdef DEBUG
      cerr << "\t" << "K_out_new: " << K_out_new << endl;
    #endif
    double K_in_new = this->total_weight_to_comm(new_comm) + k_in;
    #ifdef DEBUG
      cerr << "\t" << "K_in_new: " << K_in_new << endl;
      cerr << "\t" << "total_weight: " << total_weight << endl;
    #endif
    double diff_old = (w_to_old - this->resolution_parameter*k_out*K_in_old/total_weight) + \
               (w_from_old - this->resolution_parameter*k_in*K_out_old/total_weight);
    #ifdef DEBUG
      cerr << "\t" << "diff_old: " << diff_old << endl;
    #endif
    double diff_new = (w_to_new + self_weight - this->resolution_parameter*k_out*K_in_new/total_weight) + \
               (w_from_new + self_weight - this->resolution_parameter*k_in*K_out_new/total_weight);
    #ifdef DEBUG
      cerr << "\t" << "diff_new: " << diff_new << endl;
    #endif
    diff = diff_new - diff_old;
    #ifdef DEBUG
      cerr << "\t" << "diff: " << diff << endl;
    #endif

    if (this->pop_lambda > 0.0 || this->pop_lambda2 > 0.0)
    {
      double node_pop      = this->graph->node_pop(v);
      double pop_old_before = this->cpop(old_comm);
      double pop_old_after  = pop_old_before - node_pop;
      double pop_new_before = this->cpop(new_comm);
      double pop_new_after  = pop_new_before + node_pop;

      // pop/num_comm
      //
      // This is a local approximation: K, and therefore the target and the
      // penalty scale, are read once and held fixed across the move, and only
      // the two affected communities are re-evaluated.  A move that empties or
      // creates a community changes K, which shifts the target for every other
      // community too; that is not accounted for here.  Doing so would cost a
      // pass over all K communities per candidate move.  The approximation
      // predates the relative modes and applies to the whole adaptive-target
      // branch; the final partition is still scored exactly by quality().
      size_t K = this->n_nonempty_communities();

      double pop_target;
      if (this->pop_threshold > 0.0)
        pop_target = this->pop_threshold;
      else
        pop_target = K > 0 ? this->graph->total_pop() / K : 0.0;

      // Read once, so it is identical on both sides of the move.
      double pop_scale = pop_penalty_scale(this->pop_relative, K,
                                           this->pop_count_exponent);

      auto calc_penalty = [&](double p) -> double {
        if (p == 0.0)
          return 0.0;
        double diff_pop = pop_deviation(p, pop_target, this->pop_relative);
        return pop_scale * (this->pop_lambda * abs(diff_pop)
                            + this->pop_lambda2 * diff_pop * diff_pop);
      };

      double penalty_before = calc_penalty(pop_old_before) + calc_penalty(pop_new_before);
      double penalty_after  = calc_penalty(pop_old_after)  + calc_penalty(pop_new_after);

      diff -= (penalty_after - penalty_before);
    }

    if (this->target_communities > 0 && this->community_count_lambda > 0.0)
    {
      int count_change = 0;

      double old_comm_size = this->csize(old_comm);
      if (old_comm_size == this->graph->node_size(v))
        count_change -= 1;

      double new_comm_size = this->csize(new_comm);
      if (new_comm_size == 0.0)
        count_change += 1;

      if (count_change != 0)
      {
        int current_count = 0;
        for (size_t c = 0; c < this->n_communities(); c++)
          if (this->csize(c) > 0)
            current_count++;

        int count_before = current_count;
        int count_after = current_count + count_change;

        int deviation_before = community_count_deviation(count_before,
              this->target_communities, this->community_count_floor_only);
        int deviation_after = community_count_deviation(count_after,
              this->target_communities, this->community_count_floor_only);

        double count_penalty_before = this->community_count_lambda * deviation_before;
        double count_penalty_after = this->community_count_lambda * deviation_after;

        diff -= (count_penalty_after - count_penalty_before);
      }
    }

/*****************************************************************************
  EG Constraint

  R - Rep.
  D - Dem.
  wR - Wasted Rep. Votes
  wD - Wasted Dem. Votes
  T - Total Votes
******************************************************************************/



        
    if ((this->eg_lambda > 0.0 || this->eg_lambda2 > 0.0) && this->graph->has_votes())
    {
      double r = this->graph->rep_votes(v);
      double d = this->graph->dem_votes(v);
      double cR_old = this->crep(old_comm), cD_old = this->cdem(old_comm);
      double cR_new = this->crep(new_comm), cD_new = this->cdem(new_comm);

      double W_R = this->W_R();
      double W_D = this->W_D();
      double V   = this->V_votes();

      double wR_old_b, wD_old_b, wR_new_b, wD_new_b;
      wasted_votes(cR_old, cD_old, wR_old_b, wD_old_b);
      wasted_votes(cR_new, cD_new, wR_new_b, wD_new_b);

      double wR_old_a, wD_old_a, wR_new_a, wD_new_a;
      wasted_votes(cR_old - r, cD_old - d, wR_old_a, wD_old_a);
      wasted_votes(cR_new + r, cD_new + d, wR_new_a, wD_new_a);

      double dWR = (wR_old_a + wR_new_a) - (wR_old_b + wR_new_b);
      double dWD = (wD_old_a + wD_new_a) - (wD_old_b + wD_new_b);

      if (V > 0.0)
      {
        double EG_before = (W_D - W_R) / V;
        double EG_after  = (W_D + dWD - W_R - dWR) / V;
        double db = EG_before - this->eg_target;
        double da = EG_after  - this->eg_target;
        double penalty_before = this->eg_lambda * abs(db) + this->eg_lambda2 * db * db;
        double penalty_after  = this->eg_lambda * abs(da) + this->eg_lambda2 * da * da;
        diff -= (penalty_after - penalty_before);
      }
    }

    if ((this->cont_lambda > 0.0 || this->cont_lambda2 > 0.0) && this->graph->has_neighbors())
    {
      double frag_before = this->frag();
      double frag_after  = frag_before + this->frag_delta(v, new_comm);
      double pb = this->cont_lambda * frag_before + this->cont_lambda2 * frag_before * frag_before;
      double pa = this->cont_lambda * frag_after  + this->cont_lambda2 * frag_after  * frag_after;
      diff -= (pa - pb);
    }
  }
  #ifdef DEBUG
    cerr << "exit RBConfigurationVertexPartition::diff_move(" << v << ", " << new_comm << ")" << endl;
    cerr << "return " << diff << endl << endl;
  #endif
  return diff;
}

/*****************************************************************************
  Give the modularity of the partition.

  We here use the unscaled version of modularity, in other words, we don"t
  normalise by the number of edges.
******************************************************************************/
double RBConfigurationVertexPartition::quality(double resolution_parameter)
{
  #ifdef DEBUG
    cerr << "double ModularityVertexPartition::quality()" << endl;
  #endif
  double mod = 0.0;

  double m;
  if (this->graph->is_directed())
    m = this->graph->total_weight();
  else
    m = 2*this->graph->total_weight();

  if (m == 0)
    return 0.0;

  for (size_t c = 0; c < this->n_communities(); c++)
  {
    double w = this->total_weight_in_comm(c);
    double w_out = this->total_weight_from_comm(c);
    double w_in = this->total_weight_to_comm(c);
    #ifdef DEBUG
      double csize = this->csize(c);
      cerr << "\t" << "Comm: " << c << ", size=" << csize << ", w=" << w << ", w_out=" << w_out << ", w_in=" << w_in << "." << endl;
    #endif
    mod += w - resolution_parameter*w_out*w_in/((this->graph->is_directed() ? 1.0 : 4.0)*this->graph->total_weight());

    if (this->pop_lambda > 0.0 || this->pop_lambda2 > 0.0)
    {
      double cpop_c = this->cpop(c);
      if (cpop_c == 0.0)
        continue;
      // pop/num_comm
      size_t K = this->n_nonempty_communities();

      double pop_target;
      if (this->pop_threshold > 0.0)
        pop_target = this->pop_threshold;
      else
        pop_target = K > 0 ? this->graph->total_pop() / K : 0.0;

      double pop_scale = pop_penalty_scale(this->pop_relative, K,
                                           this->pop_count_exponent);
      double diff_pop = pop_deviation(cpop_c, pop_target, this->pop_relative);
      mod -= pop_scale * (this->pop_lambda * abs(diff_pop)
                          + this->pop_lambda2 * diff_pop * diff_pop);
    }
  }

  if (this->target_communities > 0 && this->community_count_lambda > 0.0)
  {
    int current_count = 0;
    for (size_t c = 0; c < this->n_communities(); c++)
      if (this->csize(c) > 0)
        current_count++;

    int deviation = community_count_deviation(current_count,
          this->target_communities, this->community_count_floor_only);
    mod -= this->community_count_lambda * deviation;
  }

  if ((this->eg_lambda > 0.0 || this->eg_lambda2 > 0.0) && this->graph->has_votes())
  {
    double W_R = this->W_R();
    double W_D = this->W_D();
    double V   = this->V_votes();
    if (V > 0.0)
    {
      double EG = (W_D - W_R) / V;
      double dev = EG - this->eg_target;
      mod -= this->eg_lambda * abs(dev) + this->eg_lambda2 * dev * dev;
    }
  }

  if ((this->cont_lambda > 0.0 || this->cont_lambda2 > 0.0) && this->graph->has_neighbors())
  {
    double f = this->frag();
    mod -= this->cont_lambda * f + this->cont_lambda2 * f * f;
  }

  double q = (2.0 - this->graph->is_directed())*mod;
  #ifdef DEBUG
    cerr << "exit double RBConfigurationVertexPartition::quality()" << endl;
    cerr << "return " << q << endl << endl;
  #endif
  return q;
}
