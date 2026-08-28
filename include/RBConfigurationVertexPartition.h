#ifndef RBCONFIGURATIONVERTEXPARTITION_H
#define RBCONFIGURATIONVERTEXPARTITION_H

#include "LinearResolutionParameterVertexPartition.h"

class LIBLEIDENALG_EXPORT RBConfigurationVertexPartition : public LinearResolutionParameterVertexPartition
{
  public:
    RBConfigurationVertexPartition(Graph* graph,
          vector<size_t> const& membership, double resolution_parameter);
    RBConfigurationVertexPartition(Graph* graph,
          vector<size_t> const& membership);
    RBConfigurationVertexPartition(Graph* graph,
      double resolution_parameter);
    RBConfigurationVertexPartition(Graph* graph);

    RBConfigurationVertexPartition(Graph* graph,
          vector<size_t> const& membership, double resolution_parameter,
          double pop_lambda, double pop_lambda2, double pop_threshold,
          int target_communities, double community_count_lambda,
          double eg_lambda, double eg_lambda2, double eg_target,
          double cont_lambda, double cont_lambda2,
          int pop_relative, int community_count_floor_only,
          double pop_count_exponent);
    RBConfigurationVertexPartition(Graph* graph,
          double resolution_parameter,
          double pop_lambda, double pop_lambda2, double pop_threshold,
          int target_communities, double community_count_lambda,
          double eg_lambda, double eg_lambda2, double eg_target,
          double cont_lambda, double cont_lambda2,
          int pop_relative, int community_count_floor_only,
          double pop_count_exponent);

    virtual ~RBConfigurationVertexPartition();
    virtual RBConfigurationVertexPartition* create(Graph* graph);
    virtual RBConfigurationVertexPartition* create(Graph* graph, vector<size_t> const& membership);

    virtual double diff_move(size_t v, size_t new_comm);
    virtual double quality(double resolution_parameter);

    double pop_lambda;
    double pop_lambda2;
    double pop_threshold;
    int target_communities;
    double community_count_lambda;
    double eg_lambda;
    double eg_lambda2;
    double eg_target;
    double cont_lambda;
    double cont_lambda2;
    // How the population penalty measures a community's deviation:
    //   0  absolute, p_c - target
    //   1  relative, p_c/target - 1, summed over communities
    //   2  relative and divided by the community count, i.e. the mean rather
    //      than the sum.  Equals the squared coefficient of variation of the
    //      community populations, so unlike 1 it is free of the community
    //      count and does not turn into a barrier against early merges.
    int pop_relative;
    // Exponent `a` in the mode-2 scale Lambda(K) = pop_lambda2 / K^a.  Read
    // only when pop_relative == 2; the other modes are unscaled.
    //
    // a = 1 (the default) is the mean form above, and is the only value that
    // is neutral in the community count: the same pop_lambda2 means the same
    // thing at every K, so a district count can be changed without
    // recalibrating it.  a > 1 makes the early barrier fall away faster than
    // the end-of-run one, which presses harder on the population spread at
    // the final K, but the same pop_lambda2 is then (K1/K2)^a times stronger
    // at K2 than at K1 and has to be re-tuned per district count.
    double pop_count_exponent;
    // When non-zero, the community-count penalty is one-way: it charges only
    // for having fewer communities than target_communities, never more.
    int community_count_floor_only;

  protected:
  private:
};

#endif // RBCONFIGURATIONVERTEXPARTITION_H
