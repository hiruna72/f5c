//
// Created by Hiruna Samrakoon on 2020-03-20.
// adapted from nanopolish polya_estimate algorithm ( mostly the same code)
//



#include <algorithm>
#include "f5c.h"

// ================================================================================
// Segmentation Hidden Markov Model
//   Define an HMM class `SegmentationHMM` with all relevant functions necessary for
//   segmentation of a squiggle into a series of regions.
// * struct ViterbiOutputs: contains log-prob scores and inferred state sequence
//   from a run of the viterbi algorithm.
// * struct Segmentation: contains ending sample indices for each region of a
//   squiggle's segmentation.
// * SegmentationHMM: class defining a hidden markov model for segmentation of
//   a squiggle. Contains the following members:
//   - state_transitions
//   - start_probs
//   - gaussian parameters defining emission distributions
//   - log-probabilities
//   + viterbi
//   + segment_squiggle
//   + log_probas
// ================================================================================

// The parameters of a gaussian distribution
struct GaussianParameters
{
    GaussianParameters() : mean(0.0f), stdv(1.0f), log_stdv(0.0) { }
    GaussianParameters(float m, float s) : mean(m), stdv(s) { log_stdv = log(stdv); }

    float mean;
    float stdv;
    float log_stdv; // == log(stdv), pre-computed for efficiency
};



static const float log_inv_sqrt_2pi = -0.918938f;
inline float log_normal_pdf(float x, const GaussianParameters& g)
{
    float a = (x - g.mean) / g.stdv;
    return log_inv_sqrt_2pi - g.log_stdv + (-0.5f * a * a);
}

// From SO: http://stackoverflow.com/questions/10847007/using-the-gaussian-probability-density-function-in-c
static const float inv_sqrt_2pi = 0.3989422804014327;
inline float normal_pdf(float x, const GaussianParameters& g)
{
    float a = (x - g.mean) / g.stdv;
    return inv_sqrt_2pi / g.stdv * exp(-0.5f * a * a);
}

// Segmentation struct holds endpoints of distinct regions from a segmented squiggle:
struct Segmentation {
    size_t start;   // final index of S; might not exist if skipped over
    size_t leader;  // final index of L, as indicated by 3'->5' viterbi
    size_t adapter; // final index of A, as indicated by 3'->5' viterbi
    size_t polya;   // final index of P/C, as indicated by 3'->5' viterbi
    size_t cliffs;  // number of observed 'CLIFF' samples
};

// Basic HMM struct with fixed parameters and viterbi/segmentation methods.
// (N.B.: all of the below is relative to a **scaled & shifted** set of events.)
enum HMMState
{
    HMM_START = 0,
    HMM_LEADER = 1,
    HMM_ADAPTER = 2,
    HMM_POLYA = 3,
    HMM_CLIFF = 4,
    HMM_TRANSCRIPT = 5,
    HMM_NUM_STATES = 6 // number of non-NULL states in HMM
};

// struct ViterbiOutputs composed of viterbi probs
// and a vector of integers from {0,1,2,3,4,5} == {S,L,A,P,C,T}.
struct ViterbiOutputs {
    std::vector<float> scores;
    std::vector<HMMState> labels;
};

class SegmentationHMM {
private:
    // ----- state space parameters:
    // N.B.: `state transitions` is used to compute log probabilities, as viterbi decoding is done in log-space.
    // state transition probabilities (S->L->A->[P<->C]->T):
    float state_transitions[HMM_NUM_STATES][HMM_NUM_STATES] = {
            // S -> S (10%), S -> L (90%)
            {0.10f, 0.90f, 0.00f, 0.00f, 0.00f, 0.00f},
            // L -> A (10%), L -> L (90%)
            {0.00f, 0.90f, 0.10f, 0.00f, 0.00f, 0.00f},
            // A -> P (05%), A -> A (95%)
            {0.00f, 0.00f, 0.95f, 0.05f, 0.00f, 0.00f},
            // P -> P (89%), P -> C (01%), P -> T (10%)
            {0.00f, 0.00f, 0.00f, 0.89f, 0.01f, 0.10f},
            // C -> P (99%), C -> C (01%)
            {0.00f, 0.00f, 0.00f, 0.99f, 0.01f, 0.00f},
            // T -> T (100%)
            {0.00f, 0.00f, 0.00f, 0.00f, 0.00f, 1.00f}
    };
    // All state sequences must start on S:
    float start_probs[HMM_NUM_STATES] = { 1.00f, 0.00f, 0.00f, 0.00f, 0.00f, 0.00f };

    // ----- emission parameters:
    // emission parameters, from empirical MLE on manually-flagged reads:
    // START has a mixture of Gaussian and Uniform emissions;
    // LEADER and POLYA have Gaussian emissions;
    // ADAPTER, TRANSCRIPT have Gaussian mixture emissions;
    // CLIFF has Uniform emissions.
    GaussianParameters s_emission = {70.2737f, 3.7743f};
    float s_begin = 40.0f;
    float s_end = 250.0f;
    float s_prob = 0.00476f; // == {1. / (250.0f - 40.0f)}
    float s_norm_coeff = 0.50f;
    float s_unif_coeff = 0.50f;
    GaussianParameters l_emission = {110.973f, 5.237f};
    GaussianParameters a0_emission = {79.347f, 8.3702f};
    GaussianParameters a1_emission = {63.3126f, 2.7464f};
    float a0_coeff = 0.874f;
    float a1_coeff = 0.126f;
    GaussianParameters p_emission = {108.883f, 3.257f};
    float c_begin = 70.0f;
    float c_end = 140.0f;
    float c_log_prob = -4.2485f; // natural log of [1/(140-70)]
    GaussianParameters t0_emission = {79.679f, 6.966f};
    GaussianParameters t1_emission = {105.784f, 16.022f};
    float t0_coeff = 0.346f;
    float t1_coeff = 0.654f;

    // log-probabilities are computed in the constructor:
    float log_state_transitions[HMM_NUM_STATES][HMM_NUM_STATES];
    float log_start_probs[HMM_NUM_STATES];

    // ----- inlined computation of emission log-probabilities:
    // Get the log-probability of seeing `x` given we're in state `state` of the HMM
    // N.B.: we scale the emission parameters (events are **not** scaled).
    inline float emit_log_proba(const float x, const HMMState state) const
    {
        // sometimes samples can exceed reasonable bounds due to mechanical issues;
        // in that case, we should clamp it to 100:
        float xx;
        if (x > 200.0f || x < 40.0f) {
            xx = 100.0f;
        } else {
            xx = x;
        }

        // compute on a case-by-case basis to handle heterogeneous probability distributions
        float log_probs;
        if (state == HMM_START) {
            // START state:
            float norm_term = s_norm_coeff * normal_pdf(xx, this->s_emission);
            log_probs = std::log(norm_term + s_unif_coeff * s_prob);
        }
        if (state == HMM_LEADER) {
            // LEADER state:
            log_probs = log_normal_pdf(xx, this->l_emission);
        }
        if (state == HMM_ADAPTER) {
            // ADAPTER state: compute log of gaussian mixture probability
            float mixture_proba = (this->a0_coeff*normal_pdf(xx,this->a0_emission)) + \
                (this->a1_coeff*normal_pdf(xx, this->a1_emission));
            log_probs = std::log(mixture_proba);
        }
        if (state == HMM_POLYA) {
            // POLYA state:
            log_probs = log_normal_pdf(xx, this->p_emission);
        }
        if (state == HMM_CLIFF) {
            // CLIFF state: middle-out uniform distribution
            if ((xx > this->c_begin) && (xx <  this->c_end)) {
                log_probs = this->c_log_prob;
            } else {
                log_probs = -INFINITY;
            }
        }
        if (state == HMM_TRANSCRIPT) {
            // TRANSCRIPT state: compute log of gaussian mixture probability
            float mixture_proba = (this->t0_coeff*normal_pdf(xx, this->t0_emission)) + \
                (this->t1_coeff*normal_pdf(xx, this->t1_emission));
            log_probs = std::log(mixture_proba);
        }
        return log_probs;
    }

public:
    // ----- constructor: compute logs of params & scale/shift
    SegmentationHMM(float scale, float shift, float var)
    {
        // - - - initialize log-probabilities:
        for (int i = 0; i < HMM_NUM_STATES; ++i) {
            for (int j = 0; j < HMM_NUM_STATES; ++j) {
                if (this->state_transitions[i][j] > 0.00f) {
                    this->log_state_transitions[i][j] = std::log(this->state_transitions[i][j]);
                } else {
                    this->log_state_transitions[i][j] = -INFINITY;
                }
            }
            if (this->start_probs[i] > 0.00f) {
                this->log_start_probs[i] = std::log(this->start_probs[i]);
            } else {
                this->log_start_probs[i] = -INFINITY;
            }
        }
        // - - - update all gaussian parameters by scaling/shifting:
        // START emissions:
        this->s_emission.mean = shift + scale*(this->s_emission.mean);
        this->s_emission.stdv = var * this->s_emission.stdv;
        this->s_emission.log_stdv = std::log(this->s_emission.stdv);
        // LEADER emissions:
        this->l_emission.mean = shift + scale*(this->l_emission.mean);
        this->l_emission.stdv = var * this->l_emission.stdv;
        this->l_emission.log_stdv = std::log(this->l_emission.stdv);
        // ADAPTER emissions:
        this->a0_emission.mean = shift + scale*(this->a0_emission.mean);
        this->a0_emission.stdv = var * this->a0_emission.stdv;
        this->a0_emission.log_stdv = std::log(this->a0_emission.stdv);
        this->a1_emission.mean = shift + scale*(this->a1_emission.mean);
        this->a1_emission.stdv = var * this->a1_emission.stdv;
        this->a1_emission.log_stdv = std::log(this->a1_emission.stdv);
        // POLYA emissions:
        this->p_emission.mean = shift + scale*(this->p_emission.mean);
        this->p_emission.stdv = var * this->p_emission.stdv;
        this->p_emission.log_stdv = std::log(this->p_emission.stdv);
        // TRANSCRIPT emissions:
        this->t0_emission.mean = shift + scale*(this->t0_emission.mean);
        this->t0_emission.stdv = var * this->t0_emission.stdv;
        this->t0_emission.log_stdv = std::log(this->t0_emission.stdv);
        this->t1_emission.mean = shift + scale*(this->t1_emission.mean);
        this->t1_emission.stdv = var * this->t1_emission.stdv;
        this->t1_emission.log_stdv = std::log(this->t1_emission.stdv);
    }
    // ----- destructor: nothing to clean up
    ~SegmentationHMM() { }

    // ----- for a given sample value and shift/scale parameters, return log-probs for each state:
    std::vector<float> log_probas(const float x) const
    {
        std::vector<float> log_proba(HMM_NUM_STATES);
        for (uint8_t k = 0; k < HMM_NUM_STATES; ++k) {
            log_proba[k] = this->emit_log_proba(x, static_cast<HMMState>(k));
        }
        return log_proba;
    }

    // ----- viterbi-decoding of a squiggle into region labels:
    // N.B.1: viterbi decoding happens in the 3'->5' direction.
    // N.B.2: this algorithm takes place in log-space for numerical stability;
    // the `scores` variable refers to log-prob scores.
    ViterbiOutputs viterbi(const fast5_t* f5) const
    {
        // count of raw samples:
        size_t num_samples = f5->nsample;

        // create/initialize viterbi scores and backpointers:
        std::vector<float> init_scores(HMM_NUM_STATES, -INFINITY); // log(0.0) == -INFTY
        std::vector<HMMState> init_bptrs(HMM_NUM_STATES, HMM_NUM_STATES); // HMM_NUM_STATES used as a dummy value here
        std::vector< std::vector<float> > viterbi_scores(num_samples, init_scores);
        std::vector< std::vector<HMMState> > viterbi_bptrs(num_samples, init_bptrs);

        // forward viterbi pass; fill up backpointers:
        // weight initially distributed between START and LEADER:
        viterbi_scores[0][HMM_START] = this->log_start_probs[HMM_START] + this->emit_log_proba(f5->rawptr[num_samples-1], HMM_START);
        viterbi_scores[0][HMM_LEADER] = this->log_start_probs[HMM_LEADER] + this->emit_log_proba(f5->rawptr[num_samples-1], HMM_LEADER);
        for (size_t i = 1; i < num_samples; ++i) {
            // get individual incoming state scores:
            float s_to_s = viterbi_scores.at(i-1)[HMM_START] + this->log_state_transitions[HMM_START][HMM_START];
            float s_to_l = viterbi_scores.at(i-1)[HMM_START] + this->log_state_transitions[HMM_START][HMM_LEADER];
            float l_to_l = viterbi_scores.at(i-1)[HMM_LEADER] + this->log_state_transitions[HMM_LEADER][HMM_LEADER];
            float l_to_a = viterbi_scores.at(i-1)[HMM_LEADER] + this->log_state_transitions[HMM_LEADER][HMM_ADAPTER];
            float a_to_a = viterbi_scores.at(i-1)[HMM_ADAPTER] + this->log_state_transitions[HMM_ADAPTER][HMM_ADAPTER];
            float a_to_p = viterbi_scores.at(i-1)[HMM_ADAPTER] + this->log_state_transitions[HMM_ADAPTER][HMM_POLYA];
            float p_to_p = viterbi_scores.at(i-1)[HMM_POLYA] + this->log_state_transitions[HMM_POLYA][HMM_POLYA];
            float p_to_c = viterbi_scores.at(i-1)[HMM_POLYA] + this->log_state_transitions[HMM_POLYA][HMM_CLIFF];
            float p_to_t = viterbi_scores.at(i-1)[HMM_POLYA] + this->log_state_transitions[HMM_POLYA][HMM_TRANSCRIPT];
            float c_to_c = viterbi_scores.at(i-1)[HMM_CLIFF] + this->log_state_transitions[HMM_CLIFF][HMM_CLIFF];
            float c_to_p = viterbi_scores.at(i-1)[HMM_CLIFF] + this->log_state_transitions[HMM_CLIFF][HMM_POLYA];
            float t_to_t = viterbi_scores.at(i-1)[HMM_TRANSCRIPT] + this->log_state_transitions[HMM_TRANSCRIPT][HMM_TRANSCRIPT];

            // update the viterbi scores for each state at this timestep:
            viterbi_scores.at(i)[HMM_START] = s_to_s + this->emit_log_proba(f5->rawptr[i], HMM_START);
            viterbi_scores.at(i)[HMM_LEADER] = std::max(l_to_l, s_to_l) + this->emit_log_proba(f5->rawptr[i], HMM_LEADER);
            viterbi_scores.at(i)[HMM_ADAPTER] = std::max(a_to_a, l_to_a) + this->emit_log_proba(f5->rawptr[i], HMM_ADAPTER);
            viterbi_scores.at(i)[HMM_POLYA] = std::max(p_to_p, std::max(a_to_p, c_to_p)) + this->emit_log_proba(f5->rawptr[i], HMM_POLYA);
            viterbi_scores.at(i)[HMM_CLIFF] = std::max(c_to_c, p_to_c) + this->emit_log_proba(f5->rawptr[i], HMM_CLIFF);
            viterbi_scores.at(i)[HMM_TRANSCRIPT] = std::max(p_to_t, t_to_t) + this->emit_log_proba(f5->rawptr[i], HMM_TRANSCRIPT);

            // backpointers:
            // START: S can only come from S
            viterbi_bptrs.at(i)[HMM_START] = HMM_START;
            // LEADER: L->L or S->L
            if (s_to_l < l_to_l) {
                viterbi_bptrs.at(i)[HMM_LEADER] = HMM_LEADER;
            } else {
                viterbi_bptrs.at(i)[HMM_LEADER] = HMM_START;
            }
            // ADAPTER:
            if (l_to_a < a_to_a) {
                viterbi_bptrs.at(i)[HMM_ADAPTER] = HMM_ADAPTER;
            } else {
                viterbi_bptrs.at(i)[HMM_ADAPTER] = HMM_LEADER;
            }
            // POLYA:
            if ((a_to_p < p_to_p) && (c_to_p < p_to_p)) {
                viterbi_bptrs.at(i)[HMM_POLYA] = HMM_POLYA;
            } else if ((p_to_p < a_to_p) && (c_to_p < a_to_p)) {
                viterbi_bptrs.at(i)[HMM_POLYA] = HMM_ADAPTER;
            } else {
                viterbi_bptrs.at(i)[HMM_POLYA] = HMM_CLIFF;
            }
            // CLIFF:
            if (p_to_c < c_to_c) {
                viterbi_bptrs.at(i)[HMM_CLIFF] = HMM_CLIFF;
            } else {
                viterbi_bptrs.at(i)[HMM_CLIFF] = HMM_POLYA;
            }
            // TRANSCRIPT:
            if (p_to_t < t_to_t) {
                viterbi_bptrs.at(i)[HMM_TRANSCRIPT] = HMM_TRANSCRIPT;
            } else {
                viterbi_bptrs.at(i)[HMM_TRANSCRIPT] = HMM_POLYA;
            }
        }

        // backwards viterbi pass:
        // allocate `regions` vector of same dimensions as sample sequence;
        // clamp final state to 'T' ~ transcript:
        std::vector<HMMState> regions(num_samples, HMM_START);
        std::vector<float> scores(num_samples, 0);
        regions[num_samples-1] = HMM_TRANSCRIPT;
        scores[num_samples-1] = viterbi_scores.at(num_samples-1)[HMM_TRANSCRIPT];
        // loop backwards and keep appending best states:
        for (size_t j=(num_samples-2); j > 0; --j) {
            regions[j] = viterbi_bptrs.at(j)[regions.at(j+1)];
            scores[j] = viterbi_scores.at(j)[regions.at(j+1)];
        }

        // format as ViterbiOutputs struct and return:
        ViterbiOutputs output_vectors = { scores, regions };
        return output_vectors;
    }

    // ----- parse a squiggle's viterbi labels into a regional segmentation:
    Segmentation segment_squiggle(const fast5_t* f5) const
    {
        ViterbiOutputs viterbi_outs = this->viterbi(f5);

        // compute final sample indices of each region:
        std::vector<HMMState>& region_labels = viterbi_outs.labels;

        // initial values for indices should preserve expected order:
        Segmentation ixs = { 0, 1, 2, 3, 0 };

        // loop through sequence and collect values:
        for (std::vector<uint8_t>::size_type i = 0; i < region_labels.size(); ++i) {
            // call end of START:
            if (region_labels[i] == HMM_START && region_labels[i+1] == HMM_LEADER) {
                ixs.start = static_cast<size_t>(i);
            }
            // call end of leader:
            if (region_labels[i] == HMM_LEADER && region_labels[i+1] == HMM_ADAPTER) {
                ixs.leader = static_cast<size_t>(i);
            }
            // call end of adapter:
            if (region_labels[i] == HMM_ADAPTER && region_labels[i+1] == HMM_POLYA) {
                ixs.adapter = static_cast<size_t>(i);
            }
            // call end of polya:
            if (region_labels[i] == HMM_POLYA && region_labels[i+1] == HMM_TRANSCRIPT) {
                ixs.polya = static_cast<size_t>(i);
            }
            // increment cliff counter:
            if (region_labels[i] == HMM_CLIFF) {
                ixs.cliffs++;
            }
        }

        // set sensible (easy to QC-filter) default values if not all four detected;
        // S-end is always detected (min value == 0)
        if (ixs.leader == 1 || ixs.adapter == 2 || ixs.polya == 3) {
            ixs.leader = region_labels.size() - 3;
            ixs.adapter = region_labels.size() - 2;
            ixs.polya = region_labels.size() - 1;
        }
        return ixs;
    }
};


inline float get_duration_seconds(const event_table* events, uint32_t event_idx, float sample_rate)
{
    assert(event_idx < events->n);
    return ((events->event)[event_idx].length)/sample_rate;
}

// ================================================================================
// Estimate the duration profile for a single read.
//   Estimate the underlying read rate.
// * estimate_eventalign_duration_profile : compute median read rate via collapsed-
//     duration event-alignment.
// * estimate_unaligned_duration_profile : compute median read rate via collapsed-
//     durations, without event-alignment.
// ================================================================================
// Compute a read-rate based on event-alignment, collapsed by consecutive 5mer identity
// (N.B.: deprecated; using non-eventaligned durations seems to work just as well
// while being faster to run.)
// compute a read-rate based on kmer-to-event mapping, collapsed by consecutive 5mer identity:
double estimate_unaligned_duration_profile(int32_t read_length,
                                           index_pair_t* base_to_event_map,
                                           float sample_rate,
                                           event_table* events)
{
    // get kmer stats:
    size_t num_kmers = read_length - KMER_SIZE + 1;

    // collect durations, collapsing by k-mer:
    std::vector<double> durations_per_kmer(num_kmers);
    for (size_t i = 0; i < num_kmers; ++i) {
        size_t start_idx = base_to_event_map[i].start;
        size_t end_idx = base_to_event_map[i].stop;
        // no events for this k-mer
        if (start_idx == -1) {
            continue;
        }
        assert(start_idx <= end_idx);
        for (size_t j = start_idx; j <= end_idx; ++j) {
            durations_per_kmer[i] += get_duration_seconds(events,j,sample_rate);
        }
    }

    std::sort(durations_per_kmer.begin(), durations_per_kmer.end());
    assert(durations_per_kmer.size() > 0);
    double median_duration = durations_per_kmer[durations_per_kmer.size() / 2];

    // this is our estimator of read rate, currently we use the median duration
    // per k-mer as its more robust to outliers caused by stalls
    double read_rate = 1.0 / median_duration;

    return read_rate;
}

// ================================================================================
// Poly-A Tail Length Estimation
//   Estimate the number of nucleotides in the poly-A region.
// * estimate_polya_length : return an estimate of the read rate for this read.
// ================================================================================
// Compute an estimate of the number of nucleotides in the poly-A tail
double estimate_polya_length(const float sample_rate, const Segmentation& region_indices, const double read_rate)
{
    // start and end times (sample indices) of the poly(A) tail, in original 3'->5' time-direction:
    // (n.b.: everything in 5'->3' order due to inversion in SquiggleRead constructor, but our
    // `region_indices` struct has everything in 3'->5' order)

    // calculate duration of poly(A) region (in seconds)
    double polya_duration = (region_indices.polya - (region_indices.adapter + 1)) / sample_rate;

    // Empirically determined offset to handle modal bias of the estimator:
    double estimation_error_offset = -5;

    // length of the poly(A) tail, in nucleotides:
    double polya_length = polya_duration * read_rate + estimation_error_offset;

    // ensure estimated length is non-negative:
    polya_length = std::max(0.0, polya_length);

    return polya_length;
}

void polya_estimate(bam1_t* record, int32_t read_length,
                                    index_pair_t* base_to_event_map,
                                    scalings_t scaling,
                                    fast5_t* f5,
                                    float sample_rate,
                                    event_table* events,
                                    PolyaEstimateSummary *summary){
    if((record->core.flag & BAM_FSECONDARY)){
        return;
    }
    std::string qname = bam_get_qname(record);
//    std::string ref_name(hdr->target_name[record->core.tid]);
    size_t strand_idx = 0;

    //----- get length of suffix of the read that was softclipped:
    size_t n_cigar = record->core.n_cigar;
    uint32_t prefix_cigar = bam_get_cigar(record)[0];
    uint32_t suffix_cigar = bam_get_cigar(record)[n_cigar - 1];

    uint32_t prefix_clip = bam_cigar_oplen(prefix_cigar);
    uint32_t suffix_clip = bam_cigar_oplen(suffix_cigar);


    //----- Perform pre-segmentation QC:
    // original code of pre_segmentation_qc(...)
    std::string pre_segmentation_qc_flag;
    if (suffix_clip > 200) {
        // fail if this read has a long skip at end:
        pre_segmentation_qc_flag = "SUFFCLIP";
    } else {
        // pass if none of the above fail:
        pre_segmentation_qc_flag = "PASS";
    }

    //----- perform HMM-based regional segmentation & post-segmentation QC:
    SegmentationHMM hmm(static_cast<float>(scaling.scale),
                        static_cast<float>(scaling.shift),
                        static_cast<float>(scaling.var));

    Segmentation region_indices = hmm.segment_squiggle(f5);

    //----- Perform post-segmentation QC:
    std::string post_segmentation_qc_flag;
    // fetch sizes of ADAPTER and POLYA regions:
    double num_adapter_samples = (region_indices.adapter+1) - region_indices.leader;
    double num_polya_samples = region_indices.polya - (region_indices.adapter+1);

    // check for NOREGION:
    if (num_adapter_samples < 200.0 || num_polya_samples < 200.0) {
        post_segmentation_qc_flag = "NOREGION";
    } else {
        post_segmentation_qc_flag = "PASS";
    }

    //----- compute duration profile for the read:
    double read_rate = estimate_unaligned_duration_profile(read_length,base_to_event_map,sample_rate,events);
    //    estimate_eventalign_duration_profile

    //----- estimate number of nucleotides in poly-A tail & post-estimation QC:
    double polya_length = estimate_polya_length(sample_rate, region_indices, read_rate);

    //----- Perform post_estimation QC:
    std::string post_estimation_qc_flag;
    // `adapter_qc_tol` is the (empirically-discovered) upper-tolerance for number of estimated adapter nucleotides:
    double adapter_qc_tol = 300.0f;
    // estimated adapter length, in nucleotides:
    double adapter_duration = (region_indices.adapter - (region_indices.leader - 1)) / sample_rate;
    double adapter_length = adapter_duration * read_rate;
    if (adapter_length > adapter_qc_tol) {
        post_estimation_qc_flag = "ADAPTER";
    } else {
        post_estimation_qc_flag = "PASS";
    }

    //----- Resolve QC flag based on priority:
    std::string qc_tag;
    if (post_segmentation_qc_flag.compare("PASS") != 0) {
        qc_tag = post_segmentation_qc_flag;
    } else if (post_estimation_qc_flag.compare("PASS") != 0) {
        qc_tag = post_estimation_qc_flag;
    } else if (pre_segmentation_qc_flag.compare("PASS") != 0) {
        qc_tag = pre_segmentation_qc_flag;
    } else {
        qc_tag = "PASS";
    }

    summary->leader_sample_start = region_indices.start+1;
    summary->adapter_sample_start = region_indices.leader+1;
    summary->polya_sample_start = region_indices.adapter+1;
    summary->polya_sample_end = region_indices.polya;
    summary->transcr_sample_start = region_indices.polya+1;

    summary->polya_length = polya_length;
    summary->read_rate = read_rate;
    summary->qc_tag = (char*)malloc(sizeof(char)*qc_tag.size());
    strcpy(summary->qc_tag,qc_tag.c_str());
    return;
}