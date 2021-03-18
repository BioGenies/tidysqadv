#pragma once

#include <cmath>
#include <tidysq/tidysq-typedefs.h>
#include "tidysqadv/constants/BLOSUM62_2.h"

namespace tidysq {
    inline double kernel_1(const LetterValue &codon_1, const LetterValue &codon_2, const double &exponential = 0.1) {
        return pow(internal::read_BLOSUM62_2(codon_1, codon_2), exponential);
    }

    template<typename INTERNAL>
    inline double kernels_2_3(const Sequence<INTERNAL> &sequence_1,
                              const Sequence<INTERNAL> &sequence_2,
                              const AlphSize &alph_size,
                              const LenSq &max_kmer_length,
                              const double &exponential = 0.1) {
        // Constructs a matrix of kernel^2_1 scores, but rotated by 45 degrees.
        // The reason is that k-mers for k >= 2 can be built based on 1-mers, but the shift between sequences must be preserved for all k-mer elements.
        // Such matrix ensures that each vector has all possible k-mers for given shift.
        // E.g. one of the vectors in the middle contains the following pairs: [(1,1), (2,2), (3,3), (4,4), ...],
        // while the next contains: [(1,2), (2,3), (3,4), ...].
        std::vector<std::vector<double>> kernel_2_1(sequence_1.size() + sequence_2.size() - 1);
        for (LenSq i = sequence_1.original_length() - 1; i >= 0; --i) {

        }

        return 0.0;
    }

    template<typename INTERNAL>
    inline double kernel_2(typename tidysq::Sequence<INTERNAL>::const_iterator iterator_1_start,
                           typename tidysq::Sequence<INTERNAL>::const_iterator iterator_2_start,
                           const LenSq &kmer_length,
                           const double &exponential = 0.1) {
        return std::inner_product(iterator_1_start, iterator_1_start + kmer_length, iterator_2_start, 1.0,
                std::multiplies<>(), [&exponential](auto const &codon_1, auto const &codon_2) {
            return kernel_1(codon_1, codon_2, exponential);
        });
    }

    template<typename INTERNAL>
    inline double kernel_3(const Sequence<INTERNAL> &sequence_1,
                           const Sequence<INTERNAL> &sequence_2,
                           const AlphSize &alph_size,
                           const LenSq &max_kmer_length,
                           const double &exponential = 0.1) {
        double ret = 0.0;
        for (LenSq kmer_length = 1; kmer_length <= max_kmer_length; ++kmer_length) {
            // Execute the loop only if sequence 1 can fit a k-mer of this length
            if (kmer_length <= sequence_1.original_length()) {
                auto sequence_1_end = sequence_1.end(alph_size) - kmer_length;
                for (auto iterator_1 = sequence_1.begin(alph_size); iterator_1 <= sequence_1_end; ++iterator_1) {
                    // Execute the loop only if sequence 2 can fit a k-mer of this length
                    if (kmer_length <= sequence_2.original_length()) {
                        auto sequence_2_end = sequence_2.end(alph_size) - kmer_length;
                        for (auto iterator_2 = sequence_2.begin(alph_size); iterator_2 <= sequence_2_end; ++iterator_2) {
                            ret += kernel_2<INTERNAL>(iterator_1, iterator_2, kmer_length, exponential);
                        }
                    }
                }
            }
        }
        return ret;
    }

    template<typename INTERNAL>
    inline double correlation_kernel_3(const Sequence<INTERNAL> &sequence_1,
                                       const Sequence<INTERNAL> &sequence_2,
                                       const AlphSize &alph_size,
                                       const LenSq &max_kmer_length,
                                       const double &exponential = 0.1) {
        return kernel_3(sequence_1, sequence_2, alph_size, max_kmer_length, exponential) / sqrt(
                kernel_3(sequence_1, sequence_1, alph_size, max_kmer_length, exponential) *
                kernel_3(sequence_2, sequence_2, alph_size, max_kmer_length, exponential)
        );
    }

    template<typename INTERNAL>
    inline std::vector<std::vector<double>> correlation_kernel_3(const Sq<INTERNAL> &sq,
                                       const LenSq &max_kmer_length,
                                       const double &exponential = 0.1) {
        if (sq.type() != AMI_BSC) {
            throw std::invalid_argument("sq object must be of basic amino acid alphabet");
        }

        const AlphSize alph_size = sq.alphabet().alphabet_size();
        // Initializes returned correlation matrix with 1s, because correlation of a sequence to itself is equal to 1
        std::vector<std::vector<double>> ret(sq.size(), std::vector<double>(sq.size(), 1));
        for (LenSq i = 0; i < sq.size(); ++i) {
            for (LenSq j = i + 1; j < sq.size(); ++j) {
                const double correlation_score = correlation_kernel_3(
                        sq[i].get(), sq[j].get(), alph_size, max_kmer_length, exponential);
                ret[i][j] = correlation_score;
                ret[j][i] = correlation_score;
            }
        }
        return ret;
    }
}
