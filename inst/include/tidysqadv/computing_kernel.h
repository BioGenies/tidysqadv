#pragma once

#include <cmath>
#include <tidysq/tidysq-typedefs.h>
#include "tidysqadv/constants/BLOSUM62_2.h"

namespace tidysq {
    inline double kernel_1(const LetterValue &codon_1, const LetterValue &codon_2, const double &exponential = 0.1) {
        return pow(internal::read_BLOSUM62_2(codon_1, codon_2), exponential);
    }

    template<typename INTERNAL>
    inline double kernel_2(typename tidysq::Sequence<INTERNAL>::const_iterator iterator_1_start,
                           typename tidysq::Sequence<INTERNAL>::const_iterator iterator_2_start,
                           const LenSq &kmer_length,
                           const double &exponential = 0.1) {
        return std::inner_product(iterator_1_start, iterator_1_start + kmer_length, iterator_2_start, 0.0,
                std::plus<>(), [&exponential](auto const &codon_1, auto const &codon_2) {
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
            for (auto iterator_1 = sequence_1.begin(alph_size); iterator_1 + kmer_length < sequence_1.end(alph_size); ++iterator_1) {
                for (auto iterator_2 = sequence_2.begin(alph_size); iterator_2 + kmer_length < sequence_2.end(alph_size); ++iterator_2) {
                    ret += kernel_2<INTERNAL>(iterator_1, iterator_2, kmer_length, exponential);
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
}
