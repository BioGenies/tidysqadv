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
        const LenSq sequence_1_size = sequence_1.original_length();
        const LenSq sequence_2_size = sequence_2.original_length();

        // Constructs a matrix of kernel^2_1 scores, but rotated by 45 degrees.
        // The reason is that k-mers for k >= 2 can be built based on 1-mers, but the shift between sequences must be preserved for all k-mer elements.
        // Such matrix ensures that each vector has all possible k-mers for given shift.
        // E.g. one of the vectors in the middle contains the following pairs: [(1,1), (2,2), (3,3), (4,4), ...],
        // while the next contains: [(1,2), (2,3), (3,4), ...].
        std::vector<std::vector<double>> kernel_2_1(sequence_1_size + sequence_2_size - 1);

        // Index keeps track of which vector is being filled
        LenSq index = 0;
        // Filling the first half of the matrix
        for (; index < sequence_1_size; ++index) {
            const LenSq sequence_1_start_index = sequence_1_size - index - 1;
            const LenSq sequence_2_start_index = 0;
            kernel_2_1[index] = std::vector<double>(std::min(index + 1, sequence_2_size));

            for (LenSq j = 0; j < kernel_2_1[index].size(); ++j) {
                kernel_2_1[index][j] = kernel_1(
                        sequence_1[{sequence_1_start_index + j, alph_size}],
                        sequence_2[{sequence_2_start_index + j, alph_size}],
                        exponential);
            }
        }

        // Filling the second half
        for (; index < kernel_2_1.size(); ++index) {
            const LenSq sequence_1_start_index = 0;
            // We start from (0, 1) here, because (0, 0) was already computed above
            const LenSq sequence_2_start_index = index - sequence_1_size + 1;
            kernel_2_1[index] = std::vector<double>(std::min(sequence_1_size + sequence_2_size - index - 1, sequence_1_size));

            for (LenSq j = 0; j < kernel_2_1[index].size(); ++j) {
                kernel_2_1[index][j] = kernel_1(
                        sequence_1[{sequence_1_start_index + j, alph_size}],
                        sequence_2[{sequence_2_start_index + j, alph_size}],
                        exponential);
            }
        }

        double sum = 0.0;
        index = 0;
        for (; index < kernel_2_1.size(); ++index) {
            for (LenSq i = 0; i < kernel_2_1[index].size(); ++i) {
                double kmer_product = 1.0;
                for (LenSq j = 0; j < max_kmer_length; ++j) {
                    if (i + j >= kernel_2_1[index].size()) {
                        break;
                    }
                    kmer_product *= kernel_2_1[index][i + j];
                    sum += kmer_product;
                }
            }
        }

        return sum;
    }

    template<typename INTERNAL>
    inline double correlation_kernel_3(const Sequence<INTERNAL> &sequence_1,
                                       const Sequence<INTERNAL> &sequence_2,
                                       const double &self_similarity_1,
                                       const double &self_similarity_2,
                                       const AlphSize &alph_size,
                                       const LenSq &max_kmer_length,
                                       const double &exponential = 0.1) {
        return kernels_2_3(sequence_1, sequence_2, alph_size, max_kmer_length, exponential) /
            sqrt(self_similarity_1 * self_similarity_2);
    }

    template<typename INTERNAL>
    inline std::vector<std::vector<double>> correlation_kernel_3(const Sq<INTERNAL> &sq,
                                       const LenSq &max_kmer_length,
                                       const double &exponential = 0.1) {
        if (sq.type() != AMI_BSC) {
            throw std::invalid_argument("sq object must be of basic amino acid alphabet");
        }

        const AlphSize alph_size = sq.alphabet().alphabet_size();
        // Initializes kernel_3 scores for similarity of sequences to themselves so that they aren't computed for each standardization
        std::vector<double> self_similarity(sq.size(), 0);
        for (LenSq i = 0; i < sq.size(); ++i) {
            self_similarity[i] = kernels_2_3(sq[i].get(), sq[i].get(), alph_size, max_kmer_length, exponential);
        }

        // Initializes returned correlation matrix with 1s, because correlation of a sequence to itself is equal to 1
        std::vector<std::vector<double>> ret(sq.size(), std::vector<double>(sq.size(), 1));
        for (LenSq i = 0; i < sq.size(); ++i) {
            for (LenSq j = i + 1; j < sq.size(); ++j) {
                const double correlation_score = correlation_kernel_3(
                        sq[i].get(), sq[j].get(), self_similarity[i], self_similarity[j], alph_size, max_kmer_length, exponential);
                ret[i][j] = correlation_score;
                ret[j][i] = correlation_score;
            }
        }
        return ret;
    }
}
