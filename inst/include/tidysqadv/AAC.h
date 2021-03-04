#pragma once

#include "tidysqadv/internal/AACFrame.h"

namespace tidysq {
    template<typename INTERNAL>
    std::vector<double> count_frequencies(const Sequence<INTERNAL> &sequence, const AlphSize &alph_size) {
        std::vector<double> frequencies(20);
        for (auto it = sequence.begin(alph_size); it != sequence.end(alph_size); ++it) {
            if (*it < 20u) {
                frequencies[*it]++;
            }
        }
        for (LetterValue i = 0; i < 20; ++i) {
            frequencies[i] = frequencies[i] / sequence.original_length();
        }
        return frequencies;
    }

    template<typename INTERNAL>
    internal::AACFrame<INTERNAL> AAC(const Sq<INTERNAL> &sq) {
        if (sq.type() != AMI_BSC) {
            throw std::invalid_argument("sq object must be of basic amino acid alphabet");
        }

        internal::AACFrame<INTERNAL> ret(sq.size());
        for (LenSq i = 0; i < sq.size(); ++i) {
            ret.insert(count_frequencies(sq[i].get(), sq.alphabet().alphabet_size()), i);
        }
        return ret;
    }
}
