#pragma once

#include "tidysq/tidysq-typedefs.h"

namespace tidysq::constants {
    enum SingleValueFeature {
        MOLECULAR_WEIGHT_AVG_EXPASY,
        MOLECULAR_WEIGHT_MONOISOTOPIC
    };

    template<SingleValueFeature FEATURE, LetterValue CODON>
    inline const double SINGLE_VALUE_FEATURE = 0;
}
