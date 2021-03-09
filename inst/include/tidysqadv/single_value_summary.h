#pragma once

#include <vector>
#include <tidysq.h>
#include <tidysqadv/constants/molecular_weights.h>

namespace tidysq {
    namespace internal {
        template<constants::SingleValueFeature FEATURE>
        double extract_value_for_codon(const LetterValue codon) {
            switch (codon) {
                case 0u:    return constants::SINGLE_VALUE_FEATURE<FEATURE, 0u>;
                case 1u:    return constants::SINGLE_VALUE_FEATURE<FEATURE, 1u>;
                case 2u:    return constants::SINGLE_VALUE_FEATURE<FEATURE, 2u>;
                case 3u:    return constants::SINGLE_VALUE_FEATURE<FEATURE, 3u>;
                case 4u:    return constants::SINGLE_VALUE_FEATURE<FEATURE, 4u>;
                case 5u:    return constants::SINGLE_VALUE_FEATURE<FEATURE, 5u>;
                case 6u:    return constants::SINGLE_VALUE_FEATURE<FEATURE, 6u>;
                case 7u:    return constants::SINGLE_VALUE_FEATURE<FEATURE, 7u>;
                case 8u:    return constants::SINGLE_VALUE_FEATURE<FEATURE, 8u>;
                case 9u:    return constants::SINGLE_VALUE_FEATURE<FEATURE, 9u>;
                case 10u:   return constants::SINGLE_VALUE_FEATURE<FEATURE, 10u>;
                case 11u:   return constants::SINGLE_VALUE_FEATURE<FEATURE, 11u>;
                case 12u:   return constants::SINGLE_VALUE_FEATURE<FEATURE, 12u>;
                case 13u:   return constants::SINGLE_VALUE_FEATURE<FEATURE, 13u>;
                case 14u:   return constants::SINGLE_VALUE_FEATURE<FEATURE, 14u>;
                case 15u:   return constants::SINGLE_VALUE_FEATURE<FEATURE, 15u>;
                case 16u:   return constants::SINGLE_VALUE_FEATURE<FEATURE, 16u>;
                case 17u:   return constants::SINGLE_VALUE_FEATURE<FEATURE, 17u>;
                case 18u:   return constants::SINGLE_VALUE_FEATURE<FEATURE, 18u>;
                case 19u:   return constants::SINGLE_VALUE_FEATURE<FEATURE, 19u>;
                default:    throw std::invalid_argument("translation must be made with twenty standard amino acid letters only");
            }
        }
    }

    // TODO: replace double with maximal possible precision, because values may vary wildly
    namespace ops {
        template<constants::SingleValueFeature FEATURE, typename INTERNAL>
        class OperationSingleValueSum : public OperationVectorToVector<Sq<INTERNAL>, Sequence<INTERNAL>, std::vector<double>, double> {
            const AlphSize alph_size_;
            const SqType& type_;

        public:
            explicit OperationSingleValueSum(const AlphSize alph_size, const SqType &type) :
                    alph_size_(alph_size),
                    type_(type) {};

            inline std::vector<double> initialize_vector_out(const Sq<INTERNAL> &sq, LenSq from, LenSq to) override {
                return std::vector<double>(to - from, 0);
            }

            inline double initialize_element_out(const Sequence<INTERNAL> &sequence) override {
                return 0;
            }

            void operator()(const Sequence<INTERNAL> &sequence, double &value) override {
                for (auto it = sequence.begin(alph_size_); it != sequence.end(alph_size_); ++it) {
                    value += internal::extract_value_for_codon<FEATURE>(*it);
                }
            }

            inline double operator() (const Sequence<INTERNAL> &sequence) override {
                //TODO: issue #57
                double value = OperationSingleValueSum<FEATURE, INTERNAL>::initialize_element_out(sequence);
                operator()(sequence, value);
                return value;
            }
        };

        template<constants::SingleValueFeature FEATURE, typename INTERNAL>
        class OperationSingleValueMean : public OperationVectorToVector<Sq<INTERNAL>, Sequence<INTERNAL>, std::vector<double>, double> {
            const AlphSize alph_size_;
            const SqType& type_;

        public:
            explicit OperationSingleValueMean(const AlphSize alph_size, const SqType &type) :
                    alph_size_(alph_size),
                    type_(type) {};

            void operator()(const Sequence<INTERNAL> &sequence, double &value) override {
                for (auto it = sequence.begin(alph_size_); it != sequence.end(alph_size_); ++it) {
                    value += internal::extract_value_for_codon<FEATURE>(*it);
                }
                value /= sequence.original_length();
            }

            inline double operator() (const Sequence<INTERNAL> &sequence) override {
                //TODO: issue #57
                double value = OperationSingleValueMean<FEATURE, INTERNAL>::initialize_element_out(sequence);
                operator()(sequence, value);
                return value;
            }
        };
    }

    template<constants::SingleValueFeature FEATURE, typename INTERNAL>
    std::vector<double> single_value_sum(const Sq<INTERNAL> &sq) {
        return sqapply(sq, ops::OperationSingleValueSum<FEATURE, INTERNAL>(sq.alphabet().alphabet_size(), sq.type()));
    }

    template<constants::SingleValueFeature FEATURE, typename INTERNAL>
    double single_value_sum(const Sequence<INTERNAL> &sequence, const AlphSize alph_size, const SqType &type) {
        return ops::OperationSingleValueSum<FEATURE, INTERNAL>(alph_size, type)(sequence);
    }

    template<constants::SingleValueFeature FEATURE, typename INTERNAL>
    std::vector<double> single_value_mean(const Sq<INTERNAL> &sq) {
        return sqapply(sq, ops::OperationSingleValueMean<FEATURE, INTERNAL>(sq.alphabet().alphabet_size(), sq.type()));
    }

    template<constants::SingleValueFeature FEATURE, typename INTERNAL>
    double single_value_mean(const Sequence<INTERNAL> &sequence, const AlphSize alph_size, const SqType &type) {
        return ops::OperationSingleValueMean<FEATURE, INTERNAL>(alph_size, type)(sequence);
    }

    template<typename INTERNAL>
    std::vector<double> molecular_weight(const Sq<INTERNAL> &sq, const std::string &type) {
        if (sq.type() != AMI_BSC) {
            throw std::invalid_argument("sq object must be of basic amino acid alphabet");
        }

        if (type == "expasy")       return single_value_sum<constants::MOLECULAR_WEIGHT_AVG_EXPASY, INTERNAL>(sq);
        if (type == "monoisotopic") return single_value_sum<constants::MOLECULAR_WEIGHT_MONOISOTOPIC, INTERNAL>(sq);
        else    throw std::invalid_argument("molecular_weight type must be one of: \"expasy\", \"monoisotopic\"");
    }
}
