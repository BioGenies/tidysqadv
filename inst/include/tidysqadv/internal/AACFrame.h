#pragma once

#include <vector>

namespace tidysq::internal {
    template<typename INTERNAL>
    class AACFrame {
        std::vector<double> A_frequency_;
        std::vector<double> C_frequency_;
        std::vector<double> D_frequency_;
        std::vector<double> E_frequency_;
        std::vector<double> F_frequency_;
        std::vector<double> G_frequency_;
        std::vector<double> H_frequency_;
        std::vector<double> I_frequency_;
        std::vector<double> K_frequency_;
        std::vector<double> L_frequency_;
        std::vector<double> M_frequency_;
        std::vector<double> N_frequency_;
        std::vector<double> P_frequency_;
        std::vector<double> Q_frequency_;
        std::vector<double> R_frequency_;
        std::vector<double> S_frequency_;
        std::vector<double> T_frequency_;
        std::vector<double> V_frequency_;
        std::vector<double> W_frequency_;
        std::vector<double> Y_frequency_;

    public:
        explicit AACFrame(LenSq size) :
            A_frequency_(size), C_frequency_(size), D_frequency_(size),
            E_frequency_(size), F_frequency_(size), G_frequency_(size),
            H_frequency_(size), I_frequency_(size), K_frequency_(size),
            L_frequency_(size), M_frequency_(size), N_frequency_(size),
            P_frequency_(size), Q_frequency_(size), R_frequency_(size),
            S_frequency_(size), T_frequency_(size), V_frequency_(size),
            W_frequency_(size), Y_frequency_(size) {}

        inline void insert(const std::vector<double> &frequencies_for_sequence, LenSq index) {
            A_frequency_[index] = frequencies_for_sequence[0];
            C_frequency_[index] = frequencies_for_sequence[1];
            D_frequency_[index] = frequencies_for_sequence[2];
            E_frequency_[index] = frequencies_for_sequence[3];
            F_frequency_[index] = frequencies_for_sequence[4];
            G_frequency_[index] = frequencies_for_sequence[5];
            H_frequency_[index] = frequencies_for_sequence[6];
            I_frequency_[index] = frequencies_for_sequence[7];
            K_frequency_[index] = frequencies_for_sequence[8];
            L_frequency_[index] = frequencies_for_sequence[9];
            M_frequency_[index] = frequencies_for_sequence[10];
            N_frequency_[index] = frequencies_for_sequence[11];
            P_frequency_[index] = frequencies_for_sequence[12];
            Q_frequency_[index] = frequencies_for_sequence[13];
            R_frequency_[index] = frequencies_for_sequence[14];
            S_frequency_[index] = frequencies_for_sequence[15];
            T_frequency_[index] = frequencies_for_sequence[16];
            V_frequency_[index] = frequencies_for_sequence[17];
            W_frequency_[index] = frequencies_for_sequence[18];
            Y_frequency_[index] = frequencies_for_sequence[19];
        }

        [[nodiscard]] const std::vector<double> &A_frequency() const {
            return A_frequency_;
        }

        [[nodiscard]] const std::vector<double> &C_frequency() const {
            return C_frequency_;
        }

        [[nodiscard]] const std::vector<double> &D_frequency() const {
            return D_frequency_;
        }

        [[nodiscard]] const std::vector<double> &E_frequency() const {
            return E_frequency_;
        }

        [[nodiscard]] const std::vector<double> &F_frequency() const {
            return F_frequency_;
        }

        [[nodiscard]] const std::vector<double> &G_frequency() const {
            return G_frequency_;
        }

        [[nodiscard]] const std::vector<double> &H_frequency() const {
            return H_frequency_;
        }

        [[nodiscard]] const std::vector<double> &I_frequency() const {
            return I_frequency_;
        }

        [[nodiscard]] const std::vector<double> &K_frequency() const {
            return K_frequency_;
        }

        [[nodiscard]] const std::vector<double> &L_frequency() const {
            return L_frequency_;
        }

        [[nodiscard]] const std::vector<double> &M_frequency() const {
            return M_frequency_;
        }

        [[nodiscard]] const std::vector<double> &N_frequency() const {
            return N_frequency_;
        }

        [[nodiscard]] const std::vector<double> &P_frequency() const {
            return P_frequency_;
        }

        [[nodiscard]] const std::vector<double> &Q_frequency() const {
            return Q_frequency_;
        }

        [[nodiscard]] const std::vector<double> &R_frequency() const {
            return R_frequency_;
        }

        [[nodiscard]] const std::vector<double> &S_frequency() const {
            return S_frequency_;
        }

        [[nodiscard]] const std::vector<double> &T_frequency() const {
            return T_frequency_;
        }

        [[nodiscard]] const std::vector<double> &V_frequency() const {
            return V_frequency_;
        }

        [[nodiscard]] const std::vector<double> &W_frequency() const {
            return W_frequency_;
        }

        [[nodiscard]] const std::vector<double> &Y_frequency() const {
            return Y_frequency_;
        }
    };
}
