/*
  gcc -c improved_quantum_GRS.c -lgmp && gcc -o improved_quantum_GRS ../Rank_CBC_params.o basic_operation.o init_costs.o improved_quantum_GRS.o -lgmp && ./improved_quantum_GRS
*/

#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include "include.h"


int main(){
    init_Rank_CBC_params();

    int rnd_r = 0;

    for(int index = 0; index <= 5; index++){
        mpz_init_set_ui(m, Rank_SDP_instance_list[index][0]);
        mpz_init_set_ui(n, Rank_SDP_instance_list[index][1]);
        mpz_init_set_ui(k, Rank_SDP_instance_list[index][2]);
        mpz_init_set_ui(w, Rank_SDP_instance_list[index][3]);
        mpz_init_set(r, w);

        init_GRS_computational_costs();

        if(index / 3 == 0){
            mpz_t range_end, nSk, loop_r;
            mpz_init(range_end);
            mpz_init(nSk);
            mpz_sub(nSk, n, k);
            mpz_mul(range_end, m, nSk);
            mpz_cdiv_q(range_end, range_end, n);
            mpz_sub(range_end, m, range_end);
            mpz_init_set(loop_r, r);
            while(mpz_cmp(loop_r, range_end) < 0){
                init_GRS_computational_costs();

                mpz_init_set_ui(grover_iteration, 1);
                compute_grover_iteration_n_large(grover_iteration, m, loop_r, w);

                mpz_init_set_ui(grover_processor_num, 1);

                compute_IQGRSGloop_n_large(IQGRSGloop_n_large_Gcost, IQGRSGloop_n_large_Dcost, IQGRSGloop_n_large_ancila, m, n, k, w, loop_r);
                
                compute_IQGRS_n_large_Dcost(IQGRS_n_large_Dcost_r, IQGRSGloop_n_large_Dcost, grover_iteration, grover_processor_num, m, w, loop_r);
                                
                // grover_processor_num の個数を調整
                while(mpz_cmp(IQGRS_n_large_Dcost_r, MAX_DEPTH) > 0){
                    mpz_mul_ui(grover_processor_num, grover_processor_num, 2);
                    compute_IQGRS_n_large_Dcost(IQGRS_n_large_Dcost_r, IQGRSGloop_n_large_Dcost, grover_iteration, grover_processor_num, m, w, loop_r);
                }

                compute_IQGRS_n_large_Gcost(IQGRS_n_large_Gcost_r, IQGRSGloop_n_large_Gcost, grover_iteration, grover_processor_num, m, w, loop_r);

                compute_IQGRS_n_large_Wcost(IQGRS_n_large_Wcost_r, IQGRSGloop_n_large_ancila, grover_iteration, grover_processor_num, n, m, k, w, loop_r);
                
                if((mpz_cmp_ui(IQGRS_n_large_Gcost, 0) == 0) || (mpz_cmp(IQGRS_n_large_Gcost, IQGRS_n_large_Gcost_r) > 0)){
                    mpz_set(IQGRS_n_large_Gcost, IQGRS_n_large_Gcost_r);
                    mpz_set(IQGRS_n_large_Dcost, IQGRS_n_large_Dcost_r);
                    mpz_set(IQGRS_n_large_Wcost, IQGRS_n_large_Wcost_r);
                    mpz_set(r, loop_r);
                }
                
                mpz_add_ui(loop_r, loop_r, 1);
            }
            int grover_iteration_log = compute_log(grover_iteration);
            int IQGRS_n_large_Gcost_log = compute_log(IQGRS_n_large_Gcost);
            int IQGRS_n_large_Dcost_log = compute_log(IQGRS_n_large_Dcost);
            int IQGRS_n_large_Wcost_log = compute_log(IQGRS_n_large_Wcost);
            gmp_printf("n : %Zd, m : %Zd, k : %Zd, w : %Zd, r : %Zd, grover_iteration_log : %d, Gcost_log : %d, Dcost_log : %d, Wcost_log : %d\n", n, m, k, w, r, grover_iteration_log, IQGRS_n_large_Gcost_log, IQGRS_n_large_Dcost_log, IQGRS_n_large_Wcost_log);
        }
        else{
            mpz_t range_end, loop_r;
            mpz_init(range_end);
            mpz_sub(range_end, n, k);
            mpz_init_set(loop_r, r);
            while(mpz_cmp(loop_r, range_end) < 0){
                init_GRS_computational_costs();
                
                mpz_init_set_ui(grover_iteration, 1);
                compute_grover_iteration_m_large(grover_iteration, n, loop_r, w);

                mpz_init_set_ui(grover_processor_num, 1);

                compute_IQGRSGloop_m_large(IQGRSGloop_m_large_Gcost, IQGRSGloop_m_large_Dcost, IQGRSGloop_m_large_ancila, m, n, k, w, loop_r);

                compute_IQGRS_m_large_Dcost(IQGRS_m_large_Dcost_r, IQGRSGloop_m_large_Dcost, grover_iteration, grover_processor_num, n, w, loop_r);
                
                // grover_processor_num の個数を調整
                while(mpz_cmp(IQGRS_m_large_Dcost_r, MAX_DEPTH) > 0){
                    mpz_mul_ui(grover_processor_num, grover_processor_num, 2);
                    compute_IQGRS_m_large_Dcost(IQGRS_m_large_Dcost_r, IQGRSGloop_m_large_Dcost, grover_iteration, grover_processor_num, m, w, loop_r);
                }

                compute_IQGRS_m_large_Gcost(IQGRS_m_large_Gcost_r, IQGRSGloop_m_large_Gcost, grover_iteration, grover_processor_num, n, w, loop_r);

                compute_IQGRS_m_large_Wcost(IQGRS_m_large_Wcost_r, IQGRSGloop_m_large_ancila, grover_iteration, grover_processor_num, n, m, k, w, loop_r);
                
                if((mpz_cmp_ui(IQGRS_m_large_Gcost, 0) == 0) || (mpz_cmp(IQGRS_m_large_Gcost, IQGRS_m_large_Gcost_r) > 0)){
                    mpz_set(IQGRS_m_large_Gcost, IQGRS_m_large_Gcost_r);
                    mpz_set(IQGRS_m_large_Dcost, IQGRS_m_large_Dcost_r);
                    mpz_set(IQGRS_m_large_Wcost, IQGRS_m_large_Wcost_r);
                    mpz_set(r, loop_r);
                }
                
                mpz_add_ui(loop_r, loop_r, 1);
            }
            int grover_iteration_log = compute_log(grover_iteration);
            int IQGRS_m_large_Gcost_log = compute_log(IQGRS_m_large_Gcost);
            int IQGRS_m_large_Dcost_log = compute_log(IQGRS_m_large_Dcost);
            int IQGRS_m_large_Wcost_log = compute_log(IQGRS_m_large_Wcost);
            gmp_printf("n : %Zd, m : %Zd, k : %Zd, w : %Zd, r : %Zd, grover_iteration_log : %d, Gcost_log : %d, Dcost_log : %d, Wcost_log : %d\n", n, m, k, w, r, grover_iteration_log, IQGRS_m_large_Gcost_log, IQGRS_m_large_Dcost_log, IQGRS_m_large_Wcost_log);
        }
    }
}


/*
 grover_iteration_cost, m, r, w を引数にとったとき，
 n >= m での Grover のループ回数を grover_itration に格納する関数
*/
void compute_grover_iteration_n_large(mpz_t grover_iteration, mpz_t m, mpz_t r, mpz_t w){
    mpz_t grover_iteration_n_large_denominator, grover_iteration_n_large_numerator;
    mpz_init_set_ui(grover_iteration_n_large_denominator, 1);
    mpz_init_set_ui(grover_iteration_n_large_numerator, 1);

    compute_q_binomial(grover_iteration_n_large_numerator, m, w);
    compute_q_binomial(grover_iteration_n_large_denominator, r, w);

    mpz_set(grover_iteration, grover_iteration_n_large_numerator);
    mpz_cdiv_q(grover_iteration, grover_iteration, grover_iteration_n_large_denominator);

    mpz_root(grover_iteration, grover_iteration, 2);
}


/*
 grover_iteration_cost, n, r, w を引数にとったとき，
 m > n での Grover のループ回数を grover_itration に格納する関数
*/
void compute_grover_iteration_m_large(mpz_t grover_iteration, mpz_t n, mpz_t r, mpz_t w){
    mpz_t grover_iteration_m_large_denominator, grover_iteration_m_large_numerator;
    mpz_init_set_ui(grover_iteration_m_large_denominator, 1);
    mpz_init_set_ui(grover_iteration_m_large_numerator, 1);

    compute_q_binomial(grover_iteration_m_large_numerator, n, w);
    compute_q_binomial(grover_iteration_m_large_denominator, r, w);

    mpz_set(grover_iteration, grover_iteration_m_large_numerator);
    mpz_cdiv_q(grover_iteration, grover_iteration, grover_iteration_m_large_denominator);

    mpz_root(grover_iteration, grover_iteration, 2);
}


/*
 CGRSGloop_n_large_Gcost, CGRSGloop_n_large_Dcost, CGRSGloop_n_large_ancila, m, n, k, r を引数にとったとき，
 n >= m での古典GRSアルゴリズムの for 文のループ回数の G-cost, D-cost, アンシラビット数を
 それぞれ CGRSGloop_n_large_Gcost, CGRSGloop_n_large_Dcost, CGRSGloop_n_large_ancila に格納する関数
*/
void compute_CGRSGloop_n_large(mpz_t CGRSGloop_n_large_Gcost, mpz_t CGRSGloop_n_large_Dcost, mpz_t CGRSGloop_n_large_ancila, mpz_t m, mpz_t n, mpz_t k, mpz_t r){
    mpz_t rankpow_Gcost, rankpow_Dcost, rankpow_ancila;
    mpz_t mulpow_Gcost1, mulpow_Dcost1, mulpow_ancila1;
    mpz_t GE_Gcost, GE_Dcost, GE_ancila;
    mpz_t mulpow_Gcost2, mulpow_Dcost2, mulpow_ancila2;
    mpz_t addpow_Gcost, addpow_Dcost, addpow_ancila;
    mpz_init_set_ui(rankpow_Gcost, 0);
    mpz_init_set_ui(rankpow_Dcost, 0);
    mpz_init_set_ui(rankpow_ancila, 0);
    mpz_init_set_ui(mulpow_Gcost1, 0);
    mpz_init_set_ui(mulpow_Dcost1, 0);
    mpz_init_set_ui(mulpow_ancila1, 0);
    mpz_init_set_ui(GE_Gcost, 0);
    mpz_init_set_ui(GE_Dcost, 0);
    mpz_init_set_ui(GE_ancila, 0);
    mpz_init_set_ui(mulpow_Gcost2, 0);
    mpz_init_set_ui(mulpow_Dcost2, 0);
    mpz_init_set_ui(mulpow_ancila2, 0);
    mpz_init_set_ui(addpow_Gcost, 0);
    mpz_init_set_ui(addpow_Dcost, 0);
    mpz_init_set_ui(addpow_ancila, 0);

    mpz_t nSk, nMr, nSkMnMr, mMnSk, nMrS1;
    mpz_init(nSk);
    mpz_sub(nSk, n, k);
    mpz_init(nMr);
    mpz_mul(nMr, n, r);
    mpz_init(nSkMnMr);
    mpz_mul(nSkMnMr, nSk, nMr);
    mpz_init(mMnSk);
    mpz_mul(mMnSk, m, nSk);
    mpz_init(nMrS1);
    mpz_sub(nMrS1, nMr, n);

    compute_rankpow_cost(rankpow_Gcost, rankpow_Dcost, rankpow_ancila, m, n);
    compute_mulpow_cost(mulpow_Gcost1, mulpow_Dcost1, mulpow_ancila1, m);
    compute_GE_cost(GE_Gcost, GE_Dcost, GE_ancila, nMr, mMnSk);
    compute_mulpow_cost(mulpow_Gcost2, mulpow_Dcost2, mulpow_ancila2, m);
    compute_addpow_cost(addpow_Gcost, addpow_Dcost, addpow_ancila, m);

    mpz_mul(mulpow_Gcost1, mulpow_Gcost1, nSkMnMr);
    mpz_mul(mulpow_Dcost1, mulpow_Dcost1, nSkMnMr);
    mpz_mul(mulpow_ancila1, mulpow_ancila1, nSkMnMr);
    mpz_mul(mulpow_Gcost2, mulpow_Gcost2, nMr);
    mpz_mul(mulpow_Dcost2, mulpow_Dcost2, nMr);
    mpz_mul(mulpow_ancila2, mulpow_ancila2, nMr);
    mpz_mul(addpow_Gcost, addpow_Gcost, nMr);
    mpz_mul(addpow_Dcost, addpow_Dcost, nMr);
    mpz_mul(addpow_ancila, addpow_ancila, nMr);

    // CGRSGloop_n_large_Gcost を計算する
    mpz_set(CGRSGloop_n_large_Gcost, rankpow_Gcost);
    mpz_add(CGRSGloop_n_large_Gcost, CGRSGloop_n_large_Gcost, mulpow_Gcost1);
    mpz_add(CGRSGloop_n_large_Gcost, CGRSGloop_n_large_Gcost, GE_Gcost);
    mpz_add(CGRSGloop_n_large_Gcost, CGRSGloop_n_large_Gcost, mulpow_Gcost2);
    mpz_add(CGRSGloop_n_large_Gcost, CGRSGloop_n_large_Gcost, addpow_Gcost);

    // CGRSGloop_n_large_Dcost を計算する
    mpz_set(CGRSGloop_n_large_Dcost, rankpow_Dcost);
    if(mpz_cmp(mulpow_Dcost1, CGRSGloop_n_large_Dcost) > 0){
        mpz_set(CGRSGloop_n_large_Dcost, mulpow_Dcost1);
    }
    if(mpz_cmp(GE_Dcost, CGRSGloop_n_large_Dcost) > 0){
        mpz_set(CGRSGloop_n_large_Dcost, GE_Dcost);
    }
    if(mpz_cmp(mulpow_Dcost2, CGRSGloop_n_large_Dcost) > 0){
        mpz_set(CGRSGloop_n_large_Dcost, mulpow_Dcost2);
    }
    if(mpz_cmp(addpow_Dcost, CGRSGloop_n_large_Dcost) > 0){
        mpz_set(CGRSGloop_n_large_Dcost, addpow_Dcost);
    }

    // CGRSGloop_n_large_ancila を計算する
    mpz_set(CGRSGloop_n_large_ancila, rankpow_ancila);
    mpz_add(CGRSGloop_n_large_ancila, CGRSGloop_n_large_ancila, mulpow_ancila1);
    mpz_add(CGRSGloop_n_large_ancila, CGRSGloop_n_large_ancila, GE_ancila);
    mpz_add(CGRSGloop_n_large_ancila, CGRSGloop_n_large_ancila, mulpow_ancila2);
    mpz_add(CGRSGloop_n_large_ancila, CGRSGloop_n_large_ancila, addpow_ancila);
}


/*
 IQGRSGloop_n_large_Gcost, IQGRSGloop_n_large_Dcost, IQGRSGloop_n_large_ancila, m, n, k, w, r を引数にとったとき，
 n >= m での改善版量子GRSアルゴリズムの for 文のループ回数の G-cost, D-cost, アンシラビット数を
 それぞれ IQGRSGloop_n_large_Gcost, IQGRSGloop_n_large_Dcost, IQGRSGloop_n_large_ancila に格納する関数
*/
void compute_IQGRSGloop_n_large(mpz_t IQGRSGloop_n_large_Gcost, mpz_t IQGRSGloop_n_large_Dcost, mpz_t IQGRSGloop_n_large_ancila, mpz_t m, mpz_t n, mpz_t k, mpz_t w, mpz_t r){
    mpz_t QRA_Gcost1, QRA_Dcost1, QRA_ancila1;
    mpz_t QRA_Gcost2, QRA_Dcost2, QRA_ancila2;
    mpz_t CGRSGloop_n_large_Gcost, CGRSGloop_n_large_Dcost, CGRSGloop_n_large_ancila;
    mpz_t QRA_Gcost3, QRA_Dcost3, QRA_ancila3;
    mpz_t QRA_Gcost4, QRA_Dcost4, QRA_ancila4;
    mpz_t OPF_Gcost, OPF_Dcost, OPF_ancila;
    mpz_t dif_Gcost, dif_Dcost, dif_ancila;
    mpz_init_set_ui(QRA_Gcost1, 0);
    mpz_init_set_ui(QRA_Dcost1, 0);
    mpz_init_set_ui(QRA_ancila1, 0);
    mpz_init_set_ui(QRA_Gcost2, 0);
    mpz_init_set_ui(QRA_Dcost2, 0);
    mpz_init_set_ui(QRA_ancila2, 0);
    mpz_init_set_ui(CGRSGloop_n_large_Gcost, 0);
    mpz_init_set_ui(CGRSGloop_n_large_Dcost, 0);
    mpz_init_set_ui(CGRSGloop_n_large_ancila, 0);
    mpz_init_set_ui(QRA_Gcost3, 0);
    mpz_init_set_ui(QRA_Dcost3, 0);
    mpz_init_set_ui(QRA_ancila3, 0);
    mpz_init_set_ui(QRA_Gcost4, 0);
    mpz_init_set_ui(QRA_Dcost4, 0);
    mpz_init_set_ui(QRA_ancila4, 0);
    mpz_init_set_ui(OPF_Gcost, 0);
    mpz_init_set_ui(OPF_Dcost, 0);
    mpz_init_set_ui(OPF_ancila, 0);
    mpz_init_set_ui(dif_Gcost, 0);
    mpz_init_set_ui(dif_Dcost, 0);
    mpz_init_set_ui(dif_ancila, 0);

    mpz_t nMr, nSkMnMr, mMr, Vsize, mSw, rSw, coefficient, Msize;
    mpz_init(nMr);
    mpz_mul(nMr, n, r);
    mpz_init(nSkMnMr);
    mpz_sub(nSkMnMr, n, k);
    mpz_mul(nSkMnMr, nSkMnMr, nMr);
    mpz_init(mMr);
    mpz_mul(mMr, m, r);
    mpz_init_set_ui(Vsize, 1);
    compute_q_binomial(Vsize, m, r);
    mpz_init(mSw);
    mpz_sub(mSw, m, w);
    mpz_init(rSw);
    mpz_sub(rSw, r, w);
    mpz_init_set_ui(coefficient, 1);
    compute_q_binomial(coefficient, mSw, rSw);
    mpz_cdiv_q(Vsize, Vsize, coefficient);
    mpz_init_set_ui(Msize, 1);

    compute_QRA_cost(QRA_Gcost1, QRA_Dcost1, QRA_ancila1, m, Vsize);
    compute_QRA_cost(QRA_Gcost2, QRA_Dcost2, QRA_ancila2, m, r);
    compute_CGRSGloop_n_large(CGRSGloop_n_large_Gcost, CGRSGloop_n_large_Dcost, CGRSGloop_n_large_ancila, m, n, k, r);
    compute_OPF_cost(OPF_Gcost, OPF_Dcost, OPF_ancila, Vsize, Msize);
    compute_dif_cost(dif_Gcost, dif_Dcost, dif_ancila, Vsize);

    // mpz_set(QRA_Gcost3, QRA_Gcost1);
    // mpz_set(QRA_Dcost3, QRA_Dcost1);
    // mpz_set(QRA_ancila3, QRA_ancila1);
    // mpz_set(QRA_Gcost4, QRA_Gcost2);
    // mpz_set(QRA_Dcost4, QRA_Dcost2);
    // mpz_set(QRA_ancila4, QRA_ancila2);

    mpz_mul(QRA_Gcost1, QRA_Gcost1, r);
    mpz_mul(QRA_Dcost1, QRA_Dcost1, r);
    mpz_mul(QRA_ancila1, QRA_ancila1, r);
    mpz_mul(QRA_Gcost2, QRA_Gcost2, r);
    mpz_mul(QRA_Dcost2, QRA_Dcost2, r);
    mpz_mul(QRA_ancila2, QRA_ancila2, r);
    // mpz_mul(QRA_Gcost3, QRA_Gcost3, nMr);
    // mpz_mul(QRA_Dcost3, QRA_Dcost3, nMr);
    // mpz_mul(QRA_ancila3, QRA_ancila3, nMr);
    // mpz_mul(QRA_Gcost4, QRA_Gcost4, nMr);
    // mpz_mul(QRA_Dcost4, QRA_Dcost4, nMr);
    // mpz_mul(QRA_ancila4, QRA_ancila4, nMr);

    // IQGRSGloop_n_large_Gcost の計算
    mpz_set(IQGRSGloop_n_large_Gcost, QRA_Gcost1);
    mpz_add(IQGRSGloop_n_large_Gcost, IQGRSGloop_n_large_Gcost, QRA_Gcost2);
    mpz_add(IQGRSGloop_n_large_Gcost, IQGRSGloop_n_large_Gcost, CGRSGloop_n_large_Gcost);
    // mpz_add(IQGRSGloop_n_large_Gcost, IQGRSGloop_n_large_Gcost, QRA_Gcost3);
    // mpz_add(IQGRSGloop_n_large_Gcost, IQGRSGloop_n_large_Gcost, QRA_Gcost4);
    mpz_add(IQGRSGloop_n_large_Gcost, IQGRSGloop_n_large_Gcost, OPF_Gcost);
    mpz_add(IQGRSGloop_n_large_Gcost, IQGRSGloop_n_large_Gcost, dif_Gcost);

    // IQGRSGloop_n_large_Dcost の計算
    mpz_set(IQGRSGloop_n_large_Dcost, QRA_Dcost1);
    if(mpz_cmp(QRA_Dcost2, IQGRSGloop_n_large_Dcost) > 0){
        mpz_set(IQGRSGloop_n_large_Dcost, QRA_Dcost2);
    }
    if(mpz_cmp(CGRSGloop_n_large_Dcost, IQGRSGloop_n_large_Dcost) > 0){
        mpz_set(IQGRSGloop_n_large_Dcost, CGRSGloop_n_large_Dcost);
    }
    // if(mpz_cmp(QRA_Dcost3, IQGRSGloop_n_large_Dcost) > 0){
    //     mpz_set(IQGRSGloop_n_large_Dcost, QRA_Dcost3);
    // }
    // if(mpz_cmp(QRA_Dcost4, IQGRSGloop_n_large_Dcost) > 0){
    //     mpz_set(IQGRSGloop_n_large_Dcost, QRA_Dcost4);
    // }
    if(mpz_cmp(OPF_Dcost, IQGRSGloop_n_large_Dcost) > 0){
        mpz_set(IQGRSGloop_n_large_Dcost, OPF_Dcost);
    }
    if(mpz_cmp(dif_Dcost, IQGRSGloop_n_large_Dcost) > 0){
        mpz_set(IQGRSGloop_n_large_Dcost, dif_Dcost);
    }

    // IQGRSGloop_n_large_ancila の計算
    mpz_set(IQGRSGloop_n_large_ancila, QRA_ancila1);
    mpz_add(IQGRSGloop_n_large_ancila, IQGRSGloop_n_large_ancila, QRA_ancila2);
    mpz_add(IQGRSGloop_n_large_ancila, IQGRSGloop_n_large_ancila, CGRSGloop_n_large_ancila);
    // mpz_add(IQGRSGloop_n_large_ancila, IQGRSGloop_n_large_ancila, QRA_ancila3);
    // mpz_add(IQGRSGloop_n_large_ancila, IQGRSGloop_n_large_ancila, QRA_ancila4);
    mpz_add(IQGRSGloop_n_large_ancila, IQGRSGloop_n_large_ancila, OPF_ancila);
    mpz_add(IQGRSGloop_n_large_ancila, IQGRSGloop_n_large_ancila, dif_ancila);
}


/*
 IQGRS_n_large_Gcost, IQGRSGloop_n_large_Gcost, grover_iteration, grover_processor_num, m, r を引数にとったとき，
 n >= m での改善版量子GRSアルゴリズム全体の G-cost を IQGRS_n_large_Gcost に格納する関数
*/
void compute_IQGRS_n_large_Gcost(mpz_t IQGRS_n_large_Gcost, mpz_t IQGRSGloop_n_large_Gcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t m, mpz_t w, mpz_t r){
    mpz_t QRA_Gcost, QRA_Dcost, QRA_ancila;
    mpz_init_set_ui(QRA_Gcost, 0);
    mpz_init_set_ui(QRA_Dcost, 0);
    mpz_init_set_ui(QRA_ancila, 0);
    
    mpz_t mMr, coefficient, mSw, rSw, coefficient2;
    mpz_init(mMr);
    mpz_mul(mMr, m, r);
    mpz_init_set_ui(coefficient, 1);
    compute_q_binomial(coefficient, m, r);
    mpz_init(mSw);
    mpz_sub(mSw, m, w);
    mpz_init(rSw);
    mpz_sub(rSw, r, w);
    mpz_init_set_ui(coefficient2, 1);
    compute_q_binomial(coefficient2, mSw, rSw);
    mpz_cdiv_q(coefficient, coefficient, coefficient2);

    compute_QRA_cost(QRA_Gcost, QRA_Dcost, QRA_ancila, m, coefficient);

    mpz_mul(coefficient, coefficient, mMr);
    
    mpz_set(IQGRS_n_large_Gcost, IQGRSGloop_n_large_Gcost);
    mpz_mul(IQGRS_n_large_Gcost, IQGRS_n_large_Gcost, grover_iteration);
    mpz_mul(IQGRS_n_large_Gcost, IQGRS_n_large_Gcost, grover_processor_num);
    mpz_add(IQGRS_n_large_Gcost, IQGRS_n_large_Gcost, coefficient);
    mpz_add(IQGRS_n_large_Gcost, IQGRS_n_large_Gcost, QRA_Gcost);

    mpz_add(IQGRS_n_large_Gcost, IQGRS_n_large_Gcost, IQGRSGloop_n_large_Gcost);
}


/*
 IQGRS_n_large_Dcost, IQGRSGloop_n_large_Dcost, grover_iteration, grover_processor_num, m, r を引数にとったとき，
 n >= m での改善版量子GRSアルゴリズム全体の D-cost を IQGRS_n_large_Dcost に格納する関数
*/
void compute_IQGRS_n_large_Dcost(mpz_t IQGRS_n_large_Dcost, mpz_t IQGRSGloop_n_large_Dcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t m, mpz_t w, mpz_t r){
    mpz_t QRA_Gcost, QRA_Dcost, QRA_ancila;
    mpz_init_set_ui(QRA_Gcost, 0);
    mpz_init_set_ui(QRA_Dcost, 0);
    mpz_init_set_ui(QRA_ancila, 0);
    
    mpz_t mMr, coefficient, mSw, rSw, coefficient2;
    mpz_init(mMr);
    mpz_mul(mMr, m, r);
    mpz_init_set_ui(coefficient, 1);
    compute_q_binomial(coefficient, m, r);
    mpz_init(mSw);
    mpz_sub(mSw, m, w);
    mpz_init(rSw);
    mpz_sub(rSw, r, w);
    mpz_init_set_ui(coefficient2, 1);
    compute_q_binomial(coefficient2, mSw, rSw);
    mpz_cdiv_q(coefficient, coefficient, coefficient2);

    compute_QRA_cost(QRA_Gcost, QRA_Dcost, QRA_ancila, m, coefficient);
    
    mpz_set(IQGRS_n_large_Dcost, IQGRSGloop_n_large_Dcost);
    mpz_mul(IQGRS_n_large_Dcost, IQGRS_n_large_Dcost, grover_iteration);
    mpz_cdiv_q(IQGRS_n_large_Dcost, IQGRS_n_large_Dcost, grover_processor_num);
    if(mpz_cmp(QRA_Dcost, IQGRS_n_large_Dcost) > 0){
        mpz_set(IQGRS_n_large_Dcost, QRA_Dcost);
    }

    if(mpz_cmp(IQGRSGloop_n_large_Dcost, IQGRS_n_large_Dcost) > 0){
        mpz_set(IQGRS_n_large_Dcost, IQGRSGloop_n_large_Dcost);
    }
}


/*
 IQGRS_n_large_Wcost, IQGRSGloop_n_large_ancila, grover_iteration, grover_processor_num, n, m, k, r を引数にとったとき，
 n >= m での改善版量子GRSアルゴリズム全体の W-cost を IQGRS_n_large_Wcost に格納する関数
*/
void compute_IQGRS_n_large_Wcost(mpz_t IQGRS_n_large_Wcost, mpz_t IQGRSGloop_n_large_ancila, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t m, mpz_t k, mpz_t w, mpz_t r){
    mpz_t IQGRS_n_large_input;
    mpz_init_set_ui(IQGRS_n_large_input, 0);

    mpz_t QRA_Gcost, QRA_Dcost, QRA_ancila;
    mpz_init_set_ui(QRA_Gcost, 0);
    mpz_init_set_ui(QRA_Dcost, 0);
    mpz_init_set_ui(QRA_ancila, 0);
    
    mpz_t mMr, coefficient, mSw, rSw, coefficient2;
    mpz_init(mMr);
    mpz_mul(mMr, m, r);
    mpz_init_set_ui(coefficient, 1);
    compute_q_binomial(coefficient, m, r);
    mpz_init(mSw);
    mpz_sub(mSw, m, w);
    mpz_init(rSw);
    mpz_sub(rSw, r, w);
    mpz_init_set_ui(coefficient2, 1);
    compute_q_binomial(coefficient2, mSw, rSw);
    mpz_cdiv_q(coefficient, coefficient, coefficient2);

    compute_QRA_cost(QRA_Gcost, QRA_Dcost, QRA_ancila, m, coefficient);

    mpz_t nSkMm, nSkMmMn, nMm;
    mpz_init(nSkMm);
    mpz_sub(nSkMm, n, k);
    mpz_mul(nSkMm, nSkMm, m);
    mpz_init(nSkMmMn);
    mpz_mul(nSkMmMn, nSkMm, n);
    mpz_init(nMm);
    mpz_mul(nMm, n, m);
    mpz_mul(coefficient, coefficient, mMr);

    // IQGRS_n_large_input を計算する
    mpz_add(IQGRS_n_large_input, IQGRS_n_large_input, nSkMmMn);
    mpz_add(IQGRS_n_large_input, IQGRS_n_large_input, nSkMm);
    mpz_add(IQGRS_n_large_input, IQGRS_n_large_input, coefficient);
    mpz_add(IQGRS_n_large_input, IQGRS_n_large_input, nMm);
    
    mpz_set(IQGRS_n_large_Wcost, IQGRSGloop_n_large_ancila);
    mpz_mul(IQGRS_n_large_Wcost, IQGRS_n_large_Wcost, grover_iteration);
    mpz_cdiv_q(IQGRS_n_large_Wcost, IQGRS_n_large_Wcost, grover_processor_num);
    mpz_cdiv_q(IQGRS_n_large_Wcost, IQGRS_n_large_Wcost, grover_processor_num);
    mpz_add(IQGRS_n_large_Wcost, IQGRS_n_large_Wcost, QRA_ancila);
    mpz_add(IQGRS_n_large_Wcost, IQGRS_n_large_Wcost, IQGRS_n_large_input);

    mpz_add(IQGRS_n_large_Wcost, IQGRS_n_large_Wcost, IQGRSGloop_n_large_ancila);
}


/*
 CGRSGloop_m_large_Gcost, CGRSGloop_m_large_Dcost, CGRSGloop_m_large_ancila, m, n, k, r を引数にとったとき，
 m > n での古典GRSアルゴリズムの for 文のループ回数の G-cost, D-cost, アンシラビット数を
 それぞれ CGRSGloop_m_large_Gcost, CGRSGloop_m_large_Dcost, CGRSGloop_m_large_ancila に格納する関数
*/
void compute_CGRSGloop_m_large(mpz_t CGRSGloop_m_large_Gcost, mpz_t CGRSGloop_m_large_Dcost, mpz_t CGRSGloop_m_large_ancila, mpz_t m, mpz_t n, mpz_t k, mpz_t r){
    mpz_t rankpow_Gcost, rankpow_Dcost, rankpow_ancila;
    mpz_t mulpow_Gcost1, mulpow_Dcost1, mulpow_ancila1;
    mpz_t addpow_Gcost1, addpow_Dcost1, addpow_ancila1;
    mpz_t GE_Gcost, GE_Dcost, GE_ancila;
    mpz_t mulpow_Gcost2, mulpow_Dcost2, mulpow_ancila2;
    mpz_t addpow_Gcost2, addpow_Dcost2, addpow_ancila2;
    mpz_t mulpow_Gcost3, mulpow_Dcost3, mulpow_ancila3;
    mpz_t addpow_Gcost3, addpow_Dcost3, addpow_ancila3;
    mpz_init_set_ui(rankpow_Gcost, 0);
    mpz_init_set_ui(rankpow_Dcost, 0);
    mpz_init_set_ui(rankpow_ancila, 0);
    mpz_init_set_ui(mulpow_Gcost1, 0);
    mpz_init_set_ui(mulpow_Dcost1, 0);
    mpz_init_set_ui(mulpow_ancila1, 0);
    mpz_init_set_ui(addpow_Gcost1, 0);
    mpz_init_set_ui(addpow_Dcost1, 0);
    mpz_init_set_ui(addpow_ancila1, 0);
    mpz_init_set_ui(GE_Gcost, 0);
    mpz_init_set_ui(GE_Dcost, 0);
    mpz_init_set_ui(GE_ancila, 0);
    mpz_init_set_ui(mulpow_Gcost2, 0);
    mpz_init_set_ui(mulpow_Dcost2, 0);
    mpz_init_set_ui(mulpow_ancila2, 0);
    mpz_init_set_ui(addpow_Gcost2, 0);
    mpz_init_set_ui(addpow_Dcost2, 0);
    mpz_init_set_ui(addpow_ancila2, 0);
    mpz_init_set_ui(mulpow_Gcost3, 0);
    mpz_init_set_ui(mulpow_Dcost3, 0);
    mpz_init_set_ui(mulpow_ancila3, 0);
    mpz_init_set_ui(addpow_Gcost3, 0);
    mpz_init_set_ui(addpow_Dcost3, 0);
    mpz_init_set_ui(addpow_ancila3, 0);

    mpz_t nSk, nMr, mMnSk, mMmMnSkMnMr, mMr, mMrS1, mMn;
    mpz_init(nSk);
    mpz_sub(nSk, n, k);
    mpz_init(nMr);
    mpz_mul(nMr, n, r);
    mpz_init(mMnSk);
    mpz_mul(mMnSk, m, nSk);
    mpz_init(mMmMnSkMnMr);
    mpz_mul(mMmMnSkMnMr, mMnSk, nMr);
    mpz_mul(mMmMnSkMnMr, mMmMnSkMnMr, m);
    mpz_init(mMr);
    mpz_mul(mMr, m, r);
    mpz_init(mMrS1);
    mpz_sub_ui(mMrS1, mMr, 1);
    mpz_init(mMn);
    mpz_mul(mMn, m, n);

    compute_rankpow_cost(rankpow_Gcost, rankpow_Dcost, rankpow_ancila, m, n);
    compute_mulpow_cost(mulpow_Gcost1, mulpow_Dcost1, mulpow_ancila1, m);
    compute_addpow_cost(addpow_Gcost1, addpow_Dcost1, addpow_ancila1, m);
    compute_GE_cost(GE_Gcost, GE_Dcost, GE_ancila, mMr, mMnSk);
    compute_mulpow_cost(mulpow_Gcost2, mulpow_Dcost2, mulpow_ancila2, m);
    compute_addpow_cost(addpow_Gcost2, addpow_Dcost2, addpow_ancila2, m);
    compute_mulpow_cost(mulpow_Gcost3, mulpow_Dcost3, mulpow_ancila3, m);
    compute_addpow_cost(addpow_Gcost3, addpow_Dcost3, addpow_ancila3, m);

    mpz_mul(mulpow_Gcost1, mulpow_Gcost1, mMmMnSkMnMr);
    mpz_mul(mulpow_Dcost1, mulpow_Dcost1, mMmMnSkMnMr);
    mpz_mul(mulpow_ancila1, mulpow_ancila1, mMmMnSkMnMr);
    mpz_mul(addpow_Gcost1, addpow_Gcost1, mMmMnSkMnMr);
    mpz_mul(addpow_Dcost1, addpow_Dcost1, mMmMnSkMnMr);
    mpz_mul(addpow_ancila1, addpow_ancila1, mMmMnSkMnMr);
    mpz_mul(mulpow_Gcost2, mulpow_Gcost2, mMr);
    mpz_mul(mulpow_Dcost2, mulpow_Dcost2, mMr);
    mpz_mul(mulpow_ancila2, mulpow_ancila2, mMr);
    mpz_mul(addpow_Gcost2, addpow_Gcost2, mMr);
    mpz_mul(addpow_Dcost2, addpow_Dcost2, mMr);
    mpz_mul(addpow_ancila2, addpow_ancila2, mMr);
    mpz_mul(mulpow_Gcost2, mulpow_Gcost2, mMn);
    mpz_mul(mulpow_Dcost2, mulpow_Dcost2, mMn);
    mpz_mul(mulpow_ancila2, mulpow_ancila2, mMn);
    mpz_mul(addpow_Gcost2, addpow_Gcost2, mMn);
    mpz_mul(addpow_Dcost2, addpow_Dcost2, mMn);
    mpz_mul(addpow_ancila2, addpow_ancila2, mMn);

    // CGRSGloop_m_large_Gcost を計算する
    mpz_set(CGRSGloop_m_large_Gcost, rankpow_Gcost);
    mpz_add(CGRSGloop_m_large_Gcost, CGRSGloop_m_large_Gcost, mulpow_Gcost1);
    mpz_add(CGRSGloop_m_large_Gcost, CGRSGloop_m_large_Gcost, addpow_Gcost1);
    mpz_add(CGRSGloop_m_large_Gcost, CGRSGloop_m_large_Gcost, GE_Gcost);
    mpz_add(CGRSGloop_m_large_Gcost, CGRSGloop_m_large_Gcost, mulpow_Gcost2);
    mpz_add(CGRSGloop_m_large_Gcost, CGRSGloop_m_large_Gcost, addpow_Gcost2);
    mpz_add(CGRSGloop_m_large_Gcost, CGRSGloop_m_large_Gcost, mulpow_Gcost3);
    mpz_add(CGRSGloop_m_large_Gcost, CGRSGloop_m_large_Gcost, addpow_Gcost3);

    // CGRSGloop_m_large_Dcost を計算する
    mpz_set(CGRSGloop_m_large_Dcost, rankpow_Dcost);
    if(mpz_cmp(mulpow_Dcost1, CGRSGloop_m_large_Dcost) > 0){
        mpz_set(CGRSGloop_m_large_Dcost, mulpow_Dcost1);
    }
    if(mpz_cmp(addpow_Dcost1, CGRSGloop_m_large_Dcost) > 0){
        mpz_set(CGRSGloop_m_large_Dcost, addpow_Dcost1);
    }
    if(mpz_cmp(GE_Dcost, CGRSGloop_m_large_Dcost) > 0){
        mpz_set(CGRSGloop_m_large_Dcost, GE_Dcost);
    }
    if(mpz_cmp(mulpow_Dcost2, CGRSGloop_m_large_Dcost) > 0){
        mpz_set(CGRSGloop_m_large_Dcost, mulpow_Dcost2);
    }
    if(mpz_cmp(addpow_Dcost2, CGRSGloop_m_large_Dcost) > 0){
        mpz_set(CGRSGloop_m_large_Dcost, addpow_Dcost2);
    }
    if(mpz_cmp(mulpow_Dcost3, CGRSGloop_m_large_Dcost) > 0){
        mpz_set(CGRSGloop_m_large_Dcost, mulpow_Dcost3);
    }
    if(mpz_cmp(addpow_Dcost3, CGRSGloop_m_large_Dcost) > 0){
        mpz_set(CGRSGloop_m_large_Dcost, addpow_Dcost3);
    }

    // CGRSGloop_m_large_ancila を計算する
    mpz_set(CGRSGloop_m_large_ancila, rankpow_ancila);
    mpz_add(CGRSGloop_m_large_ancila, CGRSGloop_m_large_ancila, mulpow_ancila1);
    mpz_add(CGRSGloop_m_large_ancila, CGRSGloop_m_large_ancila, addpow_ancila1);
    mpz_add(CGRSGloop_m_large_ancila, CGRSGloop_m_large_ancila, GE_ancila);
    mpz_add(CGRSGloop_m_large_ancila, CGRSGloop_m_large_ancila, mulpow_ancila2);
    mpz_add(CGRSGloop_m_large_ancila, CGRSGloop_m_large_ancila, addpow_ancila2);
    mpz_add(CGRSGloop_m_large_ancila, CGRSGloop_m_large_ancila, mulpow_ancila3);
    mpz_add(CGRSGloop_m_large_ancila, CGRSGloop_m_large_ancila, addpow_ancila3);
}


/*
 IQGRSGloop_m_large_Gcost, IQGRSGloop_m_large_Dcost, IQGRSGloop_m_large_ancila, m, n, k, w, r を引数にとったとき，
 m > n での改善版量子GRSアルゴリズムの for 文のループ回数の G-cost, D-cost, アンシラビット数を
 それぞれ IQGRSGloop_m_large_Gcost, IQGRSGloop_m_large_Dcost, IQGRSGloop_m_large_ancila に格納する関数
*/
void compute_IQGRSGloop_m_large(mpz_t IQGRSGloop_m_large_Gcost, mpz_t IQGRSGloop_m_large_Dcost, mpz_t IQGRSGloop_m_large_ancila, mpz_t m, mpz_t n, mpz_t k, mpz_t w, mpz_t r){
    mpz_t QRA_Gcost1, QRA_Dcost1, QRA_ancila1;
    mpz_t QRA_Gcost2, QRA_Dcost2, QRA_ancila2;
    mpz_t CGRSGloop_m_large_Gcost, CGRSGloop_m_large_Dcost, CGRSGloop_m_large_ancila;
    mpz_t QRA_Gcost3, QRA_Dcost3, QRA_ancila3;
    mpz_t QRA_Gcost4, QRA_Dcost4, QRA_ancila4;
    mpz_t OPF_Gcost, OPF_Dcost, OPF_ancila;
    mpz_t dif_Gcost, dif_Dcost, dif_ancila;
    mpz_init_set_ui(QRA_Gcost1, 0);
    mpz_init_set_ui(QRA_Dcost1, 0);
    mpz_init_set_ui(QRA_ancila1, 0);
    mpz_init_set_ui(QRA_Gcost2, 0);
    mpz_init_set_ui(QRA_Dcost2, 0);
    mpz_init_set_ui(QRA_ancila2, 0);
    mpz_init_set_ui(CGRSGloop_m_large_Gcost, 0);
    mpz_init_set_ui(CGRSGloop_m_large_Dcost, 0);
    mpz_init_set_ui(CGRSGloop_m_large_ancila, 0);
    mpz_init_set_ui(QRA_Gcost3, 0);
    mpz_init_set_ui(QRA_Dcost3, 0);
    mpz_init_set_ui(QRA_ancila3, 0);
    mpz_init_set_ui(QRA_Gcost4, 0);
    mpz_init_set_ui(QRA_Dcost4, 0);
    mpz_init_set_ui(QRA_ancila4, 0);
    mpz_init_set_ui(OPF_Gcost, 0);
    mpz_init_set_ui(OPF_Dcost, 0);
    mpz_init_set_ui(OPF_ancila, 0);
    mpz_init_set_ui(dif_Gcost, 0);
    mpz_init_set_ui(dif_Dcost, 0);
    mpz_init_set_ui(dif_ancila, 0);

    mpz_t nMr, nSkMmMnMr, mMr, Vsize, nSw, rSw, coefficient, Msize;
    mpz_init(nMr);
    mpz_mul(nMr, n, r);
    mpz_init(nSkMmMnMr);
    mpz_sub(nSkMmMnMr, n, k);
    mpz_mul(nSkMmMnMr, nSkMmMnMr, m);
    mpz_mul(nSkMmMnMr, nSkMmMnMr, nMr);
    mpz_init(mMr);
    mpz_mul(mMr, m, r);
    mpz_init_set_ui(Vsize, 1);
    compute_q_binomial(Vsize, n, r);
    mpz_init(nSw);
    mpz_sub(nSw, n, w);
    mpz_init(rSw);
    mpz_sub(rSw, r, w);
    mpz_init_set_ui(coefficient, 1);
    compute_q_binomial(coefficient, nSw, rSw);
    mpz_cdiv_q(Vsize, Vsize, coefficient);
    mpz_init_set_ui(Msize, 1);

    compute_QRA_cost(QRA_Gcost1, QRA_Dcost1, QRA_ancila1, nMr, Vsize);
    compute_QRA_cost(QRA_Gcost2, QRA_Dcost2, QRA_ancila2, n, r);
    compute_CGRSGloop_m_large(CGRSGloop_m_large_Gcost, CGRSGloop_m_large_Dcost, CGRSGloop_m_large_ancila, m, n, k, r);
    compute_OPF_cost(OPF_Gcost, OPF_Dcost, OPF_ancila, Vsize, Msize);
    compute_dif_cost(dif_Gcost, dif_Dcost, dif_ancila, Vsize);

    // mpz_set(QRA_Gcost3, QRA_Gcost1);
    // mpz_set(QRA_Dcost3, QRA_Dcost1);
    // mpz_set(QRA_ancila3, QRA_ancila1);
    // mpz_set(QRA_Gcost4, QRA_Gcost2);
    // mpz_set(QRA_Dcost4, QRA_Dcost2);
    // mpz_set(QRA_ancila4, QRA_ancila2);

    mpz_mul(QRA_Gcost1, QRA_Gcost1, nSkMmMnMr);
    mpz_mul(QRA_Dcost1, QRA_Dcost1, nSkMmMnMr);
    mpz_mul(QRA_ancila1, QRA_ancila1, nSkMmMnMr);
    mpz_mul(QRA_Gcost2, QRA_Gcost2, nSkMmMnMr);
    mpz_mul(QRA_Dcost2, QRA_Dcost2, nSkMmMnMr);
    mpz_mul(QRA_ancila2, QRA_ancila2, nSkMmMnMr);
    // mpz_mul(QRA_Gcost3, QRA_Gcost3, mMr);
    // mpz_mul(QRA_Dcost3, QRA_Dcost3, mMr);
    // mpz_mul(QRA_ancila3, QRA_ancila3, mMr);
    // mpz_mul(QRA_Gcost4, QRA_Gcost4, mMr);
    // mpz_mul(QRA_Dcost4, QRA_Dcost4, mMr);
    // mpz_mul(QRA_ancila4, QRA_ancila4, mMr);

    // IQGRSGloop_m_large_Gcost の計算
    mpz_set(IQGRSGloop_m_large_Gcost, QRA_Gcost1);
    mpz_add(IQGRSGloop_m_large_Gcost, IQGRSGloop_m_large_Gcost, QRA_Gcost2);
    mpz_add(IQGRSGloop_m_large_Gcost, IQGRSGloop_m_large_Gcost, CGRSGloop_m_large_Gcost);
    // mpz_add(IQGRSGloop_m_large_Gcost, IQGRSGloop_m_large_Gcost, QRA_Gcost3);
    // mpz_add(IQGRSGloop_m_large_Gcost, IQGRSGloop_m_large_Gcost, QRA_Gcost4);
    mpz_add(IQGRSGloop_m_large_Gcost, IQGRSGloop_m_large_Gcost, OPF_Gcost);
    mpz_add(IQGRSGloop_m_large_Gcost, IQGRSGloop_m_large_Gcost, dif_Gcost);

    // IQGRSGloop_m_large_Dcost の計算
    mpz_set(IQGRSGloop_m_large_Dcost, QRA_Dcost1);
    if(mpz_cmp(QRA_Dcost2, IQGRSGloop_m_large_Dcost) > 0){
        mpz_set(IQGRSGloop_m_large_Dcost, QRA_Dcost2);
    }
    if(mpz_cmp(CGRSGloop_m_large_Dcost, IQGRSGloop_m_large_Dcost) > 0){
        mpz_set(IQGRSGloop_m_large_Dcost, CGRSGloop_m_large_Dcost);
    }
    // if(mpz_cmp(QRA_Dcost3, IQGRSGloop_m_large_Dcost) > 0){
    //     mpz_set(IQGRSGloop_m_large_Dcost, QRA_Dcost3);
    // }
    // if(mpz_cmp(QRA_Dcost4, IQGRSGloop_m_large_Dcost) > 0){
    //     mpz_set(IQGRSGloop_m_large_Dcost, QRA_Dcost4);
    // }
    if(mpz_cmp(OPF_Dcost, IQGRSGloop_m_large_Dcost) > 0){
        mpz_set(IQGRSGloop_m_large_Dcost, OPF_Dcost);
    }
    if(mpz_cmp(dif_Dcost, IQGRSGloop_m_large_Dcost) > 0){
        mpz_set(IQGRSGloop_m_large_Dcost, dif_Dcost);
    }

    // IQGRSGloop_m_large_ancila の計算
    mpz_set(IQGRSGloop_m_large_ancila, QRA_ancila1);
    mpz_add(IQGRSGloop_m_large_ancila, IQGRSGloop_m_large_ancila, QRA_ancila2);
    mpz_add(IQGRSGloop_m_large_ancila, IQGRSGloop_m_large_ancila, CGRSGloop_m_large_ancila);
    // mpz_add(IQGRSGloop_m_large_ancila, IQGRSGloop_m_large_ancila, QRA_ancila3);
    // mpz_add(IQGRSGloop_m_large_ancila, IQGRSGloop_m_large_ancila, QRA_ancila4);
    mpz_add(IQGRSGloop_m_large_ancila, IQGRSGloop_m_large_ancila, OPF_ancila);
    mpz_add(IQGRSGloop_m_large_ancila, IQGRSGloop_m_large_ancila, dif_ancila);
}


/*
 IQGRS_m_large_Gcost, IQGRSGloop_m_large_Gcost, grover_iteration, grover_processor_num, n, r を引数にとったとき，
 m > n での改善版量子GRSアルゴリズム全体の G-cost を IQGRS_m_large_Gcost に格納する関数
*/
void compute_IQGRS_m_large_Gcost(mpz_t IQGRS_m_large_Gcost, mpz_t IQGRSGloop_m_large_Gcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t w, mpz_t r){
    mpz_t QRA_Gcost, QRA_Dcost, QRA_ancila;
    mpz_init_set_ui(QRA_Gcost, 0);
    mpz_init_set_ui(QRA_Dcost, 0);
    mpz_init_set_ui(QRA_ancila, 0);
    
    mpz_t nMr, coefficient, nSw, rSw, coefficient2;
    mpz_init(nMr);
    mpz_mul(nMr, n, r);
    mpz_init_set_ui(coefficient, 1);
    compute_q_binomial(coefficient, n, r);
    mpz_init(nSw);
    mpz_sub(nSw, n, w);
    mpz_init(rSw);
    mpz_sub(rSw, r, w);
    mpz_init_set_ui(coefficient2, 1);
    compute_q_binomial(coefficient2, nSw, rSw);
    mpz_cdiv_q(coefficient, coefficient, coefficient2);

    compute_QRA_cost(QRA_Gcost, QRA_Dcost, QRA_ancila, nMr, coefficient);

    mpz_mul(coefficient, coefficient, nMr);
    int log_coefficient = compute_log(coefficient);

    mpz_t num_Hgate;
    mpz_init_set_ui(num_Hgate, log_coefficient);

    
    mpz_set(IQGRS_m_large_Gcost, IQGRSGloop_m_large_Gcost);
    mpz_mul(IQGRS_m_large_Gcost, IQGRS_m_large_Gcost, grover_iteration);
    mpz_mul(IQGRS_m_large_Gcost, IQGRS_m_large_Gcost, grover_processor_num);
    mpz_add(IQGRS_m_large_Gcost, IQGRS_m_large_Gcost, num_Hgate);
    mpz_add(IQGRS_m_large_Gcost, IQGRS_m_large_Gcost, QRA_Gcost);

    mpz_add(IQGRS_m_large_Gcost, IQGRS_m_large_Gcost, IQGRSGloop_m_large_Gcost);
}


/*
 IQGRS_m_large_Dcost, IQGRSGloop_m_large_Dcost, grover_iteration, grover_processor_num, n, r を引数にとったとき，
 m > n での改善版量子GRSアルゴリズム全体の D-cost を IQGRS_m_large_Dcost に格納する関数
*/
void compute_IQGRS_m_large_Dcost(mpz_t IQGRS_m_large_Dcost, mpz_t IQGRSGloop_m_large_Dcost, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t w, mpz_t r){
    mpz_t QRA_Gcost, QRA_Dcost, QRA_ancila;
    mpz_init_set_ui(QRA_Gcost, 0);
    mpz_init_set_ui(QRA_Dcost, 0);
    mpz_init_set_ui(QRA_ancila, 0);
    
    mpz_t nMr, coefficient, nSw, rSw, coefficient2;
    mpz_init(nMr);
    mpz_mul(nMr, n, r);
    mpz_init_set_ui(coefficient, 1);
    compute_q_binomial(coefficient, n, r);
    mpz_init(nSw);
    mpz_sub(nSw, n, w);
    mpz_init(rSw);
    mpz_sub(rSw, r, w);
    mpz_init_set_ui(coefficient2, 1);
    compute_q_binomial(coefficient2, nSw, rSw);
    mpz_cdiv_q(coefficient, coefficient, coefficient2);

    compute_QRA_cost(QRA_Gcost, QRA_Dcost, QRA_ancila, nMr, coefficient);
    
    mpz_set(IQGRS_m_large_Dcost, IQGRSGloop_m_large_Dcost);
    mpz_mul(IQGRS_m_large_Dcost, IQGRS_m_large_Dcost, grover_iteration);
    mpz_cdiv_q(IQGRS_m_large_Dcost, IQGRS_m_large_Dcost, grover_processor_num);
    if(mpz_cmp(QRA_Dcost, IQGRS_m_large_Dcost) > 0){
        mpz_set(IQGRS_m_large_Dcost, QRA_Dcost);
    }

    if(mpz_cmp(IQGRSGloop_m_large_Dcost, IQGRS_m_large_Dcost) > 0){
        mpz_set(IQGRS_m_large_Dcost, IQGRSGloop_m_large_Dcost);
    }
}


/*
 IQGRS_m_large_Wcost, IQGRSGloop_m_large_ancila, grover_iteration, grover_processor_num, n, m, k, r を引数にとったとき，
 m > n での改善版量子GRSアルゴリズム全体の W-cost を IQGRS_m_large_Wcost に格納する関数
*/
void compute_IQGRS_m_large_Wcost(mpz_t IQGRS_m_large_Wcost, mpz_t IQGRSGloop_m_large_ancila, mpz_t grover_iteration, mpz_t grover_processor_num, mpz_t n, mpz_t m, mpz_t k, mpz_t w, mpz_t r){
    mpz_t IQGRS_m_large_input;
    mpz_init_set_ui(IQGRS_m_large_input, 0);

    mpz_t QRA_Gcost, QRA_Dcost, QRA_ancila;
    mpz_init_set_ui(QRA_Gcost, 0);
    mpz_init_set_ui(QRA_Dcost, 0);
    mpz_init_set_ui(QRA_ancila, 0);
    
    mpz_t nMr, coefficient, nSw, rSw, coefficient2;
    mpz_init(nMr);
    mpz_mul(nMr, n, r);
    mpz_init_set_ui(coefficient, 1);
    compute_q_binomial(coefficient, n, r);
    mpz_init(nSw);
    mpz_sub(nSw, n, w);
    mpz_init(rSw);
    mpz_sub(rSw, r, w);
    mpz_init_set_ui(coefficient2, 1);
    compute_q_binomial(coefficient2, nSw, rSw);
    mpz_cdiv_q(coefficient, coefficient, coefficient2);

    compute_QRA_cost(QRA_Gcost, QRA_Dcost, QRA_ancila, nMr, coefficient);

    mpz_t nSkMm, nSkMmMn, nMm;
    mpz_init(nSkMm);
    mpz_sub(nSkMm, n, k);
    mpz_mul(nSkMm, nSkMm, m);
    mpz_init(nSkMmMn);
    mpz_mul(nSkMmMn, nSkMm, n);
    mpz_init(nMm);
    mpz_mul(nMm, n, m);
    mpz_mul(coefficient, coefficient, nMr);

    int log_coefficient = compute_log(coefficient);

    mpz_t num_Hgate;
    mpz_init_set_ui(num_Hgate, log_coefficient);

    // IQGRS_m_large_input を計算する
    mpz_add(IQGRS_m_large_input, IQGRS_m_large_input, nSkMmMn);
    mpz_add(IQGRS_m_large_input, IQGRS_m_large_input, nSkMm);
    mpz_add(IQGRS_m_large_input, IQGRS_m_large_input, num_Hgate);
    mpz_add(IQGRS_m_large_input, IQGRS_m_large_input, nMm);
    
    mpz_set(IQGRS_m_large_Wcost, IQGRSGloop_m_large_ancila);
    mpz_mul(IQGRS_m_large_Wcost, IQGRS_m_large_Wcost, grover_iteration);
    mpz_mul(IQGRS_m_large_Wcost, IQGRS_m_large_Wcost, grover_processor_num);
    mpz_mul(IQGRS_m_large_Wcost, IQGRS_m_large_Wcost, grover_processor_num);
    mpz_add(IQGRS_m_large_Wcost, IQGRS_m_large_Wcost, QRA_ancila);
    mpz_add(IQGRS_m_large_Wcost, IQGRS_m_large_Wcost, IQGRS_m_large_input);

    mpz_add(IQGRS_m_large_Wcost, IQGRS_m_large_Wcost, IQGRSGloop_m_large_ancila);
}