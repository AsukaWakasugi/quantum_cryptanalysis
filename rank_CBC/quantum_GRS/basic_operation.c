#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include "basic_operation.h"


/*
  q_analogue, m を引数にとったとき，m の q 類似を q_analogue に格納する関数
  * q = 2 のとき，m の q 類似は 2^m - 1
*/
void compute_q_analogue(mpz_t q_analogue, mpz_t m){
    mpz_t loop_index;
    mpz_init_set_ui(loop_index, 0);

    while(mpz_cmp(loop_index, m) <= 0){
        mpz_mul_ui(q_analogue, q_analogue, 2);
        mpz_add_ui(loop_index, loop_index, 1);
    }
    mpz_sub_ui(q_analogue, q_analogue, 1);
}

/*
  q_binomial, m, r を引数にとったとき，
  m と r の q 二項係数を q_binomial に格納する関数
*/
void compute_q_binomial(mpz_t q_binomial, mpz_t m, mpz_t r){
    mpz_t q_binomial_denominator, q_binomial_numerator;
    mpz_init(q_binomial_denominator);
    mpz_init(q_binomial_numerator);

    mpz_t up_counter, down_counter;
    mpz_init_set_ui(up_counter, 1);
    mpz_init_set(down_counter, m);

    mpz_set_ui(q_binomial, 1);

    while(mpz_cmp(up_counter, r) <= 0){
        mpz_set_ui(q_binomial_numerator, 1);
        compute_q_analogue(q_binomial_numerator, down_counter);
        mpz_mul(q_binomial, q_binomial, q_binomial_numerator);

        mpz_set_ui(q_binomial_denominator, 1);
        compute_q_analogue(q_binomial_denominator, up_counter);
        mpz_cdiv_q(q_binomial, q_binomial, q_binomial_denominator);

        mpz_sub_ui(down_counter, down_counter, 1);
        mpz_add_ui(up_counter, up_counter, 1);
    }
}


/*
 GE_Gcost, GE_Dcost, GE_ancila, r, n を引数にとったとき，
 GE_Gcost に 4(r - 1)r(3 n - r + 5), GE_Dcost に 16 (r - 1)，GE_ancila に r^2 を格納する関数
*/
void compute_GE_cost(mpz_t GE_Gcost, mpz_t GE_Dcost, mpz_t GE_ancila, mpz_t r, mpz_t n){
    mpz_t rS1, c36n, c20r, c9n, c4r;
    mpz_init(rS1);
    mpz_sub_ui(rS1, r, 1);
    mpz_init_set_ui(c36n, 36);
    mpz_mul(c36n, c36n, n);
    mpz_init_set_ui(c20r, 20);
    mpz_mul(c20r, c20r, r);
    mpz_init_set_ui(c9n, 9);
    mpz_mul(c9n, c9n, n);
    mpz_init_set_ui(c4r, 4);
    mpz_mul(c4r, c4r, r);

    mpz_set(GE_Gcost, c36n);
    mpz_sub(GE_Gcost, GE_Gcost, c20r);
    mpz_add_ui(GE_Gcost, GE_Gcost, 22);
    mpz_mul(GE_Gcost, GE_Gcost, r);
    mpz_add_ui(GE_Gcost, GE_Gcost, 8);
    mpz_mul(GE_Gcost, GE_Gcost, rS1);

    mpz_set(GE_Dcost, c9n);
    mpz_sub(GE_Dcost, GE_Dcost, c4r);
    mpz_add_ui(GE_Dcost, GE_Dcost, 5);
    mpz_mul(GE_Dcost, GE_Dcost, rS1);
    mpz_mul(GE_Dcost, GE_Dcost, r);
    mpz_mul_ui(GE_Dcost, GE_Dcost, 8);
    mpz_cdiv_q_ui(GE_Dcost, GE_Dcost, 3);

    mpz_set(GE_ancila, r);
    mpz_mul(GE_ancila, GE_ancila, rS1);
    mpz_mul_ui(GE_ancila, GE_ancila, 3);
    mpz_cdiv_q_ui(GE_ancila, GE_ancila, 2);
}


/*
 addpow_Gcost, addpow_Dcost, addpow_ancila, m を引数にとったとき，
 addpow_Gcost に 26 m - 23, addpow_Dcost に 16，addpow_ancila に 2 m - 1 を格納する関数
*/
void compute_addpow_cost(mpz_t addpow_Gcost, mpz_t addpow_Dcost, mpz_t addpow_ancila, mpz_t m){
    mpz_set_ui(addpow_Gcost, 2);
    mpz_mul(addpow_Gcost, addpow_Gcost, m);

    mpz_set_ui(addpow_Dcost, 2);

    mpz_set(addpow_ancila, m);
}


/*
 mulpow_Gcost, mulpow_Dcost, mulpow_ancila, m を引数にとったとき，
 mulpow_Gcost に 12 m (m + 1), mulpow_Dcost に 16 m，mulpow_ancila に m を格納する関数
*/
void compute_mulpow_cost(mpz_t mulpow_Gcost, mpz_t mulpow_Dcost, mpz_t mulpow_ancila, mpz_t m){
    mpz_t mMm;
    mpz_init_set(mMm, m);
    mpz_mul(mMm, m, m);

    mpz_set_ui(mulpow_Gcost, 40);
    mpz_mul(mulpow_Gcost, mulpow_Gcost, mMm);
    mpz_sub_ui(mulpow_Gcost, mulpow_Gcost, 16);
    
    mpz_set_ui(mulpow_Dcost, 18);
    mpz_mul(mulpow_Dcost, mulpow_Dcost, m);
    mpz_sub_ui(mulpow_Dcost, mulpow_Dcost, 2);

    mpz_set(mulpow_ancila, m);
}


/*
 rankpow_Gcost, rankpow_Dcost, rankpow_ancila, n, m を引数にとったとき，
 rankpow_Gcost に 4 (n - 1)n(3 n - m + 5), rankpow_Dcost に 16 (m + 1)，rankpow_ancila に n m を格納する関数
*/
void compute_rankpow_cost(mpz_t rankpow_Gcost, mpz_t rankpow_Dcost, mpz_t rankpow_ancila, mpz_t n, mpz_t m){
    mpz_t nS1, c36m, c20n, nMm;
    mpz_init_set(nS1, n);
    mpz_sub_ui(nS1, nS1, 1);
    mpz_init_set_ui(c36m, 36);
    mpz_mul(c36m, c36m, m);
    mpz_init_set_ui(c20n, 20);
    mpz_mul(c20n, c20n, n);
    mpz_set(nMm, n);
    mpz_mul(nMm, nMm, m);

    mpz_set(rankpow_Gcost, c36m);
    mpz_sub(rankpow_Gcost, rankpow_Gcost, c20n);
    mpz_add_ui(rankpow_Gcost, rankpow_Gcost, 22);
    mpz_mul(rankpow_ancila, rankpow_Gcost, n);
    mpz_add_ui(rankpow_Gcost, rankpow_Gcost, 8);
    mpz_mul(rankpow_Gcost, rankpow_Gcost, nS1);
    mpz_add(rankpow_Gcost, rankpow_Gcost, nMm);
    
    mpz_set_ui(rankpow_Dcost, 16);
    mpz_addmul_ui(rankpow_Dcost, m, 16);

    mpz_set(rankpow_ancila, nMm);
}


/*
 QRA_Gcost, QRA_Dcost, QRA_ancila, n, m を引数にとったとき，
 QRA_Gcost に n m + n log_2 (n), QRA_Dcost に log_2 (m) + log_2 (n)，QRA_ancila に n m + n log_2 (n) を格納する関数
*/
void compute_QRA_cost(mpz_t QRA_Gcost, mpz_t QRA_Dcost, mpz_t QRA_ancila, mpz_t n, mpz_t m){
    mpz_t cp_m, cp_n, log_m, log_n;
    mpz_init_set(cp_m, m);
    mpz_init(log_m);
    int logm = compute_log(cp_m);
    mpz_set_ui(log_m, logm);
    mpz_init_set(cp_n, n);
    mpz_init(log_n);
    int logn = compute_log(cp_n);
    mpz_set_ui(log_n, logn);

    mpz_set(QRA_Gcost, n);
    mpz_mul(QRA_Gcost, QRA_Gcost, log_n);
    mpz_addmul(QRA_Gcost, n, m);

    mpz_set(QRA_Dcost, log_m);
    mpz_add(QRA_Dcost, QRA_Dcost, log_n);

    mpz_set(QRA_ancila, QRA_Gcost);
}


/*
 Dicke_Gcost, Dicke_Dcost, Dicke_ancila, n, r を引数にとったとき，
 Dicke_Gcost に 1221 n r - 1221 r^2 - 601 n + 4 r + 304, 
 Dicke_Dcost に (27 n r - 12 n - 27 r^2 + 3) / (r - 2), 
 Dicke_ancila に 0 を格納する関数
*/
void compute_Dicke_cost(mpz_t Dicke_Gcost, mpz_t Dicke_Dcost, mpz_t Dicke_ancila, mpz_t n, mpz_t r){
    mpz_t rS2, nSr, c601n, c12n;
    mpz_init(rS2);
    mpz_sub_ui(rS2, r, 2);
    mpz_init(nSr);
    mpz_sub(nSr, n, r);
    mpz_init_set_ui(c601n, 601);
    mpz_mul(c601n, c601n, n);
    mpz_init_set_ui(c12n, 12);
    mpz_mul(c12n, c12n, n);

    mpz_set_ui(Dicke_Gcost, 1221);
    mpz_mul(Dicke_Gcost, Dicke_Gcost, r);
    mpz_mul(Dicke_Gcost, Dicke_Gcost, nSr);
    mpz_sub(Dicke_Gcost, Dicke_Gcost, c601n);
    mpz_addmul_ui(Dicke_Gcost, r, 4);
    mpz_add_ui(Dicke_Gcost, Dicke_Gcost, 304);

    mpz_set_ui(Dicke_Dcost, 27);
    mpz_mul(Dicke_Dcost, Dicke_Dcost, r);
    mpz_mul(Dicke_Dcost, Dicke_Dcost, nSr);
    mpz_sub(Dicke_Dcost, Dicke_Dcost, c12n);
    mpz_cdiv_q(Dicke_Dcost, Dicke_Dcost, rS2);

    mpz_set_ui(Dicke_ancila, 0);
}


/*
 sp_Gcost, sp_Dcost, sp_ancila, Vsize を引数にとったとき，
 sp_Gcost に Vsize, 
 sp_Dcost に 1, 
 sp_ancila に 0 を格納する関数
*/
void compute_sp_cost(mpz_t sp_Gcost, mpz_t sp_Dcost, mpz_t sp_ancila, mpz_t Vsize){
    mpz_set(sp_Gcost, Vsize);

    mpz_set_ui(sp_Dcost, 1);

    mpz_set_ui(sp_ancila, 0);
}


/*
 OPF_Gcost, OPF_Dcost, OPF_ancila, Vsize, Msize を引数にとったとき，
 OPF_Gcost に 4 Msize, 
 OPF_Dcost に 4, 
 OPF_ancila に 0 を格納する関数
*/
void compute_OPF_cost(mpz_t OPF_Gcost, mpz_t OPF_Dcost, mpz_t OPF_ancila, mpz_t Vsize, mpz_t Msize){
    mpz_set_ui(OPF_Gcost, 4);
    mpz_mul(OPF_Gcost, OPF_Gcost, Msize);

    mpz_set_ui(OPF_Dcost, 4);

    mpz_set_ui(OPF_ancila, 0);
}


/*
 dif_Gcost, dif_Dcost, dif_ancila, Vsize を引数にとったとき，
 dif_Gcost に 48 Vsize - 94, 
 dif_Dcost に 32, 
 dif_ancila に Vsize - 1 を格納する関数
*/
void compute_dif_cost(mpz_t dif_Gcost, mpz_t dif_Dcost, mpz_t dif_ancila, mpz_t Vsize){
    mpz_set(dif_Gcost, Vsize);
    mpz_mul_ui(dif_Gcost, dif_Gcost, 48);
    mpz_sub_ui(dif_Gcost, dif_Gcost, 94);

    mpz_set_ui(dif_Dcost, 32);

    mpz_set(dif_ancila, Vsize);
    mpz_sub_ui(dif_ancila, dif_ancila, 1);
}

/*
 nfac, n を引数にとったとき，n_fac に n! を格納する関数
*/
void compute_fac(mpz_t nfac, mpz_t n){
    mpz_t fac_index;
    mpz_init_set_ui(fac_index, 1);
    mpz_set_ui(nfac, 1);
    while(mpz_cmp(fac_index, n) < 0){
        mpz_mul(nfac, nfac, fac_index);
        mpz_add_ui(fac_index, fac_index, 1);
    }
}


/*
 cost を引数にとったとき，log2 cost を返す関数
*/
int compute_log(mpz_t cost){
    int cost_log = 0;
    while(mpz_cmp_ui(cost, 1) > 0){
        mpz_cdiv_q_ui(cost, cost, 2);
        cost_log++;
    }

    return cost_log;
}