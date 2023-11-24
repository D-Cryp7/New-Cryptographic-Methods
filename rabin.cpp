#include "rabin.h"
#include <gmpxx.h>

int rabin_square_roots(mpz_class c, mpz_class n, mpz_class p, mpz_class q, mpz_t* r1, mpz_t* r2, mpz_t* r3, mpz_t* r4) {
    mpz_t mp, yp, mq, yq, r, s;
    mpz_init(mp);
    mpz_init(yp);
    mpz_init(mq);
    mpz_init(yq);
    mpz_init(r);
    mpz_init(s);

    mpz_class exp_p = (p + 1) / 4;
    mpz_class exp_q = (q + 1) / 4;

    mpz_powm(mp, c.get_mpz_t(), exp_p.get_mpz_t(), p.get_mpz_t());
    mpz_powm(mq, c.get_mpz_t(), exp_q.get_mpz_t(), q.get_mpz_t());

    mpz_invert(yp, p.get_mpz_t(), q.get_mpz_t());
    mpz_invert(yq, q.get_mpz_t(), p.get_mpz_t());

    mpz_mul(r, yp, p.get_mpz_t()); // yp * p
    mpz_mul(r, r, mq); // yp * p * mq

    mpz_mul(s, yq, q.get_mpz_t()); // yq * q
    mpz_mul(s, s, mp); // yq * q * mp

    mpz_add(*r1, r, s);
    mpz_mod(*r1, *r1, n.get_mpz_t());

    mpz_sub(*r2, n.get_mpz_t(), *r1);
    mpz_mod(*r2, *r2, n.get_mpz_t());

    mpz_sub(*r3, r, s);
    mpz_mod(*r3, *r3, n.get_mpz_t());

    mpz_sub(*r4, n.get_mpz_t(), *r3);
    mpz_mod(*r4, *r4, n.get_mpz_t());

    return 0;
}