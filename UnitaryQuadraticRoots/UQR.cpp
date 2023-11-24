#include "../rabin.h"
#include <iostream>
#include <gmpxx.h>
#include <chrono>

using namespace std;

int sign(mpz_class n, mpz_t y1, mpz_t y2, mpz_class m, mpz_t* c) {
    mpz_class exp{"3"};

    mpz_mul(*c, m.get_mpz_t(), y1);
    mpz_add(*c, *c, y2);
    mpz_mod(*c, *c, n.get_mpz_t());

    mpz_powm(*c, *c, exp.get_mpz_t(), n.get_mpz_t());

    return 0;
}

bool verify(mpz_class n, mpz_class m, mpz_t c) {
    mpz_t cc, ver;
    mpz_init(cc);
    mpz_init(ver);

    mpz_mul(cc, c, c); // c^2
    mpz_mod(cc, cc, n.get_mpz_t());

    mpz_class rhs = m*m*m*m*m*m - 6*m*m*m*m*m + 15*m*m*m*m -20*m*m*m + 15*m*m - 6*m + 1; // right hand side
    mpz_mod(ver, rhs.get_mpz_t(), n.get_mpz_t());

    if (mpz_cmp(cc, ver) == 0) return true;
    else return false;
}

int main() {
    mpz_class p{"109026040621161452594025497836680861802534281730394140002371179458964657129700176801446176507554210077052835471226105923870687559284502333509648461513725730062670760907457042040738837080427325016894830164228603220867733094527227125044245018472888028364512986853864174996319534121644398359766733914429864709961"};
    mpz_class q{"113197702091728559608979782791539933446334660335869257813634629981436788529276518451670986454045751295891386770190106142954710351522156266622709125914124859317862471342767449062799189204783429001620854503649142148491792534170018491689276955542092134515214775942340611363773930802274552448445026561852200522189"};
    mpz_class n = p * q;

    mpz_t r1, r2, r3, r4;
    mpz_init(r1);
    mpz_init(r2);
    mpz_init(r3);
    mpz_init(r4);

    mpz_class c{"1"};

    rabin_square_roots(c, n, p, q, &r1, &r2, &r3, &r4);

    cout << "p " << p << "\n";
    cout << "q " << q << "\n";
    cout << "n " << n << "\n";
    cout << "r1 " << r1 << "\n";
    cout << "r2 " << r2 << "\n";
    cout << "r3 " << r3 << "\n";
    cout << "r4 " << r4 << "\n";

    mpz_class m{"1946536587464453183016553848187676538276828860871122657249572479508861"}; // H3xTEL{Unitary_Square_Roots!}
    mpz_t s;
    mpz_init(s);

    auto start = chrono::high_resolution_clock::now();
    sign(n, r3, r4, m, &s);
    auto end = chrono::high_resolution_clock::now();
    cout << "Signature time: " << chrono::duration_cast<chrono::nanoseconds>(end - start).count() << "ns\n";
    
    start = chrono::high_resolution_clock::now();
    bool verification = verify(n, m, s);
    end = chrono::high_resolution_clock::now();
    cout << "Verification time: " << chrono::duration_cast<chrono::nanoseconds>(end - start).count() << "ns\n";

    cout << "signature " << s << "\n";
    cout << "verify " << verification << "\n";

    return 0;
}




