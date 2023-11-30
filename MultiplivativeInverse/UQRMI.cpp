#include "../rabin.h"
#include <iostream>
#include <fstream>
#include <gmpxx.h>
#include <string>
#include <chrono>

using namespace std;

bool read_file(const string& filename, string& p, string& q, string& m, string& x, string& xx, string& xx_inv) {
    ifstream file(filename);

    if (!file.is_open()) {
        cerr << "Unable to open the file " << filename << endl;
        return false;
    }

    if (!(file >> p >> q >> m >> x >> xx >> xx_inv)) {
        cerr << "Error reading numbers from the file." << endl;
        file.close();
        return false;
    }

    file.close();
    return true;
}

int sign(mpz_class n, mpz_class x, mpz_class m, mpz_t* c) {
    mpz_class exp{"3"};
    mpz_t x_inv;
    mpz_init(x_inv);

    mpz_mul(*c, m.get_mpz_t(), x.get_mpz_t());
    mpz_mod(*c, *c, n.get_mpz_t());

    mpz_invert(x_inv, x.get_mpz_t(), n.get_mpz_t());

    mpz_add(*c, *c, x_inv);
    mpz_mod(*c, *c, n.get_mpz_t());

    mpz_powm(*c, *c, exp.get_mpz_t(), n.get_mpz_t());

    return 0;
}

bool verify(mpz_class n, mpz_class m, mpz_t c, mpz_class xx, mpz_class xx_inv) {
    mpz_t cc, ver;
    mpz_init(cc);
    mpz_init(ver);

    mpz_mul(cc, c, c); // c^2
    mpz_mod(cc, cc, n.get_mpz_t());

    mpz_class rhs = m*m*m*m*m*m*xx*xx*xx + 6*m*m*m*m*m*xx*xx + 15*m*m*m*m*xx + 20*m*m*m + 15*m*m*xx_inv + 6*m*xx_inv*xx_inv + xx_inv*xx_inv*xx_inv; // right hand side
    mpz_mod(ver, rhs.get_mpz_t(), n.get_mpz_t());

    if (mpz_cmp(cc, ver) == 0) return true;
    else return false;
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <filename>" << endl;
        return 1;
    }

    string filename = argv[1];
    string pp, qq, mm, nonce, nonce_squared, nonce_squared_inverse;

    if (read_file(filename, pp, qq, mm, nonce, nonce_squared, nonce_squared_inverse)) {
        mpz_class p{pp};
        mpz_class q{qq};
        mpz_class m{mm};
        mpz_class x{nonce};
        mpz_class xx{nonce_squared};
        mpz_class xx_inv{nonce_squared_inverse};

        mpz_class n = p * q;

        mpz_t r1, r2, r3, r4;
        mpz_init(r1);
        mpz_init(r2);
        mpz_init(r3);
        mpz_init(r4);

        mpz_class c{"1"};

        rabin_square_roots(c, n, p, q, &r1, &r2, &r3, &r4);

        /*
        cout << "p " << p << "\n";
        cout << "q " << q << "\n";
        cout << "n " << n << "\n";
        cout << "r1 " << r1 << "\n";
        cout << "r2 " << r2 << "\n";
        cout << "r3 " << r3 << "\n";
        cout << "r4 " << r4 << "\n";
        */

        mpz_t s;
        mpz_init(s);

        auto start = chrono::high_resolution_clock::now();
        sign(n, x, m, &s);
        auto end = chrono::high_resolution_clock::now();
        cout << "Signature time: " << chrono::duration_cast<chrono::nanoseconds>(end - start).count() << "ns\n";

        start = chrono::high_resolution_clock::now();
        bool verification = verify(n, m, s, xx, xx_inv);
        end = chrono::high_resolution_clock::now();
        cout << "Verification time: " << chrono::duration_cast<chrono::nanoseconds>(end - start).count() << "ns\n";
        
        cout << "signature " << s << "\n";
        cout << "verify " << verification << "\n";
    }

    return 0;
}




