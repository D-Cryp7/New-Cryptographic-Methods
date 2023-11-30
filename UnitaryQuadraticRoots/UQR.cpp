#include "../rabin.h"
#include <iostream>
#include <fstream>
#include <gmpxx.h>
#include <string>
#include <chrono>

using namespace std;

bool read_file(const string& filename, string& p, string& q, string& m) {
    ifstream file(filename);

    if (!file.is_open()) {
        cerr << "Unable to open the file " << filename << endl;
        return false;
    }

    if (!(file >> p >> q >> m)) {
        cerr << "Error reading numbers from the file." << endl;
        file.close();
        return false;
    }

    file.close();
    return true;
}

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

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <filename>" << endl;
        return 1;
    }

    string filename = argv[1];
    string pp, qq, mm;

    if (read_file(filename, pp, qq, mm)) {
        mpz_class p{pp};
        mpz_class q{qq};
        mpz_class m{mm};
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
        sign(n, r3, r4, m, &s);
        auto end = chrono::high_resolution_clock::now();
        cout << "Signature time: " << chrono::duration_cast<chrono::nanoseconds>(end - start).count() << "ns\n";
        
        start = chrono::high_resolution_clock::now();
        bool verification = verify(n, m, s);
        end = chrono::high_resolution_clock::now();
        cout << "Verification time: " << chrono::duration_cast<chrono::nanoseconds>(end - start).count() << "ns\n";

        cout << "signature " << s << "\n";
        cout << "verify " << verification << "\n";
    }

    return 0;
}




