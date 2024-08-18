#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include <fstream>
using namespace std;

#define ENCODED_LENGTH 15
#define MESSAGE_LENGTH 5
#define CODE_RATE 5.0 / 15.0
#define ERROR_CORRECTION 3
#define PATTERN_NUM 1e0
#define PRIMITIVE_POLYNOMIAL 19
#define DEGREE 4
// g(x) = 11101100101
vector<int> g = {1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1};

unsigned long long SEED = 155;
unsigned long long RANV;
int RANI = 0;
int ZERO = pow(2, DEGREE) - 1;
ofstream fout("Part_3_BER.txt");

// generate U1, U2 (independent random variables uniformly distributed on (0,1))
double Ranq1()
{
    if (RANI == 0)
    {
        RANV = SEED ^ 4101842887655102017LL;
        RANV ^= RANV >> 21;
        RANV ^= RANV << 35;
        RANV ^= RANV >> 4;
        RANV *= 2685821657736338717LL;
        RANI++;
    }
    RANV ^= RANV >> 21;
    RANV ^= RANV << 35;
    RANV ^= RANV >> 4;
    return RANV * 2685821657736338717LL * 5.42101086242752217e-20;
}

// generate X
double gaussian_random()
{
    double S = 0, V1 = 0, V2 = 0;
    do
    {
        double U1 = Ranq1();
        double U2 = Ranq1();
        V1 = 2 * U1 - 1;
        V2 = 2 * U2 - 1;
        S = V1 * V1 + V2 * V2;
    } while (S >= 1);
    return V1 * sqrt(-2 * log(S) / S);
}

int modulo_negative(int a, int b)
{
    return (a % b + b) % b;
}

vector<int> decimalToBinary(int num)
{
    vector<int> binary;
    while (num > 0)
    {
        binary.push_back(num % 2);
        num /= 2;
    }
    return binary;
}

vector<int> XOR (vector<int> a, vector<int> b) {
    vector<int> c(a.size(), 0);
    for (int i = 0; i < a.size(); i++) {
        c[i] = a[i] ^ b[i];
    }
    return c;
}

vector<int> shift_right(vector<int> a, int n) {
    vector<int> b(a.size(), 0);
    for (int i = 0; i < a.size(); i++) {
        if (i + n < a.size()) {
            b[i + n] = a[i];
        }
    }
    return b;
}

vector<int> new_poly_mod(vector<int> dividend, vector<int> divisor) {
    vector<int> result(ENCODED_LENGTH, 0);
    vector<int> remainder(divisor.size() - 1, 0);
    vector<int> temp(divisor.size() - 1, 0);
    vector<vector<int>> table(dividend.size() - divisor.size() + 1, temp);
    for (int i = 0; i < table.size(); i++) {
        if (i == 0) {
            for (int j = 0; j < table[i].size(); j++) {
                table[i][j] = divisor[j];
            }
        }
        else {
            table[i] = shift_right(table[i - 1], 1);
            if (table[i - 1][divisor.size() - 2] == 1) {
                table[i] = XOR(table[i], divisor);
            }
        }
    }
    // cout << "Table: " << endl;
    // for (int i = 0; i < table.size(); i++) {
    //     for (int j = 0; j < table[i].size(); j++) {
    //         cout << table[i][j] << " ";
    //     }
    //     cout << endl;
    // }
    for (int i = 8; i < dividend.size(); i++) {
        if (dividend[i] == 1) {
            remainder = XOR(remainder, table[i - 8]);
        }
    }
    for (int i = 0; i < remainder.size(); i++) {
        result[i] = remainder[i];
    }
    return result;
}

vector<int> polynomial_mod(vector<int> dividend, vector<int> divisor) {
    vector<int> remainder(divisor.size() - 1, 0);
    vector<int> temp(divisor.size(), 0); 
    vector<int> r(ENCODED_LENGTH, 0);
    temp[divisor.size() - 1] = 0;
    for (int i = 0; i < divisor.size() - 1; i++) {
        temp[i] = divisor[i];
    }   

    for (int i = 0; i < dividend.size() - divisor.size() + 1; i++) {
        if (dividend[i + divisor.size() - 1] == 1) {
            remainder = XOR(remainder, temp);
        }
        if (i + divisor.size() - 1 == dividend.size() - 1) {
            break;
        }
        // shift temp right by 1
        for (int j = divisor.size() - 1; j > 0; j--) {
            temp[j] = temp[j - 1];
        }
        temp[0] = 0;
        if (temp[divisor.size() - 1] == 1) {
            for (int j = 0; j < divisor.size() - 1; j++) {
                temp[j] = temp[j] ^ divisor[j];
            }
        }
    }
    for (int i = 0; i < remainder.size(); i++) {
        r[i] = remainder[i];
    }
    return r;
}

vector<int> encoder(vector<int> u) {
    vector<int> c(ENCODED_LENGTH, 0);
    vector<int> remainder(ENCODED_LENGTH, 0);
    vector<int> temp(ENCODED_LENGTH, 0);
    for (int i = 0; i < temp.size(); i++) {
        if (i < ENCODED_LENGTH - MESSAGE_LENGTH) {
            temp[i] = 0;
        }
        else {
            temp[i] = u[i - (ENCODED_LENGTH - MESSAGE_LENGTH)];
        }
    }
    // r(x) = temp(x) mod g(x)
    remainder = polynomial_mod(temp, g);
    // c = temp + r
    c = XOR(temp, remainder);
    return c;
}

vector<vector<int>> construct_GF(void) {
    vector<vector<int>> galois_field;
    // initialize galois field
    for (int i = 0; i < pow(2, DEGREE); i++) {
        vector<int> temp(DEGREE, 0);
        galois_field.push_back(temp);
    }
    vector<int> primitive;
    primitive = decimalToBinary(PRIMITIVE_POLYNOMIAL);
    vector<int> d_of_x(DEGREE, 0);
    for (int i = 0; i < DEGREE; i++) d_of_x[i] = primitive[i];
    for (int i = 0; i < pow(2, DEGREE); i++)
    {
        vector<int> mod(DEGREE + 1);
        if (i == 0) mod[0] = 1;
        else
        {
            // mod = table[i - 1] >> 1
            mod[0] = 0;
            for (int j = 0; j < DEGREE; j++) mod[j + 1] = galois_field[i - 1][j];
        }
        // if mod[degree] == 1 then mod = mod ^ d_of_x
        if (mod[DEGREE] == 1)
        {
            for (int j = 0; j < DEGREE; j++) mod[j] = mod[j] ^ d_of_x[j];
        }
        for (int j = 0; j < DEGREE; j++) galois_field[i][j] = mod[j];
    }
    galois_field.pop_back();
    // push back 0
    vector<int> zero(DEGREE, 0);
    galois_field.push_back(zero);
    return galois_field;
}

int GF_multiplication(int a, int b, vector<vector<int>> galois_field) {
    // a b are the index of the galois field
    if (a == ZERO || b == ZERO) return ZERO;
    int exp = modulo_negative(a + b, (int)pow(2, DEGREE) - 1);
    return exp;
}

int GF_division(int a, int b, vector<vector<int>> galois_field) {
    // a b are the index of the galois field
    if (a == ZERO) return ZERO;
    int exp = modulo_negative(a - b, (int)pow(2, DEGREE) - 1);
    return exp;
}

int GF_addition(int a, int b, vector<vector<int>> galois_field) {
    // a b are the index of the galois field
    vector<int> result(DEGREE, 0);
    result = XOR(galois_field[a], galois_field[b]);
    // search for the index of the result in the galois field
    int exp = 0;
    for (int i = 0; i < galois_field.size(); i++) {
        if (galois_field[i] == result) {
            exp = i;
            break;
        }
    }
    return exp;
}

vector<int> syndrome_computation(vector<int> r, vector<vector<int>> galois_field) {
    vector<int> syndrome(2 * ERROR_CORRECTION + 1, 0);
    syndrome[0] = -1;
    for (int i = 1; i < 2 * ERROR_CORRECTION + 1; i++) {
        int sum = ZERO;
        for (int j = 0; j < r.size(); j++) {
            if (r[j] == 1) {
                int temp = (i * j) % ((int)pow(2, DEGREE) - 1);
                sum = GF_addition(sum, temp, galois_field);
            }
        }
        syndrome[i] = sum;
    }
    return syndrome;
}

vector<int> key_equation_slover(vector<int> syndrome, vector<vector<int>> galois_field) {
    // Berlekamp-Massey Algorithm
    vector<int> sigma(ENCODED_LENGTH, ZERO); // coefficients of sigma is in galois field
    vector<int> last_sigma(ENCODED_LENGTH, ZERO); 
    vector<int> B(ENCODED_LENGTH, ZERO); // coefficients of B is in galois field
    vector<int> last_B(ENCODED_LENGTH, ZERO);
    int L;
    int last_L;
    int delta; // delta is in galois field
    int discrepancy = ZERO; // discrepancy is in galois field
    // initialization
    sigma[0] = 0;
    B[0] = 0;
    L = 0;
    delta = 0;
    // u = 1 to 2t
    for (int u = 1; u < 2 * ERROR_CORRECTION + 1; u++) {
        // store last sigma, B, L
        for (int i = 0; i < last_sigma.size(); i++) {
            last_sigma[i] = sigma[i];
            last_B[i] = B[i];
        }
        last_L = L;
        // check last sigma, B, L
        // fout << endl << "i: " << u;
        // fout << endl << "Last Sigma: ";
        // for (int i = 0; i < last_sigma.size(); i++) {
        //     fout << last_sigma[i] << " ";
        // }
        // fout << endl;
        // fout << "Last B: ";
        // for (int i = 0; i < last_B.size(); i++) {
        //     fout << last_B[i] << " ";
        // }
        // fout << endl << "Last L: " << last_L << endl;
        // compute discrepancy
        int temp2 = ZERO;
        for (int i = 0; ; i++) {
            int temp = GF_multiplication(last_sigma[i], syndrome[u - i], galois_field);  
            temp2 = GF_addition(temp2, temp, galois_field);
            if (i == last_L) break;
        }
        discrepancy = temp2;
        // compute sigma
        int temp = GF_division(discrepancy, delta, galois_field);
        vector<int> temp_vector(ENCODED_LENGTH, ZERO); // x * last_B
        for (int j = 1; j < last_B.size(); j++) {
            temp_vector[j] = last_B[j - 1];
        }
        // temp_vector = temp_vector * temp
        for (int j = 0; j < temp_vector.size(); j++) {
            temp_vector[j] = GF_multiplication(temp_vector[j], temp, galois_field);
        }
        // sigma = sigma + temp_vector
        for (int j = 0; j < sigma.size(); j++) {
            sigma[j] = GF_addition(sigma[j], temp_vector[j], galois_field);
        }
        // if 2 * last_L <= u and discrepancy != 0
        if (2 * last_L < u && discrepancy != ZERO) {
            // B = last_sigma
            for (int j = 0; j < B.size(); j++) {
                B[j] = last_sigma[j];
            }
            // delta = discrepancy
            delta = discrepancy;
            // L = u - last_L
            L = u - last_L;
        }
        else {
            // B = x * last_B
            B[0] = ZERO;
            for (int j = 1; j < B.size(); j++) {
                B[j] = last_B[j - 1];
            }
            // L = last_L
            L = last_L;
            // delta = delta
            delta = delta;
        }
        // check
        // fout << "Discrepancy: " << discrepancy << endl;
        // fout << "Sigma: ";
        // for (int i = 0; i < sigma.size(); i++) {
        //     fout << sigma[i] << " ";
        // }
        // fout << endl;
        // fout << "B: ";
        // for (int i = 0; i < B.size(); i++) {
        //     fout << B[i] << " ";
        // }
        // fout << endl << "L: " << L << endl;
        // fout << "delta: " << delta << endl;
        // fout << endl;
    }
    return sigma;
}

vector<int> chien_search(vector<int> sigma, vector<vector<int>> galois_field) {
    vector<int> error_location;
    vector<int> error_pattern(ENCODED_LENGTH, 0);
    // if sigma(a^(-i)) = 0 then i is the error location
    for (int i = 0; i < ENCODED_LENGTH; i++) {
        int sum = ZERO;
        for (int j = 0; j < sigma.size(); j++) {
            int index = modulo_negative(-i * j, (int)pow(2, DEGREE) - 1);
            int temp = GF_multiplication(sigma[j], index, galois_field);
            sum = GF_addition(sum, temp, galois_field);
        }
        if (sum == ZERO) {
            error_location.push_back(i);
            error_pattern[i] = 1;
        }
    }
    // compute degree of sigma
    int degree = 0;
    for (int i = sigma.size() - 1; i >= 0; i--) {
        if (sigma[i] != ZERO) {
            degree = i;
            break;
        }
    }
    fout << endl << "error locator polynomial: " << degree;
    fout << endl << "error location's size: " << error_location.size();
    // compute error pattern
    if (error_location.size() == degree) {
        for (int i = 0; i < error_location.size(); i++) {
            error_pattern[error_location[i]] = 1;
        }
    }
    return error_pattern;
}

vector<int> decoder(vector<int> r) {
    vector<int> v(ENCODED_LENGTH, 0);
    vector<int> Syndrome;
    vector<int> sigma(ENCODED_LENGTH, 0);
    vector<int> e(ENCODED_LENGTH, 0);
    vector<vector<int>> galois_field;
    // costruct galois field GF(2^degree)
    galois_field = construct_GF();
    // syndrome computation
    Syndrome = syndrome_computation(r, galois_field);
    fout << endl << "Syndrome: ";
    for (int i = 1; i < Syndrome.size(); i++) {
        fout << Syndrome[i] << " ";
    }
    // key equation slover
    sigma = key_equation_slover(Syndrome, galois_field);
    // chien search
    e = chien_search(sigma, galois_field);
    // correction
    v = XOR(r, e);
    return v;
}

vector<int> random_bits(int length) {
    vector<int> bits;
    for (int i = 0; i < length; i++) {
        bits.push_back(rand() % 2);
    }
    return bits;
}

vector<int> bpsk_modulation(vector<int> c) {
    vector<int> s(ENCODED_LENGTH, 0);
    for (int i = 0; i < ENCODED_LENGTH; i++) {
        if (c[i] == 0) s[i] = 1;
        else s[i] = -1;
    }
    return s;
}

vector<int> AWGN_channel(vector<int> s, double Eb_N0_dB) {
    vector<int> r(ENCODED_LENGTH, 0);
    for (int i = 0; i < ENCODED_LENGTH; i++) {
        double sigma = sqrt(1 / (2 * CODE_RATE * pow(10, Eb_N0_dB / 10)));
        r[i] = s[i] + gaussian_random() * sigma;
        if (r[i] > 0) r[i] = 1;
        else r[i] = -1;
    }
    return r;
}

vector<int> bpsk_demodulation(vector<int> r) {
    vector<int> c(ENCODED_LENGTH, 0);
    for (int i = 0; i < ENCODED_LENGTH; i++) {
        if (r[i] > 0) c[i] = 0;
        else c[i] = 1;
    }
    return c;
}

int main() {
    srand(time(0));
    ifstream fin("Eb_N0_dB.txt");
    if (!fin.is_open() || !fout.is_open())
    {
        cout << "File not found!";
        return 0;
    }

    double Eb_N0_dB;
    while (fin >> Eb_N0_dB)
    {
        int error = 0;
        fout << "Eb/N0: " << Eb_N0_dB;
        // Generate random message
        vector<int> u;
        u = random_bits(MESSAGE_LENGTH);
        fout << endl << "Message: ";
        for (int i = 0; i < u.size(); i++) {
            fout << u[i];
        }

        // Encode message
        vector<int> c;
        c = encoder(u);
        fout << endl << "Encoded message: ";
        for (int i = 0; i < c.size(); i++) {
            fout << c[i];
        }

        // BPSK modulation
        vector<int> s;
        s = bpsk_modulation(c);
        fout << endl << "Modulated message: ";
        for (int i = 0; i < s.size(); i++) {
            fout << s[i] << " ";
        }
        
        // AWGN channel
        vector<int> noised_sig;
        noised_sig = AWGN_channel(s, Eb_N0_dB);
        fout << endl << "Noised signal: ";
        for (int i = 0; i < noised_sig.size(); i++) {
            fout << noised_sig[i] << " ";
        }

        // Demodulation
        vector<int> r;
        r = bpsk_demodulation(noised_sig);
        fout << endl << "Demodulated message: ";
        for (int i = 0; i < r.size(); i++) {
            fout << r[i];
        }

        int real_error_num = 0;
        for (int i = 0; i < ENCODED_LENGTH; i++) {
            if (c[i] != r[i]) real_error_num++;
        }

        // Decode message
        vector<int> v;
        v = decoder(r);
        fout << endl << "Real error num: " << real_error_num;
        fout << endl << "Decoded message: ";
        for (int i = 0; i < v.size(); i++) {
            fout << v[i];
        }
        // Calculate BER
        for (int i = 0; i < MESSAGE_LENGTH; i++) {
            if (u[i] != v[i + ENCODED_LENGTH - MESSAGE_LENGTH]) error++;
        }
        fout << endl << "Error: " << error << endl;
        fout << "BER: " << (double)error / MESSAGE_LENGTH << endl << endl;
        // if (Eb_N0_dB == 1) break;
    }
}