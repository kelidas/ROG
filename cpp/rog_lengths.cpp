#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <chrono>
#include <algorithm>
#include <iterator>
#include <cstdlib>
#include <climits>
#include <unordered_map>
#include "_generate.h"
#include <boost/range/adaptor/map.hpp>
#include <boost/assign.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/format.hpp>

#ifdef _WIN32
#  define DLLEXPORT extern "C" __declspec(dllexport)
#else
#  define DLLEXPORT extern "C"
#endif

class Timer
{
public:
    Timer() : beg_(clock_::now()) {}
    void reset() { beg_ = clock_::now(); }
    double elapsed() const {
        return std::chrono::duration_cast<second_>
            (clock_::now() - beg_).count(); }

private:
    typedef std::chrono::high_resolution_clock clock_;
    typedef std::chrono::duration<double, std::ratio<1> > second_;
    std::chrono::time_point<clock_> beg_;
};

double factorial(double  n);
double power2(int n);

typedef std::unordered_map<double, double> LenMap;

DLLEXPORT int ogrid_p_ae(std::unordered_map<double, double>& len_ae_map,
  std::unordered_map<double, double>& len_pae_map,
  int n, int nv, double* t, bool timeit){
    unsigned int *dt = NULL;
    int gen_dt;

    //alloc memory for vector
    dt = (unsigned int *)malloc(sizeof(unsigned int) * nv);
    if(dt == NULL)
    {
    fprintf(stderr, "error: insufficient memory\n");
    exit(EXIT_FAILURE);
    }

    //initialize
    gen_dt = gen_comb_rep_lex_init(dt, n, nv);
    //skip the first zero vector
    gen_dt = gen_comb_rep_lex_next(dt, n, nv);

    Timer tmr;

    double len_ae=0, len_pae=0;

    int nids, min_dx;
    unsigned nips;
    double nin, nid, nip;
    double ni;
    unsigned long tmp;
    double nipf = factorial(nv);

    if (timeit) tmr.reset();

    while (gen_dt == GEN_NEXT){
        len_ae = 0;
	      len_pae = 0;
        nin = 1;
        nip = 0;
        nid = 0;
        nids = -1;
        nip = nipf;
        min_dx = 0;
        tmp = -1;
        for (int j=0; j<nv; j++){
            nin *= n - dt[j];
            if (dt[j] > 0) nids += 1;

            nips = 0;
            if (tmp < 0 or tmp != dt[j]){
                tmp = dt[j];
                for (int l=0; l<nv; l++){
                    if(dt[l] - tmp == 0) nips += 1;
                };
                nip /= factorial(nips);
            };
        };

        nid = power2(nids);

	ni = (double)nin *(double)nip * (double)nid;
        //for(int x = 0; x < nv; x++)
        //    printf("%u ", dt[x]);
        //cout <<endl;
        gen_dt = gen_comb_rep_lex_next(dt, n, nv);
        };
    if (timeit) *t=tmr.elapsed();
}

int main(int argc, char* argv[], char *envp[]){
    int n = atoi(argv[1]);
    int nv = atoi(argv[2]);

    double t=0;
    LenMap len_ae_map;
    LenMap len_pae_map;
    ogrid_p_ae(len_ae_map, len_pae_map, n, nv, &t, true);

    std::cout.width(20);
    std::cout.precision(std::numeric_limits<double>::digits10 + 1);

    // Store results in binary file
    std::ofstream lens;
    lens.open(boost::str(boost::format("lens_ae_%06i_%02i.bin") % n % nv), std::ios::binary);
    std::ofstream freq;
    freq.open(boost::str(boost::format("freq_ae_%06i_%02i.bin") % n % nv), std::ios::binary);
    std::vector<double> keys_bin;
    // Retrieve all keys
    boost::copy(len_ae_map | boost::adaptors::map_keys, std::back_inserter(keys_bin));
    for(int i = 0; i < keys_bin.size(); i++)
        {lens.write((char *) &keys_bin[i], sizeof(double));
         freq.write((char *) &len_ae_map[keys_bin[i]], sizeof(double));
         //std::cout<< keys[i] << "---" << len_ae_map[keys[i]] << std::endl;
        };
    lens.close();
    freq.close();

    lens.open(boost::str(boost::format("lens_pae_%06i_%02i.bin") % n % nv), std::ios::binary);
    freq.open(boost::str(boost::format("freq_pae_%06i_%02i.bin") % n % nv), std::ios::binary);
    keys_bin.clear();
    // Retrieve all keys
    boost::copy(len_pae_map | boost::adaptors::map_keys, std::back_inserter(keys_bin));
    for(int i = 0; i < keys_bin.size(); i++)
        {lens.write((char *) &keys_bin[i], sizeof(double));
         freq.write((char *) &len_pae_map[keys_bin[i]], sizeof(double));
         //std::cout<< keys[i] << "---" << len_ae_map[keys[i]] << std::endl;
        };
    lens.close();
    freq.close();

    std::cout.fill('#');
    std::cout << std::endl;
    std::cout << "execution clock [s] = " << t << std::endl;
    std::cout.width(20);
    std::cout.fill('_');
    std::cout << std::endl;
    //std::cout << std::numeric_limits<unsigned long>::max() << std::endl;
}


double factorial(double n){
    double factorial=1;
    for(int a=1;a<=n;a++){
        factorial*=a;
    }
    return factorial;
}

double power2(int n) {
    double res = 1;
    for(int i=0;i<=n-1;i++)   res *= 2;
    return res;
}