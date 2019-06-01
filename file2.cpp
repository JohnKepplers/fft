#include <bits/stdc++.h>
#include <complex>
#include <iostream>
#include <valarray>
#include <chrono>
#include <sstream>
#include <cmath>
#include <fstream>
using namespace std;


typedef complex<double> base;
typedef std::valarray<base> my_array;
const double PI = 3.141592653589793238460;


 
void fft (my_array & a, bool invert, int k) {
	int n = (int) a.size();
	
	for (int i=1, j=0; i<n; ++i) {
		int bit = n >> 1;
		for (; j>=bit; bit>>=1)
			j -= bit;
		j += bit;
		if (i < j)
			swap (a[i], a[j]);
	}
 
	for (int len=2; len<=n; len<<=1) {
		double ang = 2*PI/len * (invert ? -1 : 1);
		base wlen (cos(ang), sin(ang));
		for (int i=0; i<n; i+=len) {
			base w (1);
			for (int j=0; j<len/2; ++j) {
				base u = a[i+j],  v = a[i+j+len/2] * w;
				a[i+j] = u + v;
				a[i+j+len/2] = u - v;
				w *= wlen;
			}
		}
	}
	if (invert)
		for (int i=0; i<n; ++i)
			a[i] /= n;
	base temp = a[n-k];
	a[n - k] = a[n-k] - a[k];
	a[k] = temp + a[k];

}

int main(int argc, char *argv[])
{
  if (argc == 0)
	{
        std::cout << "Error" << '\n';
        exit(1);
    }

  if (argc >= 2)
    {
        std::istringstream iss( argv[1] );
        int val;
        if ((iss >> val) && iss.eof()) 
        {
          
        }
		int k = 8;
        const int N = pow(2, val);
        base* test = new base[N];

        for (int i = 0; i < N; ++i)
        {
             test[i] = cos(k * 2*i* PI/N), sin(k*2*PI*i/N);
        }
        my_array data(test, N);
		ofstream fout;
		fout.open("original_data.txt");
     	for (int i = 0; i < N; ++i)
     	{
			fout << data[i] << '\n';

     	}
		fout.close();

        fft(data, false, k);

		ofstream ffout;
		fout.open("after_fft_data.txt");
		std::cout << std::endl << "Result of forward fft:" << std::endl;
		for (int i = 0; i < N; ++i)
		{
			fout << data[i] << '\n';
		}
		fout.close();
		fft(data, true, k);

		ofstream ifout;
		ifout.open("after_ifft_data.txt");
		for (int i = 0; i < N; ++i)
		{
			ifout << data[i] << '\n';
		}
		ifout.close();

		// Time measurement:

        double time_taken = 0;
        for (int i = 0; i < 10; ++i)
        {		
            auto start = std::chrono::high_resolution_clock::now();
            fft(data, false, k);
            auto end = std::chrono::high_resolution_clock::now();
            auto time_taken_1 = end - start;
            std::cout << "Running took "
            << std::chrono::duration_cast<std::chrono::milliseconds>(time_taken_1).count() << " ms\n";
            time_taken += time_taken_1 / std::chrono::milliseconds(1);
        }
        printf("Average time is %10f ms\n", time_taken / 10);
	}
    return 0;
}
