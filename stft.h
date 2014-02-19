#ifndef _STFT_H_
#define _STFT_H_

#include <vector>

struct fftw_complex;
typedef std::vector<fftw_complex*> stft_t;
typedef std::vector<double> istft_t;

stft_t stft(istft_t& data, window_fct_t window_function, uint window_size, uint step_size);
istft_t istft(stft_t& data, window_fct_t window_function, uint window_size, uint step_size);
#endif
