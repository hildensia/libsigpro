#include "stft.h"
#include <fftw3.h>

stft_t phase_adjustment(stft_t& data, int bin_num, double stretching) {
  stft_t stretched(data.size() * stretching);
  for(uint i = 0; i<stretched.size(); ++i) {
    stretched[i] = (fftw_complex*) fftw_malloc(bin_num * sizeof(fftw_complex));
  }

  for (uint bin = 0; bin < bin_num/2 + 1; ++bin) {
    stretched[0][bin][0] = data[0][bin][0];
    stretched[0][bin][1] = data[0][bin][1];
  }


  for (uint frame = 0; frame < stretched.size()-1; frame++) {
    double data_pos = frame / stretching;
    uint data_frame = (uint) data_pos;
    if (data_frame + 1 >= data.size()) break;
    double frac = data_pos - data_frame;

    for (uint bin = 0; bin < bin_num/2 + 1; ++bin) {

      std::complex<double> *curr = reinterpret_cast<std::complex<double>*>(data[data_frame][bin]);
      std::complex<double> *next = reinterpret_cast<std::complex<double>*>(data[data_frame+1][bin]);

      // amplitude interpolation
      double amp = (1-frac) * std::abs(*curr) + frac * std::abs(*next);

      // Calculate phase advance
      double delta_phi = std::arg(*next) - std::arg(*curr);
      delta_phi = fmod(delta_phi + M_PI, 2+M_PI) - M_PI;

      double phase = std::arg(*reinterpret_cast<std::complex<double>*>(stretched[frame][bin])) + delta_phi;

      std::complex<double> s = std::polar(amp, phase);
      stretched[frame + 1][bin][0] = std::real(s);
      stretched[frame + 1][bin][1] = std::imag(s);
    }
  }
  return stretched;
}
