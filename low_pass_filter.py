from numpy import *
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter, filtfilt, freqz


def low_pass_filter(x, cutoff = 12.4*10, sampling = 10000, order = 6):
      nyquist = sampling/2
      normalize_freq = (cutoff/nyquist)
      b, a = butter(order, normalize_freq,
                  btype = 'lowpass', analog = False)
      w, h = freqz(b, a)
      plt.semilogx(w/(2*pi)*sampling, 20*log10(abs(h)))
      plt.axvline(cutoff, color='k')
      plt.grid(which='both', axis='both')
      plt.title('Butterworth filter frequency response')
      plt.xlabel('Frequency [Hz]')
      plt.ylabel('Amplitude [dB]')
      plt.show()
      # return filtfilt(b, a, x)
      return filtfilt(b, a, x)
