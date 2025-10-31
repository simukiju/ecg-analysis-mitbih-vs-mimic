"""
ECG Veri Analizi: Normal Bireyler vs Yoğun Bakım Hastaları
MIT-BIH Arrhythmia Database ve MIMIC-IV-ECG karşılaştırması

"""

import wfdb
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.signal import butter, filtfilt, find_peaks
import warnings
warnings.filterwarnings('ignore')

class ECGAnalyzer:
    def __init__(self, sampling_rate=360):
      
        self.fs = sampling_rate
        
    def load_mitdb_record(self, record_name='100', duration=10, channel=0):
     
        try:
            
            record = wfdb.rdrecord(record_name, 
                                   pn_dir='mitdb',
                                   sampfrom=0, 
                                   sampto=int(duration * self.fs))
            
            
            self.fs = record.fs
            
           
            signal_data = record.p_signal[:, channel]
            
            return {
                'signal': signal_data,
                'fs': self.fs,
                'units': record.units[channel],
                'sig_name': record.sig_name[channel],
                'record_name': record_name,
                'duration': duration
            }
        except Exception as e:
            print(f"Hata: MIT-BIH kaydı yüklenemedi - {e}")
            return None
    
    def load_mimic_record(self, record_path, duration=10, channel=0):
      
        try:
           
            record = wfdb.rdrecord(record_path,
                                   sampfrom=0,
                                   sampto=int(duration * self.fs))
            
            self.fs = record.fs
            signal_data = record.p_signal[:, channel]
            
            return {
                'signal': signal_data,
                'fs': self.fs,
                'units': record.units[channel] if hasattr(record, 'units') else 'mV',
                'sig_name': record.sig_name[channel] if hasattr(record, 'sig_name') else 'ECG',
                'record_name': record_path.split('/')[-1],
                'duration': duration
            }
        except Exception as e:
            print(f"Hata: MIMIC kaydı yüklenemedi - {e}")
            print("Not: MIMIC-IV-ECG erişimi için PhysioNet credential gereklidir")
            return None
    
    def normalize_signal(self, signal_data):
       
        mean_val = np.mean(signal_data)
        std_val = np.std(signal_data)
        normalized = (signal_data - mean_val) / std_val
        return normalized
    
    def calculate_energy_power(self, signal_data):
    
        N = len(signal_data)
       
        energy = np.sum(signal_data ** 2)
        
     
        power = energy / N
        
        return {
            'energy': energy,
            'power': power,
            'rms': np.sqrt(power) 
        }
    
    def design_bandpass_filter(self, lowcut=0.5, highcut=45, order=4):
       
        nyquist = 0.5 * self.fs
        low = lowcut / nyquist
        high = highcut / nyquist
        
        b, a = butter(order, [low, high], btype='band')
        return b, a
    
    def apply_filter(self, signal_data, b, a):
       
        filtered = filtfilt(b, a, signal_data)
        return filtered
    
    def detect_noise(self, signal_data):
       
        fft_vals = np.fft.fft(signal_data)
        fft_freq = np.fft.fftfreq(len(signal_data), 1/self.fs)
      
        pos_mask = fft_freq > 0
        fft_freq_pos = fft_freq[pos_mask]
        fft_magnitude = np.abs(fft_vals[pos_mask])
        
        noise_report = {}
        
        powerline_50hz = np.sum(fft_magnitude[(fft_freq_pos > 48) & (fft_freq_pos < 52)])
        powerline_60hz = np.sum(fft_magnitude[(fft_freq_pos > 58) & (fft_freq_pos < 62)])
        
        if powerline_50hz > powerline_60hz:
            noise_report['powerline'] = {'freq': 50, 'magnitude': powerline_50hz}
        else:
            noise_report['powerline'] = {'freq': 60, 'magnitude': powerline_60hz}
        
        high_freq_noise = np.sum(fft_magnitude[fft_freq_pos > 45])
        noise_report['high_freq_emg'] = high_freq_noise
        
        baseline_wander = np.sum(fft_magnitude[fft_freq_pos < 0.5])
        noise_report['baseline_wander'] = baseline_wander
        
        signal_band = np.sum(fft_magnitude[(fft_freq_pos > 0.5) & (fft_freq_pos < 45)])
        total_power = np.sum(fft_magnitude)
        noise_report['snr_estimate'] = signal_band / (total_power - signal_band + 1e-10)
        
        return noise_report
    
    def detect_r_peaks(self, signal_data, filtered=True):
       
        if not filtered:
            b, a = self.design_bandpass_filter()
            signal_data = self.apply_filter(signal_data, b, a)
   
        diff_signal = np.diff(signal_data)
        
        squared_signal = diff_signal ** 2
        
        window_size = int(0.12 * self.fs) 
        integrated = np.convolve(squared_signal, np.ones(window_size)/window_size, mode='same')
        
        min_distance = int(0.6 * self.fs)
        
        threshold = 0.35 * np.max(integrated)
        
        peaks, properties = find_peaks(integrated, 
                                       height=threshold,
                                       distance=min_distance)
      
        r_peaks = peaks
        
        return r_peaks
    
    def calculate_heart_rate(self, r_peaks):
       
        if len(r_peaks) < 2:
            return None
        
        rr_intervals = np.diff(r_peaks) / self.fs
        
        mean_rr = np.mean(rr_intervals)
        hr_mean = 60 / mean_rr
        
        hr_std = np.std(60 / rr_intervals)
        
        return {
            'mean_hr': hr_mean,
            'std_hr': hr_std,
            'rr_intervals': rr_intervals,
            'mean_rr': mean_rr
        }
    
    def plot_comparison(self, data_normal, data_icu, save_path=None):
       
        fig, axes = plt.subplots(4, 2, figsize=(16, 14))
        fig.suptitle('ECG Analizi: Normal Birey vs Yoğun Bakım Hastası', 
                     fontsize=16, fontweight='bold')
        
        time_normal = np.arange(len(data_normal['raw'])) / self.fs
        time_icu = np.arange(len(data_icu['raw'])) / self.fs
        
        axes[0, 0].plot(time_normal, data_normal['raw'], 'b-', linewidth=0.8)
        axes[0, 0].set_title('Normal Birey - Ham ECG Sinyali')
        axes[0, 0].set_xlabel('Zaman (s)')
        axes[0, 0].set_ylabel('Genlik')
        axes[0, 0].grid(True, alpha=0.3)
        
        axes[0, 1].plot(time_icu, data_icu['raw'], 'r-', linewidth=0.8)
        axes[0, 1].set_title('Yoğun Bakım - Ham ECG Sinyali')
        axes[0, 1].set_xlabel('Zaman (s)')
        axes[0, 1].set_ylabel('Genlik')
        axes[0, 1].grid(True, alpha=0.3)
        
        axes[1, 0].plot(time_normal, data_normal['raw'], 'b-', 
                       alpha=0.5, label='Ham', linewidth=0.8)
        axes[1, 0].plot(time_normal, data_normal['filtered'], 'g-', 
                       label='Filtrelenmiş', linewidth=0.8)
        axes[1, 0].set_title('Normal - Ham vs Filtrelenmiş')
        axes[1, 0].set_xlabel('Zaman (s)')
        axes[1, 0].set_ylabel('Genlik')
        axes[1, 0].legend()
        axes[1, 0].grid(True, alpha=0.3)
        
        axes[1, 1].plot(time_icu, data_icu['raw'], 'r-', 
                       alpha=0.5, label='Ham', linewidth=0.8)
        axes[1, 1].plot(time_icu, data_icu['filtered'], 'g-', 
                       label='Filtrelenmiş', linewidth=0.8)
        axes[1, 1].set_title('Yoğun Bakım - Ham vs Filtrelenmiş')
        axes[1, 1].set_xlabel('Zaman (s)')
        axes[1, 1].set_ylabel('Genlik')
        axes[1, 1].legend()
        axes[1, 1].grid(True, alpha=0.3)
        
        axes[2, 0].plot(time_normal, data_normal['filtered'], 'b-', linewidth=0.8)
        axes[2, 0].plot(time_normal[data_normal['r_peaks']], 
                       data_normal['filtered'][data_normal['r_peaks']], 
                       'ro', markersize=8, label='R-tepeleri')
        axes[2, 0].set_title(f"Normal - R-tepe Tespiti (HR: {data_normal['hr']:.1f} bpm)")
        axes[2, 0].set_xlabel('Zaman (s)')
        axes[2, 0].set_ylabel('Genlik')
        axes[2, 0].legend()
        axes[2, 0].grid(True, alpha=0.3)
        
        axes[2, 1].plot(time_icu, data_icu['filtered'], 'r-', linewidth=0.8)
        axes[2, 1].plot(time_icu[data_icu['r_peaks']], 
                       data_icu['filtered'][data_icu['r_peaks']], 
                       'ro', markersize=8, label='R-tepeleri')
        axes[2, 1].set_title(f"Yoğun Bakım - R-tepe Tespiti (HR: {data_icu['hr']:.1f} bpm)")
        axes[2, 1].set_xlabel('Zaman (s)')
        axes[2, 1].set_ylabel('Genlik')
        axes[2, 1].legend()
        axes[2, 1].grid(True, alpha=0.3)
        
        fft_normal = np.fft.fft(data_normal['raw'])
        freq_normal = np.fft.fftfreq(len(data_normal['raw']), 1/self.fs)
        pos_mask_n = freq_normal > 0
        
        fft_icu = np.fft.fft(data_icu['raw'])
        freq_icu = np.fft.fftfreq(len(data_icu['raw']), 1/self.fs)
        pos_mask_i = freq_icu > 0
        
        axes[3, 0].plot(freq_normal[pos_mask_n], 
                       np.abs(fft_normal[pos_mask_n]), 'b-', linewidth=0.8)
        axes[3, 0].set_title('Normal - Frekans Spektrumu')
        axes[3, 0].set_xlabel('Frekans (Hz)')
        axes[3, 0].set_ylabel('Genlik')
        axes[3, 0].set_xlim([0, 60])
        axes[3, 0].grid(True, alpha=0.3)
        
        axes[3, 1].plot(freq_icu[pos_mask_i], 
                       np.abs(fft_icu[pos_mask_i]), 'r-', linewidth=0.8)
        axes[3, 1].set_title('Yoğun Bakım - Frekans Spektrumu')
        axes[3, 1].set_xlabel('Frekans (Hz)')
        axes[3, 1].set_ylabel('Genlik')
        axes[3, 1].set_xlim([0, 60])
        axes[3, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        plt.show()
    
    def generate_report(self, data_normal, data_icu):
        
        print("=" * 80)
        print("ECG ANALİZ RAPORU")
        print("=" * 80)
        
        print("\n### NORMAL BİREY ###")
        print(f"Kayıt: {data_normal['record_name']}")
        print(f"Örnekleme Frekansı: {self.fs} Hz")
        print(f"Sinyal Enerjisi: {data_normal['energy']:.2f}")
        print(f"Sinyal Gücü: {data_normal['power']:.4f}")
        print(f"RMS: {data_normal['rms']:.4f}")
        print(f"Ortalama Kalp Hızı: {data_normal['hr']:.1f} ± {data_normal['hr_std']:.1f} bpm")
        print(f"Ortalama R-R Aralığı: {data_normal['rr_mean']:.3f} s")
        print(f"Tespit Edilen R-tepe Sayısı: {len(data_normal['r_peaks'])}")
        
        print("\nGürültü Analizi:")
        print(f"  - Powerline ({data_normal['noise']['powerline']['freq']} Hz): "
              f"{data_normal['noise']['powerline']['magnitude']:.2f}")
        print(f"  - Kas Artefaktı (EMG): {data_normal['noise']['high_freq_emg']:.2f}")
        print(f"  - Baseline Wander: {data_normal['noise']['baseline_wander']:.2f}")
        print(f"  - SNR Tahmini: {data_normal['noise']['snr_estimate']:.2f}")
        
        print("\n### YOĞUN BAKIM HASTASI ###")
        print(f"Kayıt: {data_icu['record_name']}")
        print(f"Örnekleme Frekansı: {self.fs} Hz")
        print(f"Sinyal Enerjisi: {data_icu['energy']:.2f}")
        print(f"Sinyal Gücü: {data_icu['power']:.4f}")
        print(f"RMS: {data_icu['rms']:.4f}")
        print(f"Ortalama Kalp Hızı: {data_icu['hr']:.1f} ± {data_icu['hr_std']:.1f} bpm")
        print(f"Ortalama R-R Aralığı: {data_icu['rr_mean']:.3f} s")
        print(f"Tespit Edilen R-tepe Sayısı: {len(data_icu['r_peaks'])}")
        
        print("\nGürültü Analizi:")
        print(f"  - Powerline ({data_icu['noise']['powerline']['freq']} Hz): "
              f"{data_icu['noise']['powerline']['magnitude']:.2f}")
        print(f"  - Kas Artefaktı (EMG): {data_icu['noise']['high_freq_emg']:.2f}")
        print(f"  - Baseline Wander: {data_icu['noise']['baseline_wander']:.2f}")
        print(f"  - SNR Tahmini: {data_icu['noise']['snr_estimate']:.2f}")
        
        print("\n### KARŞILAŞTIRMA ###")
        hr_diff = abs(data_normal['hr'] - data_icu['hr'])
        print(f"Kalp Hızı Farkı: {hr_diff:.1f} bpm")
        
        energy_ratio = data_icu['energy'] / data_normal['energy']
        print(f"Enerji Oranı (ICU/Normal): {energy_ratio:.2f}")
        
        print("\n" + "=" * 80)


def main():
    analyzer = ECGAnalyzer(sampling_rate=360)
    
    print("ECG verilerini yüklüyor...")
    
    normal_data = analyzer.load_mitdb_record(record_name='100', duration=10, channel=0)
    
    if normal_data is None:
        print("HATA: Normal birey verisi yüklenemedi!")
        return
    
    print("\nNOT: MIMIC-IV-ECG erişimi için credential gerektiğinden,")
    print("demo amaçlı MIT-BIH'den ikinci bir kayıt kullanılıyor (kayıt 200).")
    print("Gerçek analizde MIMIC verisini kullanmalısınız!\n")
    
    icu_data = analyzer.load_mitdb_record(record_name='200', duration=10, channel=0)
    
    if icu_data is None:
        print("HATA: Yoğun bakım verisi yüklenemedi!")
        return
    
    print("Veriler başarıyla yüklendi!\n")
    
    print("Sinyaller işleniyor...\n")
    
    normal_normalized = analyzer.normalize_signal(normal_data['signal'])
    normal_energy_power = analyzer.calculate_energy_power(normal_data['signal'])
    normal_noise = analyzer.detect_noise(normal_data['signal'])
    
    b, a = analyzer.design_bandpass_filter(lowcut=0.5, highcut=45, order=4)
    normal_filtered = analyzer.apply_filter(normal_data['signal'], b, a)
    normal_r_peaks = analyzer.detect_r_peaks(normal_filtered, filtered=True)
    normal_hr_info = analyzer.calculate_heart_rate(normal_r_peaks)
    
    icu_normalized = analyzer.normalize_signal(icu_data['signal'])
    icu_energy_power = analyzer.calculate_energy_power(icu_data['signal'])
    icu_noise = analyzer.detect_noise(icu_data['signal'])
    
    icu_filtered = analyzer.apply_filter(icu_data['signal'], b, a)
    icu_r_peaks = analyzer.detect_r_peaks(icu_filtered, filtered=True)
    icu_hr_info = analyzer.calculate_heart_rate(icu_r_peaks)
    
    results_normal = {
        'raw': normal_data['signal'],
        'normalized': normal_normalized,
        'filtered': normal_filtered,
        'r_peaks': normal_r_peaks,
        'energy': normal_energy_power['energy'],
        'power': normal_energy_power['power'],
        'rms': normal_energy_power['rms'],
        'hr': normal_hr_info['mean_hr'],
        'hr_std': normal_hr_info['std_hr'],
        'rr_mean': normal_hr_info['mean_rr'],
        'noise': normal_noise,
        'record_name': normal_data['record_name']
    }
    
    results_icu = {
        'raw': icu_data['signal'],
        'normalized': icu_normalized,
        'filtered': icu_filtered,
        'r_peaks': icu_r_peaks,
        'energy': icu_energy_power['energy'],
        'power': icu_energy_power['power'],
        'rms': icu_energy_power['rms'],
        'hr': icu_hr_info['mean_hr'],
        'hr_std': icu_hr_info['std_hr'],
        'rr_mean': icu_hr_info['mean_rr'],
        'noise': icu_noise,
        'record_name': icu_data['record_name']
    }
    
    analyzer.generate_report(results_normal, results_icu)
    
    print("\nGrafikler oluşturuluyor...")
    analyzer.plot_comparison(results_normal, results_icu, save_path='ecg_comparison.png')
    print("\nAnaliz tamamlandı! Grafik 'ecg_comparison.png' olarak kaydedildi.")


if __name__ == "__main__":
    main()
