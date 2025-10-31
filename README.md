🩺 ECG Analizi: Normal Bireyler ve Yoğun Bakım Hastaları
MIT-BIH Arrhythmia Database ve MIMIC-IV-ECG Verilerinin Karşılaştırılması

Python tabanlı biyomedikal sinyal işleme projesi

📘 Genel Bakış

Bu proje, sağlıklı bireylerden (MIT-BIH Arrhythmia Database) ve yoğun bakım hastalarından (MIMIC-IV-ECG) alınan elektrokardiyogram (ECG) sinyallerini analiz eder.
Sinyalleri filtreleyerek, gürültü tespiti yaparak ve R-tepe noktalarını bularak kalp atım hızını, sinyal enerjisini ve frekans özelliklerini hesaplar.
Görsel ve istatistiksel karşılaştırmalarla normal ve kritik durumlar arasındaki fizyolojik farkları ortaya koyar.

⚙️ Özellikler

✅ MIT-BIH ve MIMIC-IV-ECG verilerini yükleme

✅ Gürültü azaltma için bant geçiren filtre (0.5–45 Hz)

✅ Otomatik R-tepe tespiti ve kalp atım hızı hesaplama

✅ Gürültü analizi: temel kayma, EMG artefaktları, ve şebeke paraziti

✅ Frekans spektrumu görselleştirmesi (0–60 Hz)

✅ Karşılaştırmalı grafikler ve detaylı rapor çıktısı

🧩 Gereksinimler

Python ≥ 3.8 sürümüne ve aşağıdaki kütüphanelere ihtiyacınız var:

    pip install wfdb numpy scipy matplotlib

⚠️ PhysioNet üzerindeki MIMIC-IV-ECG verilerine erişim için kimlik doğrulaması (credential) gereklidir.

🚀 Nasıl Çalıştırılır

1. Depoyu klonlayın:

   git clone https://github.com/kullaniciadiniz/ecg-analysis-mitbih-vs-mimic.git
   
   cd ecg-analysis-mitbih-vs-mimic

2. Main script’i (ana dosyayı) çalıştırın:

   python ecg_analysis.py

3. Program şunları yapar:

   • ECG verilerini yükler (varsayılan: MIT-BIH kayıtları 100 ve 200)

   • Sinyal analizi ve karşılaştırma gerçekleştirir

   • Terminalde detaylı bir rapor oluşturur

   • Karşılaştırma grafiğini ecg_comparison.png olarak kaydeder


