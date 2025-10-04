# 📻 SDR FM Receiver

SDR FM Receiver is a software-defined radio system built on a Raspberry Pi 4 with an RF dongle for real-time reception of FM mono/stereo audio and RDS (Radio Data System) data. It replaces traditional radio hardware components with flexible, reconfigurable software to process and decode live radio signals.

## 🚀 Features

- 📡 **FM Demodulation**  
  Extracts mono, stereo, and RDS sub-channels (0–100 kHz) from raw RF samples.

- 🔊 **Mono Audio Processing**  
  Low-pass filtering and resampling for clear 48 kHz audio playback.

- 🎵 **Stereo Audio Processing**  
  Uses PLL-based stereo pilot recovery, band-pass filtering, and resampling to separate left and right channels.

- 🛰️ **RDS Decoding**  
  Demodulates and extracts embedded digital information (station name, program ID, song metadata).

- 🧮 **Resampler & FIR Filtering**  
  Efficient block processing with linear-phase FIR filters for clean and stable audio.

- ⏱️ **Real-Time Reception**  
  Processes I/Q samples from RF hardware and outputs live FM audio in real time.

## 🛠️ Tech Stack

- 🖥️ **Hardware:** Raspberry Pi 4, RTL-SDR dongle, antenna  
- 💻 **Languages:** C++, Python  
- ⚙️ **Core DSP Components:** FIR filters, PLL, resampling, FM demodulator  
- 🧰 **Tools:** Git, custom DSP implementations in C++/Python  

## 📊 Performance Metrics

- 🔢 1111–1313 multiplications per output sample for mono processing (mode-dependent)
- ➕ +102 multiplications per sample for stereo processing
- ✔️ Accurate RDS carrier recovery and symbol resampling achieved
- ⏳ Stereo resampling and PLL runtime: ~6.4–10.2 ms per block depending on mode

## 🧰 Future Improvements

- 🖥️ Add a user interface to display signal strength, station metadata, and song information
- ⚡ Optimize runtime with multi-threading for simultaneous stereo and RDS processing
- 💾 Reduce memory usage by eliminating redundant vectors and allocations
