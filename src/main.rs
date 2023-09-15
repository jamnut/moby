use spectrum_analyzer::{samples_fft_to_spectrum, FrequencyLimit};
#[allow(unused_imports)]
use spectrum_analyzer::scaling::divide_by_N_sqrt;
use libm::cosf;
// use core::num;
use std::f32::consts::PI;

fn generate_sine_wave(freq_hz: f32, sampling_rate_hz: f32, num_samples: usize, z : usize) -> Vec<f32> {

    // Calculate the time step
    let time_step = 1.0 / sampling_rate_hz;

    // Calculate the angular frequency (in radians per second)
    let angular_freq = 2.0 * PI * freq_hz;

    // Generate the sine wave samples
    let mut sine_wave = Vec::with_capacity(num_samples);

    for i in 0..num_samples {
        let t = i as f32 * time_step;
        let sample = (angular_freq * t).sin();
        if i < z {
            sine_wave.push(sample);
        } else {
            sine_wave.push(0.0);
        }
    }

    sine_wave
}

fn main() {
    let sin_freq_hz: f32 = 256.0;
    let sampling_rate_hz: u32 = 2*1024;
    let window: usize = 256;
    let num_samples: usize = 1024;
    let bin: f32 = sampling_rate_hz as f32 / num_samples as f32;
    
    let samples: &[f32] = &generate_sine_wave(sin_freq_hz, sampling_rate_hz as f32, num_samples, window);
    
    let hann_window = samples;
    // apply hann window for smoothing; length must be a power of 2 for the FFT
    // let hann_window = hann_window_zero(&samples, window);
    // let hann_window = blackman_harris_4term(samples, window);
    // let hann_window = blackman_harris_7term(samples, window);
    
    // calc spectrum
    let spectrum_hann_window = samples_fft_to_spectrum(
        &hann_window,
        sampling_rate_hz,
        FrequencyLimit::All,
        None, //Some(&divide_by_N_sqrt),
    ).unwrap();

    let lower_bin = sin_freq_hz - 10.0 * bin;
    let upper_bin = sin_freq_hz + 10.0 * bin;

    for (fr, fr_val) in spectrum_hann_window.data().iter() {
        if *fr > lower_bin.into() && *fr < upper_bin.into() {
            println!("{fr:.1}Hz => {fr_val:.3}");
            // println!("{fr:.1},{fr_val:.3}")
        }
    }
    
}



#[must_use]
#[allow(dead_code)]
fn hann_window_zero(samples: &[f32], n: usize) -> Vec<f32> {
    let mut windowed_samples = Vec::with_capacity(samples.len());
    let samples_len_f32 = samples.len() as f32;

    for (i, sample) in samples.iter().enumerate().take(n) {
        let two_pi_i = 2.0 * PI * i as f32;
        let idontknowthename = (two_pi_i / samples_len_f32).cos();
        let multiplier = 0.5 * (1.0 - idontknowthename);
        windowed_samples.push(multiplier * sample);
    }

    // Zero out the remaining samples
    for _ in n..samples.len() {
        windowed_samples.push(0.0);
    }

    windowed_samples
}

#[must_use]
#[allow(dead_code)]
pub fn blackman_harris_4term(samples: &[f32], n: usize) -> Vec<f32> {
    // constants come from here:
    // https://en.wikipedia.org/wiki/Window_function#Blackman%E2%80%93Harris_window
    const ALPHA: [f32; 4] = [0.35875, -0.48829, 0.14128, -0.01168];

    blackman_harris_xterm(samples, &ALPHA, n)
}

/// Applies a Blackman-Harris 7-term window to an array of samples.
///
/// ## More information
/// * <https://en.wikipedia.org/wiki/Window_function#Blackman%E2%80%93Harris_window>
/// * <https://ieeexplore.ieee.org/document/940309>
/// * <https://dsp.stackexchange.com/questions/51095/seven-term-blackman-harris-window>
///
/// ## Return value
/// New vector with Blackman-Harris 7-term window applied to the values.
#[must_use]
#[allow(dead_code)]
pub fn blackman_harris_7term(samples: &[f32], n: usize) -> Vec<f32> {
    // constants come from here:
    // https://dsp.stackexchange.com/questions/51095/seven-term-blackman-harris-window
    const ALPHA: [f32; 7] = [
        0.271_051_4,
        -0.433_297_93,
        0.218_123,
        -0.065_925_45,
        0.010_811_742,
        -0.000_776_584_84,
        0.000_013_887_217,
    ];

    blackman_harris_xterm(samples, &ALPHA, n)
}

/// Applies a Blackman-Harris x-term window
/// (<https://en.wikipedia.org/wiki/Window_function#Blackman%E2%80%93Harris_window>)
/// to an array of samples. The x is specified by `alphas.len()`.
///
/// ## Return value
/// New vector with Blackman-Harris x-term window applied to the values.
#[must_use]
#[allow(dead_code)]
fn blackman_harris_xterm(samples: &[f32], alphas: &[f32], n: usize) -> Vec<f32> {
    let mut windowed_samples = Vec::with_capacity(samples.len());

    let samples_len_f32 = samples.len() as f32;

    for sample in samples.iter().take(n) {
        // Will result in something like that:
        /* ALPHA0
            + ALPHA1 * ((2.0 * PI * *samples[i])/samples_len_f32).cos()
            + ALPHA2 * ((4.0 * PI * *samples[i])/samples_len_f32).cos()
            + ALPHA3 * ((6.0 * PI * *samples[i])/samples_len_f32).cos()
        */

        let mut acc = 0.0;
        for (alpha_i, alpha) in alphas.iter().enumerate() {
            // in 1. iter. 0PI, then 2PI, then 4 PI, then 6 PI
            let two_pi_iteration = 2.0 * alpha_i as f32 * PI;
            let cos = cosf((two_pi_iteration * sample) / samples_len_f32);
            acc += alpha * cos;
        }

        windowed_samples.push(acc)
    }

    // Zero out the remaining samples
    for _ in n..samples.len() {
        windowed_samples.push(0.0);
    }
    
    
    windowed_samples
}
