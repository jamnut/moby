#![deny(clippy::all)]
#![forbid(unsafe_code)]

use std::f32::consts::PI;

// function to apply hann window to a slice of samples
pub fn hann(samples: &mut [f32]) {
    let len = samples.len();
    for (i, item) in samples.iter_mut().enumerate().take(len) {
        *item *= 0.5 * (1.0 - (2.0 * PI * i as f32 / (len - 1) as f32).cos());
    }
}

pub fn normalize(aiff: &[f64]) -> (Vec<f64>, f64) {
    let amax: Vec<f64> = aiff.iter().map(|s| s.abs()).collect();
    let amax = *amax.iter().max_by(|a, b| a.total_cmp(b)).unwrap();
    let amax = (amax * 1.1).ceil();

    (
        aiff.iter().map(|s| *s / amax).collect(),
        amax,
    )
}

pub fn signal(aiff: &[f64]) -> Vec<(usize, f64)> {
  aiff
    .iter()
    .enumerate()
    .map(|(i, s)| (i, *s as f64))
    .collect()
}

pub fn amplitudes(window: usize, aiff: &[f32]) -> Vec<f64> {

  let mut samples: [f32; 1024] = [0.0; 1024];
  samples[0..window].copy_from_slice(&aiff[..window]);
  hann(&mut samples);

  let spectrum = microfft::real::rfft_1024(&mut samples);
  spectrum[0].im = 0.0;

  spectrum.iter().map(|c| c.norm() as f64).collect()
}