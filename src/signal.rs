#![deny(clippy::all)]
#![forbid(unsafe_code)]

use std::f32::consts::PI;

use num::cast::NumCast;
use num::traits::{Num, Float};

fn to_f64<T>(s: T) -> f64
where   
T: Num + NumCast + Copy,
{
    NumCast::from(s).unwrap()
}

fn to_float<T, T2>(s: T) -> T2
where
T: Num + NumCast + Copy,
T2: Float + NumCast + Copy,
{
    NumCast::from(s).unwrap()
}

pub fn max<T>(data: &[T]) -> T
where
T: Float + NumCast + Copy,
{
    data.iter().fold(T::neg_infinity(), |a, b| a.max(*b))
}

pub fn normalize<T, T2>(data: &[T]) -> Vec<T2>
where
T: Num + NumCast + Copy,
T2: Float + NumCast + Copy,
{
    let data: Vec<T2> = data.iter().map(|s| to_float(*s)).collect();

    let amax = max(&data);

    data.iter().map(|x| *x / amax).collect()
}    

pub fn series<T>(aiff: &[T]) -> Vec<(usize, f64)> 
where
T: Num + NumCast + Copy,
{
    aiff
    .iter()
    .enumerate()
    .map(|(i, s)| (i, to_f64(*s)))
    .collect()
}

pub fn _hann(samples: &mut [f32]) {
    let len = samples.len();
    for (i, item) in samples.iter_mut().enumerate().take(len) {
        *item *= 0.5 * (1.0 - (2.0 * PI * i as f32 / (len - 1) as f32).cos());
    }
}    

pub fn amplitudes(window: usize, aiff: &[f32]) -> Vec<f32> {

  let mut samples: [f32; 1024] = [0.0; 1024];  
  samples[0..200].copy_from_slice(&aiff[window..window+200]);
//   hann(&mut samples);

  let spectrum = microfft::real::rfft_1024(&mut samples);
  spectrum[0].im = 0.0;

  spectrum.iter().map(|c| c.norm()).collect()
}  

pub fn _median(numbers: &mut Vec<i32>) -> i32 {
    let length = numbers.len();

    if length == 0 {
        return 0;
    }

    numbers.sort_unstable();

    if length % 2 == 0 {
        let left = numbers[length / 2];
        let right = numbers[length / 2 + 1];

        (left + right) / 2
    } else {
        numbers[length / 2]
    }
}
