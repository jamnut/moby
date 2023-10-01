#![deny(clippy::all)]
#![forbid(unsafe_code)]

use std::f32::consts::PI;
use std::ops::Range;

use num::cast::NumCast;
use num::traits::{Float, Num};

use crate::Model;

pub fn _hann(samples: &mut [f32]) {
    let len = samples.len();
    for (i, item) in samples.iter_mut().enumerate().take(len) {
        *item *= 0.5 * (1.0 - (2.0 * PI * i as f32 / (len - 1) as f32).cos());
    }
}

pub fn _median(numbers: &mut [f32]) -> f32 {
    let length = numbers.len();

    if length == 0 {
        return 0.0;
    }

    numbers.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());

    if length % 2 == 0 {
        let left = numbers[length / 2];
        let right = numbers[length / 2 + 1];

        (left + right) / 2.0
    } else {
        numbers[length / 2]
    }
}

fn _to_float<T, T2>(s: T) -> T2
where
    T: Num + NumCast + Copy,
    T2: Float + NumCast + Copy,
{
    NumCast::from(s).unwrap()
}

pub fn maximum<T>(data: &[T]) -> T
where
    T: Float + NumCast + Copy,
{
    data.iter().fold(T::neg_infinity(), |a, b| a.max(*b))
}

// TODO: This really should be abs max, not plain max
pub fn max_window(m: &Model, range: Range<usize>) -> f64 {
    *(m.aiff_data[range].iter().max().unwrap()) as f64
}

pub fn max_all(m: &Model, _range: Range<usize>) -> f64 {
    *(m.aiff_data.iter().max().unwrap()) as f64
}

pub fn max_none(_m: &Model, _range: Range<usize>) -> f64 {
    i16::MAX as f64
}

pub fn fft(m: &Model, range: Range<usize>) -> (Vec<f32>, f32) {
    let dmax: f32 = m.max64[0].0(m, range.clone()) as f32;
    let data: Vec<f32> = m.aiff_data.iter().map(|x| *x as f32 / dmax).collect();

    // hann(&mut samples);

    let fft: Vec<f32> = match m.fft_size[0] {
        1024 => {
            let mut samples: [f32; 1024] = [0.0; 1024];
            samples[0..200].copy_from_slice(&data[range.clone()]);
            let fft = microfft::real::rfft_1024(&mut samples);
            fft[0].im = 0.0;
            fft.iter().map(|c| c.norm()).collect()
        }
        512 => {
            let mut samples: [f32; 512] = [0.0; 512];
            samples[0..200].copy_from_slice(&data[range.clone()]);
            let fft = microfft::real::rfft_512(&mut samples);
            fft[0].im = 0.0;
            fft.iter().map(|c| c.norm()).collect()
        }
        256 => {
            let mut samples: [f32; 256] = [0.0; 256];
            samples[0..200].copy_from_slice(&data[range.clone()]);
            let fft = microfft::real::rfft_256(&mut samples);
            fft[0].im = 0.0;
            fft.iter().map(|c| c.norm()).collect()
        }
        128 => {
            let end = range.end.min(range.start + 128);
            let mut samples: [f32; 128] = [0.0; 128];
            samples[0..128].copy_from_slice(&data[range.start..end]);
            let fft = microfft::real::rfft_128(&mut samples);
            fft[0].im = 0.0;
            fft.iter().map(|c| c.norm()).collect()
        }
        _ => panic!("fft size not supported"),
    };

    let fft: Vec<f32> = fft.iter().map(|x| x * x).collect();

    let bin_size = 2000.0 / m.fft_size[0] as f32;
    let bin_range = (25.0 / bin_size) as usize..(400.0 / bin_size + 1.0) as usize;

    let fmax = maximum(&fft[bin_range]);
    let fft = fft.iter().map(|x| *x / fmax).collect();

    (fft, fmax)
}

pub fn spectrogram(m: &Model) -> (Vec<Vec<f32>>, Vec<f32>, f32) {
    let mut spec: Vec<Vec<f32>> = Vec::new();
    let mut fmax: Vec<f32> = Vec::new();
    let mut smax: f32 = 0.0;

    for i in (0..4000 - m.window).step_by(m.slide[0]) {
        let range = i..i + m.window;
        let (fft, max) = fft(m, range);

        spec.push(fft);
        fmax.push(max);
        smax = smax.max(max);
    }

    (spec, fmax, smax)
}


// coorelates two sprectrograms
pub fn _correlate(s: &Vec<Vec<f32>>, t: &Vec<Vec<f32>>) -> Option<usize> {
    if t.len() > s.len() || t[0].len() > s[0].len() {
        return None;
    }
    let mut max_corr = 0.0;
    let mut max_index = None;
    for i in 0..s.len()-t.len()+1 {
        for j in 0..s[i].len()-t[0].len()+1 {
            let mut sum = 0.0;
            for k in 0..t.len() {
                for l in 0..t[k].len() {
                    sum += s[i+k][j+l] * t[k][l];
                }
            }
            if sum > max_corr {
                max_corr = sum;
                max_index = Some(i);
            }
        }
    }
    max_index
}


// Apply the Sobel operator to detect edges in the spectrogram
pub fn _sobel_operator(image: &Vec<Vec<f32>>) -> Vec<Vec<f32>> {
    // Convert the spectrogram to a grayscale image
    let mut output_image = image.clone();
    for y in output_image.iter_mut() {
        for x in y.iter_mut() {
            *x = 0.0;
        }
    }
    // Apply the Sobel operator to detect edges in the grayscale image
    let kernel_x = [[-1, 0, 1], [-2, 0, 2], [-1, 0, 1]];
    let kernel_y = [[-1, -2, -1], [0, 0, 0], [1, 2, 1]];

    for x in 1..image.len()-1 {
        for y in 1..image[0].len()-1 {
            let mut gx = 0.0;
            let mut gy = 0.0;
            for i in 0..3 {
                for j in 0..3 {
                    let pixel = image[x+i-1][y+j-1];
                    gx += pixel * kernel_x[i][j] as f32;
                    gy += pixel * kernel_y[i][j] as f32;
                }
            }
            output_image[x][y] = ((gx.powf(2.) + gy.powf(2.)) as f32).sqrt();
        }
    }

    output_image
}


