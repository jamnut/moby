#![deny(clippy::all)]
#![forbid(unsafe_code)]

use std::f32::consts::PI;
use std::ops::Range;

use num::cast::NumCast;
use num::traits::{Float, Num};

use crate::{Model, WHALE_RANGE};



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

pub fn _maximum<T>(data: &[T]) -> T
where
    T: Float + NumCast + Copy,
{
    data.iter().fold(T::neg_infinity(), |a, b| a.max(*b))
}

// TODO: This really should be abs max, not plain max
pub fn max_window(m: &Model, range: Range<usize>) -> f32 {
    *(m.aiff_data[range].iter().max().unwrap()) as f32
}

pub fn max_all(m: &Model, _range: Range<usize>) -> f32 {
    *(m.aiff_data.iter().max().unwrap()) as f32
}

pub fn max_none(_m: &Model, _range: Range<usize>) -> f32 {
    i16::MAX as f32
}

pub fn fft(m: &Model, start: usize) -> (Vec<f32>, f32, f32, f32) {
    
    let dmax: f32 = m.max64[0].0(m, m.range()) as f32;
    let data: Vec<f32> = m.aiff_data.iter().map(|x| *x as f32 / dmax).collect();

    let fft_fn = || {
        let size = m.fft_size[0];
        let width = size.min(m.window);
        let mut samples = vec![0.0; size];
        samples[0..width].copy_from_slice(&data[start..(start + width)]);
        let mut fft = match size {
            2048 => microfft::real::rfft_2048(&mut samples.try_into().unwrap()).to_vec(),
            1024 => microfft::real::rfft_1024(&mut samples.try_into().unwrap()).to_vec(),
            512 => microfft::real::rfft_512(&mut samples.try_into().unwrap()).to_vec(),
            256 => microfft::real::rfft_256(&mut samples.try_into().unwrap()).to_vec(),
            128 => microfft::real::rfft_128(&mut samples.try_into().unwrap()).to_vec(),
            _ => panic!("Unsupported FFT size: {}", size),
        };    
        fft[0].im = 0.0;
        fft.iter().map(|c| c.norm()).collect::<Vec<f32>>()
    };    

    // hann(&mut samples);

    let fft: Vec<f32> = fft_fn();
    let fft: Vec<f32> = fft.iter().map(|x| x * x).collect();

    let bin_size = 2000.0 / fft.len() as f32;
    let whale_bins
        = (WHALE_RANGE.start / bin_size) as usize .. (WHALE_RANGE.end / bin_size) as usize;

    // Calculate peak power and frequency
    let (fmax, peak) = fft[whale_bins.clone()].iter().enumerate().fold(
        (0.0, 0),
        |(fmax, peak), (i, f)| {
            if *f > fmax {
                (*f, i)
            } else {
                (fmax, peak)
            }
        },
    );
    let fpeak = (peak + whale_bins.start) as f32 * bin_size;

    // Use median to estimate the snr
    let mut snr = fft[whale_bins].to_vec();
    snr.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let spower = snr.last().unwrap();
    let npower = snr[snr.len() / 2];
    let snr = spower / npower;

    // Normalize the fft
    let fft: Vec<f32> = fft.iter().map(|x| *x / fmax).collect();

    (fft, fmax, fpeak, snr)
}


pub fn spectrogram(m: &Model) -> (Vec<Vec<f32>>, Vec<f32>, Vec<f32>, Vec<f32>, f32) {
    let mut spec: Vec<Vec<f32>> = Vec::new();
    let mut fmax: Vec<f32> = Vec::new();
    let mut fpeak: Vec<f32> = Vec::new();
    let mut fsnr: Vec<f32> = Vec::new();
    let mut smax: f32 = 0.0;

    for start in (0..4000 - m.window).step_by(m.slide[0]) {
        let (fft, max, peak, snr) = fft(m, start);

        spec.push(fft);
        fmax.push(max);
        fpeak.push(peak);
        fsnr.push(snr);
        smax = smax.max(max);
    }

    (spec, fmax, fpeak, fsnr, smax)
}


fn find_peaks(fft: &[f32]) -> Vec<usize> {
    let bin_size = 2000.0 / fft.len() as f32;
    let whale_bins = (WHALE_RANGE.start / bin_size) as usize .. (WHALE_RANGE.end / bin_size) as usize;

    let mut peaks: Vec<usize> = whale_bins.filter(|i| {
        fft[*i] > fft[*i-1] && fft[*i] > fft[*i+1]
    }).collect();
    peaks.sort_by(|a, b| fft[*b].partial_cmp(&fft[*a]).unwrap());

    peaks[0..peaks.len().min(3)].to_vec()
}

pub fn tracking(spec: &Vec<Vec<f32>>) -> Vec<Vec<usize>> {
    spec.iter().map(|fft| { find_peaks(&fft) }).collect()
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


pub fn _normalize(data: &[f32]) -> Vec<f32> {
    let max = data.iter().fold(0.0, |a, b| a.max(*b));
    data.iter().map(|x| x / max).collect()
}

