#![deny(clippy::all)]
#![forbid(unsafe_code)]

use std::f32::consts::PI;
use std::ops::Range;

use crate::{Model, WHALE_RANGE};

use derive_more::{Add, Display};

#[derive (Add, Display, Debug, Clone, Copy, PartialEq)]
pub struct Bin(usize);

#[derive (Add, Display, Debug, Clone, Copy, PartialEq)]
pub struct Freq(usize);



pub fn hann(samples: &mut [f32]) {
    let len = samples.len();
    for (i, item) in samples.iter_mut().enumerate().take(len) {
        *item *= 0.5 * (1.0 - (2.0 * PI * i as f32 / (len - 1) as f32).cos());
    }
}

pub fn no_windowing(_: &mut [f32]) {}


// TODO: This really should be abs max, not plain max
pub fn _max_window(m: &Model, range: Range<usize>) -> f32 {
    *(m.aiff_data[range].iter().max().unwrap()) as f32
}

pub fn max_all(m: &Model, _range: Range<usize>) -> f32 {
    *(m.aiff_data.iter().max().unwrap()) as f32
}

pub fn max_none(_m: &Model, _range: Range<usize>) -> f32 {
    i16::MAX as f32
}

pub fn old_fft(m: &Model, start: usize) -> (Vec<f32>, f32, f32, (f32, f32)) {
    
    let dmax: f32 = m.max64[0].0(m, m.range()) as f32;
    let data: Vec<f32> = m.aiff_data.iter().map(|x| *x as f32 / dmax).collect();

    let fft_fn = || {
        let size = m.fft_size[0];
        let width = size.min(m.window);
        let mut samples = vec![0.0; size];
        samples[0..width].copy_from_slice(&data[start..(start + width)]);
        m.window_fn[0](&mut samples);
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


    let fft: Vec<f32> = fft_fn();
    let fft: Vec<f32> = fft.iter().map(|x| x * x).collect();

    let bin_size = 1000.0 / fft.len() as f32;
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
    snr.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
    let signal = snr[0];
    let noise = snr[snr.len() / 2];

    // Normalize the fft
    let fft: Vec<f32> = fft.iter().map(|x| *x / fmax).collect();

    (fft, fmax, fpeak, (signal, noise))
}


pub fn old_spectrogram(m: &Model) -> (Vec<Vec<f32>>, Vec<f32>, Vec<f32>, Vec<(f32, f32)>, f32) {
    let mut spec: Vec<Vec<f32>> = Vec::new();
    let mut fmax: Vec<f32> = Vec::new();
    let mut fpeak: Vec<f32> = Vec::new();
    let mut fsnr: Vec<(f32, f32)> = Vec::new();
    let mut smax: f32 = 0.0;

    for start in (0..4000 - m.window).step_by(m.slide[0]) {
        let (fft, max, peak, snr) = old_fft(m, start);

        spec.push(fft);
        fmax.push(max);
        fpeak.push(peak);
        fsnr.push(snr);
        smax = smax.max(snr.0 / snr.1);
    }

    (spec, fmax, fpeak, fsnr, smax)
}


fn detections(fft: &[f32]) -> Vec<usize> {
    let bin_size = 1000.0 / fft.len() as f32;
    let whale_bins = (WHALE_RANGE.start / bin_size) as usize .. (WHALE_RANGE.end / bin_size) as usize;

    let mut peaks = Vec::new();

    let threshold = 0.5;

    let is_col = |i| fft[i] < fft[i-1] && fft[i] < fft[i+1];
    let is_peak = |i| fft[i] > fft[i-1] && fft[i] > fft[i+1];

    let valid_col = |i, peak| { fft[i] < threshold * fft[peak] };
    let valid_peak = |i,col| { fft[col] < threshold * fft[i] };
    
    let mut col = 0;
    let mut peak = 0;

    for i in 1 .. fft.len() - 1 {
        if is_peak(i) && valid_peak(i, col) { 
            if whale_bins.contains(&i) { peaks.push(i); }
            peak = i; 
        }
        if is_col(i) && valid_col(i, peak) {
            col = i;
        }
    }

    peaks
}


pub fn tracking(spec: &Vec<Vec<f32>>) -> Tracks {
    let bin_size = 1000.0 / spec[0].len() as f32;
    let max_bin = ((WHALE_RANGE.start + (WHALE_RANGE.end - WHALE_RANGE.start) / 2.0) / bin_size) as usize;

    let mut tracks: Tracks = Tracks::new();

    for (t, fft) in spec.iter().enumerate() {
        let detects = detections(fft);
        let unused = tracks.associate(t, &detects);
        tracks.prune(t);

        // new tracks must start below the mid point, anything above is noise
        let unused = unused.iter().filter(|d| **d < max_bin).cloned().collect();
        tracks.new_tracks(t, &unused);
    }

    tracks
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct Track(Vec<(usize, usize)>);

impl Track {
    fn new(time: usize, bin: usize) -> Track {
        Track(vec![(time, bin)])
    }

    pub fn age(&self, now: usize) -> usize {
        let times: Vec<usize> = self.0.iter().map(|(t, _)| *t).collect();
        let max = times.iter().max().unwrap();
        now - max
    }

    pub fn height(&self) -> usize {
        let bins: Vec<usize> = self.0.iter().map(|(_, b)| *b).collect();
        let min = bins.iter().min().unwrap();    
        let max = bins.iter().max().unwrap();
        max - min
    }

    pub fn width(&self) -> usize {
        let times: Vec<usize> = self.0.iter().map(|(t, _)| *t).collect();
        let min = times.iter().min().unwrap();    
        let max = times.iter().max().unwrap();
        max - min
    }

    pub fn associated(&mut self, time: usize, bin: usize) -> bool {

        let diff = 2;

        if self.0.iter().any(|(t, b)| t.abs_diff(time) < diff && b.abs_diff(bin) < diff) {
            self.0.push((time, bin));
            return true;
        }
        
        return false;
    }

    pub fn iter(&self) -> std::slice::Iter<(usize, usize)> {
        self.0.iter()
    }
}


//-------------- #[derive(New)]
#[derive(Debug, Default, PartialEq)]
pub struct Tracks {
    pub current: Vec<Track>,
    pub history: Vec<Track>,
}

impl Tracks {
    fn new() -> Tracks {
        Tracks { current: Vec::new(), history: Vec::new() }
    }

    pub fn new_tracks(&mut self, time: usize, detects: &Vec<usize>) {
        for d in detects {
            self.current.push(Track::new(time, *d));
        }
    }
    
    pub fn associate(&mut self, i: usize, detects: &Vec<usize>) -> Vec<usize> {

        let mut unused = Vec::new();

        for d in detects {
            let mut taken = false;
            for t in self.current.iter_mut() {
                if t.associated(i, *d) {
                    taken = true;
                }
            }
            if !taken {
                unused.push(*d);
            }
        }

        unused
    }

    pub fn prune(&mut self, now: usize) {
        let mut i = 0;
        while i < self.current.len() {
            if self.current[i].age(now) > 5 {
                let track = self.current.swap_remove(i);
                self.history.push(track);
                continue;
            } 
            i += 1;
        }
    }

    pub fn tallest(&self) -> Option<&Track> {
        self.current.iter()
            .chain(self.history.iter())
            .max_by_key(|t| t.height())
    }

    pub fn widest(&self) -> Option<&Track> {
        self.current.iter()
            .chain(self.history.iter())
            .max_by_key(|t| t.width())
    }

}


pub fn fft(m: &Model, start: usize) -> Vec<f32> {
    
    let dmax: f32 = m.max64[0].0(m, m.range()) as f32;
    let data: Vec<f32> = m.aiff_data.iter().map(|x| *x as f32 / dmax).collect();

    let fft_fn = || {
        let size = m.fft_size[0];
        let width = size.min(m.window);
        let mut samples = vec![0.0; size];
        samples[0..width].copy_from_slice(&data[start..(start + width)]);
        m.window_fn[0](&mut samples);
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


    fft_fn().iter().map(|x| x * x).collect()

}

pub fn spectrogram(m: &Model) -> Vec<Vec<f32>> {
    (0..4000 - m.window)
        .step_by(m.slide[0])
        .map(|start| 
            fft(m, start))
        .collect()
}
