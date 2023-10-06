#![deny(clippy::all)]
#![forbid(unsafe_code)]

use std::f32::consts::PI;
use std::ops::Range;

use crate::{Model, WHALE_RANGE};


pub fn _hann(samples: &mut [f32]) {
    let len = samples.len();
    for (i, item) in samples.iter_mut().enumerate().take(len) {
        *item *= 0.5 * (1.0 - (2.0 * PI * i as f32 / (len - 1) as f32).cos());
    }
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


fn find_peaks(num: usize, fft: &[f32]) -> Vec<usize> {
    let bin_size = 2000.0 / fft.len() as f32;
    let whale_bins = (WHALE_RANGE.start / bin_size) as usize .. (WHALE_RANGE.end / bin_size) as usize;

    let mut peaks: Vec<usize> = whale_bins.filter(|i| {
        fft[*i] > fft[*i-1] && fft[*i] > fft[*i+1]
    }).collect();
    peaks.sort_by(|a, b| fft[*b].partial_cmp(&fft[*a]).unwrap());

    peaks[0..num.min(peaks.len())].to_vec()
}

const TRACK_MIN: usize = 10;
const NUM_PEAKS: usize = 5;

pub fn tracking(spec: &Vec<Vec<f32>>) -> Tracks {
    let mut tracks: Tracks = Tracks::new();

    for (i, fft) in spec.iter().enumerate() {
        let unused_peaks = tracks.extend(i, find_peaks(NUM_PEAKS, &fft));
        tracks.new_tracks(i, unused_peaks);
    }
    // filter out tracks that are too short
    tracks.0.iter().filter(|t| t.len() > TRACK_MIN).cloned().collect::<Tracks>()
}

#[derive(Debug, Clone, Default)]
pub struct Track(Vec<(usize, usize)>);

impl Track {
    fn new(bin: usize, peak: usize) -> Track {
        Track(vec![(bin, peak)])
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    fn extend(&mut self, bin: usize, peak: usize) -> bool {
        // if the peak is close to the last peak in the track, extend the track
        if let Some((last_bin, last_peak)) = self.0.last() {
            if bin - *last_bin < 2 && peak >= *last_peak && peak < *last_peak + 2 {
                self.0.push((bin, peak));
                return true
            }
        }
        false
    }
}


//-------------- #[derive(New)]
#[derive(Debug, Default)]
pub struct Tracks(Vec<Track>);

impl Tracks {
    fn new() -> Tracks {
        Tracks(Vec::new())
    }

    fn add(&mut self, track: Track) {
        self.0.push(track);
    }

    pub fn _len(&self) -> usize {
        self.0.len()
    }

    pub fn _iter(&self) -> std::slice::Iter<Track> {
        self.0.iter()
    }
    
    fn new_tracks(&mut self, bin: usize, peaks: Vec<usize>) {
        peaks.iter().for_each(|p| {
            self.add(Track::new(bin, *p));
        });
    }

    fn extend(&mut self, bin: usize, peaks: Vec<usize>) -> Vec<usize> {
        let mut unused_peaks = Vec::new();
        for peak in peaks {
            if !self.0.iter_mut().any(|t| t.extend(bin, peak)) {
                unused_peaks.push(peak);
            }
        }
        unused_peaks
    }

    pub fn as_points(&self) -> Vec<(usize, usize)> {
        self.0.iter().flat_map(|t| t.0.iter().copied()).collect()
    }

}

use std::iter::FromIterator;

impl FromIterator<Track> for Tracks {
    fn from_iter<I: IntoIterator<Item=Track>>(iter: I) -> Self {
        let mut tracks = Tracks::new();
        for track in iter {
            tracks.0.push(track);
        }
        tracks
    }
}