use libm::sin;
use std::f32::consts::PI;

const SAMPLE_RATE: u32 = 2048;
const NUM_SAMPLES: usize = 1024;
const WINDOW_SIZE: usize = 256;
const BIN_SIZE: f32 = SAMPLE_RATE as f32 / NUM_SAMPLES as f32;


fn sin32(x: f32) -> f32 {
    sin(x as f64) as f32
}

fn main() {
    let mut samples: [f32; NUM_SAMPLES] = [0.0; NUM_SAMPLES];
    generate_sine_wave(&mut samples, 256.0, 0.0, 1.0);
    generate_sine_wave(&mut samples, 250.0, 0.0, 0.5);
    samples = normalize_array(samples);

    let spectrum = microfft::real::rfft_1024(&mut samples);
    spectrum[0].im = 0.0;

    let amplitudes: [f32; NUM_SAMPLES/2] 
        = spectrum
        .iter()
        .map(|c| c.norm())
        .collect::<Vec<_>>()
        .try_into()
        .unwrap();
    let amplitudes = normalize_array(amplitudes);

    let lower_bin = 128.0 - 10.0 * BIN_SIZE as f32;
    let upper_bin = 256.0 + 15.0 * BIN_SIZE as f32;

    for (i, magnitude) in amplitudes.iter().enumerate() {
        if i as f32 * BIN_SIZE > lower_bin && i as f32 * BIN_SIZE < upper_bin {
            let stars = stars(magnitude * 100.0);
            let star_string: String = stars.iter().collect();        
            println!("{}  {magnitude:.3}  {}", i as f32 * BIN_SIZE, star_string);
        }
    }
}


fn normalize_array<T>(mut data: T) -> T
where
    T: AsMut<[f32]> + AsRef<[f32]>,
{
    let data_ref = data.as_mut();
    let max_abs_value = data_ref.iter().map(|&x| x.abs()).fold(data_ref[0].abs(), f32::max);

    print!("max_abs_value: {}\n", max_abs_value);

    if max_abs_value != 0.0 {
        for element in data_ref.as_mut() {
            *element /= max_abs_value;
        }
    }

    data
}


fn generate_sine_wave(data: &mut [f32; 1024], freq_hz: f32, _phase_rad: f32, amp: f32) {

    let time_step = 1.0 / SAMPLE_RATE as f32;
    let angular_freq = 2.0*PI * freq_hz;

    for (i, value) in data.iter_mut().enumerate().take(WINDOW_SIZE) {
        let t: f32 = i as f32 * time_step;
        *value = *value + sin32(angular_freq * t) * amp;
    }

    for i in WINDOW_SIZE..data.len() {
        data[i] = 0.0;
    }
}


fn stars(length: f32) -> Vec<char> {
    let mut stars = Vec::with_capacity(length as usize);
    for _ in 0..length as usize {
        stars.push('*');
    }
    stars
}
