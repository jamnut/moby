#![deny(clippy::all)]
#![forbid(unsafe_code)]

use error_iter::ErrorIter as _;
use log::error;
use pixels::{Error, Pixels, SurfaceTexture};
use plotters::prelude::*;
use plotters::backend::BGRXPixel;
use winit::dpi::LogicalSize;
use winit::event::{Event, VirtualKeyCode};
use winit::event_loop::{ControlFlow, EventLoop};
use winit::window::WindowBuilder;
use winit_input_helper::WinitInputHelper;

use libm::sin;
use std::f32::consts::PI;

const SAMPLE_RATE: u32 = 2048;
const NUM_SAMPLES: usize = 1024;
const WINDOW_SIZE: usize = 256;
const BIN_SIZE: f32 = SAMPLE_RATE as f32 / NUM_SAMPLES as f32;


fn sin32(x: f32) -> f32 {
    sin(x as f64) as f32
}

fn old_main() {
    let mut samples: [f32; NUM_SAMPLES] = [0.0; NUM_SAMPLES];
    generate_sine_wave(&mut samples, 256.0, 0.0, 1.0);
    // generate_sine_wave(&mut samples, 250.0, 0.0, 1.0);
    let (max_signal_value, max_signl_index) = normalize_array(&mut samples);

    let spectrum = microfft::real::rfft_1024(&mut samples);
    spectrum[0].im = 0.0;

    let mut amplitudes: [f32; NUM_SAMPLES/2] 
        = spectrum
        .iter()
        .map(|c| c.norm())
        .collect::<Vec<_>>()
        .try_into()
        .unwrap();
    let (max_amplitude, max_amplitude_index) = normalize_array(&mut amplitudes);

    let lower_bin = (max_amplitude_index - 15) as f32 * BIN_SIZE;
    let upper_bin = (max_amplitude_index + 15) as f32 * BIN_SIZE;

    for (i, magnitude) in amplitudes.iter().enumerate() {
        if i as f32 * BIN_SIZE > lower_bin && i as f32 * BIN_SIZE < upper_bin {
            let stars = stars(magnitude * 100.0);
            let star_string: String = stars.iter().collect();        
            println!("{}  {magnitude:.3}  {}", i as f32 * BIN_SIZE, star_string);
        }
    }

    let median_amplitude = find_median(&mut amplitudes);
    let db_median = db(max_amplitude, median_amplitude);

    println!("max signal value: {:.3} at index {}", max_signal_value, max_signl_index);
    println!("max amplitude: {:.3} at {}Hz", max_amplitude, max_amplitude_index as f32 * BIN_SIZE);
    println!("median amplitude: {:.3}  {:.3} dB", median_amplitude, db_median);

    let center_index = 256/BIN_SIZE as usize;
    let slice_start = center_index - 4;
    let slice_end = center_index + 5;
    let range = (slice_end - slice_start) as f32 * BIN_SIZE;

    assert!(slice_end < 512);

    let slice = &amplitudes[slice_start..=slice_end];

    let center = sinc_center_index(slice, range) + slice_start as f32;
    println!("center: {:.1}", center*BIN_SIZE);

}

fn db(value1: f32, value2: f32) -> f32 {
    // Calculate the dB difference using the formula: dB = 10 * log10(value1 / value2)
    if value1 == 0.0 || value2 == 0.0 {
        // Avoid division by zero and return negative infinity in such cases
        f32::NEG_INFINITY
    } else {
        10.0 * (value1 / value2).log10()
    }
}


fn normalize_array<T>(data: &mut T) -> (f32, usize)
where
    T: AsMut<[f32]> + AsRef<[f32]>,
{
    let data_ref = data.as_mut();

    let mut max_abs_value = 0.0;
    let mut max_index = 0;

    // Find the max absolute value and its index
    for (index, &element) in data_ref.iter().enumerate() {
        let abs_value = element.abs();
        if abs_value > max_abs_value {
            max_abs_value = abs_value;
            max_index = index;
        }
    }

    // Normalize the data
    if max_abs_value > 0.0 {
        for element in data.as_mut().iter_mut() {
            *element /= max_abs_value;
        }
    }

    (max_abs_value, max_index)
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

fn sinc_center_index(data: &[f32], frequency_range: f32) -> f32 {
    // Define the sinc function
    let sinc = |x: f32| {
        if x == 0.0 {
            1.0
        } else {
            (std::f32::consts::PI * x).sin() / (std::f32::consts::PI * x)
        }
    };

    let data_len = data.len() as isize;
    let center_index = (data_len / 2) as f32; // Initial guess for the center index

    // Calculate the width of the sinc function based on the frequency range
    let sinc_width = frequency_range / (1.0 / data_len as f32);

    // Find the fractional index with the highest convolution result
    let (best_center_index, max_convolution_result) = (0..(data_len*10))
        .map(|i| {
            let offset = (i as f32)/10.0 - center_index;
            let convolution_result = (data[(i/10) as usize] * sinc(offset / sinc_width)) as f32;
            (i as f32/10.0, convolution_result)
        })
        .max_by(|&(_, a), &(_, b)| a.partial_cmp(&b).unwrap_or(std::cmp::Ordering::Equal))
        .unwrap_or((center_index, 0.0));

    if max_convolution_result > 0.0 {
        best_center_index
    } else {
        center_index // If no better center is found, return the initial guess
    }
}

fn _binary_sinc_center_index(data: &[f32], frequency_range: f32) -> f32 {
    // Define the sinc function
    let sinc = |x: f32| {
        if x == 0.0 {
            1.0
        } else {
            (std::f32::consts::PI * x).sin() / (std::f32::consts::PI * x)
        }
    };

    let data_len = data.len() as f32;
    let center_index = (data_len / 2.0) as isize; // Initial guess for the center index

    // Calculate the width of the sinc function based on the frequency range
    let sinc_width = frequency_range / (1.0 / data_len);

    let mut left = center_index;
    let mut right = center_index;
    let mut best_center_index = center_index as f32;

    while left >= 0 || right < data_len as isize {
        // Calculate convolution results for both sides
        let left_offset = left as f32 - center_index as f32;
        let right_offset = right as f32 - center_index as f32;
        let left_result = data[left as usize] * sinc(left_offset / sinc_width);
        let right_result = data[right as usize] * sinc(right_offset / sinc_width);

        // Check if left or right side has a higher convolution result
        if left_result >= right_result && left >= 0 {
            best_center_index = left_offset;
            left -= 1;
        } else if right < data_len as isize {
            best_center_index = right_offset;
            right += 1;
        }
    }

    best_center_index
}

fn find_median(arr: &mut [f32]) -> f32 {
    arr.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let len = arr.len();

    if len % 2 == 0 {
        // If the array has an even number of elements, return the average of the two middle values
        let mid1 = arr[len / 2 - 1];
        let mid2 = arr[len / 2];
        (mid1 + mid2) / 2.0
    } else {
        // If the array has an odd number of elements, return the middle value
        arr[len / 2]
    }
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

const WIDTH: usize = 800; //320;
const HEIGHT: usize = 512; //240;
const BOX_SIZE: i16 = 10; //64;

/// Representation of the application state. In this example, a box will bounce around the screen.
fn draw(frame: &mut [u8]) {

    let root_drawing_area 
        = BitMapBackend::<BGRXPixel>::with_buffer_and_format(frame, (WIDTH as u32, HEIGHT as u32))
        .unwrap().into_drawing_area();

    root_drawing_area.fill(&WHITE).unwrap();

    let mut chart = ChartBuilder::on(&root_drawing_area)
        .build_cartesian_2d(-3.14..3.14, -1.2..1.2)
        .unwrap();

    chart.draw_series(LineSeries::new(
        (-314..314).map(|x| x as f64 / 100.0).map(|x| (x, x.sin())),
        &RED
    )).unwrap();

    root_drawing_area.present().unwrap();

}

fn main() -> Result<(), Error> {
    env_logger::init();

    let event_loop = EventLoop::new();
    let mut input = WinitInputHelper::new();
    let window = {
        let size = LogicalSize::new(WIDTH as f64, HEIGHT as f64);
        WindowBuilder::new()
        .with_title("BlueIQ")
        .with_inner_size(size)
        .with_min_inner_size(size)
        .build(&event_loop)
        .unwrap()
    };
    
    let mut pixels = {
        let window_size = window.inner_size();
        let surface_texture = SurfaceTexture::new(window_size.width, window_size.height, &window);
        Pixels::new(WIDTH as u32, HEIGHT as u32, surface_texture)?
    };
    
    event_loop.run(move |event, _, control_flow| {
        // Draw the current frame
        if let Event::RedrawRequested(_) = event {
            draw(pixels.frame_mut());
            if let Err(err) = pixels.render() {
                log_error("pixels.render", err);
                *control_flow = ControlFlow::Exit;
                return;
            }
        }

        // Handle input events
        if input.update(&event) {
            if input.key_pressed(VirtualKeyCode::Escape) || input.close_requested() {
                *control_flow = ControlFlow::Exit;
                return;
            }

            if let Some(size) = input.window_resized() {
                if let Err(err) = pixels.resize_surface(size.width, size.height) {
                    log_error("pixels.resize_surface", err);
                    *control_flow = ControlFlow::Exit;
                    return;
                }
            }

            // Update internal state and request a redraw
            window.request_redraw();
        }
    });
}

fn log_error<E: std::error::Error + 'static>(method_name: &str, err: E) {
    error!("{method_name}() failed: {err}");
    for source in err.sources().skip(1) {
        error!("  Caused by: {source}");
    }
}

