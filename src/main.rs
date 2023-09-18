#![deny(clippy::all)]
#![deny(clippy::pedantic)]
#![forbid(unsafe_code)]

// use env_logger::fmt::ColorSpec;
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

// use libm::sin;
use std::f64::consts::PI;


fn _db(value1: f32, value2: f32) -> f32 {
    // Calculate the dB difference using the formula: dB = 10 * log10(value1 / value2)
    if value1 == 0.0 || value2 == 0.0 {
        // Avoid division by zero and return negative infinity in such cases
        f32::NEG_INFINITY
    } else {
        10.0 * (value1 / value2).log10()
    }
}


fn _normalize_array<T>(data: &mut T) -> (f32, usize)
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
        for element in &mut *data.as_mut() {
            *element /= max_abs_value;
        }
    }

    (max_abs_value, max_index)
}


fn _find_median(arr: &mut [f32]) -> f32 {
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

fn sine(samples: i32, freq: f64, phase: f64) -> impl Iterator<Item = (i32, f64)> {
    (0..samples)
    .map(move |x| {
        let samples = f64::from(samples);
        let phase: f64 = (x as f64) + phase * samples / freq;
        let radians = 2.0 * PI * phase / samples;
        (x, (freq * radians).sin())
    })
}


const WIDTH: u32 = 1200; //320;
const HEIGHT: u32 = 800; //240;

/// Representation of the application state. In this example, a box will bounce around the screen.
fn draw(frame: &mut [u8]) {

    let root_drawing_area 
        = BitMapBackend::<BGRXPixel>::with_buffer_and_format(frame, (WIDTH, HEIGHT))
            .unwrap().into_drawing_area();
    root_drawing_area.fill(&WHITE).unwrap();

    let (upper, lower) = root_drawing_area.split_vertically(240);
    let upper = upper.margin(0, 0, 0, 20);
    let lower = lower.margin(0, 0, 0, 20);

    let mut upper_chart = ChartBuilder::on(&upper)
        .caption("Sin Waves", ("sans-serif", 40))
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .build_cartesian_2d(0..200, -2.1..2.1)
        .unwrap();

    let mut lower_chart = ChartBuilder::on(&lower)
        .caption("FFT", ("sans-serif", 40))
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .build_cartesian_2d(0..512, 0.0..100.0)
        .unwrap();

    upper_chart.configure_mesh().draw().unwrap();
    lower_chart.configure_mesh().draw().unwrap();

    const FREQ: f64 = 90.0;

    let s1 = sine(200, FREQ, 0.0);
    // let s2 = sine(200, 4.9, 0.1);

    // let s1a = sine(200, FREQ, 0.0);
    // let s2a = sine(200, 4.9, 0.1);
    // let s3 = s1a.zip(s2a).map(|(x, y)| (x.0, x.1 + y.1));

    upper_chart.draw_series( LineSeries::new( s1, BLACK.stroke_width(3)) ).unwrap();
    // upper_chart.draw_series( LineSeries::new( s2, GREEN.stroke_width(3)) ).unwrap();
    // upper_chart.draw_series( LineSeries::new( s3, BLUE.stroke_width(3)) ).unwrap();

    let s4: Vec<f32> = sine(200, FREQ, 0.0).map(|(_x,y)| y as f32).collect();
    let s5: Vec<f32> = sine(200, FREQ + 1.0, 0.0).map(|(_x,y)| y as f32).collect();
    let s6: Vec<f32> = s4.iter().zip(s5.iter()).map(|(x,y)| x + y).collect();

    let mut samples: [f32; 1024] = [0.0; 1024];
    samples[0..200].copy_from_slice(&s6[..200]);

    let spectrum = microfft::real::rfft_1024(&mut samples);
    spectrum[0].im = 0.0;

    let amplitudes: Vec<_> = spectrum.iter().map(|c| c.norm() as i32).collect();
    let amplitudes: Vec<(i32, f64)> = amplitudes.iter().enumerate().map(|(x, y)| (x as i32, *y as f64)).collect();
    // let max: f64 = amplitudes.iter().fold(0.0, |a, &b| a.max(b.1));
    
    lower_chart.draw_series( LineSeries::new( amplitudes, BLACK.stroke_width(3)) ).unwrap();

    root_drawing_area.present().unwrap();

}

fn main() -> Result<(), Error> {
    env_logger::init();

    let event_loop = EventLoop::new();
    let mut input = WinitInputHelper::new();
    let window = {
        let size = LogicalSize::new(WIDTH, HEIGHT);
        WindowBuilder::new()
        .with_title("Finding Moby")
        .with_inner_size(size)
        .with_min_inner_size(size)
        .build(&event_loop)
        .unwrap()
    };
    
    let mut pixels = {
        let window_size = window.inner_size();
        let surface_texture = SurfaceTexture::new(window_size.width, window_size.height, &window);
        Pixels::new(WIDTH, HEIGHT, surface_texture)?
    };
    
    event_loop.run(move |event, _, control_flow| {
        *control_flow = ControlFlow::Wait;

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

