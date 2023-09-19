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
use winit::event_loop::EventLoop;
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

fn sine(samples: i32, freq: f64, offset: f64) ->  Vec<(i32, f64)> {
    (0..samples).map(move |x| {
        let phase: f64 = 2.*PI * (x as f64 * freq / samples as f64);
        (x, (phase + offset*2.*PI).sin())
    }).collect()
}

fn vec_add(a: &Vec<(i32, f64)>, b: &Vec<(i32, f64)>) -> Vec<(i32, f64)> {
    a.iter().zip(b.iter()).map(|(x,y)| (x.0, x.1 + y.1)).collect()
}

const WIDTH: u32 = 1200; //320;
const HEIGHT: u32 = 800; //240;

/// Representation of the application state. In this example, a box will bounce around the screen.
fn draw(frame: &mut [u8], m: &Model) {
    
    let s1: Vec<(i32, f64)> = sine(m.sample_count, m.s1_freq, m.s1_offset);
    let s2: Vec<(i32, f64)> = sine(m.sample_count, m.s2_freq, m.s2_offset);
    let s3: Vec<(i32, f64)> = vec_add(&s1, &s2);
    let s4: Vec<f32> = s3.iter().map(|x| x.1 as f32).collect();

    let mut samples: [f32; 1024] = [0.0; 1024];
    samples[0..200].copy_from_slice(&s4[..200]);

    let spectrum = microfft::real::rfft_1024(&mut samples);
    spectrum[0].im = 0.0;

    let amplitudes: Vec<f64> 
            = spectrum.iter()
                .map(|c| c.norm() as f64)
                .collect();

    let max_amp: f64 = amplitudes.iter().fold(0.0, |a, &b| a.max(b));

    let root_drawing_area 
        = BitMapBackend::<BGRXPixel>::with_buffer_and_format(frame, (WIDTH, HEIGHT))
            .unwrap().into_drawing_area();
    root_drawing_area.fill(&WHITE).unwrap();

    let (upper, lower) = root_drawing_area.split_vertically(240);

    let upper = upper.margin(0, 0, 0, 20);
    let lower = lower.margin(0, 0, 0, 20);

    let upper_caption = format!("{} Hz  {} Hz", m.s1_freq*10.0, m.s2_freq*10.0);

    let mut upper_chart = ChartBuilder::on(&upper)
        .caption(upper_caption, ("sans-serif", 40))
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .build_cartesian_2d(0..200, -2.1..2.1)
        .unwrap();

    let lower_captiion = format!("Bin {:.2} Hz  -  Max {:.0}", m.sample_count as f64 / 1024.0 * 10.0, max_amp);

    let mut lower_chart = ChartBuilder::on(&lower)
        .caption(lower_captiion, ("sans-serif", 40))
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .build_cartesian_2d((0..512).into_segmented(), 0..(max_amp*1.1) as i32)
        .unwrap();

    upper_chart.configure_mesh().draw().unwrap();
    lower_chart.configure_mesh().draw().unwrap();

    upper_chart.draw_series( LineSeries::new( s1, RED.stroke_width(1)) ).unwrap();
    upper_chart.draw_series( LineSeries::new( s2, BLUE.stroke_width(2)) ).unwrap();
    upper_chart.draw_series( LineSeries::new( s3, BLACK.stroke_width(4)) ).unwrap();

    lower_chart.draw_series(amplitudes.iter().enumerate().map(|(x, y)| {
        let x0 = SegmentValue::Exact(x as i32);
        let x1 = SegmentValue::Exact(x as i32 + 1);
        let y0: i32 = *y as i32;
        let bar = Rectangle::new([(x0, 0), (x1, y0)], BLUE.filled());
        // bar.set_margin(0, 0, 5, 5);
        bar
    }))
    .unwrap();




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
    
    let mut model = Model::new(200);

    event_loop.run(move |event, _, control_flow| {
        control_flow.set_wait();

        // Draw the current frame
        if let Event::RedrawRequested(_) = event {
            draw(pixels.frame_mut(), &model);
            if let Err(err) = pixels.render() {
                log_error("pixels.render", err);
                control_flow.set_exit();
                return;
            }
        }

        // Handle input events
        if input.update(&event) {
            if input.key_pressed(VirtualKeyCode::Escape) || input.close_requested() {
                control_flow.set_exit();
                return;
            }

            if input.key_pressed(VirtualKeyCode::Equals) {
                model.s2_freq += 1.0;
            }

            if let Some(size) = input.window_resized() {
                if let Err(err) = pixels.resize_surface(size.width, size.height) {
                    log_error("pixels.resize_surface", err);
                    control_flow.set_exit();
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

struct Model {

    sample_count: i32,

    s1_freq: f64,
    s2_freq: f64,

    s1_offset: f64,
    s2_offset: f64,
}

impl Model {
    
    fn new(sample_count: i32) -> Self {
        Self {
            sample_count,
            s1_freq: 10.0,
            s2_freq: 30.0,
            s1_offset: 0.0,
            s2_offset: 0.0,
        }
    }
}