#![deny(clippy::all)]
#![forbid(unsafe_code)]

mod draw;
mod signal;

use core::ops::Range;
use std::fs::File;
use std::io::{self, Read};

use byteorder::{BigEndian, ByteOrder, ReadBytesExt};
use draw::draw_upper_signal;
use pixels::{Error, Pixels, SurfaceTexture};
use winit::dpi::LogicalSize;
use winit::event::{ElementState, Event, KeyboardInput, VirtualKeyCode, WindowEvent};
use winit::event_loop::EventLoop;
use winit::window::WindowBuilder;

use crate::draw::*;
use crate::signal::*;

const WIDTH: u32 = 1200;
const HEIGHT: u32 = 800;

/// Representation of the application state. In this example, a box will bounce around the screen.
fn main() -> Result<(), Error> {
    let event_loop = EventLoop::new();

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

    let mut m = Model::new(74);

    event_loop.run(move |event, _, control_flow| {
        control_flow.set_wait();

        // println!("Event: {:?}", event);

        match event {
            Event::WindowEvent {
                event: WindowEvent::CloseRequested,
                ..
            } => {
                println!("The close button was pressed; stopping");
                control_flow.set_exit();
            }

            Event::RedrawRequested(_) => {
                draw::draw(pixels.frame_mut(), &m).unwrap();
                if let Err(err) = pixels.render() {
                    println!("Error: {err:?}");
                    control_flow.set_exit();
                }
            }

            Event::MainEventsCleared => {
                // window.request_redraw();
            }

            Event::WindowEvent {
                event: WindowEvent::ModifiersChanged(mods),
                ..
            } => {
                m.modifiers = mods;
                println!("{:?}", mods);
            }

            Event::WindowEvent {
                event:
                    WindowEvent::KeyboardInput {
                        input:
                            KeyboardInput {
                                virtual_keycode: Some(key),
                                state: ElementState::Pressed,
                                ..
                            },
                        ..
                    },
                ..
            } => {
                use VirtualKeyCode::*;
                match key {
                    Up => m.next_file(),
                    Down => m.prev_file(),
                    Left => m.prev_window(),
                    Right => m.next_window(),
                    F => m.fft_size.rotate_left(1),
                    N => m.max64.rotate_left(1),
                    P => play_aiff(&m),
                    Q => control_flow.set_exit(),
                    S => m.slide.rotate_left(1),
                    U => m.upper.rotate_left(1),
                    _ => println!("Key pressed {key:?}"),
                }
                window.request_redraw();
            }
            _ => (),
        }
    })
}

// Takes a integer and returns the file name
fn aiff_name(file_number: usize) -> (String, String) {
    let file_name = format!("train{file_number}.aiff");
    let dir = file_number / 100;
    let full_name = format!("/Users/steveweeks/dev/whales/data/train/{dir}/{file_name}");

    (file_name, full_name)
}

fn aiff_samples(full_name: &str) -> Result<Vec<i16>, io::Error> {
    let mut file = File::open(full_name)?;

    // Read and skip the AIFF header
    let mut header_buf = [0u8; 54];
    file.read_exact(&mut header_buf)?;

    // Parse the AIFF header to determine audio data size
    let num_samples = BigEndian::read_u32(&header_buf[22..26]).min(4000);

    // Read and convert the audio samples
    let mut samples = Vec::with_capacity(num_samples as usize);
    for _ in 0..num_samples {
        let sample = file.read_i16::<BigEndian>()?;
        samples.push(sample);
    }

    Ok(samples)
}

pub struct Model<'a> {
    aiff_number: usize,
    aiff_name: String,
    full_name: String,
    aiff_data: Vec<i16>,
    fft_size: Vec<usize>,
    start: usize,
    window: usize,
    slide: Vec<usize>,
    modifiers: winit::event::ModifiersState,
    max64: Vec<(fn(&Model, Range<usize>) -> f64, &'a str)>,
    upper: Vec<fn(draw::MyDrawingArea, &Model) -> Result<(), Box<dyn std::error::Error>>>,
}

impl Model<'_> {
    fn new(aiff_number: usize) -> Self {
        let (aiff_name, aiff_path) = aiff_name(aiff_number);
        let aiff_data = aiff_samples(&aiff_path).unwrap();

        Self {
            aiff_number,
            aiff_name,
            full_name: aiff_path,
            aiff_data,
            fft_size: vec![1024, 512, 256, 128],
            start: 2000,
            window: 200,
            slide: vec![10, 50, 100, 200],
            modifiers: winit::event::ModifiersState::default(),
            max64: vec![
                (max_window, "nwindow"),
                (max_none, "nnone"),
                (max_all, "nall"),
            ],
            upper: vec![draw_upper_fft, draw_upper_signal],
        }
    }

    fn set_aiff(&mut self, aiff_number: usize) {
        if aiff_number > 30_000 {
            return;
        }
        if aiff_number < 1 {
            return;
        }

        self.aiff_number = aiff_number;
        (self.aiff_name, self.full_name) = aiff_name(self.aiff_number);
        self.aiff_data = aiff_samples(&self.full_name).unwrap();
    }

    fn next_file(&mut self) {
        self.set_aiff(self.aiff_number + 1);
    }

    fn prev_file(&mut self) {
        self.set_aiff(self.aiff_number - 1);
    }

    fn range(&self) -> std::ops::Range<usize> {
        self.start..self.start + self.window
    }

    fn next_window(&mut self) {
        use winit::event::*;
        self.start = (4000 - self.window).min(
            self.start
                + match self.modifiers {
                    ModifiersState::SHIFT => 10,
                    _ => 100.max(self.slide[0]),
                },
        );
    }

    fn prev_window(&mut self) {
        use winit::event::*;
        self.start = self.start.saturating_sub(match self.modifiers {
            ModifiersState::SHIFT => 10,
            _ => 100.max(self.slide[0]),
        });
    }
}

fn play_aiff(m: &Model) {
    use rodio::buffer::SamplesBuffer;
    use rodio::{OutputStream, Sink};

    let max = 0.7 * m.max64[0].0(m, m.range()) as f32;
    let signal: Vec<f32> = m.aiff_data.iter().map(|x| *x as f32 / max).collect();

    let (_stream, stream_handle) = OutputStream::try_default().unwrap();
    let sink = Sink::try_new(&stream_handle).unwrap();

    let source = SamplesBuffer::new(1, 2000, signal);
    sink.append(source);
    sink.sleep_until_end();
}
