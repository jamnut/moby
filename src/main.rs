#![deny(clippy::all)]
#![forbid(unsafe_code)]

mod signal;
mod draw;

use std::fs::File;
use std::io::{self, Read};


use byteorder::{BigEndian, ByteOrder, ReadBytesExt};
use pixels::{Error, Pixels, SurfaceTexture};
use winit::dpi::LogicalSize;
use winit::event::{ElementState, Event, KeyboardInput, VirtualKeyCode, WindowEvent};
use winit::event_loop::EventLoop;
use winit::window::WindowBuilder;

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

    let mut m = Model::new(28);

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
                draw::draw(pixels.frame_mut(), &m);
                if let Err(err) = pixels.render() {
                    println!("Error: {err:?}");
                    control_flow.set_exit();
                }
            }

            Event::MainEventsCleared => {
                // window.request_redraw();
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
                match key {
                    VirtualKeyCode::Up => m.next(),
                    VirtualKeyCode::Down => m.prev(),
                    VirtualKeyCode::Left => m.window_start = m.window_start.saturating_sub(100),
                    VirtualKeyCode::Right => m.window_start = 3800.min(m.window_start + 100),
                    VirtualKeyCode::Q => control_flow.set_exit(),
                    _ => println!("Key pressed {key:?}"),
                };
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
    let num_samples = BigEndian::read_u32(&header_buf[22..26]);
    
    // Calculate the number of samples to read (16-bit samples)
    let num_samples_to_read = std::cmp::min(num_samples, 4000);
    
    // Read and convert the audio samples
    let mut samples = Vec::with_capacity(num_samples_to_read as usize);
    for _ in 0..num_samples_to_read {
        let sample = file.read_i16::<BigEndian>()?;
        samples.push(sample);
    }
    
    Ok(samples)
}

pub struct Model {
    aiff_number: usize,
    aiff_name: String,
    aiff_path: String,
    aiff_data: Vec<i16>,
    window_start: usize,
}

impl Model {
    fn new(aiff_number: usize) -> Self {
        let (aiff_name, aiff_path) = aiff_name(aiff_number);
        let aiff_data = aiff_samples(&aiff_path).unwrap();

        Self { aiff_number, aiff_name, aiff_path, aiff_data, window_start: 0 }
    }

    fn set_aiff(&mut self, aiff_number: usize) {
        if aiff_number > 30_000 { return; }
        if aiff_number < 1 { return; }

        self.aiff_number = aiff_number;
        (self.aiff_name, self.aiff_path) = aiff_name(self.aiff_number);
        self.aiff_data = aiff_samples(&self.aiff_path).unwrap();
    }

    fn next(&mut self) {
        self.set_aiff(self.aiff_number + 1);
    }

    fn prev(&mut self) {
        self.set_aiff(self.aiff_number - 1);
    }

}
