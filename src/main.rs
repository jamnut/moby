#![deny(clippy::all)]
#![forbid(unsafe_code)]

mod draw;
mod signal;

use core::ops::Range;
use std::fs::File;
use std::io::{self, Read};

use byteorder::{BigEndian, ByteOrder, ReadBytesExt};
use draw::draw_signal;
use pixels::{Error, Pixels, SurfaceTexture};
use winit::dpi::LogicalSize;
use winit::event::{ElementState, Event, KeyboardInput, VirtualKeyCode, WindowEvent};
use winit::event_loop::EventLoop;
use winit::window::WindowBuilder;

use crate::draw::*;
use crate::signal::*;

const WHALE_RANGE: Range<f32> = 175.0 .. 600f32;
const WHALE_VIEW: Range<f32> = 100.0 .. 700f32;


/// Representation of the application state. In this example, a box will bounce around the screen.
fn main() -> Result<(), Error> {
    
    let event_loop = EventLoop::new();
    
    let window = {
        WindowBuilder::new()
        .with_maximized(true)
        .with_title("Finding Moby")
        .build(&event_loop)
        .unwrap()
    };
    
    let psize= window.inner_size();
    let lsize: LogicalSize<u32> = psize.to_logical(window.scale_factor());

    let mut pixels = {
        let surface_texture = SurfaceTexture::new(psize.width, psize.height, &window);
        Pixels::new(lsize.width, lsize.height, surface_texture)?
    };
    
    let mut m = Model::new(5);

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
                draw::draw((lsize.width, lsize.height), pixels.frame_mut(), &m).unwrap();
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
                    W => m.aiff_type.rotate_left(1),
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

#[derive(Debug, Clone, Copy)]
enum AiffType {
    All,
    Whales,
}

type DrawingFn = fn(&draw::Drawing, &Model) -> Result<(), Box<dyn std::error::Error>>;

pub struct Model<'a> {
    is_whale: Vec<bool>,
    aiff_type: Vec<AiffType>,
    aiff_number: usize,
    aiff_name: String,
    full_name: String,
    aiff_data: Vec<i16>,
    fft_size: Vec<usize>,
    start: usize,
    window: usize,
    slide: Vec<usize>,
    modifiers: winit::event::ModifiersState,
    max64: Vec<(fn(&Model, Range<usize>) -> f32, &'a str)>,
    upper: Vec<DrawingFn>,
}

impl Model<'_> {
    fn new(aiff_number: usize) -> Self {
        let (aiff_name, aiff_path) = aiff_name(aiff_number);
        let aiff_data = aiff_samples(&aiff_path).unwrap();

        Self {
            is_whale: load_train_csv(),
            aiff_type: vec![AiffType::Whales, AiffType::All],
            aiff_number,
            aiff_name,
            full_name: aiff_path,
            aiff_data,
            fft_size: vec![1024, 512, 256, 128, 2048],
            start: 2000,
            window: 200,
            slide: vec![10, 50, 100, 200],
            modifiers: winit::event::ModifiersState::default(),
            max64: vec![
                (max_window, "nwindow"),
                (max_none, "nnone"),
                (max_all, "nall"),
            ],
            upper: vec![draw_tracking, draw_signal],
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
        match self.aiff_type[0] {
            AiffType::All => {
                self.set_aiff(self.aiff_number + 1)
            },
            AiffType::Whales => {
                for i in (self.aiff_number + 1)..30_000 {
                    if self.is_whale[i] {
                        self.set_aiff(i);
                        break;
                    }
                }
            },
        };
    }

    fn prev_file(&mut self) {
        match self.aiff_type[0] {
            AiffType::All => {
                self.set_aiff(self.aiff_number - 1)
            },
            AiffType::Whales => {
                for i in (2..self.aiff_number).rev() {
                    if self.is_whale[i] {
                        self.set_aiff(i);
                        break;
                    }
                }
            },
        };
    }

    fn range(&self) -> std::ops::Range<usize> {
        self.start..self.start + self.window
    }

    fn next_window(&mut self) {
        use winit::event::*;
        self.start = (4000 - self.window).min(
            self.start
                + match self.modifiers {
                    ModifiersState::SHIFT => self.slide[0].max(10),
                    _ => 100.max(self.slide[0]),
                },
        );
    }

    fn prev_window(&mut self) {
        use winit::event::*;
        self.start = self.start.saturating_sub(match self.modifiers {
            ModifiersState::SHIFT => self.slide[0].max(10),
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

fn load_train_csv() -> Vec<bool>{
    use std::io::{BufRead, BufReader};

    let mut is_whale = vec![false; 30_000];
    let file = File::open("/Users/steveweeks/dev/whales/data/train.csv").unwrap();
    let reader = BufReader::new(file);

    for line in reader.lines().skip(1) {
        if let Ok(line) = line {
            let mut iter = line.split(',');
            let file_number = match iter.next().unwrap().chars()
                .skip_while(|c| !c.is_digit(10))
                .take_while(|c| c.is_digit(10))
                .collect::<String>()
                .parse::<usize>() {
                Ok(file_number) => file_number,
                Err(err) => {
                    println!("Error parsing file number {err:?}  {line}");
                    continue;
                }
            };
            if iter.next().unwrap() == "1" {
                is_whale[file_number] = true;
            }
        } 
    }
    is_whale
}

