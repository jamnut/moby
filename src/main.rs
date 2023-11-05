#![deny(clippy::all)]
#![forbid(unsafe_code)]

mod draw;
mod signal;

use core::ops::Range;
use std::fs::File;
use std::io::{self, Read};
use std::vec;

use byteorder::{BigEndian, ByteOrder, ReadBytesExt};
use draw::draw_signal;
use pixels::{Error, Pixels, SurfaceTexture};
use winit::dpi::LogicalSize;
use winit::event::*;
use winit::event_loop::EventLoop;
use winit::window::WindowBuilder;

use crate::draw::*;
use crate::signal::*;

const WHALE_RANGE: Range<f32> = 75.0 .. 375f32;
const WHALE_VIEW: Range<f32> = 00.0 .. 500f32;


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
                    Up => m.up(),
                    Down => m.down(),
                    Left => m.prev_window(),
                    Right => m.next_window(),
                    C => chirp(&mut m),
                    D => rotate(&mut m.fft_scale, m.modifiers),
                    F => rotate(&mut m.fft_size, m.modifiers),
                    H => rotate(&mut m.window_fn, m.modifiers),
                    P => play_aiff(&m),
                    Q => control_flow.set_exit(),
                    R => record(&mut m),
                    S => rotate(&mut m.slide, m.modifiers),
                    U => rotate(&mut m.upper, m.modifiers),
                    W => rotate(&mut m.aiff_type, m.modifiers),
                    Key1 => m.digits(1),
                    Key2 => m.digits(2),
                    Key3 => m.digits(3),
                    Key4 => m.digits(4),
                    Key5 => m.digits(5),
                    Key6 => m.digits(6),
                    Key7 => m.digits(7),
                    Key8 => m.digits(8),
                    Key9 => m.digits(9),
                    Key0 => m.digits(0),
                    Escape => m.digits = String::new(),
                    Return => {
                        m.set_aiff(m.digits.parse().unwrap_or(1));
                        m.digits = String::new();
                    },
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
    AllFiles,
    WhalesOnly,
}

type DrawingFn = fn(&draw::Drawing, &Model) -> Result<(), Box<dyn std::error::Error>>;

pub struct Model {
    is_whale: Vec<bool>,
    digits: String,
    aiff_type: Vec<AiffType>,
    aiff_number: usize,
    aiff_name: String,
    full_name: String,
    aiff_raw: Vec<i16>,
    aiff_max: i16,
    aiff_data: Vec<f32>,
    spectrogram: Vec<Vec<f32>>,
    noise: Vec<(f32, f32)>,
    fft_size: Vec<usize>,
    fft_scale: Vec<(f32, f32, fn (f32) -> f32)>,
    freq: f32, // horizonl cursur on spectrogram
    start: usize,  // vertical cursor on spectrogram
    window: usize,
    slide: Vec<usize>,
    window_fn: Vec<fn(&mut [f32])>,
    modifiers: winit::event::ModifiersState,
    upper: Vec<DrawingFn>,
}

impl Model {
    fn new(aiff_number: usize) -> Self {

        let mut new = Self {
            is_whale: load_train_csv(),
            digits: String::new(),
            aiff_type: vec![AiffType::WhalesOnly, AiffType::AllFiles],
            aiff_number,
            aiff_name: String::new(),
            full_name: String::new(),
            aiff_raw: vec![0; 0],
            aiff_max: 0,
            aiff_data: vec![0.0; 0],
            spectrogram: vec![], // vec![vec![0.0; 0]; 0],
            noise: vec![], // vec![(0.0, 0.0); 0],
            fft_size: vec![256, 128, 2048, 1024, 512],
            fft_scale: vec![(0.0, 1.0, |x| x), (-90.0, 0.0, |x| f32::log10(x) * 20.0)],
            freq: 200.0,
            start: 2000,
            window: 200,
            slide: vec![20, 40, 80, 140, 200],
            window_fn: vec![hann, no_windowing],
            modifiers: winit::event::ModifiersState::default(),
            upper: vec![draw_noise, draw_tracking, draw_signal, draw_raw],
        };

        new.set_aiff(aiff_number);
        
        new
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
        self.aiff_name.push(if self.is_whale[self.aiff_number] { '*' } else { ' ' });
        self.aiff_raw = aiff_samples(&self.full_name).unwrap();
        self.aiff_max = self.aiff_raw.iter().max_by(|x, y| x.abs().cmp(&y.abs())).unwrap().abs();
        self.aiff_data = self.aiff_raw.iter().map(|x| *x as f32 / self.aiff_max as f32).collect();

        self.recalculate();
    }

    fn recalculate(&mut self) {
        self.spectrogram = spectrogram(&self);
        self.noise = vec![(0.0, 0.0); self.fft_size[0]/2];
        for bin in 0 .. self.noise.len() {
            let mut row: Vec<f32> = vec![0.0; self.spectrogram.len()];
            self.spectrogram.iter()
                .enumerate()
                .for_each(|(i, fft)| row[i] = fft[bin]);
            self.noise[bin] = snr(&row);
        }
    }

    fn next_file(&mut self) {
        match self.aiff_type[0] {
            AiffType::AllFiles => {
                self.set_aiff(self.aiff_number + 1)
            },
            AiffType::WhalesOnly => {
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
            AiffType::AllFiles => {
                self.set_aiff(self.aiff_number - 1)
            },
            AiffType::WhalesOnly => {
                for i in (2..self.aiff_number).rev() {
                    if self.is_whale[i] {
                        self.set_aiff(i);
                        break;
                    }
                }
            },
        };
    }

    fn up(&mut self) {
        match self.modifiers {
            ModifiersState::LOGO => self.next_file(),
            ModifiersState::SHIFT => self.freq = (self.freq + 10.*self.bin_size()).min(WHALE_VIEW.end - self.bin_size()),
            _ => self.freq = (self.freq + self.bin_size()).min(WHALE_VIEW.end - self.bin_size()),
        }
    }

    fn down(&mut self) {
        match self.modifiers {
            ModifiersState::LOGO => self.prev_file(),
            ModifiersState::SHIFT => self.freq = (self.freq - 10.*self.bin_size()).max(WHALE_VIEW.start),
            _ => self.freq = (self.freq - self.bin_size()).max(WHALE_VIEW.start),
        }
    }

    fn range(&self) -> std::ops::Range<usize> {
        self.start..self.start + self.window
    }

    fn next_window(&mut self) {
        let inc = match self.modifiers {
            ModifiersState::SHIFT => self.slide[0] * 10,
            _ => self.slide[0],
        };
        self.start = (self.start + inc).min(4000 - self.window - self.slide[0]);
    }

    fn prev_window(&mut self) {
        let inc = match self.modifiers {
            ModifiersState::SHIFT => self.slide[0] * 10,
            _ => self.slide[0],
        };
        self.start = self.start.saturating_sub(inc);
    }

    fn digits(&mut self, digit: usize) {
        if self.modifiers == ModifiersState::empty() {
            self.digits.push_str(&digit.to_string());
        } else if self.modifiers == ModifiersState::SHIFT {
            match digit {
                1 => self.set_aiff(32),
                2 => self.set_aiff(55),
                3 => self.set_aiff(673),
                4 => self.set_aiff(19536),
                _ => (),
            }
        }
    }

    fn bin_size(&self) -> f32 {
        2000.0 / self.fft_size[0] as f32
    }

}

fn rotate<T>(vec: &mut Vec<T>, modifiers: ModifiersState) {
    if modifiers == ModifiersState::empty() {
        vec.rotate_left(1);
    } else if modifiers == ModifiersState::SHIFT {
        vec.rotate_right(1);
    }
}

fn play_aiff(m: &Model) {
    use rodio::buffer::SamplesBuffer;
    use rodio::{OutputStream, Sink};

    let (_stream, stream_handle) = OutputStream::try_default().unwrap();
    let sink = Sink::try_new(&stream_handle).unwrap();

    let source = SamplesBuffer::new(1, 2000, m.aiff_data.clone());
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

fn chirp(m: &mut Model) {
    use std::f64::consts::PI;

    const RATE: f64 = 2000.0; // samples per second
    let f0 = 100f64; // start frequency, Hz
    let f1 = 200f64; // end frequency, Hz
    let flen = f1 - f0;
    
    let s0 = 1000; // starting sample
    let s1 = 3000;
    let slen = s1 - s0;
    
    let mut df = flen / slen as f64; // change in freq per sample

    let mut f = f0; // frequency in Hz
    let mut r = 2.0 * PI * f / RATE;  // radians per sample
    let mut phase = 0.0f64; // phase in radians

    for i in 0 .. m.aiff_data.len() {

        let s = phase.sin();
        m.aiff_data[i] = (s * (i16::MAX as f64)) as f32;
        
        if (s0 .. s1).contains(&i) && f < 400.0 {
            f += df;
            r = 2.0 * PI * f / RATE;
            df *= 1.002;
        }
        
        phase += r;
    }

    m.recalculate();

}

fn record(m: &mut Model) {
    use std::time::Duration;
    use std::sync::Arc;
    use std::sync::Mutex;
    use cpal::*;
    use cpal::traits::*;


    let host = cpal::default_host();
    println!("Default host: {:?}", host.id());

    let device = host.default_input_device().unwrap();
    println!("Default input device: {:?}", device.name());

    let config = device.default_input_config().unwrap();
    println!("Default input config: {:?}", config);

    let config = StreamConfig {
        channels: 1,
        sample_rate: cpal::SampleRate(48000),
        buffer_size: BufferSize::Default, // 512 samples
    };

    let stream_data: Arc<Mutex<Vec<f32>>> = Arc::new(Mutex::new(Vec::new()));
    let stream_data_clone = stream_data.clone();
    let stream = device.build_input_stream_raw
        ( &config
        , SampleFormat::F32
        , move |data, _info| {
            if let Some(slice) = data.as_slice() {
                let mut stream_data = stream_data_clone.lock().unwrap();
                stream_data.extend_from_slice(slice);
            }
            // println!("recording {}", data.len());
        } 
        , |_| {}
        , Some(Duration::new(2, 0))
        ).unwrap();

    stream.play().unwrap();

    println!("Recording for ~2 seconds ...");
    std::thread::sleep(Duration::from_secs(2));
    stream.pause().unwrap();

    let stream_data = stream_data.lock().unwrap();

    println!("Recorded {} samples", stream_data.len());
    stream_data.iter().step_by(24).enumerate().for_each(|(i, x)| {
        if i < m.aiff_data.len() {
            m.aiff_data[i] = (*x * (i16::MAX as f32)) as f32;
        }
    });

    m.recalculate();

}

