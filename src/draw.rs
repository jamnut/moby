#![deny(clippy::all)]
#![forbid(unsafe_code)]

use num::Integer;
use plotters::backend::RGBXPixel;
use plotters::coord::Shift;
use plotters::prelude::*;

use crate::signal::{self, old_spectrogram, tracking, spectrogram};
use crate::{Model, WHALE_RANGE, WHALE_VIEW};

pub type Drawing<'a> = DrawingArea<BitMapBackend<'a, RGBXPixel>, Shift>;

pub fn draw((width, height): (u32, u32), frame: &mut [u8], m: &Model) -> Result<(), Box<dyn std::error::Error>> {

    let root =
        BitMapBackend::<RGBXPixel>::with_buffer_and_format(frame, (width, height))
            .unwrap()
            .into_drawing_area();

    root.fill(&WHITE)?;

    let areas = root.split_by_breakpoints([width*3/4], [height/4]);
    let [ref tl, ref tr, ref bl, ref br] 
        = areas[0..4] else {panic!("Not enough areas")};
        
    m.upper[0](tl, m)?;
    draw_spectrogram(bl, m)?;
    draw_fft(br, m)?;
    draw_info(tr, m)?;

    root.present()?;
    Ok(())
}

fn draw_info(area: &Drawing, m: &Model) -> Result<(), Box<dyn std::error::Error>> {
    let font_size = 20;
    let text_style = TextStyle::from(("serif", font_size).into_font()).color(&BLACK);

    // Area in pixels (378, 229)
    let text = 
        [ &format!("{}", if m.digits.len() == 0 { "[DIGITS]" } else { &m.digits })
        , &format!("W: {:?}",m.aiff_type[0])
        , "Cmd-Up/Down: file"
        , "Left/Right: slice"
        , "F: FFT size"
        , "D: dB scale"
        , &format!("S: Slide ({})", m.slide[0])
        , "H: Hann window"
        , "U: Upper display"
        , "P: Play audio"
        , "R: Record audo"
        , "Q: Quit"
        ];

    for (i, label) in text.iter().enumerate() {
        let row = (i / 2) as i32 * font_size + 6; 
        let col = if i.is_even() { 5 } else { 185 };
        area.draw_text(*label, &text_style , (col, row+25))?;
    }

    Ok(())
}

pub fn draw_raw( area: &Drawing, m: &Model, ) -> Result<(), Box<dyn std::error::Error>> {

    let upper = area.margin(0, 0, 0, 20);
    let upper_caption = format!(
        "{}  {:?}",
        m.aiff_name,
        m.range(),
    );

    let yrange = -32_768 .. 32_768;
    let mut upper_chart = ChartBuilder::on(&upper)
        .caption(upper_caption, ("sans-serif", 30))
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .build_cartesian_2d(m.range(), yrange)
        .unwrap();

    upper_chart
        .configure_mesh()
        .max_light_lines(5)
        .x_labels(5)
        .y_labels(3)
        .draw()?;

    upper_chart
        .draw_series(
            LineSeries::new(
                m.aiff_raw[m.range()].iter().enumerate().map(|(i, f)| (i + m.start, *f as i32)),
                BLUE.stroke_width(3)))?;

    upper.present()?;
    Ok(())
}
    

pub fn draw_signal( area: &Drawing, m: &Model, ) -> Result<(), Box<dyn std::error::Error>> {
    let signal: Vec<(usize, f32)> = m.aiff_data[m.range()]
        .iter()
        .enumerate()
        .map(|(i, f)| (i + m.start, *f))
        .collect();

    let upper = area.margin(0, 0, 0, 20);
    let upper_caption = format!(
        "{}  {:?}",
        m.aiff_name,
        m.range(),
    );

    let mut upper_chart = ChartBuilder::on(&upper)
        .caption(upper_caption, ("sans-serif", 30))
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .build_cartesian_2d(m.range(), -1.0..1.0f32)
        .unwrap();

    upper_chart
        .configure_mesh()
        .max_light_lines(5)
        .x_labels(5)
        .y_labels(3)
        .draw()?;

    upper_chart
        .draw_series(
            LineSeries::new(signal, BLUE.stroke_width(3)))?;

    upper.present()?;
    Ok(())
}
    

pub fn draw_tracking(drawing: &Drawing, m: &Model) -> Result<(), Box<dyn std::error::Error>> {
    let (spec, _fmax, _fpeak, snr, _smax) = old_spectrogram(m);
    let bin_size = 1000.0 / spec[0].len() as f32;

    let tracks = tracking(&spec);

    let drawing = drawing.margin(0, 0, 0, 20);
    let caption = format!( "{}" , m.aiff_name );

    // let whale_view = (WHALE_VIEW.start / bin_size) as usize .. (WHALE_VIEW.end / bin_size + 1.0) as usize;
    // let whale_range = (WHALE_RANGE.start / bin_size) as usize .. (WHALE_RANGE.end / bin_size) as usize;

    let chart = ChartBuilder::on(&drawing)
        .x_label_area_size(35)
        .y_label_area_size(40)
        .margin(5)
        .caption(caption, ("sans-serif", 30))
        .build_cartesian_2d(0..4000usize, WHALE_RANGE)?;
    
    let mut chart 
        = chart.set_secondary_coord(0..4000usize, 0.0..0.5f32);

    chart
        .configure_mesh()
        .disable_x_mesh()
        .disable_y_mesh()
        .max_light_lines(5)
        .x_labels(5)
        .y_labels(3)
        .draw()?;


    chart.draw_secondary_series(
        AreaSeries::new(
            snr.iter().enumerate().map(|(i, (s, n))| (i * m.slide[0], n / s)),
            0.0,
            RED.mix(0.1).filled(),
    ))?;

    // let to_ms_hz = |(w, h)| (w * m.slide[0]/2, h * bin_size as usize);

    let draw_height = match m.fft_size[0] {
        2048 => 4,
        1024 => 2,
        _ => 1,
    };

    chart.draw_series(
        tracks.current.iter().chain(tracks.history.iter()).flat_map(|t| {
            let mut mix = 0.1;

            if let Some(track) = tracks.tallest() {
                if *track == *t {
                    mix = 1.0;
                }
            }

            if let Some(track) = tracks.widest() {
                if *track == *t {
                    mix = 1.0;
                }
            }

            let style = BLUE.mix(mix).filled();
        
            t.iter().map(move |(i, bin)| {
                let x0 = i * m.slide[0];
                let x1 = x0 + m.slide[0];
                let y0 = *bin as f32 * bin_size;
                let y1 = (*bin + draw_height) as f32 * bin_size;
                Rectangle::new([(x0, y0), (x1, y1)], style)
            })
        })
    )?;

  drawing.present()?;
  Ok(())
}
      
      
pub fn draw_fft(drawing: &Drawing, m: &Model) -> Result<(), Box<dyn std::error::Error>> {

    let fft = signal::fft(m, m.start);
    let snr = signal::snr(&fft);

    let snr_db = f32::log10(snr.1 / snr.0) * 20.0;

    let caption = format!("{:?}  SNR ({:.0}/{:.0}) {:.0}  {:.0}dB", 
        m.range(), snr.0, snr.1, snr.0 / snr.1, snr_db);

    let (x0, x1, fscale) = m.fft_scale[0];

    let mut chart = ChartBuilder::on(&drawing)
        .caption(caption, ("sans-serif", 20))
        .set_label_area_size(LabelAreaPosition::Right, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .build_cartesian_2d(x0 .. x1, 0.0 .. 500f32)
        .unwrap();

    chart
        .configure_mesh()
        .axis_desc_style(("sans-serif", 30))
        .x_labels(5)
        .y_labels(5)
        .draw()?;

    // let gradient = colorous::VIRIDIS;

    // chart
    //     .draw_series(fft.iter().enumerate()
    //     // .filter(|(bin, _mag)| *bin as f32 * m.bin_size() >= WHALE_RANGE.start && *bin as f32 * m.bin_size() <= WHALE_RANGE.end)
    //     .map(|(bin, mag)| {
    //         let x1 = fscale(*mag);
    //         let y0 = bin as f32 * m.bin_size();
    //         let y1 = y0 + m.bin_size();
    //         let color = gradient.eval_continuous(*mag as f64);
    //         let style = RGBColor(color.r, color.g, color.b).filled();
    //         Rectangle::new([(x0, y0), (x1, y1)], style)
    //     }))?;

    chart.plotting_area().draw(&Rectangle::new(
        [(x0, 0.0), (fscale(snr.1/snr.0), WHALE_VIEW.end)],
        RED.mix(0.1).filled(),
    ))?;

    drawing.present()?;
    Ok(())
}


pub fn draw_noise(drawing: &Drawing, m: &Model) -> Result<(), Box<dyn std::error::Error>> {

    let spec = spectrogram(m);

    let bin = (m.freq / m.bin_size()) as usize;
    let vs = spec.iter().map(|f| f[bin]).collect::<Vec<f32>>();
    let vmax: f32 = vs.iter().fold(0.0, |a, b| a.max(*b));

    let mut noise = vs.clone();
    noise.sort_by(|a, b| a.total_cmp(b));
    let noise: f32 = noise[noise.len()/2] / vmax;

    let noise = vec![ noise; spec.len() ];

    let drawing = drawing.margin(0, 0, 0, 20);
    let caption = format!(
        "{}  --  noise max {:.3}  median {:.3} ratio {:.3}"
        , m.aiff_name
        , vmax
        , noise[0] * vmax
        , 1.0 / noise[0]
    );

    let chart = ChartBuilder::on(&drawing)
        .x_label_area_size(35)
        .y_label_area_size(40)
        .margin(5)
        .caption(caption, ("sans-serif", 30))
        .build_cartesian_2d(0..4000usize, 0.0 .. 1.0f32)?;
    
    let mut chart 
        = chart.set_secondary_coord(0..4000usize, 0.0 .. 1.0f32);

    chart
        .configure_mesh()
        .disable_x_mesh()
        .disable_y_mesh()
        .max_light_lines(5)
        .x_labels(5)
        .y_labels(3)
        .draw()?;


    chart
        .draw_series(vs.iter().enumerate().map(|(i, v)| {
            let x0 = i * m.slide[0];
            let x1 = x0 + m.slide[0];
            let style = BLUE.stroke_width(1);
            Rectangle::new([(x0, 0.0), (x1, *v/vmax)], style)
        }))?;

    chart.draw_secondary_series(
        AreaSeries::new(
            noise.iter().enumerate().map(|(i, n)| (i * m.slide[0], *n)),
            0.0,
            RED.mix(0.1).filled(),
    ))?;

    drawing.present()?;
    Ok(())      
}


fn draw_spectrogram(drawing: &Drawing, m: &Model) -> Result<(), Box<dyn std::error::Error>> {
            
    let gradient = colorous::VIRIDIS;
            
    let drawing = drawing.margin(0, 0, 0, 20);
    let caption = format!(
        "FFT {}  Bin {:.2}Hz  Window {}  Slide {} ({}ms)",
        m.fft_size[0], m.bin_size(), m.window, m.slide[0], m.slide[0] as f32 / 2.0
    );

    let mut chart = ChartBuilder::on(&drawing)
        .caption(caption, ("sans-serif", 20))
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .build_cartesian_2d(0..4000usize, WHALE_VIEW)
        .unwrap();

    chart
        .configure_mesh()
        .x_labels(5)
        .y_labels(5)
        .x_desc("Sample")
        .axis_desc_style(("sans-serif", 20))
        .draw()?;

    for (i, fft) in m.spectrogram.iter().enumerate() {
        let x0 = i * m.slide[0];
        let x1 = x0 + m.slide[0];
        chart
            .draw_series(fft.iter().enumerate().map(| (bin, power)| {
                let y0 = bin as f32 * m.bin_size();
                let y1 = y0 + m.bin_size();
                let color = gradient.eval_continuous((*power/m.noise[bin].0) as f64);
                let style = RGBColor(color.r, color.g, color.b).filled();
                Rectangle::new([(x0, y0), (x1, y1)], style)
            }))?;
    }

    // draw cursors
    let y0 = m.freq;
    let y1 = y0 + m.bin_size();
    let right = 4000 - m.window;
    let x0 = m.start;
    let x1 = x0 + m.slide[0];
    let style = WHITE.filled();

    let cursors = vec![
        [(0, y0), (200, y1)],  // left
        [(right-200, y0), (right, y1)],  // right
        [(x0, 0.0), (x1, 50.0)],  // top
        [(x0, 450.0), (x1, 500.0)],  // bottom
    ];

    chart
        .draw_series(
            cursors.iter().map(|c| 
                Rectangle::new(*c, style)
            )
        )?;


    drawing.present()?;
    Ok(())
}
