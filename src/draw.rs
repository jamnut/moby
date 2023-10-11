#![deny(clippy::all)]
#![forbid(unsafe_code)]

use plotters::backend::RGBXPixel;
use plotters::coord::Shift;
use plotters::prelude::*;

use crate::signal::{self, spectrogram, tracking};
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
    draw_spec(bl, m)?;
    draw_fft(br, m)?;
    draw_info(tr, m)?;

    let text = format!("{:?}", m.aiff_type[0]);
    let text_style = TextStyle::from(("serif", 30).into_font()).color(&BLACK);
    root.draw_text(&text, &text_style , (5,5))?;

    root.present()?;
    Ok(())
}

fn draw_info(area: &Drawing, m: &Model) -> Result<(), Box<dyn std::error::Error>> {
    let text = format!("{}", m.digits);
    let text_style = TextStyle::from(("serif", 30).into_font()).color(&BLACK);
    area.draw_text(&text, &text_style , (5,5))?;
    Ok(())
}

pub fn draw_signal( area: &Drawing, m: &Model, ) -> Result<(), Box<dyn std::error::Error>> {
    let max = m.max64[0].0(m, m.range()) as f32;
    let signal: Vec<(usize, f32)> = m.aiff_data[m.range()]
        .iter()
        .enumerate()
        .map(|(i, f)| (i + m.start, *f as f32 / max))
        .collect();

    let upper = area.margin(0, 0, 0, 20);
    let upper_caption = format!(
        "{}  {:?}  {} {}  {:.3}",
        m.aiff_name,
        m.range(),
        m.max64[0].1,
        max,
        max / i16::MAX as f32
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
    let (spec, _fmax, _fpeak, snr, _smax) = spectrogram(m);
    let bin_size = 2000.0 / spec[0].len() as f32;

    let tracks = tracking(&spec);

    let drawing = drawing.margin(0, 0, 0, 20);
    let caption = format!(
        "{}{}"
        , m.aiff_name
        , if m.is_whale[m.aiff_number] { "*" } else { " " }
    );

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
    let (fft, _fmax, _fpeak, snr) = signal::fft(m, m.start);
    let bin_size = 2000.0 / fft.len() as f32;

    let snr_db = f32::log10(snr.1 / snr.0) * 20.0;

    let caption = format!("{:?}  SNR ({:.0}/{:.0}) {:.0}  {:.0}dB", 
        m.range(), snr.0, snr.1, snr.0 / snr.1, snr_db);

    let (x0, x1, fscale) = m.fft_scale[0];

    let mut chart = ChartBuilder::on(&drawing)
        .caption(caption, ("sans-serif", 20))
        .set_label_area_size(LabelAreaPosition::Right, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .build_cartesian_2d(x0 .. x1, WHALE_VIEW)
        .unwrap();

    chart
        .configure_mesh()
        .axis_desc_style(("sans-serif", 30))
        .x_labels(5)
        .y_labels(5)
        .draw()?;

    let gradient = colorous::VIRIDIS;

    chart
        .draw_series(fft.iter().enumerate()
        .filter(|(bin, _mag)| *bin as f32 * bin_size >= WHALE_RANGE.start && *bin as f32 * bin_size <= WHALE_RANGE.end)
        .map(|(bin, mag)| {
            let x1 = fscale(*mag);
            let y0 = bin as f32 * bin_size;
            let y1 = y0 + bin_size;
            let color = gradient.eval_continuous(*mag as f64);
            let style = RGBColor(color.r, color.g, color.b).filled();
            Rectangle::new([(x0, y0), (x1, y1)], style)
        }))?;

    chart.plotting_area().draw(&Rectangle::new(
        [(x0, 0.0), (fscale(snr.1/snr.0), WHALE_VIEW.end)],
        RED.mix(0.1).filled(),
    ))?;

    drawing.present()?;
    Ok(())
}


fn draw_spec(drawing: &Drawing, m: &Model) -> Result<(), Box<dyn std::error::Error>> {
    let (spec, _fmax, _fpeak, _snr, _smax) = spectrogram(m);
    let bin_size = 2000.0 / spec[0].len() as f32;

    let drawing = drawing.margin(0, 0, 0, 20);
    let caption = format!(
        "FFT {}  Bin {:.2}Hz  Window {}  Slide {} ({}ms)",
        m.fft_size[0], bin_size, m.window, m.slide[0], m.slide[0] as f32 / 2 as f32
    );

    let mut chart = ChartBuilder::on(&drawing)
        .caption(caption, ("sans-serif", 20))
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        // .build_cartesian_2d((0..spec.len()).into_segmented(), crate::WHALE_VIEW)
        .build_cartesian_2d(0..4000usize, crate::WHALE_VIEW)
        .unwrap();

    let gradient = colorous::VIRIDIS;

    chart
        .configure_mesh()
        .x_labels(5)
        .y_labels(5)
        .x_desc("Sample")
        // .y_desc("HZ")
        .axis_desc_style(("sans-serif", 20))
        .draw()?;

    for (i, mags) in spec.iter().enumerate() {
        let x0 = i * m.slide[0];
        let x1 = (i + 1) * m.slide[0];
        let show_cursor = (m.range().start .. m.range().start + m.slide[0]).contains(&x0);
        chart
            .draw_series(mags.iter().enumerate().map(| (bin, mag)| {
                let y0 = bin as f32 * bin_size;
                let y1 = y0 + bin_size;
                let mut color = gradient.eval_continuous(*mag as f64);
                // Make bottom and top  "cursors" white
                if show_cursor && (y0 < WHALE_RANGE.start || y1 > WHALE_RANGE.end) {
                    color = colorous::Color { r: 255, g: 255, b: 255, };
                }
                let style = RGBColor(color.r, color.g, color.b).filled();
                Rectangle::new([(x0, y0), (x1, y1)], style)
            }))
            .unwrap();
    }

    drawing.present()?;
    Ok(())
}
