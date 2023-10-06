#![deny(clippy::all)]
#![forbid(unsafe_code)]

use plotters::backend::RGBXPixel;
use plotters::coord::Shift;
use plotters::prelude::*;

use crate::signal::{self, spectrogram, tracking};
use crate::{Model, WHALE_RANGE};

pub type Drawing<'a> = DrawingArea<BitMapBackend<'a, RGBXPixel>, Shift>;

pub fn draw((width, height): (u32, u32), frame: &mut [u8], m: &Model) -> Result<(), Box<dyn std::error::Error>> {

    let root_drawing_area =
        BitMapBackend::<RGBXPixel>::with_buffer_and_format(frame, (width, height))
            .unwrap()
            .into_drawing_area();

    root_drawing_area.fill(&WHITE)?;

    let (upper, lower) = root_drawing_area.split_vertically(240);

    m.upper[0](upper, m)?;
    draw_lower(lower, m)?;

    let text = format!("{:?}", m.aiff_type[0]);
    let text_style = TextStyle::from(("serif", 30).into_font()).color(&BLACK);
    root_drawing_area.draw_text(&text, &text_style , (5,5))?;

    root_drawing_area.present()?;
    Ok(())
}

pub fn draw_upper_signal( upper: Drawing, m: &Model, ) -> Result<(), Box<dyn std::error::Error>> {
    let max = m.max64[0].0(m, m.range()) as f32;
    let signal: Vec<(usize, f32)> = m.aiff_data[m.range()]
        .iter()
        .enumerate()
        .map(|(i, f)| (i + m.start, *f as f32 / max))
        .collect();

    let upper = upper.margin(0, 0, 0, 20);
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
    
    
pub fn draw_upper_tracking(upper: Drawing, m: &Model) -> Result<(), Box<dyn std::error::Error>> {
    let (spec, _fmax, _fpeak, _snr, _smax) = spectrogram(m);
    let bin_size = 2000.0 / spec[0].len() as f32;

    let tracks = tracking(&spec);

    let upper = upper.margin(0, 0, 0, 20);
    let upper_caption = format!(
        "{}{}  {:?}"
        , m.aiff_name
        , if m.is_whale[m.aiff_number] { "*" } else { " " }
        , m.range()
    );

    // let whale_view = (WHALE_VIEW.start / bin_size) as usize .. (WHALE_VIEW.end / bin_size + 1.0) as usize;
    let whale_range = (WHALE_RANGE.start / bin_size) as usize .. (WHALE_RANGE.end / bin_size) as usize;

    let mut upper_chart = ChartBuilder::on(&upper)
        .x_label_area_size(35)
        .y_label_area_size(40)
        // .right_y_label_area_size(20)
        .margin(5)
        .caption(upper_caption, ("sans-serif", 30))
        .build_cartesian_2d(0..spec.len(), whale_range)?;
  
    upper_chart
        .configure_mesh()
        .disable_x_mesh()
        .disable_y_mesh()
        .max_light_lines(5)
        .x_labels(5)
        .y_labels(3)
        .draw()?;

    let points: Vec<(usize, usize)> = tracks.as_points();

    upper_chart.draw_series(points.iter().map(|(bin, y)| {
        // Circle::new((*bin, *y), 2, BLUE.filled())
        let height = match m.fft_size[0] {
            2048 => 4,
            1024 => 2,
            _ => 1,
        };
        Rectangle::new([(*bin,*y),(*bin+1,*y+height)], BLUE.filled())
    }))?;

  upper.present()?;
  Ok(())
}
      
      
pub fn draw_upper_fft(upper: Drawing, m: &Model) -> Result<(), Box<dyn std::error::Error>> {
    let (fft, fmax, _fpeak, snr) = signal::fft(m, m.start);
    let bin_size = 2000.0 / fft.len() as f32;

    let upper = upper.margin(0, 0, 0, 20);
    let caption = format!("{}{}  {:?}  fmax {fmax:.1}  snr {snr:.2}"
      , m.aiff_name
      , if m.is_whale[m.aiff_number] { "*" } else { " " }
      , m.range()
      );

    let mut chart = ChartBuilder::on(&upper)
        .caption(caption, ("sans-serif", 30))
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .build_cartesian_2d((0..fft.len()).into_segmented(), 0.0..1.0)
        .unwrap();

    chart
        .configure_mesh()
        .axis_desc_style(("sans-serif", 30))
        .x_labels(5)
        .y_labels(5)
        // .x_desc("Bin")
        .draw()?;

    let gradient = colorous::VIRIDIS;

    chart
        .draw_series(fft.iter().enumerate().map(|(i, y)| {
            let x0 = SegmentValue::Exact(i);
            let x1 = SegmentValue::Exact(i + 1);
            let color = gradient.eval_continuous(*y as f64);
            let style = RGBColor(color.r, color.g, color.b).filled();
            Rectangle::new([(x0, 0.0), (x1, *y as f64)], style)
        }))
        .unwrap();

    // place a red rectangle at the top of the whale bins
    let whale_bins = WHALE_RANGE.start / bin_size .. WHALE_RANGE.end / bin_size;
    chart.draw_series((0..1).into_iter().map(|_| {
        let x0 = SegmentValue::Exact(whale_bins.start as usize);
        let x1 = SegmentValue::Exact(whale_bins.end as usize);
        let style = RED.filled();
        Rectangle::new([(x0, 0.98), (x1, 1.0)], style)
    }))?;



    upper.present()?;
    Ok(())
}

fn draw_lower(lower: Drawing, m: &Model) -> Result<(), Box<dyn std::error::Error>> {
    let (spec, _fmax, _fpeak, _snr, _smax) = spectrogram(m);
    let bin_size = 2000.0 / spec[0].len() as f32;

    let lower = lower.margin(0, 0, 0, 20);
    let lower_captiion = format!(
        "FFT {}  Bin size {:.2}Hz  Window {}  Slide {}",
        m.fft_size[0], bin_size, m.window, m.slide[0],
    );

    let mut lower_chart = ChartBuilder::on(&lower)
        .caption(lower_captiion, ("sans-serif", 30))
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .build_cartesian_2d((0..spec.len()).into_segmented(), crate::WHALE_VIEW)
        .unwrap();

    let gradient = colorous::VIRIDIS;

    lower_chart
        .configure_mesh()
        .x_labels(5)
        .y_labels(5)
        .x_desc("Sample/Slide")
        .y_desc("HZ")
        .axis_desc_style(("sans-serif", 20))
        .draw()?;

    for (i, amplitudes) in spec.iter().enumerate() {
        lower_chart
            .draw_series(amplitudes.iter().enumerate().map(| (bin, y)| {
                let x0 = SegmentValue::Exact(i);
                let x1 = SegmentValue::Exact(i + 1);
                let flo = bin as f32 * bin_size;
                let fhi = flo + bin_size;
                let mut color = gradient.eval_continuous(*y as f64);
                if (i * m.slide[0] .. (i + 1) * m.slide[0]).contains(&m.range().start) 
                && flo < crate::WHALE_RANGE.start {
                    color = colorous::Color { r: 255, g: 255, b: 255, };
                }
                let style = RGBColor(color.r, color.g, color.b).filled();
                Rectangle::new([(x0, flo), (x1, fhi)], style)
            }))
            .unwrap();
    }

    lower.present()?;
    Ok(())
}
