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
    let [ref tl, ref _tr, ref bl, ref br] 
        = areas[0..4] else {panic!("Not enough areas")};
        
    m.upper[0](tl, m)?;
    draw_spec(bl, m)?;
    draw_fft(br, m)?;

    let text = format!("{:?}", m.aiff_type[0]);
    let text_style = TextStyle::from(("serif", 30).into_font()).color(&BLACK);
    root.draw_text(&text, &text_style , (5,5))?;

    root.present()?;
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
    let whale_range = (WHALE_RANGE.start / bin_size) as usize .. (WHALE_RANGE.end / bin_size) as usize;

    let chart = ChartBuilder::on(&drawing)
        .x_label_area_size(35)
        .y_label_area_size(40)
        .margin(5)
        .caption(caption, ("sans-serif", 30))
        .build_cartesian_2d(0..spec.len(), whale_range)?;
    
    let mut chart 
        = chart.set_secondary_coord(0..spec.len(), 0.0..0.5f32);

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
            snr.iter().enumerate().map(|(i, (s, n))| (i, n / s)),
            0.0,
            RED.mix(0.1).filled(),
    ))?;


    let points: Vec<(usize, usize)> = tracks.as_points();

    chart.draw_series(points.iter().map(|(bin, y)| {
        // Circle::new((*bin, *y), 2, BLUE.filled())
        let height = match m.fft_size[0] {
            2048 => 4,
            1024 => 2,
            _ => 1,
        };
        Rectangle::new([(*bin,*y),(*bin+1,*y+height)], BLUE.filled())
    }))?;

  drawing.present()?;
  Ok(())
}
      
      
pub fn draw_fft(drawing: &Drawing, m: &Model) -> Result<(), Box<dyn std::error::Error>> {
    let (fft, _fmax, _fpeak, snr) = signal::fft(m, m.start);
    let bin_size = 2000.0 / fft.len() as f32;

    let caption = format!("{:?}  SNR ({:.0}/{:.0}) {:.0}", 
        m.range(), snr.0, snr.1, snr.0 / snr.1);

    let mut chart = ChartBuilder::on(&drawing)
        .caption(caption, ("sans-serif", 20))
        .set_label_area_size(LabelAreaPosition::Right, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .build_cartesian_2d(0.0..1f32, WHALE_VIEW)
        .unwrap();

    chart
        .configure_mesh()
        .axis_desc_style(("sans-serif", 30))
        .x_labels(5)
        .y_labels(5)
        .draw()?;

    let gradient = colorous::VIRIDIS;

    chart
    .draw_series(fft.iter().enumerate().map(|(bin, mag)| {
        let y0 = bin as f32 * bin_size;
        let y1 = y0 + bin_size;
        let color = gradient.eval_continuous(*mag as f64);
        let style = RGBColor(color.r, color.g, color.b).filled();
        Rectangle::new([(0.0, y0), (*mag as f32, y1)], style)
    }))?;
    
    chart.plotting_area().draw(&Rectangle::new(
        [(0.0, 0.0), (snr.1/snr.0, WHALE_VIEW.end)],
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
        "FFT {}  Bin {:.2}Hz  Window {}  Slide {}",
        m.fft_size[0], bin_size, m.window, m.slide[0],
    );

    let mut chart = ChartBuilder::on(&drawing)
        .caption(caption, ("sans-serif", 30))
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .build_cartesian_2d((0..spec.len()).into_segmented(), crate::WHALE_VIEW)
        .unwrap();

    let gradient = colorous::VIRIDIS;

    chart
        .configure_mesh()
        .x_labels(5)
        .y_labels(5)
        .x_desc("Sample/Slide")
        .y_desc("HZ")
        .axis_desc_style(("sans-serif", 20))
        .draw()?;

    for (i, mags) in spec.iter().enumerate() {
        chart
            .draw_series(mags.iter().enumerate().map(| (bin, y)| {
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

    drawing.present()?;
    Ok(())
}
