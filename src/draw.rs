#![deny(clippy::all)]
#![forbid(unsafe_code)]

use plotters::backend::RGBXPixel;
use plotters::coord::Shift;
use plotters::prelude::*;

use crate::signal;
use crate::Model;
use crate::signal::spectrogram;

pub type MyDrawingArea<'a> = DrawingArea<BitMapBackend<'a, RGBXPixel>, Shift>;

pub fn draw(frame: &mut [u8], m: &Model) -> Result<(), Box<dyn std::error::Error>> {
    use crate::HEIGHT;
    use crate::WIDTH;

    let root_drawing_area =
        BitMapBackend::<RGBXPixel>::with_buffer_and_format(frame, (WIDTH, HEIGHT))
            .unwrap()
            .into_drawing_area();

    root_drawing_area.fill(&WHITE)?;
        
    let (upper, lower) = root_drawing_area.split_vertically(240);

    m.upper[0](upper, m)?;
    draw_lower(lower, m)?;

    root_drawing_area.present()?;
    Ok(())
}

pub fn draw_upper_signal(upper: MyDrawingArea, m: &Model) -> Result<(), Box<dyn std::error::Error>> {

  let max = m.max64[0].0(m, m.range());
  let signal: Vec<(usize, f64)> = m.aiff_data[m.range()].iter().enumerate()
    .map(|(i, f)| ((i + m.start) as usize, *f as f64 / max)).collect();

  let upper = upper.margin(0, 0, 0, 20);
  let upper_caption = format!(
      "{}  {}..{}  {}",
      m.aiff_name,
      m.start,
      m.start + m.window,
      m.max64[0].1,
  );

  let mut upper_chart = ChartBuilder::on(&upper)
      .caption(upper_caption, ("sans-serif", 40))
      .set_label_area_size(LabelAreaPosition::Left, 40)
      .set_label_area_size(LabelAreaPosition::Bottom, 40)
      .build_cartesian_2d(m.range(), -1.0..1.0)
      .unwrap();

  upper_chart
      .configure_mesh()
      .max_light_lines(5)
      .x_labels(5)
      .y_labels(3)
      .draw()?;

  upper_chart
      .draw_series(LineSeries::new(signal, BLUE.stroke_width(3)))
      .unwrap();

  upper.present()?;
  Ok(())
}


pub fn draw_upper_fft(upper: MyDrawingArea, m: &Model) -> Result<(), Box<dyn std::error::Error>> {

  let (fft, fmax) = signal::fft(m, m.range());
  let bin_size = 2000.0 / m.fft_size[0] as f32;
  let _bin_range = (25.0 / bin_size) as usize .. (400.0 / bin_size + 1.0) as usize;

  let upper = upper.margin(0, 0, 0, 20);
  let caption = format!(
    "{}  {:?}  Max {:.1}",
      m.aiff_name,
      m.range(),
      fmax,
  );

  let mut chart = ChartBuilder::on(&upper)
      .caption(caption, ("sans-serif", 30))
      .set_label_area_size(LabelAreaPosition::Left, 40)
      .set_label_area_size(LabelAreaPosition::Bottom, 40)
      .build_cartesian_2d((0..fft.len()).into_segmented(), 0.0..1.0)
      .unwrap();

  chart.configure_mesh()
      .axis_desc_style(("sans-serif", 30))
      .x_labels(5)
      .y_labels(5)
      .x_desc("Bin")
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

  upper.present()?;
  Ok(())
}


fn draw_lower(lower: MyDrawingArea, m: &Model) -> Result<(), Box<dyn std::error::Error>> {
  
  let (spec, _fmax, smax) = spectrogram(m);
  let bin_size = 2000.0 / m.fft_size[0] as f32;
  let bin_range = (25.0 / bin_size) as usize .. (400.0 / bin_size + 1.0) as usize;
  let bin_start = bin_range.start;

  let lower = lower.margin(0, 0, 0, 20);
  let lower_captiion = format!(
    "FFT {}  Bin size {:.2}Hz  Window {}  Slide {}  Max {:.1}",
    m.fft_size[0],
    bin_size,
    m.window,
    m.slide[0],
    smax,
  );
  
  let mut lower_chart = ChartBuilder::on(&lower)
  .caption(lower_captiion, ("sans-serif", 30))
  .set_label_area_size(LabelAreaPosition::Left, 40)
  .set_label_area_size(LabelAreaPosition::Bottom, 40)
  .build_cartesian_2d((0..spec.len()).into_segmented(), bin_range)
  .unwrap();

  let gradient = colorous::VIRIDIS;

  lower_chart
      .configure_mesh()
      .x_labels(5)
      .y_labels(5)
      .x_desc("Samples/10")
      .y_desc("bin")
      .draw()?;

  for (i, amplitudes) in spec.iter().enumerate() {
      lower_chart
          .draw_series(amplitudes.iter().enumerate().map(|(x, y)| {
              let x0 = SegmentValue::Exact(i);
              let x1 = SegmentValue::Exact(i + 1);
              let mut color = gradient.eval_continuous(*y as f64);
              if x == bin_start && i == m.range().start / m.slide[0] {
                color = colorous::Color { r: 255, g: 255, b: 255};
              }
              let style = RGBColor(color.r, color.g, color.b).filled();
              Rectangle::new([(x0, x), (x1, x + 1)], style)
          }))
          .unwrap();
  }

  lower.present()?;
  Ok(())
}
