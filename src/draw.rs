#![deny(clippy::all)]
#![forbid(unsafe_code)]

use plotters::backend::BGRXPixel;
use plotters::prelude::*;

use crate::signal;
use crate::Model;
use crate::WIDTH;
use crate::HEIGHT;

pub fn draw(frame: &mut [u8], m: &Model) {
  let window = 200;
  let aiff = m.aiff_data.iter().map(|s| *s as f64).collect::<Vec<f64>>();
  let (aiff, aiff_max) = signal::normalize(&aiff);
  let signal = signal::signal(&aiff);
  let amplitudes = signal::amplitudes(window, &aiff);
  let (amp, amp_max) = signal::normalize(&amplitudes);

  let root_drawing_area =
      BitMapBackend::<BGRXPixel>::with_buffer_and_format(frame, (WIDTH, HEIGHT))
          .unwrap()
          .into_drawing_area();

  root_drawing_area.fill(&WHITE).unwrap();

  let (upper, lower) = root_drawing_area.split_vertically(240);

  let upper = upper.margin(0, 0, 0, 20);
  let upper_caption = format!(
      "{} [{}..{}]  max {}",
      m.aiff_name,
      0,
      window,
      aiff_max
  );

  let mut upper_chart = ChartBuilder::on(&upper)
      .caption(upper_caption, ("sans-serif", 40))
      .set_label_area_size(LabelAreaPosition::Left, 40)
      .set_label_area_size(LabelAreaPosition::Bottom, 40)
      .build_cartesian_2d(0..window, -1.1..1.1)
      .unwrap();

  upper_chart.configure_mesh().draw().unwrap();

  upper_chart
      .draw_series(LineSeries::new(signal, RED.stroke_width(3)))
      .unwrap();


  let lower = lower.margin(0, 0, 0, 20);
  let lower_captiion = format!("Bin {:.2} Hz  -  Max {:.0}", window as f64 / 1024.0 * 10.0, amp_max,);

  let mut lower_chart = ChartBuilder::on(&lower)
      .caption(lower_captiion, ("sans-serif", 40))
      .set_label_area_size(LabelAreaPosition::Left, 40)
      .set_label_area_size(LabelAreaPosition::Bottom, 40)
      .build_cartesian_2d((0..512).into_segmented(), 0..30)
      .unwrap();

  lower_chart
      .draw_series(amplitudes.iter().enumerate().map(|(x, y)| {
          let x0 = SegmentValue::Exact(x as i32);
          let x1 = SegmentValue::Exact(x as i32 + 1);
          let y0: i32 = *y as i32;
          Rectangle::new([(x0, 0), (x1, y0)], BLUE.filled())
      }))
      .unwrap();

  root_drawing_area.present().unwrap();
}

