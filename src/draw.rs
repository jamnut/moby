#![deny(clippy::all)]
#![forbid(unsafe_code)]

use plotters::prelude::*;
use plotters::backend::RGBXPixel;

use crate::WIDTH;
use crate::HEIGHT;
use crate::signal;
use crate::Model;

pub fn draw(frame: &mut [u8], m: &Model) {
  let window = 200;

  let aiff = signal::normalize(&m.aiff_data);
  let signal = signal::series(&aiff);
 
  let amplitudes = signal::amplitudes(m.window_start, &aiff);
  let amp_max = signal::max(&amplitudes);
  let amplitudes: Vec<f64> = signal::normalize(&amplitudes);
  let amplitudes: Vec<i32> = amplitudes.iter().map(|x| (*x * *x * 1000.0) as i32).collect();
  
  let root_drawing_area =
      BitMapBackend::<RGBXPixel>::with_buffer_and_format(frame, (WIDTH, HEIGHT))
          .unwrap()
          .into_drawing_area();

  root_drawing_area.fill(&WHITE).unwrap();

  let (upper, lower) = root_drawing_area.split_vertically(240);

  let upper = upper.margin(0, 0, 0, 20);
  let upper_caption = format!(
      "{} [{}..{}]",
      m.aiff_name,
      m.window_start,
      m.window_start + window,
  );

  let mut upper_chart = ChartBuilder::on(&upper)
      .caption(upper_caption, ("sans-serif", 40))
      .set_label_area_size(LabelAreaPosition::Left, 40)
      .set_label_area_size(LabelAreaPosition::Bottom, 40)
      .build_cartesian_2d(m.window_start..m.window_start+window, -1.1..1.1)
      .unwrap();

  upper_chart.configure_mesh().draw().unwrap();

  upper_chart
      .draw_series(LineSeries::new(signal, BLUE.stroke_width(3)))
      .unwrap();


  let lower = lower.margin(0, 0, 0, 20);
  let lower_captiion = format!("Bin {:.2} Hz  -  Max {:.0}", window as f64 / 1024.0 * 10.0, amp_max,);

  let mut lower_chart = ChartBuilder::on(&lower)
      .caption(lower_captiion, ("sans-serif", 40))
      .set_label_area_size(LabelAreaPosition::Left, 40)
      .set_label_area_size(LabelAreaPosition::Bottom, 40)
      .build_cartesian_2d((0..512).into_segmented(), 0..1000)
      .unwrap();

  lower_chart
      .draw_series(amplitudes.iter().enumerate().map(|(x, y)| {
          let x0 = SegmentValue::Exact(x as i32);
          let x1 = SegmentValue::Exact(x as i32 + 1);
          Rectangle::new([(x0, 0), (x1, *y)], RED.filled())
      }))
      .unwrap();

  root_drawing_area.present().unwrap();
}

