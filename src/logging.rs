//! A simple logging based progress counter using count and timing information
//!
//! Intended as a light-weight progress counter that easily integrate with existing logging output
//! in the vein of proglog but with timing and customised for my use case.
use log::info;
use std::{sync::Arc, time::Instant};

pub type LogFn = dyn Fn(&str) + Send + Sync;

/// Progress bar
///
/// Generic container for different progress bar options. Currently just a
/// logging progress bar and a NoOp dummy.
pub enum Progress<'a> {
    /// Logging progress bar
    Log(LogProgress<'a>),

    /// Dummy progress bar
    None,
}

impl<'a> Progress<'a> {
    /// Create a logging progress bar
    pub fn log(
        message: &'a str,
        final_message: &'a str,
        use_thread_id: bool,
        total: Option<u64>,
        log_interval: u64,
        log_fn: Arc<dyn Fn(&str) + Send + Sync>,
    ) -> Self {
        Self::Log(LogProgress::new(
            message,
            final_message,
            use_thread_id,
            total,
            log_interval,
            log_fn,
        ))
    }

    /// Create a dummy progress bar
    pub fn none() -> Self {
        Self::None
    }

    /// Create a progress bar from a ProgressStyle object
    ///
    /// Using this approach makes it easier to initiate a single style for all
    /// progress bars in a script
    pub fn from_style(
        style: &ProgressStyle,
        message: &'a str,
        final_message: &'a str,
        total: Option<u64>,
        log_interval: u64,
    ) -> Self {
        match &style.log_fn {
            None => Self::none(),
            Some(log_fn) => Self::log(
                message,
                final_message,
                style.use_thread_id,
                total,
                log_interval,
                log_fn.clone(),
            ),
        }
    }

    /// Increment the progress bar
    pub fn inc(&mut self, amount: u64) {
        match self {
            Self::Log(x) => x.inc(amount),
            Self::None => {}
        }
    }

    /// Complete the progress bar
    pub fn finish(&self) {
        match self {
            Self::Log(x) => x.finish(),
            Self::None => {}
        }
    }
}

/// Progress monitor outputing via logging
pub struct LogProgress<'a> {
    /// Message to output before each update
    message: &'a str,

    /// Message to output before final update
    final_message: &'a str,

    /// Message ID associated with the thread being logged
    thread_id: String,

    /// Total number of iterations expected. None means unknown
    total: Option<u64>,

    /// Current iteration count
    current: u64,

    /// Number of iterations between logging output
    log_interval: u64,

    /// When the opperation initially started
    start_time: Instant,

    /// When the last log update occured
    last_log_time: Instant,

    /// Count at last log instant
    last_log_count: u64,

    /// Logging function to use. For instance info!()
    log_fn: Arc<dyn Fn(&str) + Send + Sync>,
}

impl<'a> LogProgress<'a> {
    /// Create a new progress tracker.
    pub fn new(
        message: &'a str,
        final_message: &'a str,
        use_thread_id: bool,
        total: Option<u64>,
        log_interval: u64,
        log_fn: Arc<dyn Fn(&str) + Send + Sync>,
    ) -> Self {
        assert_ne!(log_interval, 0);

        let thread_id = if use_thread_id {
            format!(" [{:?}]", std::thread::current().id())
        } else {
            "".to_string()
        };

        let now = Instant::now();
        Self {
            message,
            final_message,
            thread_id,
            total,
            current: 0,
            start_time: now,
            last_log_time: now,
            last_log_count: 0,
            log_interval,
            log_fn,
        }
    }

    /// Increment progress by a specific amount.
    pub fn inc(&mut self, amount: u64) {
        self.current += amount;
        if self.current % self.log_interval == 0 {
            self.log_progress();
            self.last_log_time = Instant::now();
            self.last_log_count = self.current;
        }
    }

    /// Finish the progress and log the final message.
    pub fn finish(&self) {
        let elapsed = self.start_time.elapsed();
        let avg_rate = self.current as f64 / elapsed.as_secs_f64();

        match self.total {
            None => {
                (self.log_fn)(&format!(
                    "{} {} in {:.2?} | avg. rate: {:.2} items/s{}",
                    self.final_message, self.current, elapsed, avg_rate, self.thread_id
                ));
            }
            Some(total) => {
                (self.log_fn)(&format!(
                    "{} {}/{} 100% in {:.2?} | avg. rate: {:.2} items/s{}",
                    self.final_message, self.current, total, elapsed, avg_rate, self.thread_id
                ));
            }
        }
    }

    /// Log the current progress.
    fn log_progress(&self) {
        let elapsed = self.start_time.elapsed();
        let since_last = self.last_log_time.elapsed();

        let current_rate = if since_last.as_secs_f64() > 0.0 {
            (self.current - self.last_log_count) as f64 / since_last.as_secs_f64()
        } else {
            0.0
        };

        let avg_rate = if elapsed.as_secs_f64() > 0.0 {
            self.current as f64 / elapsed.as_secs_f64()
        } else {
            0.0
        };

        match self.total {
            None => {
                (self.log_fn)(&format!(
                    "{} {} in {:.2?} | current rate: {:.2} items/s | avg. rate: {:.2} items/s{}",
                    self.message, self.current, elapsed, current_rate, avg_rate, self.thread_id
                ));
            }
            Some(total) => {
                let percent: f64 = (self.current as f64) / (total as f64) * 100.0;
                let remaining: f64 = ((total - self.current) as f64) / avg_rate;

                (self.log_fn)(&format!(
                    "{} {}/{:.1} {:.0}% in {:.2?} | current rate: {:.2} items/s | avg. rate: {:.2} items/s | est {:.0}s remaining{}",
                    self.message,
                    self.current,
                    total,
                    percent,
                    elapsed,
                    current_rate,
                    avg_rate,
                    remaining,
                    self.thread_id
                ));
            }
        }
    }
}

/// Manager for multiple progress trackers with shared settings.
#[derive(Clone)]
pub struct ProgressStyle {
    log_fn: Option<Arc<LogFn>>,
    pub use_thread_id: bool,
}

impl ProgressStyle {
    /// Create a new progress manager.
    pub fn new(log_fn: Option<Arc<LogFn>>, use_thread_id: bool) -> Self {
        Self {
            log_fn,
            use_thread_id,
        }
    }
}

impl Default for ProgressStyle {
    fn default() -> Self {
        ProgressStyle::new(Some(Arc::new(|msg| info!("{}", msg))), false)
    }
}
