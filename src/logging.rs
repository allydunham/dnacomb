//! Lightweight progress reporting integrated with logging.
//!
//! This module provides a simple progress-reporting abstraction that emits
//! periodic updates through a logging-style callback. It is designed for
//! long-running counting and comparison steps where full terminal progress bars
//! are unnecessary or awkward.
use log::info;
use std::{sync::Arc, time::Instant};

/// Logging callback used by progress reporters.
///
/// Typically this will wrap a logging macro such as `info!`.
pub type LogFn = dyn Fn(&str) + Send + Sync;

/// Generic progress reporter.
///
/// This wraps either:
/// - a real logging-based progress reporter,
/// - or a no-op implementation when progress output is disabled.
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

    /// Construct a progress reporter from shared style settings.
    ///
    /// This makes it easy to apply the same logging behaviour consistently across
    /// multiple stages of a workflow.
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

/// Progress reporter that periodically emits updates through a logging callback.
///
/// Progress is reported in terms of processed item count, elapsed time, current
/// rate, and average rate. If a total is known, percentage completion and an
/// estimated remaining time are also reported.
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

    /// When the operation initially started
    start_time: Instant,

    /// When the last log update occured
    last_log_time: Instant,

    /// Count at last log instant
    last_log_count: u64,

    /// Logging function to use. For instance info!()
    log_fn: Arc<dyn Fn(&str) + Send + Sync>,
}

impl<'a> LogProgress<'a> {
    /// Create a new logging progress reporter.
    ///
    /// `message` is used for intermediate updates, `final_message` for the final
    /// completion line, `total` optionally sets the expected total item count,
    /// and `log_interval` determines how many processed items occur between log
    /// updates.
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

    /// Increment the processed item count and emit an update if the logging
    /// interval has been reached.
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
                    "{} {}/{} {:.0}% in {:.2?} | current rate: {:.2} items/s | avg. rate: {:.2} items/s | est {:.0}s remaining{}",
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

/// Shared configuration for constructing progress reporters.
///
/// This lets different stages of a workflow share the same logging function and
/// the same choice of whether thread IDs should be included in progress output.
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    // ---------- PROGRESS STYLE ----------
    #[test]
    fn progress_style_new() {
        let style = ProgressStyle::new(None, false);
        assert!(!style.use_thread_id);
        assert!(style.log_fn.is_none());
    }

    #[test]
    fn progress_style_new_with_log_fn() {
        let log_fn = Arc::new(|_: &str| {});
        let style = ProgressStyle::new(Some(log_fn), true);
        assert!(style.use_thread_id);
        assert!(style.log_fn.is_some());
    }

    #[test]
    fn progress_style_default() {
        let style = ProgressStyle::default();
        assert!(!style.use_thread_id);
        assert!(style.log_fn.is_some());
    }

    #[test]
    fn progress_style_clone() {
        let log_fn = Arc::new(|_: &str| {});
        let style1 = ProgressStyle::new(Some(log_fn), true);
        let style2 = style1.clone();
        assert_eq!(style1.use_thread_id, style2.use_thread_id);
    }

    // Progress - None
    #[test]
    fn progress_none_inc() {
        let mut progress = Progress::none();
        progress.inc(100); // Should not panic
    }

    #[test]
    fn progress_none_finish() {
        let progress = Progress::none();
        progress.finish(); // Should not panic
    }

    // Progress - Log
    #[test]
    fn progress_log_new() {
        let messages = Arc::new(std::sync::Mutex::new(Vec::new()));
        let messages_clone = Arc::clone(&messages);

        let log_fn = Arc::new(move |msg: &str| {
            messages_clone.lock().unwrap().push(msg.to_string());
        });

        let mut progress = Progress::log("Processing", "Done", false, Some(100), 10, log_fn);

        progress.inc(10);

        let logged = messages.lock().unwrap();
        assert!(logged.len() > 0);
    }

    #[test]
    fn progress_log_increments() {
        let messages = Arc::new(std::sync::Mutex::new(Vec::new()));
        let messages_clone = Arc::clone(&messages);

        let log_fn = Arc::new(move |msg: &str| {
            messages_clone.lock().unwrap().push(msg.to_string());
        });

        let mut progress = Progress::log("Processing", "Done", false, Some(100), 10, log_fn);

        for _ in 0..5 {
            progress.inc(10);
        }

        let logged = messages.lock().unwrap();
        // Should have logged at least once (at 50 items)
        assert!(logged.len() > 0);
    }

    #[test]
    fn progress_log_finish() {
        let messages = Arc::new(std::sync::Mutex::new(Vec::new()));
        let messages_clone = Arc::clone(&messages);

        let log_fn = Arc::new(move |msg: &str| {
            messages_clone.lock().unwrap().push(msg.to_string());
        });

        let mut progress = Progress::log("Processing", "Complete", false, Some(100), 50, log_fn);

        progress.inc(100);
        progress.finish();

        let logged = messages.lock().unwrap();
        let final_msg = logged.last().unwrap();
        assert!(final_msg.contains("Complete"));
        assert!(final_msg.contains("100"));
    }

    #[test]
    fn progress_from_style_with_log_fn() {
        let messages = Arc::new(std::sync::Mutex::new(Vec::new()));
        let messages_clone = Arc::clone(&messages);

        let log_fn = Arc::new(move |msg: &str| {
            messages_clone.lock().unwrap().push(msg.to_string());
        });

        let style = ProgressStyle::new(Some(log_fn), false);
        let mut progress = Progress::from_style(&style, "Processing", "Done", Some(100), 20);

        progress.inc(20);

        let logged = messages.lock().unwrap();
        assert!(logged.len() > 0);
    }

    #[test]
    fn progress_from_style_without_log_fn() {
        let style = ProgressStyle::new(None, false);
        let mut progress = Progress::from_style(&style, "Processing", "Done", Some(100), 20);

        progress.inc(100); // Should not panic
        progress.finish(); // Should not panic
    }

    // Logprogress
    #[test]
    fn logprogress_new() {
        let log_fn = Arc::new(|_: &str| {});
        let log_progress = LogProgress::new("msg", "final", false, Some(100), 10, log_fn);

        assert_eq!(log_progress.message, "msg");
        assert_eq!(log_progress.final_message, "final");
        assert_eq!(log_progress.current, 0);
        assert_eq!(log_progress.total, Some(100));
        assert_eq!(log_progress.log_interval, 10);
    }

    #[test]
    fn logprogress_new_with_thread_id() {
        let log_fn = Arc::new(|_: &str| {});
        let log_progress = LogProgress::new("msg", "final", true, None, 5, log_fn);

        assert!(!log_progress.thread_id.is_empty());
        assert!(log_progress.thread_id.contains("["));
        assert!(log_progress.thread_id.contains("]"));
    }

    #[test]
    fn logprogress_new_without_thread_id() {
        let log_fn = Arc::new(|_: &str| {});
        let log_progress = LogProgress::new("msg", "final", false, None, 5, log_fn);

        assert_eq!(log_progress.thread_id, "");
    }

    #[test]
    #[should_panic(expected = "assertion")]
    fn logprogress_zero_log_interval_panics() {
        let log_fn = Arc::new(|_: &str| {});
        LogProgress::new("msg", "final", false, None, 0, log_fn);
    }

    #[test]
    fn logprogress_inc_single() {
        let log_fn = Arc::new(|_: &str| {});
        let mut log_progress = LogProgress::new("msg", "final", false, Some(100), 10, log_fn);

        log_progress.inc(1);
        assert_eq!(log_progress.current, 1);
    }

    #[test]
    fn logprogress_inc_multiple() {
        let log_fn = Arc::new(|_: &str| {});
        let mut log_progress = LogProgress::new("msg", "final", false, Some(100), 10, log_fn);

        log_progress.inc(5);
        assert_eq!(log_progress.current, 5);
        log_progress.inc(3);
        assert_eq!(log_progress.current, 8);
    }

    #[test]
    fn logprogress_inc_triggers_log() {
        let messages = Arc::new(std::sync::Mutex::new(Vec::new()));
        let messages_clone = Arc::clone(&messages);

        let log_fn = Arc::new(move |msg: &str| {
            messages_clone.lock().unwrap().push(msg.to_string());
        });

        let mut log_progress = LogProgress::new("msg", "final", false, Some(100), 10, log_fn);

        log_progress.inc(10);

        let logged = messages.lock().unwrap();
        assert_eq!(logged.len(), 1);
        assert!(logged[0].contains("msg"));
    }

    #[test]
    fn logprogress_inc_at_interval_boundary() {
        let messages = Arc::new(std::sync::Mutex::new(Vec::new()));
        let messages_clone = Arc::clone(&messages);

        let log_fn = Arc::new(move |msg: &str| {
            messages_clone.lock().unwrap().push(msg.to_string());
        });

        let mut log_progress = LogProgress::new("msg", "final", false, Some(100), 10, log_fn);

        log_progress.inc(9);
        let logged_count_9 = messages.lock().unwrap().len();
        assert_eq!(logged_count_9, 0);

        log_progress.inc(1);
        let logged_count_10 = messages.lock().unwrap().len();
        assert_eq!(logged_count_10, 1);
    }

    #[test]
    fn logprogress_inc_multiple_boundaries() {
        let messages = Arc::new(std::sync::Mutex::new(Vec::new()));
        let messages_clone = Arc::clone(&messages);

        let log_fn = Arc::new(move |msg: &str| {
            messages_clone.lock().unwrap().push(msg.to_string());
        });

        let mut log_progress = LogProgress::new("msg", "final", false, Some(100), 5, log_fn);

        for _ in 0..4 {
            log_progress.inc(5);
        }

        let logged = messages.lock().unwrap();
        assert_eq!(logged.len(), 4);
    }

    #[test]
    fn logprogress_inc_large_amount() {
        let log_fn = Arc::new(|_: &str| {});
        let mut log_progress = LogProgress::new("msg", "final", false, Some(1000), 10, log_fn);

        log_progress.inc(500);
        assert_eq!(log_progress.current, 500);
    }

    #[test]
    fn logprogress_output_contains_message() {
        let messages = Arc::new(std::sync::Mutex::new(Vec::new()));
        let messages_clone = Arc::clone(&messages);

        let log_fn = Arc::new(move |msg: &str| {
            messages_clone.lock().unwrap().push(msg.to_string());
        });

        let mut log_progress =
            LogProgress::new("Processing reads", "final", false, Some(100), 10, log_fn);
        log_progress.inc(10);

        let logged = messages.lock().unwrap();
        assert!(logged[0].contains("Processing reads"));
    }

    // Edge cases
    #[test]
    fn logprogress_zero_current() {
        let messages = Arc::new(std::sync::Mutex::new(Vec::new()));
        let messages_clone = Arc::clone(&messages);

        let log_fn = Arc::new(move |msg: &str| {
            messages_clone.lock().unwrap().push(msg.to_string());
        });

        let log_progress = LogProgress::new("msg", "final", false, Some(100), 10, log_fn);
        log_progress.finish();

        let logged = messages.lock().unwrap();
        assert_eq!(logged.len(), 1);
        assert!(logged[0].contains("0/100"));
    }

    #[test]
    #[should_panic]
    fn logprogress_exceed_total() {
        let messages = Arc::new(std::sync::Mutex::new(Vec::new()));
        let messages_clone = Arc::clone(&messages);

        let log_fn = Arc::new(move |msg: &str| {
            messages_clone.lock().unwrap().push(msg.to_string());
        });

        let mut log_progress = LogProgress::new("msg", "final", false, Some(100), 10, log_fn);
        log_progress.inc(50);
        log_progress.inc(60);
    }

    #[test]
    fn logprogress_exact_total() {
        let messages = Arc::new(std::sync::Mutex::new(Vec::new()));
        let messages_clone = Arc::clone(&messages);

        let log_fn = Arc::new(move |msg: &str| {
            messages_clone.lock().unwrap().push(msg.to_string());
        });

        let mut log_progress = LogProgress::new("msg", "final", false, Some(100), 100, log_fn);
        log_progress.inc(100);
        log_progress.finish();

        let logged = messages.lock().unwrap();
        let final_msg = logged.last().unwrap();
        assert!(final_msg.contains("100/100"));
        assert!(final_msg.contains("100%"));
    }

    // Thread ID
    #[test]
    fn logprogress_with_thread_id_in_output() {
        let messages = Arc::new(std::sync::Mutex::new(Vec::new()));
        let messages_clone = Arc::clone(&messages);

        let log_fn = Arc::new(move |msg: &str| {
            messages_clone.lock().unwrap().push(msg.to_string());
        });

        let mut log_progress = LogProgress::new("msg", "final", true, Some(100), 10, log_fn);
        log_progress.inc(10);

        let logged = messages.lock().unwrap();
        assert!(logged[0].contains("["));
        assert!(logged[0].contains("]"));
    }
}
