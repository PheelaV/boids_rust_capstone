/// Selenium-based integration tests for the boids WASM web demo
///
/// These tests validate that the web demo:
/// - Loads without errors
/// - Properly initializes the simulation
/// - Responds to user interactions (buttons, keyboard shortcuts)
/// - Handles dynamic boid count changes without crashing
/// - Properly toggles behaviors
/// - Updates statistics in real-time
///
/// ## Prerequisites:
///
/// 1. WebDriver must be installed and running:
///    - ChromeDriver for Chrome: https://chromedriver.chromium.org/
///    - GeckoDriver for Firefox: https://github.com/mozilla/geckodriver
///
///    Start the driver before running tests:
///    ```bash
///    # For Chrome
///    chromedriver --port=4444
///
///    # For Firefox
///    geckodriver --port=4444
///    ```
///
/// 2. The WASM module must be built:
///    ```bash
///    cd ../boids_wasm
///    wasm-pack build --target web --out-dir ../web/pkg
///    ```
///
/// ## Running the tests:
///
/// ```bash
/// cargo test -p web_tests
/// ```
///
/// ## Test Structure:
///
/// Each test:
/// 1. Starts a local web server serving the web demo
/// 2. Launches a browser via Selenium
/// 3. Navigates to the demo
/// 4. Performs interactions
/// 5. Validates expected behavior
/// 6. Cleans up (closes browser, stops server)

use std::time::Duration;
use thirtyfour::prelude::*;
use tokio::time::sleep;
use warp::Filter;

/// Start a local web server serving the web demo
/// Returns the server handle and the URL to access it
async fn start_test_server() -> (tokio::task::JoinHandle<()>, String) {
    let port = 8765; // Use a different port than the dev server
    let url = format!("http://localhost:{}", port);

    // Serve static files from the web directory
    let web_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .join("web");

    let routes = warp::fs::dir(web_dir);

    let server = tokio::spawn(async move {
        warp::serve(routes).run(([127, 0, 0, 1], port)).await;
    });

    // Give the server a moment to start
    sleep(Duration::from_millis(500)).await;

    (server, url)
}

/// Create a new WebDriver instance
/// Tries Chrome first, falls back to Firefox
async fn create_driver() -> WebDriverResult<WebDriver> {
    let caps = DesiredCapabilities::chrome();

    // Try Chrome first
    match WebDriver::new("http://localhost:4444", caps).await {
        Ok(driver) => Ok(driver),
        Err(_) => {
            // Fall back to Firefox
            let caps = DesiredCapabilities::firefox();
            WebDriver::new("http://localhost:4444", caps).await
        }
    }
}

/// Helper to wait for an element and get its text
async fn get_element_text(driver: &WebDriver, selector: &str) -> WebDriverResult<String> {
    let element = driver.query(By::Css(selector)).first().await?;
    element.text().await
}

/// Helper to click a button by ID
async fn click_button(driver: &WebDriver, button_id: &str) -> WebDriverResult<()> {
    let button = driver.find(By::Id(button_id)).await?;
    button.click().await
}

/// Test: Web demo loads without errors
#[tokio::test]
async fn test_demo_loads() -> WebDriverResult<()> {
    let (_server, url) = start_test_server().await;
    let driver = create_driver().await?;

    // Navigate to the demo
    driver.goto(&url).await?;

    // Wait for WASM to load (loading message should disappear)
    sleep(Duration::from_secs(3)).await;

    // Check that the loading message is hidden
    let loading = driver.find(By::Id("loading")).await?;
    let display_style = loading.css_value("display").await?;
    assert_eq!(display_style, "none", "Loading message should be hidden");

    // Check that the canvas is visible
    let canvas = driver.find(By::Id("canvas")).await?;
    let canvas_display = canvas.css_value("display").await?;
    assert_eq!(canvas_display, "block", "Canvas should be visible");

    // Check that initial boid count is displayed
    let boid_count_text = get_element_text(&driver, "#boid-count").await?;
    let boid_count: i32 = boid_count_text.parse().unwrap();
    assert_eq!(boid_count, 200, "Initial boid count should be 200");

    driver.quit().await?;
    Ok(())
}

/// Test: Doubling boids via button
#[tokio::test]
async fn test_double_boids_button() -> WebDriverResult<()> {
    let (_server, url) = start_test_server().await;
    let driver = create_driver().await?;

    driver.goto(&url).await?;
    sleep(Duration::from_secs(3)).await;

    // Get initial count
    let initial_count_text = get_element_text(&driver, "#boid-count").await?;
    let initial_count: i32 = initial_count_text.parse().unwrap();

    // Click double button
    click_button(&driver, "double-boids").await?;

    // Wait for update
    sleep(Duration::from_millis(500)).await;

    // Check new count
    let new_count_text = get_element_text(&driver, "#boid-count").await?;
    let new_count: i32 = new_count_text.parse().unwrap();

    assert_eq!(new_count, initial_count * 2, "Boid count should double");

    driver.quit().await?;
    Ok(())
}

/// Test: Halving boids via button (tests the bug fix)
#[tokio::test]
async fn test_halve_boids_button() -> WebDriverResult<()> {
    let (_server, url) = start_test_server().await;
    let driver = create_driver().await?;

    driver.goto(&url).await?;
    sleep(Duration::from_secs(3)).await;

    // Get initial count
    let initial_count_text = get_element_text(&driver, "#boid-count").await?;
    let initial_count: i32 = initial_count_text.parse().unwrap();

    // Click halve button
    click_button(&driver, "halve-boids").await?;

    // Wait for update
    sleep(Duration::from_millis(500)).await;

    // Check new count
    let new_count_text = get_element_text(&driver, "#boid-count").await?;
    let new_count: i32 = new_count_text.parse().unwrap();

    assert_eq!(new_count, initial_count / 2, "Boid count should halve");

    // Critical: Ensure simulation continues running without crashing
    // Get FPS to confirm it's still updating
    sleep(Duration::from_secs(1)).await;
    let fps_text = get_element_text(&driver, "#fps").await?;
    let fps: i32 = fps_text.parse().unwrap_or(0);
    assert!(fps > 0, "Simulation should still be running after halving boids");

    driver.quit().await?;
    Ok(())
}

/// Test: Setting exact boid count via input field
#[tokio::test]
async fn test_set_boid_count() -> WebDriverResult<()> {
    let (_server, url) = start_test_server().await;
    let driver = create_driver().await?;

    driver.goto(&url).await?;
    sleep(Duration::from_secs(3)).await;

    // Find the input field
    let input = driver.find(By::Id("boid-count-input")).await?;

    // Clear and set new value
    input.clear().await?;
    input.send_keys("150").await?;

    // Click set button
    click_button(&driver, "set-boid-count").await?;

    // Wait for update
    sleep(Duration::from_millis(500)).await;

    // Check count
    let count_text = get_element_text(&driver, "#boid-count").await?;
    let count: i32 = count_text.parse().unwrap();

    assert_eq!(count, 150, "Boid count should be set to 150");

    driver.quit().await?;
    Ok(())
}

/// Test: Keyboard shortcut for doubling boids (I key)
#[tokio::test]
async fn test_keyboard_double_boids() -> WebDriverResult<()> {
    let (_server, url) = start_test_server().await;
    let driver = create_driver().await?;

    driver.goto(&url).await?;
    sleep(Duration::from_secs(3)).await;

    // Get initial count
    let initial_count_text = get_element_text(&driver, "#boid-count").await?;
    let initial_count: i32 = initial_count_text.parse().unwrap();

    // Send 'i' key to double boids
    let body = driver.find(By::Tag("body")).await?;
    body.send_keys("i").await?;

    // Wait for update
    sleep(Duration::from_millis(500)).await;

    // Check new count
    let new_count_text = get_element_text(&driver, "#boid-count").await?;
    let new_count: i32 = new_count_text.parse().unwrap();

    assert_eq!(new_count, initial_count * 2, "Boid count should double with 'i' key");

    driver.quit().await?;
    Ok(())
}

/// Test: Keyboard shortcut for halving boids (D key)
#[tokio::test]
async fn test_keyboard_halve_boids() -> WebDriverResult<()> {
    let (_server, url) = start_test_server().await;
    let driver = create_driver().await?;

    driver.goto(&url).await?;
    sleep(Duration::from_secs(3)).await;

    // Get initial count
    let initial_count_text = get_element_text(&driver, "#boid-count").await?;
    let initial_count: i32 = initial_count_text.parse().unwrap();

    // Send 'd' key to halve boids
    let body = driver.find(By::Tag("body")).await?;
    body.send_keys("d").await?;

    // Wait for update
    sleep(Duration::from_millis(500)).await;

    // Check new count
    let new_count_text = get_element_text(&driver, "#boid-count").await?;
    let new_count: i32 = new_count_text.parse().unwrap();

    assert_eq!(new_count, initial_count / 2, "Boid count should halve with 'd' key");

    driver.quit().await?;
    Ok(())
}

/// Test: Behavior toggle buttons work
#[tokio::test]
async fn test_behavior_toggles() -> WebDriverResult<()> {
    let (_server, url) = start_test_server().await;
    let driver = create_driver().await?;

    driver.goto(&url).await?;
    sleep(Duration::from_secs(3)).await;

    // All behaviors should start as ON (active class)
    let separation_btn = driver.find(By::Id("toggle-separation")).await?;
    let class_attr = separation_btn.attr("class").await?.unwrap_or_default();
    assert!(class_attr.contains("active"), "Separation should start active");

    // Click to toggle off
    separation_btn.click().await?;
    sleep(Duration::from_millis(100)).await;

    // Should no longer have active class
    let class_attr_after = separation_btn.attr("class").await?.unwrap_or_default();
    assert!(!class_attr_after.contains("active"), "Separation should be inactive after click");

    driver.quit().await?;
    Ok(())
}

/// Test: Reset button works
#[tokio::test]
async fn test_reset_simulation() -> WebDriverResult<()> {
    let (_server, url) = start_test_server().await;
    let driver = create_driver().await?;

    driver.goto(&url).await?;
    sleep(Duration::from_secs(3)).await;

    // Modify the simulation first
    click_button(&driver, "double-boids").await?;
    sleep(Duration::from_millis(500)).await;

    let doubled_count_text = get_element_text(&driver, "#boid-count").await?;
    let doubled_count: i32 = doubled_count_text.parse().unwrap();
    assert_eq!(doubled_count, 400, "Count should be doubled");

    // Reset
    click_button(&driver, "reset").await?;
    sleep(Duration::from_millis(500)).await;

    // Should be back to initial count
    let reset_count_text = get_element_text(&driver, "#boid-count").await?;
    let reset_count: i32 = reset_count_text.parse().unwrap();
    assert_eq!(reset_count, 200, "Count should reset to initial 200");

    driver.quit().await?;
    Ok(())
}

/// Test: Pause/Resume functionality
#[tokio::test]
async fn test_pause_resume() -> WebDriverResult<()> {
    let (_server, url) = start_test_server().await;
    let driver = create_driver().await?;

    driver.goto(&url).await?;
    sleep(Duration::from_secs(3)).await;

    // Get initial frame count
    let initial_frame_text = get_element_text(&driver, "#frame").await?;
    let initial_frame: i32 = initial_frame_text.parse().unwrap();

    // Click pause
    let pause_btn = driver.find(By::Id("pause")).await?;
    pause_btn.click().await?;

    // Button text should change to "Resume"
    let button_text = pause_btn.text().await?;
    assert_eq!(button_text, "Resume", "Button should say 'Resume' when paused");

    // Wait and verify frame count doesn't change
    sleep(Duration::from_secs(1)).await;
    let paused_frame_text = get_element_text(&driver, "#frame").await?;
    let paused_frame: i32 = paused_frame_text.parse().unwrap();
    assert_eq!(paused_frame, initial_frame, "Frame count should not increase when paused");

    // Resume
    pause_btn.click().await?;
    let resumed_button_text = pause_btn.text().await?;
    assert_eq!(resumed_button_text, "Pause", "Button should say 'Pause' when running");

    // Frame count should increase
    sleep(Duration::from_secs(1)).await;
    let resumed_frame_text = get_element_text(&driver, "#frame").await?;
    let resumed_frame: i32 = resumed_frame_text.parse().unwrap();
    assert!(resumed_frame > paused_frame, "Frame count should increase after resume");

    driver.quit().await?;
    Ok(())
}

/// Test: Multiple rapid boid count changes don't crash (regression test)
#[tokio::test]
async fn test_rapid_boid_count_changes() -> WebDriverResult<()> {
    let (_server, url) = start_test_server().await;
    let driver = create_driver().await?;

    driver.goto(&url).await?;
    sleep(Duration::from_secs(3)).await;

    // Rapidly double and halve multiple times
    for _ in 0..5 {
        click_button(&driver, "double-boids").await?;
        sleep(Duration::from_millis(100)).await;
        click_button(&driver, "halve-boids").await?;
        sleep(Duration::from_millis(100)).await;
    }

    // Verify simulation is still running
    sleep(Duration::from_secs(1)).await;
    let fps_text = get_element_text(&driver, "#fps").await?;
    let fps: i32 = fps_text.parse().unwrap_or(0);
    assert!(fps > 0, "Simulation should still be running after rapid changes");

    driver.quit().await?;
    Ok(())
}
