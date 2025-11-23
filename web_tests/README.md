# Web Tests - Selenium-based Integration Tests

This directory contains Selenium WebDriver tests for the boids WebAssembly web demo.

## Overview

These integration tests validate that the web demo:

- ✅ Loads without errors
- ✅ Properly initializes the simulation
- ✅ Responds to user interactions (buttons, keyboard shortcuts)
- ✅ Handles dynamic boid count changes without crashing
- ✅ Properly toggles behaviors
- ✅ Updates statistics in real-time
- ✅ Pauses and resumes correctly
- ✅ Resets to initial state

## Prerequisites

### 1. Build the WASM module

Before running tests, ensure the WASM module is built:

```bash
cd ../boids_wasm
wasm-pack build --target web --out-dir ../web/pkg
```

### 2. Install WebDriver

You need a WebDriver server running. Choose one:

#### Option A: ChromeDriver (Recommended)

1. Install ChromeDriver:
   ```bash
   # macOS (with Homebrew)
   brew install chromedriver

   # Linux (Ubuntu/Debian)
   sudo apt-get install chromium-chromedriver

   # Or download from: https://chromedriver.chromium.org/
   ```

2. Start ChromeDriver:
   ```bash
   chromedriver --port=4444
   ```

#### Option B: GeckoDriver (Firefox)

1. Install GeckoDriver:
   ```bash
   # macOS (with Homebrew)
   brew install geckodriver

   # Linux (Ubuntu/Debian)
   sudo apt-get install firefox-geckodriver

   # Or download from: https://github.com/mozilla/geckodriver
   ```

2. Start GeckoDriver:
   ```bash
   geckodriver --port=4444
   ```

## Running the Tests

### Run all tests

```bash
# From the project root
cargo test -p web_tests

# Or from this directory
cargo test
```

### Run a specific test

```bash
cargo test -p web_tests test_halve_boids_button
```

### Run tests with output

```bash
cargo test -p web_tests -- --nocapture
```

## Test Coverage

### Basic Functionality

- **test_demo_loads**: Verifies the web demo loads without errors and displays initial state
- **test_pause_resume**: Tests pause/resume functionality and frame count tracking

### Boid Count Management

- **test_double_boids_button**: Tests doubling boids via UI button
- **test_halve_boids_button**: Tests halving boids via UI button (critical regression test)
- **test_set_boid_count**: Tests setting exact boid count via input field
- **test_keyboard_double_boids**: Tests doubling boids with 'I' keyboard shortcut
- **test_keyboard_halve_boids**: Tests halving boids with 'D' keyboard shortcut

### Behavior Controls

- **test_behavior_toggles**: Tests behavior toggle buttons (separation, cohesion, alignment, wander)
- **test_reset_simulation**: Tests reset button functionality

### Stress Tests

- **test_rapid_boid_count_changes**: Regression test for spatial hash tracker bug - rapidly doubles and halves boids to ensure no crashes

## Test Architecture

Each test follows this pattern:

1. **Start local server**: Serves the web demo on `http://localhost:8765`
2. **Launch browser**: Creates a WebDriver instance (Chrome or Firefox)
3. **Navigate**: Goes to the local demo URL
4. **Wait for load**: Waits for WASM initialization (3 seconds)
5. **Perform actions**: Clicks buttons, sends keyboard input, etc.
6. **Validate**: Checks DOM elements, CSS, text content
7. **Cleanup**: Quits browser and stops server

## Debugging Tests

If tests fail:

1. **Check WebDriver is running**:
   ```bash
   # Should see output from chromedriver/geckodriver
   ```

2. **Verify WASM is built**:
   ```bash
   ls ../web/pkg/
   # Should see: boids_wasm_bg.wasm, boids_wasm.js, etc.
   ```

3. **Run with verbose output**:
   ```bash
   RUST_LOG=debug cargo test -p web_tests -- --nocapture
   ```

4. **Check browser manually**:
   ```bash
   # Start the demo manually and test in browser
   cd ../web
   ./serve.sh
   # Open http://localhost:8080 in browser
   ```

5. **Increase wait times**: If tests fail on CI or slower machines, increase sleep durations in tests

## Continuous Integration

For CI environments, use headless mode:

```bash
# ChromeDriver headless
chromedriver --port=4444 --headless

# GeckoDriver headless
geckodriver --port=4444 --headless
```

Or use a containerized approach with Selenium Grid:

```bash
docker run -d -p 4444:4444 selenium/standalone-chrome
```

## Known Issues

- **Port conflicts**: If port 8765 is in use, tests will fail. Change the port in `boids_web_demo.rs`
- **Timing sensitivity**: Tests use fixed sleep durations which may need adjustment for slower systems
- **Browser dependencies**: Tests require Chrome or Firefox to be installed

## Future Improvements

- [ ] Add screenshot capture on test failure
- [ ] Test parameter sliders (separation, cohesion, alignment coefficients)
- [ ] Test visual rendering (canvas content verification)
- [ ] Add performance benchmarks (FPS tracking)
- [ ] Test keyboard shortcuts for all behaviors (1-4 keys)
- [ ] Test mobile viewport responsiveness
- [ ] Add tests for error states (invalid input values)
