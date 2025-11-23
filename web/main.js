// Import the WebAssembly module
// Note: This assumes wasm-pack has been run and output to ./pkg/
import init, { WasmSimulation } from './pkg/boids_wasm.js';

let simulation = null;
let canvas = null;
let ctx = null;
let animationId = null;
let isPaused = false;
let lastTime = performance.now();
let fps = 0;
let frameCount = 0;

// Initialize the WebAssembly module and start the simulation
async function run() {
    try {
        // Initialize the WASM module
        await init();

        // Hide loading message
        document.getElementById('loading').style.display = 'none';

        // Show canvas
        canvas = document.getElementById('canvas');
        canvas.style.display = 'block';
        ctx = canvas.getContext('2d');

        // Create simulation with configuration
        const config = {
            init_boids: 200,
            window_width: 800,
            window_height: 600,
            separation_coefficient: 1.5,
            cohesion_coefficient: 1.0,
            alignment_coefficient: 1.0,
            max_speed: 4.0,
            min_speed: 2.0,
            sensory_distance: 60.0,
            wander_on: true,
            wander_coefficient: 0.5
        };

        simulation = WasmSimulation.from_config(config);
        console.log('Simulation initialized with', simulation.get_boid_count(), 'boids');

        // Set up controls
        setupControls();

        // Set up keyboard controls
        setupKeyboardControls();

        // Start animation loop
        animate();

    } catch (error) {
        console.error('Failed to initialize:', error);
        showError(`Failed to initialize WebAssembly: ${error.message}`);
    }
}

function showError(message) {
    const loading = document.getElementById('loading');
    loading.className = 'error';
    loading.innerHTML = `<strong>Error:</strong> ${message}`;
}

function setupControls() {
    // Behavior toggles
    document.getElementById('toggle-separation').addEventListener('click', function() {
        simulation.toggle_separation();
        this.classList.toggle('active');
        updateToggleText(this, 'Separation');
    });

    document.getElementById('toggle-cohesion').addEventListener('click', function() {
        simulation.toggle_cohesion();
        this.classList.toggle('active');
        updateToggleText(this, 'Cohesion');
    });

    document.getElementById('toggle-alignment').addEventListener('click', function() {
        simulation.toggle_alignment();
        this.classList.toggle('active');
        updateToggleText(this, 'Alignment');
    });

    document.getElementById('toggle-wander').addEventListener('click', function() {
        simulation.toggle_wander();
        this.classList.toggle('active');
        updateToggleText(this, 'Wander');
    });

    // Parameter sliders
    document.getElementById('separation').addEventListener('input', function() {
        const value = parseFloat(this.value);
        simulation.set_separation_coefficient(value);
        document.getElementById('separation-value').textContent = value.toFixed(1);
    });

    document.getElementById('cohesion').addEventListener('input', function() {
        const value = parseFloat(this.value);
        simulation.set_cohesion_coefficient(value);
        document.getElementById('cohesion-value').textContent = value.toFixed(1);
    });

    document.getElementById('alignment').addEventListener('input', function() {
        const value = parseFloat(this.value);
        simulation.set_alignment_coefficient(value);
        document.getElementById('alignment-value').textContent = value.toFixed(1);
    });

    document.getElementById('max-speed').addEventListener('input', function() {
        const value = parseFloat(this.value);
        simulation.set_max_speed(value);
        document.getElementById('max-speed-value').textContent = value.toFixed(1);
    });

    // Action buttons
    document.getElementById('reset').addEventListener('click', () => {
        simulation.reset();
        console.log('Simulation reset');
    });

    document.getElementById('pause').addEventListener('click', function() {
        isPaused = !isPaused;
        this.textContent = isPaused ? 'Resume' : 'Pause';
        if (!isPaused && !animationId) {
            animate();
        }
    });
}

function setupKeyboardControls() {
    document.addEventListener('keydown', (event) => {
        const key = event.key.toLowerCase();

        switch(key) {
            case ' ': // Spacebar - pause/resume
                event.preventDefault();
                const pauseBtn = document.getElementById('pause');
                pauseBtn.click();
                break;

            case 'r': // Reset simulation
                event.preventDefault();
                simulation.reset();
                console.log('Simulation reset (keyboard)');
                break;

            case 'c': // Toggle controls panel visibility
                event.preventDefault();
                const controls = document.getElementById('controls');
                controls.style.display = controls.style.display === 'none' ? 'block' : 'none';
                break;

            case '1': // Toggle alignment
                event.preventDefault();
                document.getElementById('toggle-alignment').click();
                break;

            case '2': // Toggle cohesion
                event.preventDefault();
                document.getElementById('toggle-cohesion').click();
                break;

            case '3': // Toggle separation
                event.preventDefault();
                document.getElementById('toggle-separation').click();
                break;

            case '4': // Toggle wander
                event.preventDefault();
                document.getElementById('toggle-wander').click();
                break;

            case '+':
            case '=': // Increase max speed
                event.preventDefault();
                const speedSlider = document.getElementById('max-speed');
                const newSpeed = Math.min(10, parseFloat(speedSlider.value) + 0.5);
                speedSlider.value = newSpeed;
                speedSlider.dispatchEvent(new Event('input'));
                break;

            case '-':
            case '_': // Decrease max speed
                event.preventDefault();
                const speedSlider2 = document.getElementById('max-speed');
                const newSpeed2 = Math.max(1, parseFloat(speedSlider2.value) - 0.5);
                speedSlider2.value = newSpeed2;
                speedSlider2.dispatchEvent(new Event('input'));
                break;

            case 'h':
            case '?': // Show keyboard help
                event.preventDefault();
                showKeyboardHelp();
                break;
        }
    });

    console.log('Keyboard controls enabled. Press H or ? for help.');
}

function showKeyboardHelp() {
    const helpText = `
KEYBOARD SHORTCUTS:
  Space    - Pause/Resume
  R        - Reset simulation
  C        - Toggle controls panel
  1        - Toggle Alignment
  2        - Toggle Cohesion
  3        - Toggle Separation
  4        - Toggle Wander
  +/=      - Increase max speed
  -/_      - Decrease max speed
  H/?      - Show this help
`;
    console.log(helpText);
    alert(helpText);
}

function updateToggleText(button, name) {
    const isActive = button.classList.contains('active');
    button.textContent = `${name}: ${isActive ? 'ON' : 'OFF'}`;
}

function animate() {
    if (!isPaused) {
        // Update simulation
        simulation.update();

        // Get boids and render
        const boids = simulation.get_boids();
        render(boids);

        // Update stats
        updateStats();
    }

    animationId = requestAnimationFrame(animate);
}

function render(boids) {
    // Clear canvas
    ctx.fillStyle = '#0a1929';
    ctx.fillRect(0, 0, canvas.width, canvas.height);

    // Center the coordinate system
    ctx.save();
    ctx.translate(canvas.width / 2, canvas.height / 2);

    // Draw each boid
    for (const boid of boids) {
        const x = boid.x;
        const y = boid.y;
        const vx = boid.vx;
        const vy = boid.vy;

        // Calculate angle from velocity
        const angle = Math.atan2(vy, vx);

        // Draw boid as a triangle pointing in direction of velocity
        ctx.save();
        ctx.translate(x, y);
        ctx.rotate(angle);

        // Boid body (triangle)
        ctx.beginPath();
        ctx.moveTo(8, 0);  // nose
        ctx.lineTo(-6, -4); // top wing
        ctx.lineTo(-6, 4);  // bottom wing
        ctx.closePath();

        // Color based on velocity magnitude
        const speed = Math.sqrt(vx * vx + vy * vy);
        const hue = (speed / 6) * 120; // 0 (red) to 120 (green)
        ctx.fillStyle = `hsl(${hue}, 70%, 60%)`;
        ctx.fill();

        // Outline
        ctx.strokeStyle = '#ffffff';
        ctx.lineWidth = 1;
        ctx.stroke();

        ctx.restore();
    }

    ctx.restore();
}

function updateStats() {
    const now = performance.now();
    const delta = now - lastTime;

    frameCount++;
    if (frameCount % 10 === 0) {
        fps = Math.round(1000 / delta);
        document.getElementById('fps').textContent = fps;
    }

    lastTime = now;

    const stats = simulation.get_stats();
    document.getElementById('boid-count').textContent = stats.boid_count;
    document.getElementById('frame').textContent = stats.frame_count;
}

// Start the application
run();
