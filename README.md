# Double Pendulum Simulation

## Description
This program simulates the motion of a double pendulum using numerical methods. It implements a Runge-Kutta integration technique to solve the equations of motion for a frictionless double pendulum system. The program also uses the `cpgplot` library to visualize the pendulum's movement.

## Features
- Simulates the motion of a double pendulum using physics-based equations.
- Implements numerical integration (Runge-Kutta) for accurate simulation.
- Allows user input for mass, rod length, and initial angles.
- Provides a graphical representation using `cpgplot`.

## How to run
This program requires the following:
- `cpgplot` (PGPLOT graphics library)
- C++ libraries (`iostream`, `cstdlib`, `cmath`, `algorithm`)

## Installation
1. Ensure `cpgplot` is installed on your system.
   - On Linux: Install with your package manager (e.g., `sudo apt install pgplot5`)
   - On macOS: Use Homebrew or build from source.
2. Compile the program using a C++ compiler:
   ```sh
   g++ -o double_pendulum dpend.cpp -lcpgplot -lpgplot -lm
   ```
3. Run the executable:
   ```sh
   ./double_pendulum
   ```

## Usage
- Run the program and choose between default or custom initial conditions.
- If using custom parameters, enter:
  - Mass values for both pendulums.
  - Rod lengths for both pendulums.
  - Initial angles in degrees.
- The simulation will display the motion of the pendulum.

## Numerical Methods
The program calculates the derivatives of the angular velocities using the equations of motion for a double pendulum and integrates them using the Runge-Kutta method. The equations are broken into smaller components for debugging and accuracy.

## License
This project is open-source and available for modification and redistribution.

## Author
Original source: AR (12/2020)

## Notes
- Ensure `cpgplot` is configured properly; otherwise, the program may fail to launch the graphical interface.
- The simulation assumes an idealized system with massless rods and frictionless pivots.
- The motion of a double pendulum is highly sensitive to initial conditions, demonstrating chaotic behavior in many cases.

