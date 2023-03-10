# NACA Aerofoil Generator
This is a MATLAB script that analyzes the lift and velocity of an airfoil. The script is intended to take in user inputs, generate a velocity field vector and streamlines for the inputted airfoil, and then calculate the lift coefficient, C_L. The script also has a special case when the airfoil is the NACA 2412. In this case, the script reads XFOIL data and plots C_L vs AoA for a range of AoA values and N values.

### Table of Contents
- [Usage](#usage)
- [Outputs](#outputs)
- [License](#license)

### Usage
To run the script, open MATLAB and navigate to the directory where the script is located. Then, simply run the script by typing airfoil_analysis into the command window.

The user is prompted for the following inputs:

```
airfoilcode - A four-digit code that refers to the NACA airfoil, e.g., 2412
U_infinity - The free-stream velocity in m/s, e.g., 15 m/s
AoA - The angle of attack in degrees, e.g., 10°
N - The number of panels to be used, e.g., 100

```

### Outputs
The script outputs the following:

```
C_L - The airfoil's lift coefficient
Plots with the velocity field vectors and streamlines
C_L vs AoA plot for the special case when the airfoil is a NACA 2412
```

### License
This project is licensed under the MIT License. See the [License](https://opensource.org/licenses/MIT) file for details.
