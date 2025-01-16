# PRESim: Pressure Simulation Tool

## Overview

PRESim is a specialized program designed for calculating external sea 
pressures $P_{ex}$ on the hull and internal pressures $P_{in}$ on 
the cargo hold of a Bulk Carrier. The calculations are based on 
the **Common Structural Rules (CSR) 2024** and use the concept 
of Equivalent Design Waves (EDWs), **only for strength assessment** 
of the vessel. The resulting pressure distributions can be exported
and imported into Finite Element Analysis (FEA) software such as
ANSYS, ABAQUS, or open-source alternatives.

## Features

- **<ins>External Pressure Calculation:</ins>** Calculates sea 
pressures acting on the hull using the `ExternalSeaPressureCalc` class.
- **<ins>Internal Pressure Calculation:</ins>** Calculates pressures
due to bulk cargo within the cargo hold using the 
`IntCargoPressureCalc` class.

## Equivalent Design Waves (EDWs)

The program incorporates the following EDWs:

- **<ins>HSM-1 & HSM-2:</ins>** Head Sea EDWs to minimize and 
maximize the vertical wave bending moment amidships, respectively.

- **<ins>HSA-1 & HSA-2:</ins>** Head Sea EDWs to maximize and
minimize the head sea vertical acceleration at the forward
perpendicular (FP), respectively.

- **<ins>FSM-1 & FSM-2:</ins>** Following Sea EDWs to minimize and
maximize the vertical wave bending moment amidships, respectively.

- **<ins>BSR-1 & BSR-2:</ins>** Beam Sea EDWs to minimize and
maximize roll motion downward and upward on the port side 
(or starboard side), respectively, with waves from the port 
side (or starboard side).

- **<ins>BSP-1 & BSP-2:</ins>** Beam Sea EDWs to maximize and
minimize the hydrodynamic pressure at the waterline amidships on
the port side (or starboard side), respectively.

- **<ins>OST-1 & OST-2:</ins>** Oblique Sea EDWs to minimize and
maximize the torsional moment at 0.25L from the aft end, with waves
from the port side (or starboard side), respectively.

- **<ins>OSA-1 & OSA-2:</ins>** Oblique Sea EDWs to maximize and
minimize the pitch acceleration, with waves from the port side
(or starboard side), respectively.

<div align="center">
    <img src="images/EDW.png" alt="Equivalent Design Waves">
</div>

## Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/Lalis98/DynLoadCalc
   ```

2. Navigate to the project directory:
   ```bash
   cd DynLoadCalc
   ```

3. Install required dependencies:
   ```bash
   pip install -r requirements.txt
   ```
   Alternatively, you can install the dependencies individually:

   ```bash
   pip install numpy matplotlib pandas
   ```
## Parameters for Classes and Methods

This section describes the parameters for each class and method in PRESim.

### ExternalSeaPressureCalc

#### Initialization Parameters:

- `coordinates (numpy.ndarray)`: 2D array of calculation points' coordinates `np.array([[x1, y1, z1],..., [xN, yN, zN]])`.
- `L (float)`: Rule Length (m).
- `B (float)`: Breadth Moulded (m).
- `TLC (float)`: Loading Condition Draught (m).
- `TSC (float)`: Scantling Draught (m).
- `Cb (float)`: Block Coefficient (-).
- `Bx_distribution (numpy.ndarray | float)`: 2D array of x-coordinates and corresponding breadth values at waterline.
- `kr (float)`: Roll Radius of Gyration (m).
- `GM (float)`: Metacentric Height (m).
- `bilge_keel (bool)`: Whether to include a bilge keel (`True`/`False`).

#### Methods and Parameters:

1. `calculate_load_case`:
Calculates the external pressure load case based on 
different ship scenarios and parameters, and saves in
the class the last pressure results.

    **<ins>Parameters:</ins>**
   - `base_case (str)`: Load case type (e.g.,`"HSM"`, `"FSM"`, `"HSA"`, etc.).
   - `sub_case (str)`: Sub-case (e.g., `"1"`, `"2"`, `"1P"`, etc.).
   - `fps (float)`: Strength assessment factor.
   - `load_scenario (str)`: Type of load scenario 
   (e.g.,`"Extreme Sea Loads"`, `"Ballast Water Exchange"`).
   - `fb (float)`: Head correction factor ($f_{β}$).


2. `plot_last_pressure_data`:
Plots the last calculated external pressure data.

    **<ins>Parameters:</ins>**
   - `color (str, optional)`: Color map for the plot 
   (e.g., `'coolwarm'`).
   - `size (int, optional)`: Size of the markers.
   - `file_path (str, optional)`: Path to save the plot. If not specified, the plot is not saved.
   - `file_name (str, optional)`: Name for the saved plot file. If not specified, the default name is used.
   - `show_plot (bool, optional)`: Whether to show the plot.


No parameters required.

Returns: Pressure distribution data.

calculate_wave_effects():

No parameters required.

Returns: Wave effects data.

IntCargoPressureCalc

Initialization Parameters:

param1 (type): Description of param1.

param2 (type): Description of param2.

param3 (type): Description of param3.

Methods and Parameters:

calculate_pressures():

No parameters required.

Returns: Pressure distribution data for the cargo.

calculate_cargo_effects():

No parameters required.

Returns: Cargo effects data.


## Usage

### Methods Overview

For both classes, the following attributes for methods are used for
`base_case` and `sub_case` calculations:

- For `base_case` of `"HSM"` or`"FSM"` or `"HSA"`:
  - `sub_case="1"`
  - `sub_case="2"`


- For `base_case` of `"BSR"` or`"BSP"` or `"OSA"` or `"OST"`:
  - `sub_case="1P"`
  - `sub_case="2P"`
  - `sub_case="1S"`
  - `sub_case="2S"`

### ExternalSeaPressureCalc

#### Example
Initializes with parameters for external pressure calculation.

```python
from calculations.external_pressure import *
import numpy as np

# Initialize the ExternalSeaPressureCalc class
external_analysis = ExternalSeaPressureCalc(
    coordinates=np.array(  # Coordinates of calculation points (example)
        [0, 0, 0],  # Point 1
        #   ...
        [170.0, 12.0, 1.0]  # Point N
    ),
    L=218.372,  # Rule Length (according to CSR) (m)
    B=32.240,  # Breadth Moulded (m)
    TLC=14.555,  # Loading Condition Draught (m)
    TSC=14.555,  # Scantling Draught (m)
    Cb=0.797,  # Block Coefficient (-)
    Bx_distribution=np.array(  # Distribution of Breadth at waterline relative to x-coordinate
        [[82.27, 158.71],  # x-coordinates along the ship's length (m)
         [32.24, 32.24]]  # Corresponding breadth values at those positions at the waterline (m)
    ),
    kr=0.4 * 32.240,  # Roll Radius of Gyration (m)
    GM=0.2 * 32.240,  # Metacentric height (m)
    bilge_keel=True  # Including Bilge Keel
)

# Calculate the external pressure
ext_pressure = external_analysis.calculate_load_case(
    base_case="HSM",                    # or "FSM", "HSA", "BSP", "BSR", "OST", "OSA"
    sub_case="1",                       # "1", "2", "1P", "2P", "1S", "2S"
    fps=1.0,                            # Strength Assessment factor
    load_scenario="Extreme Sea Loads",  # "Ballast Water Exchange" or "Extreme Sea Loads"
    fb=1.05,                            # Head Correction Factor fβ
    design_load="S+D"                   # Design load ("S+D", "S", "D)"
)

# Plot
external_analysis.plot_last_pressure_data(
    color='coolwarm',  # Color of colorbar
    size=5,  # Size of Markers
    file_path=None,  # Otherwise set path to save i.e. "C:/Users/user/Downloads/"
    file_name=None,  # Otherwise set name to save i.e. "HSM-1.png"
    show_plot=True  # Show Plot
)

# Save the data into a csv
external_analysis.save_to_csv(
    file_path="C:/Users/user/Downloads/",
    file_name="HSM-1.csv"
)

external_analysis.print_data_summary()  # Get data of the analysis

```

#### Output

The following is a sample plot generated by running
the above code for the **HSM-1** scenario:

<div align="center">
    <img src="images/hsm_1.png" alt="Dynamic Load Case Result">
</div>

### Internal Cargo Pressure Calculation

#### Example
The internal pressure calculation specifically models the effect
of bulk cargo within the cargo hold. It evaluates the static and
dynamic loads imparted on the internal structure, 
adhering to CSR rules for Bulk Carriers.

```python
from calculations.internal_pressure import *
import numpy as np

# Initialize the IntCargoPressureCalc class
internal_analysis = IntCargoPressureCalc(
    L=218.372,                  # Rule Length (according to CSR) (m)
    B=32.240,                   # Breadth Moulded (m)
    D=20.200,                   # Depth Moulded
    TLC=14.555,                 # Loading Condition Draught (m)
    TSC=14.555,                 # Scantling Draught (m)
    Cb=0.797,                   # Block Coefficient (-)
    kr=0.4 * 32.240,            # Roll Radius of Gyration (m)
    GM=0.2 * 32.240,            # Metacentric height (m)
    hdb=1.74,                   # Double Hull height (m)
    cargo_type="General",       # Allowed values are:
                                # 'General', 'Iron Ore', or 'Cement'
    cargo_density=1.535,        # Cargo Density (t/m³)
    cargo_density_effective=1.600,  # Effective Cargo Density (t/m³)
    load_scenario="Extreme Sea Loads",  # Allowed values are:
                                # "Extreme Sea Loads", "Water Exchange",
                                # "Accidental Flooded", "Harbour / Sheltered Water"
    bilge_keel=True             # Including Bilge Keel
)

# Calculate the internal pressure for Full Load
internal_analysis.calculate_load_case_full(
    base_case="HSM",        # or "FSM", "HSA", "BSP", "BSR", "OST", "OSA"
    sub_case="1",           # "1", "2", "1P", "2P", "1S", "2S"
    coordinates=np.array(   # Coordinates of calculation points (example)
        [80.0, 0, 20.0],    # Point 1
        #   ...
        [90.0, 12.0, 1.0]   # Point N
    ),
    fb=1.05,                # Head Correction factor fβ
    center_of_gravity=np.array(  # Center of Gravity of Cargo Hold (x, y, z) (m)
        [95.62, 0.0, 11.19]
    ),
    hHPU=13.46,             # Height hHPU (m)
    S0=125.840,             # Area S0 (m²)
    Bh=32.240,              # Breadth of Cargo Hold (m)
    Vhc=360.78,             # Volume of Hatch Coaming (m³)
    lh=25.480,              # Length of Cargo (m)
    angle_alpha=np.array([  # Angle α distribution related to z
        [5.68, 15.20, 90.0],  # 5.68 m < z < 15.20 m, angle is 90.0 deg
        [15.20, 23.90, 147.688]  # 15.20 m < z < 23.90 m, angle is 147.0 deg
    ]),
    shear_load_hopper=False,  # Geometry of hopper tank or Lower stool? (include Shear Load)
    shear_load_inner_bottom=False,  # Geometry of inner bottom? (include Shear Load)
    fdc=1.0,                # Dry Cargo Factor 
    design_load="S+D"       # Design Load
)

# Plot
internal_analysis.plot_last_pressure_data(
    color='coolwarm',   # Color of colorbar
    size=5,             # Size of Markers
    file_path=None,     # Otherwise set path to save i.e. "C:/Users/user/Downloads/"
    file_name=None,     # Otherwise set name to save i.e. "HSM-1.png"
    show_plot=True      # Show Plot
)

# Save the data into a csv
internal_analysis.save_to_csv(
    file_path="C:/Users/user/Downloads/",
    file_name="HSM-1.csv"
)

internal_analysis.print_data_summary()  # Get data of the analysis

```
#### Output

The following is a sample plot generated by running
the above code for the **HSM-1 Full Load** scenario:

<div align="center">
    <img src="images/cargo_hold_full_load.png" alt="Dynamic Load Case Result">
</div>

## Applications

Structural assessment of Bulk Carriers.

Importing pressure distributions into FEA software for further analysis.

## Contribution

Contributions and bug reports are welcome! Please submit issues or pull requests via the project repository.

## License

PRESim is licensed under the MIT License. See the LICENSE file for details.

## Contact

For inquiries or support, contact [michalis.leukioug1@gmail.com]().
