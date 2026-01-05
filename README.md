# The Ecophysiometer: An Open-Source Multipurpose Environmental Sensor

The **Ecophysiometer** is a low-cost (\~US\$128), modular environmental monitoring platform designed to acquire rich optical data and interpret it through chemometric calibration. Developed at the Universidad Nacional de Córdoba (UNC), CONICET and the Technical University of Rheinland-Pfalz, Landau (RPTU) this project shifts analytical complexity from hardware to software, enabling multi-analyte quantification with flexible, 3D-printable instrumentation.

## 📁 Repository Structure

The project is organized to support both hardware fabrication and data analysis:

```         
/Design/Hardware: Engineering files for the flow cell, modular housing, and sensor mounts in FreeCAD and STL formats. Includes KiCAD project files for the constant-current PCB.
```

/Design/Software: Python source code including the core Measurement.py, ImageProcessor.py, and SpectralData.py layers .

/Instructions: Step-by-step assembly guides, including the critical IR filter removal tutorial for the Raspberry Pi Camera .

/Characterisation: Raw data and R scripts used to validate signal stability, linear range, and manufacturing variability.

/Implementation examples: Datasets and R models for the primary use cases: Turbidity and Phytoplankton differential quantification.

## 🛠️ Key Technical Features

```         
Software-Defined Specificity: Analyte detection is determined by chemometric calibration (PLSR, PCA) rather than fixed hardware filters.
```

High-Certainty Signal Processing: The pipeline aligns 450 vertical rows and applies 5-pixel binning, deriving each data point from the average of 2,250 measurements to maximize the signal-to-noise ratio .

Reproducible Hardware: Designed for manufacturing with minimal technical knowledge using locally available electronic materials and 3D printing in PETG .

## 🚀 Getting Started

### 1. Hardware Fabrication

```         
3D Printing: Print parts found in /Design/Hardware. Black PETG is recommended to minimize internal reflections.
```

[Circuitry]{.underline}: Fabricate the PCB using the provided KiCAD files or the CNC-ready .nc files in the Circuit Board directory.

[Sensor Modification]{.underline}: Follow Instructions/IR_REMOVAL.md to remove the IR cut-off filter from the Raspberry Pi Camera V1.3.

### 2. Software Installation

Clone this repository and install the dependencies: Bash

git clone <https://github.com/australScientist/PEF1.git> cd PEF1 pip install -r Design/Software/PEF0/requirements.txt

### 3. Basic Operation

[Calibration]{.underline}: Run test.py to identify the linear integration time for your specific hardware assembly .

[Measurement]{.underline}: Launch the graphical interface: Bash

```         
python3 Design/Software/PEF0/src/GUI2.py
```

