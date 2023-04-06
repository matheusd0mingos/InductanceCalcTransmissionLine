# Power Transmission Line Calculator

This program calculates various electrical properties of power transmission lines based on the parameters of the individual wires and the distances between them. The following functions are available:

## Wire Class
- `Wire`: Initializes a wire object with the following parameters:
  - `h`: wire height from ground (m)
  - `r`: wire radius (m)
  - `rfict`: wire radius of a fictitious conductor used in the inductance calculation (m)
  - `f`: frequency of operation (Hz)
  - `flecha`: wire sag (m)
  - `hreal`: effective wire height from ground (m)

## WireMatrix Class
- `WireMatrix`: Initializes a wire matrix object with the following parameters:
  - `wires`: a list of wire objects
  - `dist_matrix`: a matrix containing the distances between each wire

## GroundDistance Method
- `GroundDistance`: Calculates the ground distance between each wire in the wire matrix

## MonophasicInductance Method
- `MonophasicInductance`: Calculates the inductance between two wires in the wire matrix

## ThreeConductorsMonophaseInductance Method
- `ThreeConductorsMonophaseInductance`: Calculates the inductance between three wires in the wire matrix

## ThreePhaseInductanceWithoutGround Method
- `ThreePhaseInductanceWithoutGround`: Calculates the inductance between three wires in the wire matrix without considering ground effects

## ThreePhaseInductanceWithoutGroundTrancado Method
- `ThreePhaseInductanceWithoutGroundTrancado`: Calculates the inductance between three wires in the wire matrix without considering ground effects, using a simplified formula

## ThreePhaseInductandeGround Method
- `ThreePhaseInductandeGround`: Calculates the inductance between three wires in the wire matrix, considering ground effects

## ThreePhaseInductanceGroundTrancado Method
- `ThreePhaseInductanceGroundTrancado`: Calculates the inductance between three wires in the wire matrix, considering ground effects, using a simplified formula

## ThreePhaseInductanceforLighting Method
- `ThreePhaseInductanceforLighting`: Calculates the inductance between three wires in the wire matrix, considering ground effects, using a simplified formula for lighting systems

## ThreePhaseInductanceforLightingTransposto Method
- `ThreePhaseInductanceforLightingTransposto`: Calculates the inductance between three wires in the wire matrix, considering ground effects, using a simplified formula for transposed configurations

## ThreePhaseInductanceforLightingWithoutGround Method
- `ThreePhaseInductanceforLightingWithoutGround`: Calculates the inductance between three wires in the wire matrix, without considering ground effects, using a simplified formula for lighting systems

## PlotBizurado Method
- `PlotBizurado`: Plots the impedance variation of a three-phase line as a function of conductor diameter

## PlotBizuradoGuarni Method
- `PlotBizuradoGuarni`: Plots the impedance variation of a three-phase line as a function of conductor diameter, using the simplified Guarnieri formula

### Note

This program requires the `numpy` and `matplotlib` libraries.

