# Generator of exact spectrum and Green's function

## Prerequists

- Python >= 3.6
- Python packages
    - numpy
    - scipy
    - toml

## Usage

`--help` option shows the usage message:

``` bash
$ python3 ./gen_Gtau_peaks.py --help
usage: gen_Gtau_peaks.py [-h] [-n NSAMPLES] [-q] input

To generate Gtau with white noise from a spectrum

positional arguments:
  input        Input TOML file

optional arguments:
  -h, --help   show this help message and exit
  -n NSAMPLES  the number of samples to be generated (default=1)
  -q           do quietly
```

After preparing an input parameter file, say `input.toml`,

``` bash
python3 ./gen_Gtau_peaks.py input.toml
```

generates the exact spectrum `Gtau.in.dos` and the Green's function `Gtau.in` under the output directory, `work` (the name can be changed by the input file).

When `NSAMPLES` is greater than 1, independent `NSAMPLES` Green's functions with random noise, `Gtau_0.in`, `Gtau_1.in`, and so on, are generated.

## Files

### Input file

- `[parameter]` section
    - `stat`: string
        - Statistics of the system
        - "fermion" or "boson"
        - When "boson", the anti-symmetry of the spectrum, `ρ(-ω) = -ρ(ω)`, is guaranteed.
    - `beta`: float
        - Inverse temperature
    - `Ntau`: int
        - The number of points of imaginary-time
    - `omega_min`: float
        - The minimum value of real-frequency
    - `omega_max`: float
        - The maximum value of real-frequency
    - `Nomega`: int
        - The number of points of real frequency
    - `noise`: float
        - The amplitude of the random noise in the Green's function
    - `seed_base`: int
        - The origin of the seed of the random number generator
        - Default: 314159
    - `seed_step`: int
        - The step of the seed of the random number generator
        - Default: 11
        - `seed = seed_base + seed_step * i_sample`
    - `output_dir`: string
        - The name of the directory where output files are saved
        - Default: "work"
- `[[peaks]]` sections
    - `peak_type`: string
        - Type of a peak
        - "G": Gaussian
        - "L": Lorentzian
        - "S": Semi-circle
    - `position`: float
    - `width`: float
    - `intensity`: float

### Output files

- `Gtau.in.dos`
    - The exact spectrum
    - The first column means frequency, ω
    - The second column means spectrum, ρ(ω)
- `Gtau.in`
    - A Green's function with noise
    - The first column means imaginary-time, τ/β
    - The second column means Green's function with noise, G(τ)
    - The third column means Green's function without noise, G^{exact}(τ)
