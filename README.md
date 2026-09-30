## interface_stability 

This package and associated scripts are designed to analyze the interface stability. 

## Citing

If you use this library, please consider citing the following papers:

    Zhu, Yizhou, Xingfeng He, and Yifei Mo. "Origin of outstanding stability in the lithium solid electrolyte materials: insights from thermodynamic analyses based on first-principles calculations." ACS applied materials & interfaces 7.42 (2015): 23685-23693.
    
DOI: 10.1021/acsami.5b07517

    Zhu, Yizhou, Xingfeng He, and Yifei Mo. "First principles study on electrochemical and chemical stability of solid electrolyte–electrode interfaces in all-solid-state Li-ion batteries." Journal of Materials Chemistry A 4.9 (2016): 3253-3266.

DOI: 10.1039/C5TA08574H

    Han, Fudong, et al. "Electrochemical stability of Li10GeP2S12 and Li7La3Zr2O12 solid electrolytes." Advanced Energy Materials 6.8 (2016).
    
DOI: 10.1002/aenm.201501590 

    Ping Ong, Shyue, et al. "Li− Fe− P− O2 phase diagram from first principles calculations." Chemistry of Materials 20.5 (2008): 1798-1807.

DOI: 10.1021/cm702327g

    Mo, Yifei, Shyue Ping Ong, and Gerbrand Ceder. "First principles study of the Li10GeP2S12 lithium super ionic conductor material." Chemistry of Materials 24.1 (2011): 15-17.

DOI: 10.1021/cm203303y

    Jain, Anubhav, et al. "Commentary: The Materials Project: A materials genome approach to accelerating materials innovation." Apl Materials 1.1 (2013): 011002.

DOI:10.1063/1.4812323

## Install

This package needs Python 3.10 or newer. It uses the current [pymatgen](https://pymatgen.org/) and the
next-gen Materials Project API client ([mp-api](https://github.com/materialsproject/api)).

1. Clone this package from GitHub and install it. The dependencies are installed automatically.

    ```bash
    $ git clone https://github.com/mogroupumd/interface_stability.git
    $ cd interface_stability
    $ python -m venv .venv && source .venv/bin/activate
    $ pip install -e .
    ```

2. Set up your Materials Project API key. This enables you to fetch data from the Materials Project database.
   Your API key is on your dashboard at https://next-gen.materialsproject.org/api (login required).
   Legacy (pre-2022) keys do not work.

   Either set an environment variable:
   ```bash
   $ export MP_API_KEY="your API key"
   ```
   or put the key in your pymatgen settings file, ~/.pmgrc.yaml (create it if it does not exist):
   ```bash
   PMG_MAPI_KEY: your API key
   ```

3. Optional: cache downloaded entries on disk, so repeated runs are faster and can work offline.
   Set `IFS_CACHE_DIR` to a folder (or `PMG_PD_PRELOAD_PATH` in ~/.pmgrc.yaml). Delete the cached
   files to pick up changes to the Materials Project database.

4. The install creates two commands, phase_stability and pseudo_binary. Try them and read the documentation:

    ```bash
    $ phase_stability -h
    $ pseudo_binary -h
    ```

5. Run the tests. The offline tests use a small made-up data set; the live tests run only when an API key is set.

    ```bash
    $ pip install pytest
    $ pytest interface_stability/tests
    ```

### Notes on the Materials Project data

* Energies come from the `GGA_GGA+U` thermo type (GGA/GGA+U with MaterialsProject2020Compatibility corrections),
  the same method as the papers above. Use `--thermo-type GGA_GGA+U_R2SCAN` (before the sub-command) to use
  the mixed GGA/GGA+U/r2SCAN data shown on the Materials Project website, e.g.
  `phase_stability --thermo-type GGA_GGA+U_R2SCAN stability Li3PS4`.
* The Materials Project data has been updated since the examples below were made, so your numbers will differ somewhat.

## Usage

There are two executable python scripts, both work with a few sub-commands options. 
You can always use -h to see the help information. 
For example,

```bash
$ phase_stability -h
$ phase_stability evolution -h
```

The `evolution` and `plotvc` sub-commands make a figure. Add `--save plot.png` to save it to a file instead of
showing it, or `--noplot` to skip it.

### 1. scripts/phase_stability.py

**phase_stability stability composition**

This gives the phase equilibria of a given composition.

```bash
$ phase_stability stability Li10GeP2S12
------------------------------------------------------------
Reduced formula of the given composition: Li10Ge(PS6)2
Calculated phase equilibria: Li4GeS4    Li3PS4
Li10Ge(PS6)2 -> 2 Li3PS4 + Li4GeS4
------------------------------------------------------------
```

**phase_stability mu composition open_element chemical_potential**

This gives the phase equilibria of a given composition under given chemical potential.

Note: Chemical potential is always referenced to elementary phases and in eV units.

```bash
$ phase_stability mu Li3PS4 Li -5
------------------------------------------------------------
Reduced formula of the given composition: Li3PS4
Open element : Li
Chemical potential: -5 eV referenced to pure phase
------------------------------------------------------------
Reaction:Li3PS4 -> 3 Li + 0.5 S + 0.5 P2S7
Reaction energy: -13.465 eV per Li3PS4
------------------------------------------------------------
```
**phase_stability evolution [-posmu] composition open_element**

This gives the evolution profile with changing chemical potential of an open element.

A figure of reaction energy will also be generated

```bash
$ phase_stability evolution Li3PS4 Li
------------------------------------------------------------
Reduced formula of the given composition: Li3PS4

 === Evolution Profile ===
mu_high (eV) mu_low (eV)  d(n_Li) Phase equilibria                   Reaction                 
    0.00        -0.87      8.00       Li2S, Li3P                Li3PS4 + 8 Li -> Li3P + 4 Li2S
   -0.87        -0.93      6.00        Li2S, LiP                 Li3PS4 + 6 Li -> LiP + 4 Li2S
   -0.93        -1.17      5.43      Li2S, Li3P7    Li3PS4 + 5.429 Li -> 0.1429 Li3P7 + 4 Li2S
   -1.17        -1.30      5.14       Li2S, LiP7     Li3PS4 + 5.143 Li -> 0.1429 LiP7 + 4 Li2S
   -1.30        -1.72      5.00          Li2S, P                   Li3PS4 + 5 Li -> 4 Li2S + P
   -1.72        -2.36      0.00           Li3PS4                              Li3PS4 -> Li3PS4
   -2.36        -3.74     -2.88       LiS4, P2S7    Li3PS4 -> 2.875 Li + 0.125 LiS4 + 0.5 P2S7
   -3.74         -inf     -3.00          P2S7, S             Li3PS4 -> 3 Li + 0.5 S + 0.5 P2S7

 === Reaction energy ===
miu_Li (eV)  Rxn energy (eV/atom)
   0.00             -1.42        
  -0.87             -0.55        
  -0.93             -0.50        
  -1.17             -0.35        
  -1.30             -0.26        
  -1.72              0.00        
  -2.36             -0.00        
  -3.74             -0.50        
  -3.94             -0.57
Note:
Chemical potential referenced to element phase.
Reaction energy is normalized to per atom of the given composition.
------------------------------------------------------------
```
**phase_stability plotvc [-posmu] [-v VALENCE] composition open_element**

Generate a figure of voltage profile (and display all raw data). 

```bash
$ phase_stability plotvc Li3PS4 Li
d n(Li)  Voltage ref. to Li (V)
-3.00             3.74         
-2.88             3.74         
-2.88             2.36         
 0.00             2.36         
 0.00             1.72         
 5.00             1.72         
 5.00             1.30         
 5.14             1.30         
 5.14             1.17         
 5.43             1.17         
 5.43             0.93         
 6.00             0.93         
 6.00             0.87         
 8.00             0.87
```

### 2. scripts/pseudo_binary.py

**pseudo_binary pd composition_1 composition_2**

This is used to calculate the chemical stability of two phases. 

The minimum point is marked in the comment column

```bash
$ pseudo_binary pd LiCoO2 Li3PS4
---------------------------------------------------------------------------------------------------- 
The starting phases compositions are  LiCoO2 and Li3PS4
All mixing ratio based on all formula already normalized to ONE atom per fu!

 ===  Pseudo-binary evolution profile  === 
x(Li3PS4)  x(LiCoO2)  Rxn. E. (meV/atom)  Mutual Rxn. E. (meV/atom)        Phase Equilibria       Comment 
  1.00       0.00            -0.00                   0.00                                  Li3PS4         
  0.50       0.50          -402.55                -402.55               CoS2, Co3S4, Li2S, Li3PO4         
  0.48       0.52          -405.94                -405.94             Co3S4, Li2S, Li2SO4, Li3PO4         
  0.41       0.59          -406.03                -406.03             Li2S, Li3PO4, Li2SO4, Co9S8  Minimum
  0.34       0.66          -368.44                -368.44             Li3PO4, Li2O, Li2SO4, Co9S8         
  0.16       0.84          -237.56                -237.56                Co, Li2O, Li2SO4, Li3PO4         
  0.15       0.85          -233.60                -233.60             Co, Li2SO4, Li6CoO4, Li3PO4         
  0.06       0.94           -90.55                 -90.55            CoO, Li3PO4, Li2SO4, Li6CoO4         
  0.00       1.00             0.00                   0.00                                  LiCoO2
```

**pseudo_binary gppd composition_1 composition_2 open_element chemical_potential**

This is used to calculate the electrochemical stability of two phases, in a system open to an element
held at the given chemical potential (in eV, referenced to the pure element).
The mixing ratios count all atoms of each phase, as in `pd`; the energies are per atom of the elements
other than the open element.

The minima of the reaction energy and of the mutual reaction energy are marked in the comment column.

Note: the example below was made with an older version, which weighted the two phases wrongly in the mutual
reaction energy when they contain different fractions of the open element. The current version gives different
values in that column (and may mark a different minimum).

```bash
$ pseudo_binary gppd LiCoO2 Li3PS4 Li -5
---------------------------------------------------------------------------------------------------- 
The starting phases compositions are  LiCoO2 and Li3PS4
All mixing ratio based on all formula already normalized to ONE atom per fu!
Chemical potential is miu_Li = -5.0, using elementary phase as reference.
------------------------------------------------------------

 ===  Pseudo-binary evolution profile  === 
x(Li3PS4)  x(LiCoO2)  Rxn. E. (meV/atom)  Mutual Rxn. E. (meV/atom)     Phase Equilibria           Comment       
  1.00      -0.00         -1,547.69                  0.00                            P2S7, S                     
  0.99       0.01         -1,601.10                -68.74                    CoS2, P2S7, S8O                     
  0.58       0.42         -1,679.16               -606.19                 CoS2, CoP4O11, S8O                     
  0.55       0.45         -1,679.82               -631.65                Co(PO3)2, CoS2, S8O         Rxn. E. Min.
  0.41       0.59         -1,568.49               -677.13              Co(PO3)2, CoS2, CoSO4                     
  0.33       0.67         -1,496.84               -695.56             CoS2, CoSO4, Co3(PO4)2                     
  0.29       0.71         -1,458.93               -704.30            Co3S4, CoSO4, Co3(PO4)2  Mutual Rxn. E. Min.
  0.26       0.74         -1,420.64               -704.18            CoSO4, Co3(PO4)2, Co9S8                     
  0.12       0.88         -1,182.46               -618.68              CoO, CoSO4, Co3(PO4)2                     
  0.10       0.90         -1,104.35               -569.65            Co3(PO4)2, CoSO4, Co3O4                     
  0.09       0.91         -1,075.28               -545.43                CoSO4, Co3O4, CoPO4                     
  0.00       1.00           -428.07                  0.00                               CoO2
```

**pseudo_binary gppd_screen composition_1 composition_2 open_element miu_low miu_high**

This scans a chemical potential range (in eV, referenced to the pure element) and reports the
electrochemical stability of the two phases over it.

```bash
$ pseudo_binary gppd_screen Li2S P2S5 Li -4 0
```

## License


Python library interface_stability is released under the MIT License. The terms of the license are as
follows:

    The MIT License (MIT) Copyright (c) 2018 UMD 
     
    Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated 
    documentation files (the "Software"), to deal in the Software without restriction, including without limitation 
    the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and 
    to permit persons to whom the Software is furnished to do so, subject to the following conditions:
     
    The above copyright notice and this permission notice shall be included in all copies or substantial portions of 
    the Software.
     
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO 
    THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
    TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
    SOFTWARE.
