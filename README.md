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

## Install

This package works with both Python 2.7+ and Python 3.x. However, it is suggested to use Python 3.x, as the dependent package Pymatgen will be py3k only in the future.


1. Install all dependency packages 

    Use your favorite way to install [pymatgen](http://pymatgen.org/) first.
    
2. install this interface_stability package.

    clone this package from github
    
    ```bash
    $ git clone https://github.com/mogroupumd/interface_stability.git
    ```
    
    install the package
    ```bash
    $ python setup.py install --user
    ``` 
    
3. Try to import python classes in your python console

    ```bash
    $ python
    >>> from interface_stability import singlephase
    ```

4. The setup.py will automatically create an executable file phase_stability and pseudo_binary
 into your PATH. Try to call it from terminal and read the documentations:

    ```bash
    $ phase_stability -h
    $ pseudo_binary -h
    ```
    
5. Setup Materials Project API key. This enables you to fetch data from the Materials Project database.
   
   The documentation is at 
   https://materialsproject.org/open
   
   Your API key is at (login required)
   https://materialsproject.org/dashboard
   
   You need to put the API key in ~/.pmgrc.yaml file. (Create this file if not exist)
   ```bash
   PMG_MAPI_KEY: [Your API key goes here]
   ```
   
## Usage

There are two executable python scripts, both work with a few sub-commands options. 
You can always use -h to see the help information. 
For example,

```bash
$ phase_stability -h
$ phase_stability evolution -h
```

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

**pseudo_binary gppd composition_1 composition_2**

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
(py3k) Yizhous-MBP:interface_stability yizhou$ pseudo_binary gppd LiCoO2 Li3PS4 Li -5
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
