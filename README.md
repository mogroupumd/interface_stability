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

    Use your favorite way to install pymatgen and argparse first.
    
2. install this interface_stability package.

    clone this package from github
    
    ```bash
    git clone https://github.com/mogroupumd/interface_stability.git
    ```
    
    install the package
    
    ```bash
    python setup.py install --user
    ``` 
    
3. Try to import python classes in your python console

    ```bash
    python
    >>> from interface_stability import singlephase
    ```

4. The setup.py will automatically create an executable file phase_stability and pseudo_binary
 into your PATH. Try to call it from terminal and read the documentations:

    ```bash
    phase_stability -h
    pseudo_binary -h
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

There are two executable python scripts:

1. scripts/phase_stability.py

2. scripts/pseudo_binary.py

## License


Python library aimd is released under the MIT License. The terms of the license are as
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
