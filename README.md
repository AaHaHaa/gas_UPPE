# gas_UPPE
This is the shared package to simulate pulse propagation in a gas-filled hollow-core fiber with MM-UPPE with MATLAB.

It solves the pulse propagation with RK4IP if single mode and MPA if multimode. Both scalar and polarized scenarios can be simulated. The gas encompasses inert (He, Ne, Ar, Kr, Xe) and Raman-active gases (H<sub>2</sub>, N<sub>2</sub>, O<sub>2</sub>, air, and CH<sub>4</sub>). Besides, it is implemented with an adaptive step-size control for both methods, which improves the performance and allows users to be free from worrying the reliability of a simulation. Photoionization is included as well.

For multimode, GPU computations (with Nvidia cuda) is highly recommended. I have written a lot of cuda files to speed up simulations. It is controlled by `sim.gpu_yes=true or false`.

For details, please read the supplement of our paper: https://doi.org/10.1063/5.0189749.  
Please don't forget to cite our paper if you find this code useful in your work. I, the young and early-career researcher, need your support. Similarly, if you need help or have questions about the code, please feel free to send me an email.

There is a readme.pdf in the Documentations/ folder. Please find details of how to use this code in it. However, the fastest way to learn how to use this package is to learn from the examples in the Examples/ folder.

The structure of this code is developed similar to our solid-core counterpart (https://github.com/AaHaHaa/MMTools). For optimization details of multimode (transverse modes and polarization modes), please see the supplement of our paper on multimode gain fiber (https://doi.org/10.1364/JOSAB.500586).

I'm Yi-Hao Chen, the author of the code and from Frank Wise's group at Cornell Applied Physics.  
Feel free to ask questions here or by sending me an email (email address is in my paper).
