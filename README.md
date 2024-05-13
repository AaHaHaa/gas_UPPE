# gas_UPPE
This is the shared package to simulate pulse propagation in a gas-filled hollow-core fiber with MM-UPPE with MATLAB.

It solves the pulse propagation with RK4IP if single mode and MPA if multimode. Both scalar and polarized scenarios can be simulated. The gas encompasses inert (He, Ne, Ar, Kr, Xe) and Raman-active gases (H<sub>2</sub>, N<sub>2</sub>, O<sub>2</sub>, air, and CH<sub>4</sub>). Besides, it is implemented with an adaptive step-size control for both methods, which improves the performance and allows users to be free from worrying the reliability of a simulation. Photoionization is included as well.
> Note (copied from https://github.com/AaHaHaa/MMTools):<br>
Although adaptive-step-size control for RK4IP isn't new with published papers, adaptive-step-size control for MPA is new. I didn't publish a separate paper discussing this numerical scheme, which is perhaps the fastest and the most convenient numerical scheme for general multimode situations by far (written on 2/14/2024). The implementation detail is described in the supplement of https://doi.org/10.1364/JOSAB.500586.

This package should be the world's first correct first-principle implementation of polarized Raman simulations with both vibrational and rotational Raman scattering (scalar modeling has been known for years already). Due to the connection of angular momentum for the rotational Raman scattering, it's long been unclear how rotational Raman scattering affects nonlinear processes with polarization coupling. Although there are a few prior studies, I would define them as qualitative (or not quite fully quantitative) investigations. This package is able to solve quantitatively all the nonlinear interactions, electronic and both types (vibrational and rotational) of Raman scattering. In cases other than a single linearly polarized light involved in a decently nonlinear process, polarization coupling resulting from rotational Raman scattering is strong and requires significant attention. For details of the underlying physics, please read [our open-access paper](https://doi.org/10.1063/5.0189749).

For multimode, GPU computations (with Nvidia cuda) is highly recommended. I have written a lot of cuda files to speed up simulations. It is controlled by `sim.gpu_yes=true or false`.

For details, please read the supplement of our paper: https://doi.org/10.1063/5.0189749.  
Please don't forget to cite our paper if you find this code useful in your work. I, the young and early-career researcher, need your support. Similarly, if you need help or have questions about the code, please feel free to send me an email.

There is a readme.pdf in the Documentations/ folder. Please find details of how to use this code in it. However, the fastest way to learn how to use this package is to learn from the examples in the Examples/ folder.

The structure of this code is developed similar to our solid-core counterpart (https://github.com/AaHaHaa/MMTools). For optimization details of multimode (transverse modes and polarization modes), please see the supplement of our paper on multimode gain fiber (https://doi.org/10.1364/JOSAB.500586).

I'm Yi-Hao Chen, the author of the code and from Frank Wise's group at Cornell Applied Physics.  
Feel free to ask questions here or by sending me an email (email address is in my paper).

## Important notice:<br>
* 5/11/2024:<br>
I added some examples. APL Photonics data files are updated so that they can be run correctly.  
A bug regarding SRS under gradient pressure is also fixed.  
8pm (GMT-4): Extensive comments are added to examples.
