# gas_UPPE
This is the shared package to simulate pulse propagation in a gas-filled hollow-core fiber with MM-UPPE with MATLAB.

It is useful for simulating pulse compression in gases, Raman physics (e.g. soliton self-frequency shift, Raman generation, or multidimensional solitary state generation), and photoionization-induced blueshift, etc.

## Capabilities:<br>
1. It solves the pulse propagation with
   - RK4IP (Runge-Kutta under the interaction picture) if single-mode (http://www.sciencedirect.com/science/article/pii/S0010465512004262).

> [!NOTE]
> I know that split-step algorithm is common, but I'd like to advocate people to switch to RK4IP since RK4IP has a higher-order truncation error, which allows higher precision or larger step size (and faster simulation).

   - MPA (massively parallel algorithm) if multimode (https://ieeexplore.ieee.org/document/8141863)

2. Adaptive step-size control are implemented for both RK4IP and MPA, which improves the performance and allows users to be free from worrying the reliability of a simulation.

> [!NOTE]
> Although adaptive-step-size control for RK4IP isn't new with published papers, adaptive-step-size control for MPA is new. I didn't publish a separate paper discussing this numerical scheme, which is perhaps the fastest and the most convenient numerical scheme for general multimode situations (as in step-index, graded-index, or hollow-core fibers, etc.) by far (written on 2/14/2024). The implementation detail is described in the supplement of https://doi.org/10.1364/JOSAB.500586.

3. Support both scalar and polarized scenarios, controlled with `sim.scalar=true/false`.
4. The gas encompasses noble (`He`, `Ne`, `Ar`, `Kr`, `Xe`) and Raman-active gases (`H<sub>2</sub>`, `N<sub>2</sub>`, `O<sub>2</sub>`, `air`, and `CH<sub>4</sub>`). 
5. Support photoionization with the PPT model. Check the supplement of https://opg.optica.org/josab/abstract.cfm?URI=josab-40-4-796.
6. Support Raman computations with both scalar and polarized scenarios, as well as with both vibrational and rotational Raman scattering.

> [!NOTE]
> This package should be the world's first correct implementation of polarized Raman simulations with both vibrational and rotational Raman scattering (scalar modeling has been known for years already). Due to the connection of angular momentum for the rotational Raman scattering, it's long been unclear how rotational Raman scattering affects nonlinear processes with polarization coupling. Although there are a few prior studies, I would define them as qualitative (or not quite fully quantitative) investigations. This package is able to solve quantitatively all the nonlinear interactions, electronic and both types (vibrational and rotational) of Raman scattering. In cases other than a single linearly polarized light involved in a decently nonlinear process, polarization coupling resulting from rotational Raman scattering is strong and requires significant attention. For details of the underlying physics, please read [our open-access paper](https://doi.org/10.1063/5.0189749).

7. Support the addition of spontaneous Raman scattering and input-pulse shot noise.
8. For multimode, GPU computations (with Nvidia CUDA) is highly recommended. I have written a lot of CUDA files to speed up simulations. It is controlled by `sim.gpu_yes=true/false`.

## Notes:<br>
For details, please read the supplement of our paper: https://doi.org/10.1063/5.0189749.  
Please don't forget to cite our paper if you find this code useful in your work. I, the young and early-career researcher, need your support. Similarly, if you need help or have questions about the code, please feel free to send me an email.

There is a `readme.pdf` in the `Documentations/` folder. Please find details of how to use this code in it. However, the fastest way to learn how to use this package is to learn from the examples in the `Examples/` folder.

The structure of this code is developed similar to our solid-core counterpart (https://github.com/AaHaHaa/MMTools). For optimization details of multimode (transverse modes and polarization modes), please see the supplement of our paper on multimode gain fiber (https://doi.org/10.1364/JOSAB.500586).

I'm Yi-Hao Chen, the author of the code and from Frank Wise's group at Cornell Applied Physics.  
Feel free to ask questions here or by sending me an email (email address is in my paper).

## Upgrade this code together:<br>
If you have any other function that you think important, please point it out in Github's discussions or send me an email. For example, perhaps you would like to add more gas species. I implement with current gases just due to my own research interest.

## Important notice:<br>
* 5/11/2024:<br>
I added some examples. APL Photonics data files are updated so that they can be run correctly.  
A bug regarding SRS under gradient pressure is also fixed.  
8pm (GMT-4): Extensive comments are added to examples.
* 5/14/2024:<br>
I extended the photoionization model to gases other than H2 and He.
* 7/24/2024:<br>
Fixed Ar n2 to the should-be correct value. See the comment in `gas_info()` for details.
