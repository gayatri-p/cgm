## files

- `final_grid.pkl`: $\log N$ of each ion species.
- `final_flat_logT.pkl`: (flattened) average $\log T$ values.
- `final_grid_logf.pkl`: $\log f$ value for each species.

    Let $$N_H = 10^{\log N_\text{HI}} + 10^{\log N_\text{HII}}$$
    - for H, $$f_\text{HI} = \frac{10^{\log N_\text{HI}}}{N_\text{H}}$$
    and $$f_\text{HII} = \frac{10^{\log N_\text{HII}}}{N_\text{H}}$$
    - for He, $$f_\text{HeIII} = \frac{10^{\log N_\text{HeIII}}}{10^{\log N_\text{HeI}} + 10^{\log N_\text{HeII}} + 10^{\log N_\text{HeIII}}}$$
    - starting from Li for a metallicity $M$, $$f_i = \frac{10^{\log N_i}}{N_H 10^{M}}$$

## Calculate $\log U$ from $n_\text{HI}$

Given $J_\nu$ and $n_\text{HI}$.

First calculate the ionizing photon flux in units photon/s/cm$^2$.
Consider only $J_\nu$ energies corresponding to ionizing H-flux $\ge 1$ Ryd, i.e. from $\nu_0$ to $\infty$.

- $J_\nu$ is the average intensity over $4\pi$ sr in ergs/s/cm$^2$/Hz/sr, first convert to J/s/cm$^2$/Hz/sr.
- Next, to go from energy to photon count, divide by $h\nu$, the photon energy, this puts units as photon/s/cm$^2$/Hz/sr.
- Multiply by $4\pi$ sr to account for photons arriving per cm$^2$ from all $4\pi$ sr, this puts units as photon/s/cm$^2$/Hz.

$$\text{(photon density)}_\nu = \frac{4\pi J_\nu 10^{-7}}{h\nu}$$

Calculate the ionizing flux of photons >= 1 Ryd, units in photon/s/cm^2.
$$\phi = \int_{\nu_0}^\infty \frac{4\pi J_\nu 10^{-7}}{h\nu} \,\, d\nu $$

Assuming a ionizing photon flux in units photon/s/cm$^2$.
Divide by c in cm/s to get photon number density in photon/cm^3.

$$n_\gamma = \frac{\phi}{c}$$

Calculate the ionizing parameter. Divide by $n_\text{HI}$ in cm$^{-3}$.
$$U = \frac{n_\gamma}{n_\text{HI}}$$