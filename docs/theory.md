# Theory Notes

This package reproduces the `f_zy` geometric factor from equation (11) of
Harris, Bryce, and Wasylishen, Can. J. Chem. 87, 1338-1351 (2009).

## Orientation Convention

Each crystallite orientation is represented by two angles, `alpha` and `beta`,
in the rotor frame. The active rotation is

```text
R(alpha, beta) = R_y(beta) R_z(alpha)
```

with columns equal to the crystallite `x`, `y`, and `z` axes expressed in the
rotor-frame basis:

```text
x = ( cos(alpha) cos(beta),  sin(alpha), -cos(alpha) sin(beta))
y = (-sin(alpha) cos(beta),  cos(alpha),  sin(alpha) sin(beta))
z = ( sin(beta),             0,           cos(beta))
```

A third Euler angle is omitted because a common rotation about the rotor axis
does not change the scalar `f_zy` factor.

## Equation (11)

The original expression is

```text
f_zy = 1/sqrt(3) [
    sin(xi_y) sin(xi_x) cos(epsilon_yx') c_yx
  + sin(xi_z) sin(xi_x) cos(epsilon_zx') c_zx
  + sin(xi_z) sin(xi_y) cos(epsilon_zy')
]
```

where `xi_a` is the angle between crystallite axis `a` and rotor-frame `z`.
The epsilon terms are angles between projected axes in the rotor-frame `xy`
plane, and a primed axis means the projection has been rotated by `-pi/2` about
rotor-frame `z`.

The code evaluates the same expression as projected-vector dot products:

```text
f_zy = [c_yx (p_y . p_x') + c_zx (p_z . p_x') + (p_z . p_y')] / sqrt(3)
```

Here `p_a` is the projection of crystallite axis `a` into the rotor-frame `xy`
plane. The projected vectors are intentionally not normalized. Their lengths
are `sin(xi_a)`, so each dot product already contains the
`sin(xi_a) sin(xi_b) cos(epsilon_ab')` factor. This also keeps the calculation
finite at the poles where a projected vector can vanish.

## Validation Formula

The tests also compare the geometric implementation to the independently
simplified expression

```text
f_zy = (
    -c_yx cos(beta)
    + c_zx sin(beta) sin(alpha)
    + sin(beta) cos(alpha)
) / sqrt(3)
```

This closed form is used as a validation oracle only. The public function keeps
the original geometric construction visible in code.

## Figure 3-Style AB Spectrum

For each crystallite orientation, the spectrum calculation first converts
`f_zy` into the effective antisymmetric coupling:

```text
A(Theta) = J_zy_anti f_zy(Theta; c_yx, c_zx)
```

The shielding-frequency difference is supplied in Hz:

```text
Delta S = sigma_I_hz - sigma_S_hz
```

The Table 2 quantity `C` is evaluated independently for every crystallite:

```text
C(Theta) = sqrt(J_iso^2 + A(Theta)^2 + Delta S^2)
```

The common carrier/center term in Table 2 is omitted, so the four transition
frequencies are

```text
v_34 = +1/2 C + 1/2 J_iso
v_12 = +1/2 C - 1/2 J_iso
v_24 = -1/2 C + 1/2 J_iso
v_13 = -1/2 C - 1/2 J_iso
```

with relative intensities

```text
I_34 = 1 - J_iso / C
I_12 = 1 + J_iso / C
I_24 = 1 + J_iso / C
I_13 = 1 - J_iso / C
```

Each transition intensity is multiplied by the crystallite orientation weight,
binned over `[-(|J_iso| + |J_zy_anti|), +(|J_iso| + |J_zy_anti|)]`, and
Gaussian-broadened using the requested FWHM in Hz.
