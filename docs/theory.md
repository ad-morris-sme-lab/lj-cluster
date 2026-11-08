# Theory and model assumptions

For particles at Cartesian positions **r** in metres, the model sums the pair potential

\[U(r)=4\epsilon\left[(\sigma/r)^{12}-(\sigma/r)^6\right].\]

Here \(\epsilon\) is the well depth in joules and \(\sigma\) is the distance in metres where the pair potential is zero. The minimum occurs at \(r=2^{1/6}\sigma\), where \(U=-\epsilon\). Total energy is the sum over unique pairs \(i<j\); the implementation uses no cutoff or periodic boundary conditions.

The gradient for a pair, with \(\mathbf d=\mathbf r_i-\mathbf r_j\), is

\[\nabla_{\mathbf r_i} U_{ij}=24\epsilon\left[(\sigma/r)^6-2(\sigma/r)^{12}\right]\mathbf d/r^2.\]

The partner receives the negative gradient, so total gradient is zero under translation. Units are SI: coordinates and sigma in metres, energy and epsilon in joules, gradient in joules per metre. This is a classical, isotropic, homogeneous-pair model. It omits electrostatics, many-body effects, quantum effects, chemical bonding, species mixing rules, and thermal sampling; parameters must be supplied for the system of interest. It is a teaching and small-cluster calculation, not a general molecular force field.
