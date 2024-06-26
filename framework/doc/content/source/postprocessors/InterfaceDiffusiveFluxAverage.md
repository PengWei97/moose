# InterfaceDiffusiveFluxAverage

!syntax description /Postprocessors/InterfaceDiffusiveFluxAverage

The diffusive flux average is defined as
\begin{equation}
  \dfrac{\int_{\partial \Omega} -D \vec{\nabla} u \cdot \vec{n} d\Omega}{A}
\end{equation}
with $\partial \Omega$ the interface, $D$ the diffusion coefficient, $u$ the variable,
$\vec{n}$ the normal to the interface and $A$ the surface area of the interface.


!alert note
For finite volume variables, this postprocessor computes the diffusive flux using a two
point gradient. This is only accurate if the (interface) kernel is also computing gradients
this way. Currently, only [FVDiffusionInterface](/fviks/FVDiffusionInterface.md)
is computing gradients this way.

!alert warning
The expression of the diffusive flux in this object is generic, as described, and may differ from the diffusive flux in your specific physics implementation. If so, you may not use this object to compute the diffusive flux.

!syntax parameters /Postprocessors/InterfaceDiffusiveFluxAverage

!syntax inputs /Postprocessors/InterfaceDiffusiveFluxAverage

!syntax children /Postprocessors/InterfaceDiffusiveFluxAverage
