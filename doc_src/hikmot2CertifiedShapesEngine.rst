CertifiedShapesEngine: An engine to generate certified shape intervals
======================================================================

The recommeded way to obtain certified intervals for the shapes is via manifold.tetrahedra_shapes(intervals=True) as described earlier. Here we document the engine used internally to generate these intervals. It is of interest for those users who want to understand the underlying interval math and experiment with the Newton interval method.


..   automodule:: snappy.hikmot2
..   autoclass:: CertifiedShapesEngine
     :members:
     :inherited-members:
