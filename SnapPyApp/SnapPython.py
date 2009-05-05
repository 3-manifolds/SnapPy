import SnapPy

m = SnapPy.Manifold("m004")
d = SnapPy.DirichletDomain(m)
d.view()
d.viewer.window.mainloop()

#LE = SnapPy.plink.LinkEditor()
#LE.window.mainloop()
