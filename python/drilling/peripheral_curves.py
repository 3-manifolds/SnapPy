from ..snap.t3mlite import Mcomplex
from ..snap.t3mlite import simplex

from collections import deque

def install_all_peripheral_curves(mcomplex : Mcomplex) -> None:
    for tet in mcomplex.Tetrahedra:
        if tet.base_for_peripherals:
            install_peripheral_curves(tet)

def install_peripheral_curves(tet):
    install_meridian(tet)
    install_longitude(tet)

def install_meridian(tet0):
    tet1 = tet0
    while True:
        tet1.PeripheralCurves[0][0][simplex.V0][simplex.F3] = +1
        tet1.PeripheralCurves[0][0][simplex.V0][simplex.F2] = -1
        tet2 = tet1.Neighbor[simplex.F2]
        tet2.PeripheralCurves[0][1][simplex.V0][simplex.F2] = +1
        tet2.PeripheralCurves[0][1][simplex.V0][simplex.F1] = -1
        tet3 = tet2.Neighbor[simplex.F1]
        tet3.PeripheralCurves[0][0][simplex.V0][simplex.F1] = +1
        tet3.PeripheralCurves[0][0][simplex.V0][simplex.F2] = -1
        tet4 = tet3.Neighbor[simplex.F2]
        tet4.PeripheralCurves[0][1][simplex.V0][simplex.F2] = +1
        tet4.PeripheralCurves[0][1][simplex.V0][simplex.F3] = -1
        tet1 = tet4.Neighbor[simplex.F3]
        if tet1 is tet0:
            break

def has_meridian(tet):
    for sheet in tet.PeripheralCurves[0]:
        for v in sheet[simplex.V0].values():
            if v != 0:
                return True
    return False

def install_longitude(tet0):

    if not has_meridian(tet0):
        raise Exception("tet0 not meridian")
    if not has_meridian(tet0.Neighbor[simplex.F2]):
        raise Exception("tet0->2 not meridian")
    if not has_meridian(tet0.Neighbor[simplex.F3]):
        raise Exception("tet0->1 not meridian")

    tet1 = tet0.Neighbor[simplex.F1]
    tet_finish = tet0.Neighbor[simplex.F2]

    pending_tets = deque([(tet1, simplex.F1)])
    visited_tets = { }
    while True:
        tet, old_f = pending_tets.popleft()
        if tet in visited_tets:
            continue
        visited_tets[tet] = old_f
        for f in [ simplex.F1, simplex.F2, simplex.F3 ]:
            if f != old_f:
                neighbor = tet.Neighbor[f]
                if neighbor is tet_finish:
                    l = [ (tet0, simplex.F2), (neighbor, f) ]
                    while tet in visited_tets:
                        f = visited_tets[tet]
                        l.append((tet, f))
                        tet = tet.Neighbor[f]

                    for tet, f in l:
                        tet.PeripheralCurves[1][tet.Index % 2][simplex.V0][f] -= 1
                        tet.Neighbor[f].PeripheralCurves[1][(tet.Index+1) % 2][simplex.V0][f] += 1
                    return
                if not has_meridian(neighbor):
                    pending_tets.append((neighbor, f))
