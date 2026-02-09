def compute_epsilon(RF):
    return RF(0.5) ** (RF.prec() // 2)


def compute_tube_injectivity_radius_epsilon(RF):
    return RF(0.5) ** (RF.prec() // 2 - 8)
