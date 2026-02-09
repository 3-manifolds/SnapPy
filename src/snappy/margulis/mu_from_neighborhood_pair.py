from .neighborhood import Neighborhood
from .geodesic_neighborhood import GeodesicNeighborhood
from .cusp_neighborhood_neighborhood import CuspNeighborhoodNeighborhood

from ..math_basics import correct_min, correct_max, upper

def mu_from_neighborhood_pair(
        neighborhood1 : Neighborhood,
        neighborhood2 : Neighborhood,
        object_distance, verified : bool):

    is_cusp1 = isinstance(neighborhood1, CuspNeighborhoodNeighborhood)
    is_cusp2 = isinstance(neighborhood2, CuspNeighborhoodNeighborhood)
    if is_cusp1 and is_cusp2:
        return mu_from_cusp_neighborhood_pair(
            neighborhood1, neighborhood2, object_distance)

    if neighborhood1.index == neighborhood2.index:
        return neighborhood1.epsilon_from_radius(object_distance / 2)

    return mu_from_neighborhood_pair_generic(
        neighborhood1, neighborhood2, object_distance, verified=verified)

def _my_max(values):
    n = len(values)
    if n == 0:
        raise RuntimeError("No values to take max")
    if n == 1:
        return values[0]
    return correct_max(values)

def mu_from_neighborhood_pair_generic(
        neighborhood1, neighborhood2, object_distance, verified):

    lower_bound_epsilon = _my_max([
        neighborhood.geodesic_info.length.real()
        for neighborhood in [ neighborhood1, neighborhood2 ]
        if isinstance(neighborhood, GeodesicNeighborhood) ])

    lower_bound_radius1 = neighborhood1.radius_from_epsilon(lower_bound_epsilon)
    lower_bound_radius2 = neighborhood2.radius_from_epsilon(lower_bound_epsilon)
    lower_bound_total_radius = lower_bound_radius1 + lower_bound_radius2

    if lower_bound_total_radius > object_distance:
        return lower_bound_epsilon

    if lower_bound_total_radius < object_distance:
        return mu_from_neighborhood_pair_implicit(
            neighborhood1, neighborhood2,
            object_distance,
            lower_bound_epsilon,
            verified)

    if verified:
        upper_bound_epsilon1 = (
            neighborhood1.epsilon_from_radius(object_distance))
        upper_bound_epsilon2 = (
            neighborhood2.epsilon_from_radius(object_distance))
        upper_bound_epsilon = correct_min([
            upper_bound_epsilon1, upper_bound_epsilon2])
        return lower_bound_epsilon.union(upper_bound_epsilon)
    else:
        return lower_bound_epsilon

def mu_from_neighborhood_pair_implicit(
        neighborhood1, neighborhood2,
        object_distance,
        lower_bound_mu,
        verified):

    initial_mu = correct_max([neighborhood1.epsilon,neighborhood2.epsilon])

    mu = mu_from_neighborhood_pair_newton(
        neighborhood1, neighborhood2,
        object_distance,
        initial_mu,
        lower_bound_mu,
        verified=verified)

    if verified:
        return mu_from_neighborhood_pair_newton_interval(
            neighborhood1, neighborhood2,
            object_distance,
            mu)
    else:
        return mu

_max_iterations_newton = 500

def mu_from_neighborhood_pair_newton(
        neighborhood1, neighborhood2,
        object_distance,
        initial_mu,
        lower_bound_mu,
        verified):

    RF = initial_mu.parent()
    mu = initial_mu

    pos_diff = None
    neg_diff = None

    abs_diff = None
    for i in range(_max_iterations_newton):
        radius1, d1 = neighborhood1.radius_and_derivative_from_epsilon(mu)
        radius2, d2 = neighborhood2.radius_and_derivative_from_epsilon(mu)

        radius = radius1 + radius2
        d = d1 + d2

        diff = object_distance - radius

        if diff >= 0:
            if pos_diff is not None:
                if not diff < pos_diff:
                    return mu
            pos_diff = diff
        else:
            if neg_diff is not None:
                if not -diff < -neg_diff:
                    return mu
            neg_diff = diff

        new_mu = mu + diff / d
        safe_mu = (mu - lower_bound_mu) / 8 + lower_bound_mu
        mu = correct_max([new_mu, safe_mu])

        if verified:
            mu = RF(mu.center())

    raise RuntimeError("Newton method did not converge")

def mu_from_neighborhood_pair_newton_interval(
        neighborhood1 : Neighborhood,
        neighborhood2 : Neighborhood,
        object_distance,
        mu_sample):

    radius1 = neighborhood1.radius_from_epsilon(mu_sample)
    radius2 = neighborhood2.radius_from_epsilon(mu_sample)
    radius = radius1 + radius2

    diff = object_distance - radius

    mu_interval = mu_sample

    for i in range(50):
        d1 = neighborhood1.radius_derivative_from_epsilon(mu_interval)
        d2 = neighborhood2.radius_derivative_from_epsilon(mu_interval)
        d = d1 + d2

        new_mu_interval = mu_sample + diff / d

        if new_mu_interval in mu_interval:
            return new_mu_interval

        mu_interval = mu_interval.union(new_mu_interval)

    raise RuntimeError(
        "Interval Newton method did not converge for neighborhood pair %r %r, object_distance = %r, mu_sample = %r, mu_interval = %r" % (
            neighborhood1, neighborhood2, object_distance, mu_sample, mu_interval))


def mu_from_cusp_neighborhood_pair(
        neighborhood1 : CuspNeighborhoodNeighborhood,
        neighborhood2 : CuspNeighborhoodNeighborhood,
        object_distance):
    # Radius1 + radius2 = minimum
    # log(h / _factor1) + log(h / _factor2) = minimum
    # log(h^2 / (_factor1 * _factor2)) = minimum
    # h^2 = exp(minimum) / (_factor1 * _factor2)

    len_prod = neighborhood1.euclidean_length * neighborhood2.euclidean_length
    h = (object_distance.exp() * len_prod).sqrt()
    return 2 * (h/2).arcsinh()
