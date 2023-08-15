import random
from datetime import datetime
import snappy
import unittest
import spherogram
from tqdm import tqdm

start = 0
end = 1
scale = 1000
num_tests = 1000
test = "random"

if len(snappy.HTLinkExteriors(crossings=15)) == 0:
    file_name = "small-db"
else:
    file_name = "large-db"


def is_symplectic(M):
    """
    Test if the matrix M is symplectic
    :param M: square matrix
    :return: true or false
    """
    n = M.dimensions()

    for i in range(n[0]):
        for j in range(i, n[1]):
            omega = abs(symplectic_form(M.data[i], M.data[j]))

            if i % 2 == 0 and j % 2 == 1 and j == i + 1:
                if omega != 2:
                    return False
            elif omega:
                return False

    return True


def symplectic_form(u, v):
    return sum([u[2 * i] * v[2 * i + 1] - u[2 * i + 1] * v[2 * i]
                for i in range(len(u) // 2)])


def save_manifold(index: int):
    M = snappy.HTLinkExteriors[index]
    M.save(f"CuspedCensusData/link-{index}.tri")


def process_manifold(index: int):
    M = snappy.HTLinkExteriors[index]
    label = M.identify()[0] if len(M.identify()) > 0 else ""

    if index == 0:
        return True

    basis = M.symplectic_basis()
    result = is_symplectic(basis)

    if result:
        string = "Passed"
    else:
        string = "Failed"

    print(f"Testing: {str(index)} {(20 - len(str(index))) * ' '} {str(label)} {(40 - len(str(label))) * ' '} {string}")

    return result


def random_link_exteriors(n: int, n_tet: int, n_cusps: int):
    for i in range(n):
        L = spherogram.random_link(n_tet, n_cusps, alternating=True)
        M = spherogram.Link.exterior(L)
        print(M.num_cusps())
        M.save(f"CuspedCensusData/link-{n_tet}-{n_cusps}-{i}.tri")


def generate_tests(output: bool):
    print(f"[{datetime.now().strftime('%d-%m-%y %H:%M:%S')}]  Generating symplectic basis tests")
    for _ in range(2000):
        index = random.randint(1, len(snappy.HTLinkExteriors) - 1)
        process_manifold(index)

        if output:
            with open(file_name, "a") as file:
                file.write(f"{index}\n")


def testing_string(num: int = 0):
    time = datetime.now().strftime('%d-%m-%y %H:%M:%S')

    if test == "random":
        return f"[{time}]    Testing ({test}): {num} manifolds\n"
    elif test == "sequence":
        return f"[{time}]    Testing ({test}): {scale * start} - {scale * end}\n"


def test_link_complements():
    """ Pick 'num_test' random manifolds from HTLinkExteriors and test symplectic_basis() """
    # with open(file_name, "r") as file:
    #     lst = file.readlines()

    # manifolds = list(set([int(x[:-1]) for x in lst]))

    # with open("logs/total.log", "a") as file:
    #     file.write(testing_string(len(manifolds)))
    print(testing_string(num_tests))

    # if test == "random":
    #     result = [process_manifold(i) for i in range(len(manifolds))]
    # else:
    #     result = [process_manifold(i) for i in range(scale * start, scale * end)]
    manifolds = [random.randint(1, len(snappy.HTLinkExteriors)) for _ in range(num_tests)]
    result = [process_manifold(i) for i in manifolds]

    # with open("logs/total.log", "a") as file:
    #     file.write(f"[{datetime.now().strftime('%d-%m-%y %H:%M:%S')}]    Passed: {sum(result)} / {len(result)}\n")
    print(f"[{datetime.now().strftime('%d-%m-%y %H:%M:%S')}]    Passed: {sum(result)} / {len(result)}")


class TestSymplecticBasis(unittest.TestCase):
    def test_knot_complements(self):
        i = 0
        for M in tqdm(snappy.CensusKnots, desc="Knots...", ncols=120):
            with self.subTest(i=i):
                # print(M.identify()[0])
                basis = M.symplectic_basis()
                self.assertTrue(is_symplectic(basis), str(M.identify()[0]))
                i += 1

    # @unittest.skip
    def test_link_complements(self):
        i = 0
        for M in tqdm(snappy.HTLinkExteriors[1:1000], desc="Links...", ncols=120):
            with self.subTest(i=i):
                # print(M.identify()[0])
                basis = M.symplectic_basis()
                self.assertTrue(is_symplectic(basis))
                i += 1

    @unittest.skip
    def test_random_links(self):
        iterations = 10

        for i in tqdm(range(iterations), desc="Random Links...", ncols=120):
            with self.subTest(i=i):
                L = spherogram.random_link(100, num_components=random.randint(3, 10), alternating=True)
                M = spherogram.Link.exterior(L)
                basis = M.symplectic_basis()
                self.assertTrue(is_symplectic(basis))


if __name__ == "__main__":
    test_link_complements()
    # generate_tests(True)
    # unittest.main()
