from datetime import datetime
import snappy
from multiprocessing import Pool
import itertools
from test_base import is_symplectic, testing_string
import random

start = 0
end = 1
scale = 1000
test = "random"

if len(snappy.HTLinkExteriors(crossings=15)) == 0:
    file_name = "small-db"
else:
    file_name = "large-db"


def process_manifold(index: int, output: bool = True):
    if test == "sequence":
        index = random.randint(1, len(snappy.HTLinkExteriors) - 1)

    M = snappy.HTLinkExteriors[index]
    label = M.identify()[0] if len(M.identify()) > 0 else ""
    print(label)

    if index == 0:
        return True

    basis = M.symplectic_basis()
    result = is_symplectic(basis)

    if result:
        string = "Passed"
    else:
        string = "Failed"

    if output and result is False:
        with open("logs/links-0.log", "a") as file:
            file.write(f"Testing: {str(index)} {(20 - len(str(index))) * ' '} {str(label)} {(40 - len(str(label))) * ' '} {string}\n")

    return result


def test_link_complements_pool(manifolds):
    with open("logs/total.log", "a") as file:
        if test == "random":
            length = len(manifolds)
        else:
            length = scale * (end - start)

        file.write(testing_string(length))
        print(testing_string(length))

    with Pool(maxtasksperchild=25) as pool:
        if test == "random":
            result = pool.imap(process_manifold, manifolds)
        else:
            result = pool.imap(process_manifold, range(start * scale, end * scale))

        for _ in range(start, end):
            lst = list(itertools.islice(result, scale))

            time = datetime.now().strftime('%d-%m-%y %H:%M:%S')
            print(f"[{time}]    Passed: {sum(lst)} / {len(lst)}")

            with open("logs/total.log", "a") as file:
                file.write(f"[{time}]    Passed: {sum(lst)} / {len(lst)}\n")


if __name__ == "__main__":
    with open(file_name, "r") as file:
        lst = file.readlines()

    manifolds = list(set([int(x[:-1]) for x in lst]))
    test_link_complements_pool(manifolds)
    # test_link_complements()
    # generate_tests()
    # unittest.main()
