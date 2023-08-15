from datetime import datetime
import snappy
from test_base import is_symplectic, testing_string
from multiprocessing import Pool
import itertools

start = 0
end = 1
scale = 1000
test = "random"

if len(snappy.HTLinkExteriors(crossings=15)) == 0:
    file_name = "links-linux"
else:
    file_name = "links"

print(f"[{datetime.now().strftime('%d-%m-%y %H:%M:%S')}]    Building test set")

with open(file_name, "r") as file:
    lst = file.readlines()
    manifolds = list(set([int(x[:-1]) for x in lst]))
    manifolds_tri = [snappy.HTLinkExteriors[i] for i in manifolds]
    manifolds_labels = [M.identify()[0] for M in manifolds_tri if len(M.identify()) > 0]


def process_manifold(i: int, output: bool = True):
    M = manifolds_tri[i]
    index = manifolds[i]
    label = manifolds_labels[i]

    if index == 0:
        return True

    basis = M.symplectic_basis()
    result = is_symplectic(basis)

    if result:
        string = "Passed"
    else:
        string = "Failed"

    if output:
        with open("logs/links-0.log", "a") as file:
            file.write(f"Testing: {str(index)} {(20 - len(str(index))) * ' '} {str(label)} {(40 - len(str(label))) * ' '} {string}\n")

    return result


def test_link_complements_pool():
    with open("logs/total.log", "a") as file:
        file.write(testing_string(len(manifolds)))
        print(testing_string(len(manifolds)))

    with Pool(maxtasksperchild=25) as pool:
        if test == "random":
            result = pool.imap(process_manifold, range(len(manifolds)))
        else:
            result = pool.imap(process_manifold, range(start * scale, end * scale))

        for _ in range(start, end):
            lst = list(itertools.islice(result, scale))

            time = datetime.now().strftime('%d-%m-%y %H:%M:%S')
            print(f"[{time}]    Passed: {sum(lst)} / {len(lst)}")

            with open("logs/total.log", "a") as file:
                file.write(f"[{time}]    Passed: {sum(lst)} / {len(lst)}\n")


if __name__ == "__main__":
    test_link_complements_pool()
    # unittest.main()
