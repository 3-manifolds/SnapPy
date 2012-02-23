from surface_generator import *
from Twister import determine_info

for i in range(3):
	for j in range(3):
		write_surface_to_file(construct_file_contents(i,j), ".\\Surfaces\\Test\\LP_S_" + str(i) + '_' + str(j) + ".sur")
		genus, boundary, Euler_characteristic = determine_info(".\\Surfaces\\Test\\LP_S_" + str(i) + '_' + str(j) + ".sur")
		assert genus == i and boundary == j

