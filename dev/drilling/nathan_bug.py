from snappy import *

# A bug discovered by Nathan

M = ManifoldHP('K11n34(0,1)')

words = ['iFcdbEiFJ', 'iFJ']

#print(M.drill_words(words).filled_triangulation().isometry_signature())
print(M.drill_words(words[::-1]).filled_triangulation().isometry_signature())
