"""
M2 (flint, approx solve, tetshapes only, *with* incremental prec increase):

16 tets: Found 100 bit shapes in 0.002s
16 tets: Found 1000 bit shapes in 0.007s
16 tets: Found 10000 bit shapes in 0.084s
16 tets: overall time 0.096s
20 tets: Found 100 bit shapes in 0.003s
20 tets: Found 1000 bit shapes in 0.014s
20 tets: Found 10000 bit shapes in 0.181s
20 tets: overall time 0.200s
48 tets: Found 100 bit shapes in 0.011s
48 tets: Found 1000 bit shapes in 0.042s
48 tets: Found 10000 bit shapes in 0.387s
48 tets: overall time 0.443s
66 tets: Found 100 bit shapes in 0.045s
66 tets: Found 1000 bit shapes in 0.151s
66 tets: Found 10000 bit shapes in 1.889s
66 tets: overall time 2.089s
90 tets: Found 100 bit shapes in 0.359s
90 tets: Found 1000 bit shapes in 2.742s
90 tets: Found 10000 bit shapes in 11.054s
90 tets: overall time 14.163s
Overall 0: 16.991517066955566


M2 (pari, tetshapes only):

16 tets: Found 100 bit shapes in 0.012s
16 tets: Found 1000 bit shapes in 0.020s
16 tets: Found 10000 bit shapes in 0.317s
16 tets: overall time 0.354s
20 tets: Found 100 bit shapes in 0.010s
20 tets: Found 1000 bit shapes in 0.033s
20 tets: Found 10000 bit shapes in 0.576s
20 tets: overall time 0.622s
48 tets: Found 100 bit shapes in 0.057s
48 tets: Found 1000 bit shapes in 0.213s
48 tets: Found 10000 bit shapes in 3.666s
48 tets: overall time 3.939s
66 tets: Found 100 bit shapes in 0.105s
66 tets: Found 1000 bit shapes in 0.364s
66 tets: Found 10000 bit shapes in 6.211s
66 tets: overall time 6.684s
90 tets: Found 100 bit shapes in 0.264s
90 tets: Found 1000 bit shapes in 1.398s
90 tets: Found 10000 bit shapes in 34.435s
90 tets: overall time 36.105s
Overall 0: 47.703080892562866


M2 (Sage, approx solve, tetshapes only):

16 tets: Found 100 bit shapes in 0.008s
16 tets: Found 1000 bit shapes in 0.016s
16 tets: Found 10000 bit shapes in 0.239s
16 tets: overall time 0.266s
20 tets: Found 100 bit shapes in 0.011s
20 tets: Found 1000 bit shapes in 0.026s
20 tets: Found 10000 bit shapes in 0.444s
20 tets: overall time 0.484s
48 tets: Found 100 bit shapes in 0.092s
48 tets: Found 1000 bit shapes in 0.182s
48 tets: Found 10000 bit shapes in 2.976s
48 tets: overall time 3.255s
66 tets: Found 100 bit shapes in 0.139s
66 tets: Found 1000 bit shapes in 0.332s
66 tets: Found 10000 bit shapes in 5.072s
66 tets: overall time 5.548s
90 tets: Found 100 bit shapes in 0.312s
90 tets: Found 1000 bit shapes in 1.053s
90 tets: Found 10000 bit shapes in 28.238s
90 tets: overall time 29.611s
Overall 0: 39.16382122039795


M2 (flint, verified, no 90-tet manifold):

16 tets: Found 100 bit shapes in 0.011s
16 tets: Found 1000 bit shapes in 0.014s
16 tets: Found 10000 bit shapes in 0.210s
16 tets: overall time 0.240s
20 tets: Found 100 bit shapes in 0.009s
20 tets: Found 1000 bit shapes in 0.027s
20 tets: Found 10000 bit shapes in 0.457s
20 tets: overall time 0.495s
34 tets: Found 100 bit shapes in 0.022s
34 tets: Found 1000 bit shapes in 0.068s
34 tets: Found 10000 bit shapes in 1.183s
34 tets: overall time 1.275s
48 tets: Found 100 bit shapes in 0.039s
48 tets: Found 1000 bit shapes in 0.114s
48 tets: Found 10000 bit shapes in 1.688s
48 tets: overall time 1.845s
66 tets: Found 100 bit shapes in 0.114s
66 tets: Found 1000 bit shapes in 0.364s
66 tets: Found 10000 bit shapes in 7.478s
66 tets: overall time 7.959s
Overall: 11.813881158828735


M2 (Sage, verified, no 90-tet manifold):

16 tets: verified 100 bit shapes in 0.175s
16 tets: verified 1000 bit shapes in 0.081s
16 tets: verified 10000 bit shapes in 0.418s
16 tets: overall time 0.678s
20 tets: verified 100 bit shapes in 0.027s
20 tets: verified 1000 bit shapes in 0.058s
20 tets: verified 10000 bit shapes in 0.680s
20 tets: overall time 0.768s
34 tets: verified 100 bit shapes in 0.110s
34 tets: verified 1000 bit shapes in 0.149s
34 tets: verified 10000 bit shapes in 1.564s
34 tets: overall time 1.826s
48 tets: verified 100 bit shapes in 0.173s
48 tets: verified 1000 bit shapes in 0.355s
48 tets: verified 10000 bit shapes in 3.843s
48 tets: overall time 4.374s
66 tets: verified 100 bit shapes in 0.278s
66 tets: verified 1000 bit shapes in 0.635s
66 tets: verified 10000 bit shapes in 6.285s
66 tets: overall time 7.203s
Overall: 14.848s



-------

Intel is keeling-i18 with AMD EPYC 9115 @ 2.6 GHZ (circa 2024)


Intel (flint, approx solve, tetshapes only, *with* incremental prec increase):

16 tets: Found 100 bit shapes in 0.002s
16 tets: Found 1000 bit shapes in 0.008s
16 tets: Found 10000 bit shapes in 0.086s
16 tets: overall time 0.098s
20 tets: Found 100 bit shapes in 0.003s
20 tets: Found 1000 bit shapes in 0.015s
20 tets: Found 10000 bit shapes in 0.186s
20 tets: overall time 0.207s
48 tets: Found 100 bit shapes in 0.011s
48 tets: Found 1000 bit shapes in 0.043s
48 tets: Found 10000 bit shapes in 0.393s
48 tets: overall time 0.450s
66 tets: Found 100 bit shapes in 0.047s
66 tets: Found 1000 bit shapes in 0.161s
66 tets: Found 10000 bit shapes in 1.922s
66 tets: overall time 2.133s
90 tets: Found 100 bit shapes in 0.422s
90 tets: Found 1000 bit shapes in 0.695s
90 tets: Found 10000 bit shapes in 11.157s
90 tets: overall time 12.281s
Overall 0: 15.16890549659729


Intel (pari, approx solve, tetshapes only):

16 tets: Found 100 bit shapes in 0.010s
16 tets: Found 1000 bit shapes in 0.016s
16 tets: Found 10000 bit shapes in 0.256s
16 tets: overall time 0.288s
20 tets: Found 100 bit shapes in 0.033s
20 tets: Found 1000 bit shapes in 0.026s
20 tets: Found 10000 bit shapes in 0.471s
20 tets: overall time 0.533s
48 tets: Found 100 bit shapes in 0.058s
48 tets: Found 1000 bit shapes in 0.167s
48 tets: Found 10000 bit shapes in 3.052s
48 tets: overall time 3.284s
66 tets: Found 100 bit shapes in 0.133s
66 tets: Found 1000 bit shapes in 0.290s
66 tets: Found 10000 bit shapes in 5.184s
66 tets: overall time 5.618s
90 tets: Found 100 bit shapes in 0.270s
90 tets: Found 1000 bit shapes in 0.920s
90 tets: Found 10000 bit shapes in 29.905s
90 tets: overall time 31.113s
Overall 0: 40.835227489471436


Intel (Sage, approx solve, tetshapes only):

16 tets: Found 100 bit shapes in 0.007s
16 tets: Found 1000 bit shapes in 0.015s
16 tets: Found 10000 bit shapes in 0.316s
16 tets: overall time 0.341s
20 tets: Found 100 bit shapes in 0.010s
20 tets: Found 1000 bit shapes in 0.023s
20 tets: Found 10000 bit shapes in 0.476s
20 tets: overall time 0.513s
48 tets: Found 100 bit shapes in 0.101s
48 tets: Found 1000 bit shapes in 0.181s
48 tets: Found 10000 bit shapes in 3.134s
48 tets: overall time 3.420s
66 tets: Found 100 bit shapes in 0.142s
66 tets: Found 1000 bit shapes in 0.328s
66 tets: Found 10000 bit shapes in 5.256s
66 tets: overall time 5.732s
90 tets: Found 100 bit shapes in 0.330s
90 tets: Found 1000 bit shapes in 0.986s
90 tets: Found 10000 bit shapes in 30.359s
90 tets: overall time 31.685s
Overall: 41.691020011901855


Intel (flint, verified, no 90-tet manifold):

16 tets: Found 100 bit shapes in 0.007s
16 tets: Found 1000 bit shapes in 0.015s
16 tets: Found 10000 bit shapes in 0.219s
16 tets: overall time 0.242s
20 tets: Found 100 bit shapes in 0.009s
20 tets: Found 1000 bit shapes in 0.029s
20 tets: Found 10000 bit shapes in 0.474s
20 tets: overall time 0.514s
34 tets: Found 100 bit shapes in 0.023s
34 tets: Found 1000 bit shapes in 0.072s
34 tets: Found 10000 bit shapes in 1.247s
34 tets: overall time 1.344s
48 tets: Found 100 bit shapes in 0.041s
48 tets: Found 1000 bit shapes in 0.132s
48 tets: Found 10000 bit shapes in 1.787s
48 tets: overall time 1.962s
66 tets: Found 100 bit shapes in 0.117s
66 tets: Found 1000 bit shapes in 0.413s
66 tets: Found 10000 bit shapes in 7.712s
66 tets: overall time 8.247s
Overall: 12.308403015136719

M2 (Sage, verified, no 90-tet manifold):

16 tets: verified 100 bit shapes in 0.159s
16 tets: verified 1000 bit shapes in 0.093s
16 tets: verified 10000 bit shapes in 0.442s
16 tets: overall time 0.700s
20 tets: verified 100 bit shapes in 0.026s
20 tets: verified 1000 bit shapes in 0.049s
20 tets: verified 10000 bit shapes in 0.720s
20 tets: overall time 0.797s
34 tets: verified 100 bit shapes in 0.124s
34 tets: verified 1000 bit shapes in 0.169s
34 tets: verified 10000 bit shapes in 1.589s
34 tets: overall time 1.886s
48 tets: verified 100 bit shapes in 0.190s
48 tets: verified 1000 bit shapes in 0.358s
48 tets: verified 10000 bit shapes in 3.895s
48 tets: overall time 4.448s
66 tets: verified 100 bit shapes in 0.339s
66 tets: verified 1000 bit shapes in 0.547s
66 tets: verified 10000 bit shapes in 6.484s
66 tets: overall time 7.376s
Overall: 15.208s

"""


import snappy
import time
snappy.pari.allocatemem(2**25, 2**27, silent=True)

iso16 = 'qLLvLLzAPQQkcehjlnhnopommpponiiaciijpxlanfggfk_baBaDBaBbB'

iso20 = 'uLLLLvPzMLAPQQccefemllnkmqorptstsrtsiitditpatgevpcmppuggn_abBaBbHg'

iso34 = 'ILLALLMvvAzLLPPMMvzAQQQkbcdeghiklqoqstwxzvAuxCGBGAEHCDHGFFHtsfxjxajmdobxadrfcdgoaamlloocnkkoks_BaCB'

iso48 = 'WLLLwvzvAAQALzPMQLzLvAPzQMzQzAPLQcdfemopklnqmnntsxwuywxzxACIDGHIHLHKMOMPPPQRTSVUVViceaiocvfvvwfnaamofdvpaikdcrvfbbhopcafaacaaaaaagb_baDB'

iso66 = '-ccbLLLvPAvvQQvvLPPzPLMMLMwvAzAzQvPwLQPLPPAzQPQQcaeafaiajaiakanalalamaqaqasauaxazaAayaEaBaFazaEaHaBaGaJaIaEaOaLaLaKaRaVaUaVaVaQaZaYa3a2a4aWaXaXa6a1a0a1a9a8a5a7a-a-a+abbab+a9a9aabbbbbiiegvgutnfksktaafffaaabaaabffgaqevadwqucaapaaaxdifocaavpiqqdjrcrgea_bBcB'

iso90 = '-cAbLLLLLLLvzvvvPAvLPvAvAAMvzwMwzvLLQzvQAQQQQzwwLPQvMPQQMPQQLQQQcadafahaiamaiapaqaDaravauawaxaKaCaRaBaBaBaMaMaWaIaYa2a0aZa2a9a4aPaQaQaSaVaUaVaabhbgbcb1a1aYa1aYaabbbbbab5a8a6akblbqbebfbdbibkbvbobsbrbsbrbrbtbtbubxbmbqbnbmbobsbpbqbzbybxbxbvbwbwbzbzbicegvivfiuatboftfmsgnqqovadaiqalbkbbcvcmapadqmthskgsvfglmajbfqimiqtqdplxmdvggvjrgsaaeefcwfo_bBRqp'



def main(num_loops=10):
   for i in range(num_loops):
      overall_start = time.time()
      for iso in [iso16, iso20, iso48, iso66, iso90]:
         mfld_start = time.time()
         for prec in [100, 1000, 10000]:
            M = snappy.Manifold(iso)
            start = time.time()
            shapes = M.tetrahedra_shapes('rect', bits_prec=prec)
            #shapes = M.tetrahedra_shapes('rect')
            tets = M.num_tetrahedra()
            curr_time = time.time()
            print(f'{tets} tets: Found {prec} bit shapes in {time.time() - start:.3f}s')
         print(f'{tets} tets: overall time {time.time() - mfld_start:.3f}s')

      print(f'Overall {i}: {time.time() - overall_start}')

def alt():
   for iso in 1*[iso66]:
      for prec in [10000]:
         M = snappy.Manifold(iso)
         start = time.time()
         shapes = M.tetrahedra_shapes('rect', bits_prec=prec)
         # vol = M.volume(bits_prec=prec)
         tets = M.num_tetrahedra()
         print(f'{tets} tets: Found {prec} bit shapes + vol in {time.time() - start:.3f}s')

def increment():
   for iso in 1*[iso66]:
      start = time.time()
      M = snappy.Manifold(iso)
      for prec in [200, 500, 1000, 2000, 5000, 10000]:
         shapes = M.tetrahedra_shapes('rect', bits_prec=prec)

      tets = M.num_tetrahedra()
      # vol = M.volume(bits_prec=prec)
      print(f'{tets} tets: Found {prec} bit shapes + vol in {time.time() - start:.3f}s')


def verified_flint():
   from snappy.verify.shape_certifiers import KrawczykShapeCertifier
   overall_start = time.time()
   for iso in [iso16, iso20, iso34, iso48, iso66]:
      mfld_start = time.time()
      for prec in [100, 1000, 10000]:
         M = snappy.Manifold(iso)
         start = time.time()
         vfr = KrawczykShapeCertifier(M, bits_prec=prec)
         result = 'Found' if vfr.certify() else 'FAILED'
         tets = M.num_tetrahedra()
         print(f'{tets} tets: {result} {prec} bit shapes in {time.time() - start:.3f}s')
      print(f'{tets} tets: overall time {time.time() - mfld_start:.3f}s')
   print(f'Overall: {time.time() - overall_start}')


def verified_sage():
   overall_start = time.time()
   for iso in [iso16, iso20, iso34, iso48, iso66]:
      mfld_start = time.time()
      for prec in [100, 1000, 10000]:
         M = snappy.Manifold(iso)
         start = time.time()
         success, shapes = M.verify_hyperbolicity(bits_prec=prec)
         tets = M.num_tetrahedra()
         curr_time = time.time()
         print(f'{tets} tets: verified {prec} bit shapes in {time.time() - start:.3f}s')
      print(f'{tets} tets: overall time {time.time() - mfld_start:.3f}s')
   print(f'Overall: {time.time() - overall_start:.3f}s')

      

def high_precision():
   """
   ManifoldHP:

   * tets only: M2: 0.4s Intel: 0.35s
   * with vol:  M2: 0.4s Intel: 0.35s

   PARI:

   * tets only: M2: 0.57s Intel: 0.64s
   * with vol:  M2: 0.62s Intel: 0.66s

   FLINT:

   * tets only: M2: 0.35s Intel: 0.24s
   * with vol:  M2: 0.36s Intel: 0.26s

   Sage:

   * tets only: M2: 0.67s   Intel: 0.73s
   * with vol:  M2: 0.71s   Intel: 0.77s 
   """
   for i in range(5):
      overall_start = time.time()
      for iso in [iso16, iso20, iso48, iso66, iso90]:
         for prec in [200]:
            mfld_start = time.time()
            # M = snappy.ManifoldHP(iso)
            M = snappy.Manifold(iso)
            #shapes = M.tetrahedra_shapes('rect', bits_prec=prec)
            vol = M.volume(bits_prec=prec)
            #M.tetrahedra_shapes()
            # M.volume()
            tets = M.num_tetrahedra()
            print(f'{tets} tets: Found {prec} bit shapes in {time.time() - mfld_start:.3f}s')

         print(f'{tets} tets: overall time {time.time() - mfld_start:.3f}s')

      print(f'Overall {time.time() - overall_start}')
         
# main(1)
# alt()
# increment()
high_precision()
# verified_flint()
# verified_sage()
