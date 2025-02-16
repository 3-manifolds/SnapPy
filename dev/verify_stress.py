import snappy
import time

iso16 = 'qLLvLLzAPQQkcehjlnhnopommpponiiaciijpxlanfggfk_baBaDBaBbB'

iso20 = 'uLLLLvPzMLAPQQccefemllnkmqorptstsrtsiitditpatgevpcmppuggn_abBaBbHg'

iso34 = 'ILLALLMvvAzLLPPMMvzAQQQkbcdeghiklqoqstwxzvAuxCGBGAEHCDHGFFHtsfxjxajmdobxadrfcdgoaamlloocnkkoks_BaCB'

iso48 = 'WLLLwvzvAAQALzPMQLzLvAPzQMzQzAPLQcdfemopklnqmnntsxwuywxzxACIDGHIHLHKMOMPPPQRTSVUVViceaiocvfvvwfnaamofdvpaikdcrvfbbhopcafaacaaaaaagb_baDB'

iso66 = '-ccbLLLvPAvvQQvvLPPzPLMMLMwvAzAzQvPwLQPLPPAzQPQQcaeafaiajaiakanalalamaqaqasauaxazaAayaEaBaFazaEaHaBaGaJaIaEaOaLaLaKaRaVaUaVaVaQaZaYa3a2a4aWaXaXa6a1a0a1a9a8a5a7a-a-a+abbab+a9a9aabbbbbiiegvgutnfksktaafffaaabaaabffgaqevadwqucaapaaaxdifocaavpiqqdjrcrgea_bBcB'

iso90 = '-cAbLLLLLLLvzvvvPAvLPvAvAAMvzwMwzvLLQzvQAQQQQzwwLPQvMPQQMPQQLQQQcadafahaiamaiapaqaDaravauawaxaKaCaRaBaBaBaMaMaWaIaYa2a0aZa2a9a4aPaQaQaSaVaUaVaabhbgbcb1a1aYa1aYaabbbbbab5a8a6akblbqbebfbdbibkbvbobsbrbsbrbrbtbtbubxbmbqbnbmbobsbpbqbzbybxbxbvbwbwbzbzbicegvivfiuatboftfmsgnqqovadaiqalbkbbcvcmapadqmthskgsvfglmajbfqimiqtqdplxmdvggvjrgsaaeefcwfo_bBRqp'


def verified_flint():
   """
   Nathan's M2 with python-flint 0.7a
   
   FLINT
   16 tets: verified 100 bit shapes in 0.008s
   16 tets: verified 1000 bit shapes in 0.018s
   16 tets: FAILED 10000 bit shapes in 0.231s
   16 tets: overall time 0.259s
   20 tets: verified 100 bit shapes in 0.013s
   20 tets: verified 1000 bit shapes in 0.033s
   20 tets: FAILED 10000 bit shapes in 0.491s
   20 tets: overall time 0.539s
   34 tets: verified 100 bit shapes in 0.034s
   34 tets: verified 1000 bit shapes in 0.085s
   34 tets: FAILED 10000 bit shapes in 1.275s
   34 tets: overall time 1.396s
   48 tets: verified 100 bit shapes in 0.063s
   48 tets: verified 1000 bit shapes in 0.152s
   48 tets: FAILED 10000 bit shapes in 2.012s
   48 tets: overall time 2.230s
   66 tets: verified 100 bit shapes in 0.157s
   66 tets: verified 1000 bit shapes in 0.436s
   66 tets: FAILED 10000 bit shapes in 7.864s
   66 tets: overall time 8.462s
   90 tets: FAILED 100 bit shapes in 0.696s
   90 tets: FAILED 1000 bit shapes in 1.322s
   90 tets: FAILED 10000 bit shapes in 38.476s
   90 tets: overall time 40.501s
   Overall: 53.386476039886475
   """
   from snappy.verify.shape_certifiers import KrawczykShapeCertifier
   overall_start = time.time()
   for iso in [iso16, iso20, iso34, iso48, iso66, iso90]:
      mfld_start = time.time()
      for prec in [100, 1000, 10000]:
         M = snappy.Manifold(iso)
         start = time.time()
         vfr = KrawczykShapeCertifier(M, bits_prec=prec)
         result = 'verified' if vfr.certify() else 'FAILED'
         tets = M.num_tetrahedra()
         print(f'{tets} tets: {result} {prec} bit shapes in {time.time() - start:.3f}s')
      print(f'{tets} tets: overall time {time.time() - mfld_start:.3f}s')
   print(f'Overall: {time.time() - overall_start}')


def verified_sage():
   """
   SageMath 10.5 on Nathan's M2:

   SAGEMATH
   16 tets: verified 100 bit shapes in 0.167s
   16 tets: verified 1000 bit shapes in 0.087s
   16 tets: verified 10000 bit shapes in 0.445s
   16 tets: overall time 0.703s
   20 tets: verified 100 bit shapes in 0.027s
   20 tets: verified 1000 bit shapes in 0.057s
   20 tets: verified 10000 bit shapes in 0.665s
   20 tets: overall time 0.751s
   34 tets: verified 100 bit shapes in 0.067s
   34 tets: verified 1000 bit shapes in 0.200s
   34 tets: verified 10000 bit shapes in 1.548s
   34 tets: overall time 1.817s
   48 tets: verified 100 bit shapes in 0.166s
   48 tets: verified 1000 bit shapes in 0.341s
   48 tets: verified 10000 bit shapes in 3.718s
   48 tets: overall time 4.229s
   66 tets: verified 100 bit shapes in 0.303s
   66 tets: verified 1000 bit shapes in 0.592s
   66 tets: verified 10000 bit shapes in 6.193s
   66 tets: overall time 7.092s
   90 tets: verified 100 bit shapes in 0.648s
   90 tets: verified 1000 bit shapes in 1.710s
   90 tets: verified 10000 bit shapes in 29.886s
   90 tets: overall time 32.252s
   Overall: 46.843s
   """
   overall_start = time.time()
   for iso in [iso16, iso20, iso34, iso48, iso66, iso90]:
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

      
if snappy.sage_helper._within_sage:
   print('SAGEMATH')
   verified_sage()
else:
   print('FLINT')
   verified_flint()
