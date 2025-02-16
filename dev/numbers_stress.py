import snappy
import time

iso16 = 'qLLvLLzAPQQkcehjlnhnopommpponiiaciijpxlanfggfk_baBaDBaBbB'

iso20 = 'uLLLLvPzMLAPQQccefemllnkmqorptstsrtsiitditpatgevpcmppuggn_abBaBbHg'

iso34 = 'ILLALLMvvAzLLPPMMvzAQQQkbcdeghiklqoqstwxzvAuxCGBGAEHCDHGFFHtsfxjxajmdobxadrfcdgoaamlloocnkkoks_BaCB'

iso48 = 'WLLLwvzvAAQALzPMQLzLvAPzQMzQzAPLQcdfemopklnqmnntsxwuywxzxACIDGHIHLHKMOMPPPQRTSVUVViceaiocvfvvwfnaamofdvpaikdcrvfbbhopcafaacaaaaaagb_baDB'

iso66 = '-ccbLLLvPAvvQQvvLPPzPLMMLMwvAzAzQvPwLQPLPPAzQPQQcaeafaiajaiakanalalamaqaqasauaxazaAayaEaBaFazaEaHaBaGaJaIaEaOaLaLaKaRaVaUaVaVaQaZaYa3a2a4aWaXaXa6a1a0a1a9a8a5a7a-a-a+abbab+a9a9aabbbbbiiegvgutnfksktaafffaaabaaabffgaqevadwqucaapaaaxdifocaavpiqqdjrcrgea_bBcB'

iso90 = '-cAbLLLLLLLvzvvvPAvLPvAvAAMvzwMwzvLLQzvQAQQQQzwwLPQvMPQQMPQQLQQQcadafahaiamaiapaqaDaravauawaxaKaCaRaBaBaBaMaMaWaIaYa2a0aZa2a9a4aPaQaQaSaVaUaVaabhbgbcb1a1aYa1aYaabbbbbab5a8a6akblbqbebfbdbibkbvbobsbrbsbrbrbtbtbubxbmbqbnbmbobsbpbqbzbybxbxbvbwbwbzbzbicegvivfiuatboftfmsgnqqovadaiqalbkbbcvcmapadqmthskgsvfglmajbfqimiqtqdplxmdvggvjrgsaaeefcwfo_bBRqp'

def eval_gluing_equation(eqn, acb_shapes):
    """
    Evaluate the product of cross ratios in an edge equation.
    The result will be 1 if the equation is satisfied.
    """
    a, b, c = eqn
    ans = c
    for i, z in enumerate(acb_shapes):
        ans = ans * (z**int(a[i]) * (1 - z)**int(b[i]))
    return ans


def gluing_equation_errors(eqns, acb_shapes):
    """ A list containing the difference from 1 for each equation."""
    return [eval_gluing_equation(eqn, acb_shapes) - 1 for eqn in eqns]

def eval_gluing(manifold):
   shapes = manifold.tetrahedra_shapes('rect')
   eqns = manifold.gluing_equations('rect')
   gluing_equation_errors(eqns, shapes) 

def eval_gluing_complex(manifold):
   shapes = manifold.tetrahedra_shapes('rect')
   eqns = manifold.gluing_equations('rect')
   gluing_equation_errors(eqns, [complex(z) for z in shapes]) 
   

def main():
   for iso in [iso16, iso20, iso34, iso48, iso66, iso90]:
      M = snappy.Manifold(iso)
      for i in range(100):
         eval_gluing(M)

def main_complex():
   for iso in [iso16, iso20, iso34, iso48, iso66, iso90]:
      M = snappy.Manifold(iso)
      for i in range(100):
         eval_gluing_complex(M)

main()
# main_complex()
