def main():
    from snappy import Manifold

    M = Manifold("m015")

    #word0 = 'aBBAAAAAAbaaaaaabbA'
    
    #word = 'bAAAbaaaB'

    word0 = 'aaab'

    word1 = 'b'

    gui = M.inside_view(geodesics = [ word0, word1 ])
    gui.mainloop()

if __name__ == '__main__':

    main()
