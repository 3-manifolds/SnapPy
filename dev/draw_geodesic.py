def main():
    from snappy import Manifold

    M = Manifold("m015")

    word = 'aBBAAAAAAbaaaaaabbA'
    
    #word = 'bAAAbaaaB'

    # word = 'aaab'

    #word = 'b'

    gui = M.inside_view(geodesics = [ word ])
    gui.mainloop()

if __name__ == '__main__':

    main()
