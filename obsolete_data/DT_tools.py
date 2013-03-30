import os
import gzip
linkdir = '../snappy/manifolds/MTLinks'
linkfiles = os.listdir(linkdir)

def letter2int(letter):
    if letter.isupper():
        return 64 - ord(letter)
    elif letter.islower():
        return ord(letter) - 96

def linkname(DTcode, count):
    linkorknot = 'K' if DTcode[1] == 'a' else 'L'
    crossings = letter2int(DTcode[0])
    alternation = 'a' if DTcode.islower() else 'n'
    return '%s%d%s%d'%(linkorknot, crossings, alternation, count)

def sortkey(DTcode):
    components = letter2int(DTcode[1])
    crossings = letter2int(DTcode[0])
    alternation = 0 if DTcode.islower() else 1
    return (components, crossings, alternation, DTcode.lower(), DTcode)

def numeric_DT(DTcode):
    preamble = 2 + letter2int(DTcode[1])
    return [ letter2int(x)<<1 for x in DTcode[preamble:] ]

def all_links():
    linkfiles.sort()
    lines = []
    for linkfile in linkfiles:
        file_obj = gzip.open(os.path.join(linkdir, linkfile))
        lines += file_obj.readlines()
        file_obj.close()
    DTcodes = [line[:-1] for line in lines]
    DTcodes.sort(key=sortkey)
    links = []
    crossings = 4
    isknot = True
    alternating = True
    count = 0
    for code in DTcodes:
        xcrossings = letter2int(code[0])
        xalternating = code.islower()
        xisknot = (code[1] == 'a')
        rewind = False
        if xcrossings != crossings:
            crossings = xcrossings
            rewind = True
        if xalternating != alternating:
            alternating = xalternating
            rewind = True
        if xisknot != isknot:
            isknot = xisknot
            rewind = True
        if rewind:
            count = 1
        else:
            count += 1
        links.append( (code, linkname(code,count)) )
    return links
