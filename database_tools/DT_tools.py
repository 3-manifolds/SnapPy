import os
import gzip, snappy
linkdir = 'MTLinks'
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
    # first all the knots, then the links
    components = 0 if letter2int(DTcode[1]) == 1 else 1
    # order by number of crossings
    crossings = letter2int(DTcode[0])
    # alternating before nonalternating
    alternation = 0 if DTcode.islower() else 1
    # sort first by the DT code of the diagram, then the actual DT code
    return (components, crossings, alternation, DTcode.lower(), DTcode)

def numeric_DT(DTcode):
    preamble = 2 + letter2int(DTcode[1])
    return [ letter2int(x)<<1 for x in DTcode[preamble:] ]

def basic_all_links():
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
        xisknot = (code[1] == 'a')
        xalternating = code.islower()
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

def ten_and_eleven_order():
    ten, eleven = [], []
    for M in snappy.NonalternatingKnotExteriors():
        if M.name().startswith('12'):
            break
        if M.name().startswith('10'):
            ten.append(M.DT_code(alpha=True))
        if M.name().startswith('11'):
            eleven.append(M.DT_code(alpha=True))
    return ten, eleven

def all_links():
    """
    The function basic_all_links sorts the DT codes in a natural order.
    Unfortunately, by some historical accident the 10 and 11 crossing
    nonalternating knots are *not* sorted in this or any other consistent
    way.  Therefore, we have to order those according to the order in
    NonalternatingKnotExteriors.
    """
    ten, eleven = ten_and_eleven_order()
    links = basic_all_links()

    for prefix, order in [('jaj', ten), ('kak', eleven)]:
        for n, L in enumerate(links):
            dt = L[0]
            if dt.startswith(prefix) and dt.lower() != dt:
                break

        assert {L[0] for L in links[n:n+len(order)]} == set(order)
        for j in range(len(order)):
            links[n+j] = (order[j], links[n+j][1])

    return links

    
# Getting DT codes of Rolfsen links from the Christy table

joes_links = '/Users/dunfield/work/work/joes_links/LinkTables/links/'

def munge_name(M):
    cross, index = M.name().split('_')
    return 'l%d%02d%03d' % (M.num_cusps(), int(cross.split('^')[0]), int(index))

def compute_DTs_of_Christy():
    import snappy, plink
    LE = plink.LinkEditor()
    file = gzip.open('ChristyDT.gz', 'w')
    for M in snappy.LinkExteriors:
        LE.load(joes_links + munge_name(M))
        file.write(M.name() + '   ' + LE.DT_code(alpha=True)[0] + '\n')
