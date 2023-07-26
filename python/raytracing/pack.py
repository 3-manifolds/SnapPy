def pack_tet_data(prefix, tets_to_data):
    result = {}

    offset = 0
    offsets = [ offset ]

    for data in tets_to_data:
        for datum in data:
            if result:
                for k, (t, v) in datum.items():
                    result[prefix + k][1].append(v)
            else:
                for k, (t, v) in datum.items():
                    result[prefix + k] = (t + '[]', [ v ])
        offset += len(data)
        offsets.append(offset)

    result[prefix + 'Offsets'] = ('int[]', offsets)

    return result, offset

                
