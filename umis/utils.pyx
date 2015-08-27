
def weigh_evidence(aux_list):
    for aux_tag in aux_list:
            if aux_tag[0] == 'NH':
                return 1. / aux_tag[1]