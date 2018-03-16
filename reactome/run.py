from reactome_api import generate_network

if __name__ == '__main__':
    # large example
    # generate_network(109581, 'all_apoptosis')

    generate_network(5693548, 'double_stranded_break')
    # generate_network(5693532, 'double_stranded_break')
    # generate_network(109581, 'intrinsic')
    quit()
    generate_network(114294, 'bid_activates_bax')

    # pathway 111457 shows an example of needing to include upstream events
    # generate_network(111457, 'release_smac_cycs_from_mito')

    # pathway 111453 shows an example of needing to use GeneSet for BH3 proteins
    generate_network(111453, 'bcl2_interactions')
