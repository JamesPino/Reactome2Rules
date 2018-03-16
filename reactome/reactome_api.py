import textwrap
import networkx as nx
import pandas as pd
import pygraphviz as pyg
from bioservices.services import REST


pd.set_option('expand_frame_repr', False)
pd.set_option('display.width', 10000)
pd.set_option('display.max_rows', 5000)
pd.set_option('display.max_columns', 5000)
pd.options.display.max_colwidth = 5000


class Reactome(REST):
    _url = 'http://reactome.org/ContentService/data/'

    def __init__(self, cache=False):
        super(Reactome, self).__init__("Reactome(URL)", url=Reactome._url,
                                       verbose="ERROR", cache=cache)
        # http://reactome.org/ContentService/
        self._content_url = 'http://reactome.org/ContentService/data/'

    def get_reaction_info(self, reaction):
        q = 'query/{}'.format(reaction)
        return self._get(q, 'json')

    def get_reaction_detail(self, reaction):
        q = 'query/{}'.format(reaction)
        return self._get(q, 'json')

    def get_entity_info(self, entity):
        q = 'query/enhanced/{}/'.format(entity)
        return self._get(q, 'json')

    def get_entity_attribute(self, entity, attribute):
        q = 'query/{}/{}'.format(entity, attribute)
        return self._get(q, 'string')

    def get_event_participants(self, event):
        q = 'participants/{}/participants'.format(event)
        return self._get(q, 'json')

    def get_event_participating_phys_entities(self, event):
        q = 'participants/{}/participatingPhysicalEntities'.format(event)
        return self._get(q, 'json')

    def get_pathway(self, pathway_id):
        q = 'pathway/{}/Complex'.format(pathway_id)
        return self._get(q, 'json')

    def get_pathway_events(self, pathway_id):
        q = 'pathway/{}/containedEvents'.format(pathway_id)
        return self._get(q, 'json')

    def get_schema(self, schema):
        q = '/schema/{}?species=9606'.format(schema)

        return self._get(q, 'json')
    def _get(self, q, fmt):
        return self.http_get(q, frmt=fmt)


_reactome = Reactome()

verbose = True


name_dict = {
    'plasma membrane':                   'PlasmaMem',
    'cytosol':                           'CYTO',
    'mitochondrial outer membrane':      'MOM',
    'nucleoplasm':                       'NucPlasm',
    'endosome membrane':                 'EndosomeMem',
    'mitochondrial intermembrane space': 'MIM',
    'nuclear envelope':                  'NucEnv',
    'endoplasmic reticulum membrane':    'ERM',
    '[':                                 r'\n['
}


shapes = {
    'Protein':               'box',
    'Chemical Compound':     'oval',
    'Complex':               'box',
    'Set':                   'rectangle',
    'OtherEntity':           'note',
    'Genes and Transcripts': 'note',
    'DNA Sequence':          'note',
    'RNA Sequence':          'note'
}


def shorten_name(string):
    for s in name_dict:
        if s in string:
            string = string.replace(s, name_dict[s])
    return string


def extract_list_from_key(keyword, y):
    input_list = []
    for i in y[keyword]:
        if isinstance(i, int):
            input_list.append(i)
        else:
            if 'dbId' in i.keys():
                input_list.append(i['dbId'])
            else:
                if verbose:
                    print("No dbID and not an int")
    return input_list


def get_entity_info(species):
    """
    Gather information about species

    Parameters
    ----------
    species : str

    Returns
    -------
    information_dict : dict
        Information about species

    """

    entity = _reactome.get_entity_info(species)
    c_name = entity['className']

    if c_name == 'Reaction':
        dn = entity['displayName']
        return {'displayName': dn,
                'species_type': c_name}

    if isinstance(entity, int):
        print(entity)
        print("Here lies an error {}".format(species))
        return dict()
    info_needed = ['referenceType', 'className', 'startCoordinate',
                   'endCoordinate', 'referenceEntity', 'compartment',
                   'referenceEntity', 'displayName', 'hasModifiedResidue']

    dont_need = ['inDisease', 'name', 'dbId', 'stId', 'inferredTo',
                 'speciesName', 'species', 'schemaClass', 'isChimeric',
                 'literatureReference']

    need_to_investigate = ['hasComponent', 'literatureReference',
                           'hasMember', 'summation', 'hasCandidate']

    entity_info = dict()
    entity_info['species_type'] = c_name
    for n in entity:
        if n in info_needed:
            entity_info[n] = entity[n]
    if 'referenceEntity' in entity:
        info = entity['referenceEntity']
        if 'dbId' in info:
            entity_info['parent_dbid'] = info['dbId']
        if 'databaseName' in info:
            entity_info['parent_db_name'] = info['databaseName']
        if 'identifier' in info:
            entity_info['parent_identifier'] = info['identifier']

    if 'hasModifiedResidue' in entity:
        info = entity['hasModifiedResidue']
        mod_residues = []
        mods = []
        for r in info:
            if 'coordinate' in r:
                mod_residues.append(r['coordinate'])
            if 'psiMod' in r:
                ms = r['psiMod']
                if isinstance(ms, int):
                    mods.append(ms)
                elif isinstance(ms, dict):
                    if 'dbId' in ms:
                        mods.append(ms['dbId'])
                    else:
                        for m_type in ms:
                            print("hasModifiedResidue", m_type, ms[m_type])

    return entity_info


def parse_reaction():
    # get reactants, products
    pass


def parse_entity():
    # get name, compartment, parent molecule, modified residue,
    # start and end res
    pass


def add_to_graph(sample, list_of_species, graph):
    if type(sample) == int:
        return
    if sample['className'] not in shapes:
        if verbose:
            print('Is a {}'.format(sample['className']))
        return
    name = sample['dbId']
    display = shorten_name(sample['displayName'])
    graph.add_node(name,
                   label=display,
                   displayName=display,
                   dbId=name,
                   shape=shapes[sample['className']])

    list_of_species[name] = display


def create_graph(reaction_info,  graph):
    inputs = reaction_info['inputs']
    outputs = reaction_info['outputs']
    catalyst = reaction_info['catalyst']
    rxn_name = reaction_info['id']
    prev_events = reaction_info['prev_events']

    for name in inputs:
        current_input = inputs[name]
        label = current_input['displayname']
        shape = shapes[current_input['species_type']]
        graph.add_node(name,
                       label=label,
                       displayName=label,
                       dbId=name,
                       shape=shape)
        if catalyst:
            if name == catalyst:
                pass
            else:
                graph.add_edge(catalyst, rxn_name)
        graph.add_edge(name, rxn_name)
        # for i in prev_events:
        #     graph.add_edge(i, name)

    for each in outputs:
        current_output = inputs[name]
        label = current_output['displayname']
        shape = shapes[current_output['species_type']]
        graph.add_node(name,
                       label=label,
                       displayName=label,
                       dbId=name,
                       shape=shape)
        graph.add_edge(rxn_name, each)


def _extract_info_from_reactants_or_products(r_dict, in_out):

    return_info = dict()
    for each in r_dict[in_out]:
        if isinstance(each, int):
            if each in return_info:
                return_info[each]['counter'] += 1
                continue
            else:
                if verbose:
                    print(each, "Not in inputs")
                continue
        rxn_id = each['dbId']
        dict_of_info = dict()
        dict_of_info['id'] = rxn_id
        dict_of_info['counter'] = 1
        dict_of_info['species_type'] = each['className']
        dict_of_info['displayname'] = shorten_name(each['displayName'])
        return_info[rxn_id] = dict_of_info
    return return_info


def get_reaction_info(reaction):

    if 'className' not in reaction:
        return None

    # ensure class is a reaction
    if reaction['className'] != 'Reaction':
        if verbose:
            print("Not a reaction")
            print(reaction)
        return None

    # get all reaction info
    rxn_name = reaction['dbId']
    y = _reactome.get_reaction_info(rxn_name)

    # check to see if compartment exists
    if 'compartment' not in y:
        if verbose:
            print("No compartment")
        return None

    # check to make sure there is input and output
    if ('input' or 'output') not in y:
        if verbose:
            print("No input or output")
        return None

    # get reaction id and name
    rxn_display_name = reaction['displayName']

    if verbose:
        print("Reaction {} : {} ".format(rxn_name, rxn_display_name))
        for i in y:
            print("\t{} : {}".format(i, y[i]))

    # get compartment info
    comp = y['compartment']
    compartments = []

    for c in comp:
        if 'name' in c:
            compartments.append(c['name'])
    prev_events = []
    if 'precedingEvent' in y:
        prec_event = y['precedingEvent']
        for i in prec_event:
            if 'dbId' in i:
                prev_event_id = i['dbId']
                prev_events.append(prev_event_id)

    catalyst = False
    cat_all = []
    if 'catalystActivity' in y:
        catalyst = y['catalystActivity'][0]['dbId']
        if verbose:
            print(catalyst, 1)
        test = _reactome.get_event_participating_phys_entities(catalyst)
        for n in test:
            if verbose:
                for j in n:
                    print('\t\t {} : {}'.format(j, n[j]))
            if 'displayName' in n:
                catalyst = n['dbId']
                if verbose:
                    print(catalyst, 2)
                cat_all.append(catalyst)

    inputs = _extract_info_from_reactants_or_products(y, 'input')
    outputs = _extract_info_from_reactants_or_products(y, 'output')

    if verbose:
        print('\tInputs : {}'.format(inputs))
        print('\tOutputs : {}'.format(outputs))

    reaction_info = dict()
    reaction_info['inputs'] = inputs
    reaction_info['outputs'] = outputs
    reaction_info['name'] = rxn_display_name
    reaction_info['id'] = rxn_name
    reaction_info['compartment'] = compartments
    reaction_info['catalyst'] = catalyst
    reaction_info['prev_events'] = prev_events

    if len(cat_all) > 1:
        if verbose:
            print(cat_all)
        quit()

    return reaction_info


def extract_pathways(pathway_events):
    pathways = []
    for i in pathway_events:
        if verbose:
            print(i)
        if i['className'] == 'Pathway':
            pathways.append(i['dbId'])
    return pathways


def entity_info(entity_id):
    entity_dict = _reactome.get_entity_info(entity_id)

    for key, value in entity_dict.items():
        print(key, value)

    if 'schemaClass' in entity_dict:
        if entity_dict['schemaClass'] == 'EntityWithAccessionedSequence':
            reference = entity_dict['referenceEntity']
            if 'geneName' in reference:
                ref_gene_name = reference['geneName']
                if len(ref_gene_name) != 1:
                    print("More than one reference gene")
                else:
                    ref_gene_name = ref_gene_name[0]
    if 'startCoordinate' in entity_dict:
        start_aa = entity_dict['startCoordinate']
    if 'endCoordinate' in entity_dict:
        end_aa = entity_dict['endCoordinate']

    if 'compartment' in entity_dict:
        compartment = entity_dict['compartment']
        if len(compartment) > 1:
            print("More than one compartment for this species")
            all_compartments = []
            for c in compartment:
                compartment_name = c['displayName']
                all_compartments.append(compartment_name)
        else:
            compartment_name = compartment[0]['displayName']
    print("Compartment: {}".format(compartment_name))
    print("Reference gene: {}".format(ref_gene_name))
    print("Protein AA range = {}-{}".format(start_aa, end_aa))


def generate_network(ref_id, save_name):

    def extract_list_of_events(ref_id):
        found_events = _reactome.get_pathway_events(ref_id)
        for i in found_events:
            print(i['displayName'], i['dbId'], i['className'])
        # quit()
        print("{} has {} events".format(save_name, len(found_events)))
        print("Extracting events")
        all_events = []
        for i in found_events:
            info = get_reaction_info(i)
            if info is None:
                continue
            all_events.append(info)
        return all_events

    events = extract_list_of_events(ref_id)
    g = pyg.AGraph(directed=True)
    g.graph_attr['splines'] = 'true'

    for e in events:
        if isinstance(e, dict):
            create_graph(reaction_info=e, graph=g)

    for i in g.nodes():
        node = g.get_node(i)
        # print("Node and attributes = {} : {} ".format(i, node.attr))

        if 'displayName' not in node.attr or len(dict(node.attr)) == 0:
            ent_dict = get_entity_info(i)
            # print(ent_dict)
            if len(ent_dict) == 0:
                continue
            lab = ent_dict['displayName']

            if ent_dict['species_type'] != 'Reaction':
                node.attr['shape'] = shapes[ent_dict['species_type']]
                node.attr['label'] = shorten_name(lab)
            else:
                node.attr['label'] = '\n'.join(textwrap.wrap(str(lab), 20))
        else:
            lab = node.attr['label']
            node.attr['label'] = '\n'.join(textwrap.wrap(str(lab), 20))

    gml_g = nx.nx_agraph.from_agraph(g)
    nx.write_gml(gml_g, '{}.gml'.format(save_name))
    print("Created GML file!")
    g.write('{}.dot'.format(save_name), )
    g.draw('{}.png'.format(save_name), prog='dot')

    # df = pd.DataFrame(events)
    # df.to_csv('reactome_df_{}.csv'.format(save_name))
    # df = pd.read_csv('reactome_df.csv')
    # if 'compartment' in df.columns:
    #     translocations = df[df['compartment'].map(len) > 1]
    #     if len(translocations) > 0:
    #         print("Translocations found")
    #         print(translocations)
    #         print(translocations.shape)