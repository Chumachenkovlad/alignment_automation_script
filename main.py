from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import os
import time

"""
qblast(program, database, sequence, url_base='https://blast.ncbi.nlm.nih.gov/Blast.cgi', auto_format=None, composition_based_statistics=None, db_genetic_code=None, endpoints=None, entrez_query='(none)', expect=10.0, filter=None, gapcosts=None, genetic_code=None, hitlist_size=50, i_thresh=None, layout=None, lcase_mask=None, matrix_name=None, nucl_penalty=None, nucl_reward=None, other_advanced=None, perc_ident=None, phi_pattern=None, query_file=None, query_believe_defline=None, query_from=None, query_to=None, searchsp_eff=None, service=None, threshold=None, ungapped_alignment=None, word_size=None, alignments=500, alignment_view=None, descriptions=500, entrez_links_new_window=None, expect_low=None, expect_high=None, format_entrez_query=None, format_object=None, format_type='XML', ncbi_gi=None, results_file=None, show_overview=None, megablast=None)
"""

organism_pattern = r'.*?\[(.*)].*'
INDATA = {
    'gastro': {
        'dir_saved': 'saved_data_gastro',
        'result_file_name': 'result_gastro',
        'orgaism_file': 'organisms_gastro'
    },
    'skin': {
        'dir_saved': 'saved_data_skin',
        'result_file_name': 'result_skin',
        'orgaism_file': 'organisms_skin'
    },
    'eye': {
        'dir_saved': 'saved_data_vagina_eye',
        'result_file_name': 'result_eye',
        'orgaism_file': 'organisms_eye'
    },
    'vagina': {
        'dir_saved': 'saved_data_vagina_eye',
        'result_file_name': 'result_vagina',
        'orgaism_file': 'organisms_vagina'
    }
}
part = 'skin'
DIR_SAVED = INDATA[part]['dir_saved']
RESULT_FILE_NAME = INDATA[part]['result_file_name']
ORGANISMS_FILE = INDATA[part]['orgaism_file']
CONDITIONS = {
    'identities': 25,
    'e-value': 0.0001,
    'block-word': 'hypothetical',
    'min-length': 100
}
DZ_CONDITIONS = {
    'identities': 18,
    'e-value': 0.05,
    'block-word': 'hypothetical',
    'min-length': 100
}

MAM_GIDS = {
    "MamA": "AAL09996.1",
    "MamB": "AAL09999.1",
    "MamM": "CDK99590.1",
    "MamO": "CDK99588.1",
    "MamE": "CDK99594.1",
    "MamN": "CDK99589.1",
    "MamK": "CDK99592.1",
    "MamH": "CDK99596.1"
}

parsedData = {
    "MamA": [],
    "MamB": [],
    "MamM": [],
    "MamO": [],
    "MamE": [],
    "MamN": [],
    "MamK": [],
    "MamH": [],
}

goalOrganisms = ['Acinetobacter spp']
all_unique_organisms = set()
all_organisms = []
alignments = []


def get_goal_organisms():
    with open(ORGANISMS_FILE) as organisms_source:
        # print(organisms_source.read().split('\n'))
        return organisms_source.read().split('\n')
        # return goalOrganisms


def get_data(gid, goal):
    print('Start getting data for {} {}'.format(gid, goal))

    while True:
        try:
            return NCBIWWW.qblast(
                'blastp',
                'nr',
                gid,
                entrez_query=goal)
        except Exception:
            time.sleep(10)
            print('Error: Cannot get data from NCBI, trying again after 10 seconds...')


def save_data(data, filename):
    print('Start saving data to file {}.xml'.format(filename))
    with open("{}/{}.xml".format(DIR_SAVED, filename), 'w') as savedFile:
        savedFile.write(data.read())


def sort_list_of_alignments(alignment_list):
    sorted_list = []
    unique_organisms = set()
    for alg in alignment_list:
        unique_organisms.add(alg['organism'])
    for org in unique_organisms:
        all_organisms.append(org)
    for org in unique_organisms:
        best_alignment = False
        for alg in alignment_list:
            if alg['organism'] == org:
                if not best_alignment:
                    best_alignment = alg
                elif best_alignment['e-value'] < alg['e-value']:
                    best_alignment = alg
        sorted_list.append(best_alignment)
    return sorted_list


def parse_data(mam, goal, filename):
    openfile = open("{}/{}.xml".format(DIR_SAVED, filename))
    blastrecords = NCBIXML.parse(openfile)

    for element in blastrecords:
        for alignment in element.alignments:
            alignment_seqs = []

            for hsp in alignment.hsps:

                if not hsp.identities > CONDITIONS['identities'] \
                        and (CONDITIONS['block-word'] in alignment.title):
                    continue
                if hsp.align_length > DZ_CONDITIONS['min-length'] \
                        and hsp.expect < DZ_CONDITIONS['e-value']:
                    statistic = 'DZ'
                else:
                    continue
                if hsp.align_length > CONDITIONS['min-length'] \
                        and hsp.expect < CONDITIONS['e-value']:
                    statistic = 'H'

                organism, protein = parse_title(alignment.title)

                if not (goal.lower() in organism):
                    continue
                alignment_seqs.append({
                    'statistic': statistic,
                    'protein': protein,
                    'organism': organism,
                    'identities': hsp.identities,
                    'e-value': hsp.expect
                })

            parsedData[mam].append(sort_list_of_alignments(alignment_seqs))
            # parsedData[mam].append(alignment_seqs)
    openfile.close()


def filter_alignments(alignments):
    if len(alignments) == 0:
        return []
    best_alignment = alignments[0]
    for i in range(len(alignments) - 1):
        if alignments[i]['e-value'] <= best_alignment['e-value']:
            best_alignment = alignments[i]['e-value']

    return [best_alignment]


def parse_title(title):
    title = title.split('>')
    title = title[len(title) - 1]
    title = title.split('|')
    title = title[len(title) - 1]
    organism = title[title.find("[") + 1: title.find("]")]
    organism = organism.lower()
    protein = title.replace(title[title.find("["): title.find("]") + 1], '')
    return organism, protein


def get_filename(mam, goal):
    goal_str = '_'.join(goal.split(' '))
    return "{}-{}".format(goal_str, mam)


def get_path_to_file(filename):
    return "{}/{}.xml".format(DIR_SAVED, filename)


def make_mams_aligns(goal):
    for mam in MAM_GIDS:

        print('Start parsing data {} {}'.format(mam, goal))
        filename = get_filename(mam, goal)

        if not os.path.exists(get_path_to_file(filename)):
            time.sleep(5)
            save_data(get_data(MAM_GIDS[mam], goal), filename)

        parse_data(mam, goal, filename)
        print('End parsing data {} {}'.format(mam, goal))


def print_parsed_data(data):
    for mam in data:
        print(mam)
        for item in data[mam]:
            print('\t{}\t{}\t{} ({}%)'.format(item['organism'], item['protein'], item['e-value'], item['identities']))


def string_to_line(string):
    if string:
        return '{}\r\n'.format(str(string))
    return False


def save_final_alignments(data):
    with open('{}.txt'.format(RESULT_FILE_NAME), 'w') as outfile:
        # for mam in data:
        #     outfile.write(string_to_line(mam))
        #     for item in data[mam]:
        #         line = '\t{}\t{}\t{} ({}%)'.format(
        #             item['organism'],
        #             item['protein'],
        #             item['e-value'],
        #             item['identities'])
        #         outfile.write(string_to_line(line))
        for alg in data:
            algn_string = string_to_line(print_alignment(alg))
            if bool(algn_string):
                outfile.write(algn_string)


def make_mam_table(data):

    alignment = {
        "organism_name": '',
        # "statistic": '',
        "MamA": '-',
        "MamB": '-',
        "MamM": '-',
        "MamO": '-',
        "MamE": '-',
        "MamN": '-',
        "MamK": '-',
        "MamH": '-',
    }
    for organism in all_unique_organisms:
        alg = alignment.copy()
        alg['organism'] = organism
        # alignments.append(alg)
        for mam in data:
            for cur_alg in data[mam]:
                for data_alg in cur_alg:
                    if data_alg['organism'] == organism:
                        # alg['statistic'] = data_alg['statistic']
                        alg[mam] = data_alg['statistic'] #'{} ({}%) {}'.format(data_alg['e-value'], data_alg['identities'], data_alg['protein'])
        alignments.append(alg)


def old_print_alignment(a):
    table_string = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
        a["organism"],
        a["MamA"],
        a["MamB"],
        a["MamM"],
        a["MamO"],
        a["MamE"],
        a["MamN"],
        a["MamK"],
        a["MamH"])
    print(table_string)
    return table_string


def print_alignment(a):
    # table_string = '{}\t\t\t{}\t{}\t{}\t{}\t'.format(
    #     a["organism"],
    #     # a["MamA"],
    #     a["MamB"],
    #     a["MamM"],
    #     a["MamO"],
    #     a["MamE"],
    #     # a["MamN"],
    #     # a["MamK"],
    #     # a["MamH"]
    # )
    production_type = False
    statistic = False
    if not [a['MamB'], a['MamM'], a['MamO'], a['MamE']].__contains__('-'):
        production_type = 'A'
        if all(item == 'H' for item in [a['MamB'], a['MamM'], a['MamO'], a['MamE']]):
            statistic = 'H'
        else:
            statistic = 'DZ'
    else:
        return False
    if a['MamA'] != '-':
        production_type = 'C'
        if a['MamA'] == 'H' and not statistic == 'DZ':
            statistic = 'H'
        else:
            statistic = 'DZ'
    if a['MamK'] != '-' and a['MamA'] != '-':
        production_type = 'CC'
        if all(item == 'H' for item in [a['MamK'], a['MamA']]) and not statistic == 'DZ':
            statistic = 'H'
        else:
            statistic = 'DZ'
    table_string = '{}\t{}\t{}'.format(a["organism"], production_type, statistic)
    print(table_string)
    return table_string


def main():
    goal_organisms = get_goal_organisms()
    for goal in goal_organisms:
        print('\n{}/{}\n' \
            .format(
                goal_organisms.index(goal) + 1,
                len(goal_organisms)))
        make_mams_aligns(goal)


if __name__ == '__main__':
    main()
    for organism in all_organisms:
        all_unique_organisms.add(organism)
    make_mam_table(parsedData)

    # print_parsed_data(parsedData)
    save_final_alignments(alignments)
    # get_goal_organisms()
