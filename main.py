from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import os
import time

"""
qblast(program, database, sequence, url_base='https://blast.ncbi.nlm.nih.gov/Blast.cgi', auto_format=None, composition_based_statistics=None, db_genetic_code=None, endpoints=None, entrez_query='(none)', expect=10.0, filter=None, gapcosts=None, genetic_code=None, hitlist_size=50, i_thresh=None, layout=None, lcase_mask=None, matrix_name=None, nucl_penalty=None, nucl_reward=None, other_advanced=None, perc_ident=None, phi_pattern=None, query_file=None, query_believe_defline=None, query_from=None, query_to=None, searchsp_eff=None, service=None, threshold=None, ungapped_alignment=None, word_size=None, alignments=500, alignment_view=None, descriptions=500, entrez_links_new_window=None, expect_low=None, expect_high=None, format_entrez_query=None, format_object=None, format_type='XML', ncbi_gi=None, results_file=None, show_overview=None, megablast=None)
"""

organism_pattern = r'.*?\[(.*)].*'
DIR_SAVED = 'saved_data'
CONDITIONS = {
    'identities': 18,
    'e-value': 1,
    'block-word': 'hypothetical'
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


def get_goal_organisms():
    with open('organisms') as organisms_source:
        # print(organisms_source.read().split('\n'))
        return organisms_source.read().split('\n')
    # return goalOrganisms


def get_data(gid, goal):
    print('Start getting data for {} {}'.format(gid, goal))
    return NCBIWWW.qblast(
        'blastp',
        'nr',
        gid,
        entrez_query=goal)


def save_data(data, filename):
    print('Start saving data to file {}.xml'.format(filename))
    with open("{}/{}.xml".format(DIR_SAVED, filename), 'w') as savedFile:
        savedFile.write(data.read())


def parse_data(mam, goal, filename):
    openfile = open("{}/{}.xml".format(DIR_SAVED, filename))
    blastrecords = NCBIXML.parse(openfile)

    for element in blastrecords:
        for alignment in element.alignments:
            alignment_seqs = []

            for hsp in alignment.hsps:

                if hsp.identities > CONDITIONS['identities'] and not (CONDITIONS['block-word'] in alignment.title):
                    organism, protein = parse_title(alignment.title)

                    if not (goal.lower() in organism):
                        continue
                    alignment_seqs.append({
                        'protein': protein,
                        'organism': organism,
                        'identities': hsp.identities,
                        'e-value': hsp.expect
                    })
            parsedData[mam] = alignment_seqs[:10]
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
    title = title[len(title)-1]
    title = title.split('|')
    title = title[len(title)-1]
    organism = title[title.find("[")+1: title.find("]")]
    organism = organism.lower()
    protein = title.replace(title[title.find("["): title.find("]")+1], '')
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
    return '{}\r\n'.format(str(string))


def save_final_alignments(data):
    with open('result_v1.txt', 'w') as outfile:
        for mam in data:
            outfile.write(string_to_line(mam))
            for item in data[mam]:
                line = '\t{}\t{}\t{} ({}%)'.format(
                    item['organism'],
                    item['protein'],
                    item['e-value'],
                    item['identities'])
                outfile.write(string_to_line(line))


def main():
    for goal in get_goal_organisms():
        make_mams_aligns(goal)


if __name__ == '__main__':
    main()
    print_parsed_data(parsedData)
    save_final_alignments(parsedData)
    # get_goal_organisms()




