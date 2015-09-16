import argparse

def argparser():
    parser = argparse.ArgumentParser(version=0.1)
    parser.add_argument('input', help='The input shapefile')
    parser.add_argument('output', help='The results json file')
    parser.add_argument('-a', '--adjacency', default='QUEEN', dest='adjacency', help='Adjacency Criteria')
    parser.add_argument('dependent_variable', help='The dependent variable field from the shapefile')
    parser.add_argument('independent_variable(s)', nargs='*', help='The independent variable(s) from the shapefile')

    return vars(parser.parse_args())
